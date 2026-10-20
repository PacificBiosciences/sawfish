// Copyright (c) 2020 Robert Vaser

#include "sisd_alignment_engine.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>

#include "spoa/graph.hpp"

#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#endif

#if defined(__aarch64__)
#include <arm_neon.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#if defined(__SSE4_1__)
#include <smmintrin.h>
#endif
#endif

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define SPOA_LOW_MEMORY_AVX2
#endif

namespace spoa {

constexpr std::int32_t kNegativeInfinity =
    std::numeric_limits<std::int32_t>::min() + 1024;

std::unique_ptr<AlignmentEngine> SisdAlignmentEngine::Create(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
  return std::unique_ptr<AlignmentEngine>(
      new SisdAlignmentEngine(type, subtype, m, n, g, e, q, c));
}

// Buffers for LinearPrunedT, kept between alignments to avoid reallocating them
template <typename T>
struct PrunedBuffers {
  std::vector<T> vals;
  std::vector<T> row0;
  // Full width rows, with -infinity outside unpruned cells, and the free list of row indexes
  std::vector<std::vector<T>> rows;
  std::vector<std::uint32_t> free_rows;
};

struct SisdAlignmentEngine::Implementation {
  PrunedBuffers<std::int16_t> pruned16;
  PrunedBuffers<std::int32_t> pruned32;
  std::vector<std::uint32_t> node_id_to_rank;
  std::vector<std::int32_t> sequence_profile;
  std::vector<std::int16_t> sequence_profile16;
  std::vector<std::int32_t> M;
  std::int32_t* H;
  std::int32_t* F;
  std::int32_t* E;
  std::int32_t* O;
  std::int32_t* Q;

  Implementation()
      : pruned16(),
        pruned32(),
        node_id_to_rank(),
        sequence_profile(),
        sequence_profile16(),
        M(),
        H(nullptr),
        F(nullptr),
        E(nullptr),
        O(nullptr),
        Q(nullptr) {
  }
};

SisdAlignmentEngine::SisdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c)
    : AlignmentEngine(type, subtype, m, n, g, e, q, c),
      pimpl_(new Implementation()) {
}

void SisdAlignmentEngine::Prealloc(
    std::uint32_t max_sequence_len,
    std::uint8_t alphabet_size) {
  if (max_sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Prealloc] error: too large sequence!");
  }
  try {
    Realloc(
        static_cast<std::uint64_t>(max_sequence_len) + 1,
        static_cast<std::uint64_t>(max_sequence_len) * alphabet_size + alphabet_size,  // NOLINT
        alphabet_size);
  } catch (std::bad_alloc& ba) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Prealloc] error: insufficient memory!");
  }
}

void SisdAlignmentEngine::Realloc(
    std::uint64_t matrix_width,
    std::uint64_t matrix_height,
    std::uint8_t num_codes) {
  if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
    pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
  }
  if (pimpl_->sequence_profile.size() < num_codes * matrix_width) {
    pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
  }
  if (subtype_ == AlignmentSubtype::kLinear) {
    if (pimpl_->M.size() < matrix_height * matrix_width) {
      pimpl_->M.resize(matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = nullptr;
      pimpl_->E = nullptr;
    }
  } else if (subtype_ == AlignmentSubtype::kAffine) {
    if (pimpl_->M.size() < 3 * matrix_height * matrix_width) {
      pimpl_->M.resize(3 * matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = pimpl_->H + matrix_width * matrix_height;
      pimpl_->E = pimpl_->F + matrix_width * matrix_height;
    }
  } else if (subtype_ == AlignmentSubtype::kConvex) {
    if (pimpl_->M.size() < 5 * matrix_height * matrix_width) {
      pimpl_->M.resize(5 * matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = pimpl_->H + matrix_width * matrix_height;
      pimpl_->E = pimpl_->F + matrix_width * matrix_height;
      pimpl_->O = pimpl_->E + matrix_width * matrix_height;
      pimpl_->Q = pimpl_->O + matrix_width * matrix_height;
    }
  }
}

void SisdAlignmentEngine::Initialize(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph) noexcept {
  std::uint32_t matrix_width = sequence_len + 1;
  std::uint32_t matrix_height = graph.nodes().size() + 1;

  for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
    char c = graph.decoder(i);
    pimpl_->sequence_profile[i * matrix_width] = 0;
    for (std::uint32_t j = 0; j < sequence_len; ++j) {
      pimpl_->sequence_profile[i * matrix_width + (j + 1)] =
          (c == sequence[j] ? m_ : n_);
    }
  }

  const auto& rank_to_node = graph.rank_to_node();
  for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
    pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
  }

  // initialize secondary matrices
  switch (subtype_) {
    case AlignmentSubtype::kConvex:
      pimpl_->O[0] = 0;
      pimpl_->Q[0] = 0;
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->O[j] = kNegativeInfinity;
        pimpl_->Q[j] = q_ + (j - 1) * c_;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? q_ - c_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->O[pred_i * matrix_width]);
        }
        pimpl_->O[i * matrix_width] = penalty + c_;
        pimpl_->Q[i * matrix_width] = kNegativeInfinity;
      }
      // fall through
    case AlignmentSubtype::kAffine:
      pimpl_->F[0] = 0;
      pimpl_->E[0] = 0;
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->F[j] = kNegativeInfinity;
        pimpl_->E[j] = g_ + (j - 1) * e_;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->F[pred_i * matrix_width]);
        }
        pimpl_->F[i * matrix_width] = penalty + e_;
        pimpl_->E[i * matrix_width] = kNegativeInfinity;
      }
      // fall through
    case AlignmentSubtype::kLinear:
      pimpl_->H[0] = 0;
      break;
    default:
      break;
  }

  // initialize primary matrix
  switch (type_) {
    case AlignmentType::kSW:
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->H[j] = 0;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        pimpl_->H[i * matrix_width] = 0;
      }
      break;
    case AlignmentType::kNW:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = std::max(pimpl_->Q[j], pimpl_->E[j]);
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = std::max(
                pimpl_->O[i * matrix_width],
                pimpl_->F[i * matrix_width]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = pimpl_->F[i * matrix_width];
          }
          break;
        case AlignmentSubtype::kLinear:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = j * g_;
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            const auto& edges = rank_to_node[i - 1]->inedges;
            std::int32_t penalty = edges.empty() ? 0 : kNegativeInfinity;
            for (const auto& it : edges) {
              std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
              penalty = std::max(penalty, pimpl_->H[pred_i * matrix_width]);
            }
            pimpl_->H[i * matrix_width] = penalty + g_;
          }
        default:
          break;
      }
      break;
    case AlignmentType::kOV:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = std::max(pimpl_->Q[j], pimpl_->E[j]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          break;
        case AlignmentSubtype::kLinear:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = j * g_;
          }
          break;
        default:
          break;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        pimpl_->H[i * matrix_width] = 0;
      }
      break;
    default:
      break;
  }
}

Alignment SisdAlignmentEngine::Align(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) {
  if (sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: too large sequence!");
  }

  if (graph.nodes().empty() || sequence_len == 0) {
    return Alignment();
  }

  if (WorstCaseAlignmentScore(sequence_len, graph.nodes().size()) < kNegativeInfinity) {  // NOLINT
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: possible overflow!");
  }

  if (subtype_ == AlignmentSubtype::kLinear &&
      (type_ == AlignmentType::kSW || type_ == AlignmentType::kOV)) {
    Alignment alignment;
    std::int32_t low_memory_score = 0;
    if (LinearLowMemory(sequence, sequence_len, graph, &low_memory_score, &alignment)) {  // NOLINT
      static const bool verify =
          std::getenv("SPOA_VERIFY_LOW_MEMORY_LINEAR") != nullptr;
      if (verify) {
        std::int32_t full_score = 0;
        Alignment full_alignment =
            AlignFullMatrix(sequence, sequence_len, graph, &full_score);
        if (full_alignment != alignment ||
            (!alignment.empty() && full_score != low_memory_score)) {
          std::fprintf(stderr,
              "[spoa::SisdAlignmentEngine::Align] error: low memory linear "
              "alignment mismatch (score %d vs %d, length %zu vs %zu)\n",
              low_memory_score, full_score, alignment.size(),
              full_alignment.size());
          std::abort();
        }
      }
      if (score && !alignment.empty()) {
        *score = low_memory_score;
      }
      return alignment;
    }
  }

  return AlignFullMatrix(sequence, sequence_len, graph, score);
}

Alignment SisdAlignmentEngine::AlignFullMatrix(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) {
  try {
    Realloc(sequence_len + 1, graph.nodes().size() + 1, graph.num_codes());
  } catch (std::bad_alloc& ba) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: insufficient memory!");
  }
  Initialize(sequence, sequence_len, graph);

  if (subtype_ == AlignmentSubtype::kLinear) {
    return Linear(sequence_len, graph, score);
  } else if (subtype_ == AlignmentSubtype::kAffine) {
    return Affine(sequence_len, graph, score);
  } else if (subtype_ == AlignmentSubtype::kConvex) {
    return Convex(sequence_len, graph, score);
  }
  return Alignment();
}

Alignment SisdAlignmentEngine::Linear(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ?
        0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    // update H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      H_row[j] = std::max(
          H_pred_row[j - 1] + char_profile[j],
          H_pred_row[j] + g_);
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        H_row[j] = std::max(
            H_pred_row[j - 1] + char_profile[j],
            std::max(
                H_row[j],
                H_pred_row[j] + g_));
      }
    }

    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      H_row[j] = std::max(H_row[j - 1] + g_, H_row[j]);

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          it->outedges.empty() && j == matrix_width - 1) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && it->outedges.empty()) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;
  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
          pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          std::uint32_t pred_i =
              pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ? 0 :
          pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          std::uint32_t pred_i =
              pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && H_ij == pimpl_->H[i * matrix_width + j - 1] + g_) {  // NOLINT
      prev_i = i;
      prev_j = j - 1;
      predecessor_found = true;
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

// Low memory version of Linear for kSW and kOV alignment
//
// This gives exactly the same alignment and score as Linear, without storing
// the full score matrix:
//
// 1. The forward pass keeps each score row only until all of its successor
//    rows are computed. To support the backtrack, it also stores the scores of
//    every kColumnStride-th column of every row, and the full score row of any
//    row with a successor in a different block of kRankBlockSize rows.
//
// 2. The backtrack follows the same rules as Linear, on scores recomputed for
//    one tile at a time, covering the rows of one block between two stored
//    columns. The predecessors of the tile rows are either in the tile or
//    stored in full, so the tile scores are the same as in the forward pass.
//
namespace {

constexpr std::uint32_t kRankBlockSize = 64;
constexpr std::uint32_t kColumnStride = 64;

// Score buffer which is returned to the OS as soon as it is freed, if large
//
// Memory allocators may keep freed memory resident for reuse. The low memory alignment allocates
// and frees many score rows and other large buffers, so to keep these from accumulating across
// alignments on many threads, large buffers are mapped directly where possible.
//
template <typename T>
class ScoreBuffer {
 public:
  ScoreBuffer() = default;

  explicit ScoreBuffer(std::uint64_t size)
      : bytes_(size * sizeof(T)) {
#if defined(__unix__) || defined(__APPLE__)
    constexpr std::uint64_t kMinMappedBytes = 1 << 16;
    if (bytes_ >= kMinMappedBytes) {
      void* data = mmap(nullptr, bytes_, PROT_READ | PROT_WRITE,
          MAP_PRIVATE | MAP_ANON, -1, 0);
      if (data == MAP_FAILED) {
        throw std::bad_alloc();
      }
      data_ = static_cast<T*>(data);
      is_mapped_ = true;
      return;
    }
#endif
    data_ = new T[size];
  }

  ~ScoreBuffer() {
    Release();
  }

  ScoreBuffer(const ScoreBuffer&) = delete;
  ScoreBuffer& operator=(const ScoreBuffer&) = delete;

  ScoreBuffer(ScoreBuffer&& other) noexcept {
    *this = std::move(other);
  }

  ScoreBuffer& operator=(ScoreBuffer&& other) noexcept {
    if (this != &other) {
      Release();
      data_ = other.data_;
      bytes_ = other.bytes_;
      is_mapped_ = other.is_mapped_;
      other.data_ = nullptr;
      other.bytes_ = 0;
      other.is_mapped_ = false;
    }
    return *this;
  }

  T* get() const { return data_; }
  T& operator[](std::uint64_t i) const { return data_[i]; }
  explicit operator bool() const { return data_ != nullptr; }

 private:
  void Release() {
    if (data_ == nullptr) {
      return;
    }
#if defined(__unix__) || defined(__APPLE__)
    if (is_mapped_) {
      munmap(data_, bytes_);
      data_ = nullptr;
      return;
    }
#endif
    delete[] data_;
    data_ = nullptr;
  }

  T* data_ = nullptr;
  std::uint64_t bytes_ = 0;
  bool is_mapped_ = false;
};

// Value below every score in the low memory alignment, which is only ever compared, never added to
template <typename T>
inline std::int32_t ScanNegativeInfinity();
template <>
inline std::int32_t ScanNegativeInfinity<std::int32_t>() {
  return std::numeric_limits<std::int32_t>::min() / 2;
}
template <>
inline std::int32_t ScanNegativeInfinity<std::int16_t>() {
  return -30000;
}

// Check if all scores of a linear gap alignment of a sequence of this length fit in int16
//
// Every cell can be reached with left moves from column zero, so scores are at least
// sequence_len * g. Each cell consumes at most one more sequence base than its predecessor, so
// scores are at most sequence_len * m. The row computation only offsets scores by up to the vector
// width before taking prefix maximums. A margin is added for the vector width.
//
bool IsInt16ScoreRange(
    std::uint32_t sequence_len,
    std::int32_t m,
    std::int32_t n,
    std::int32_t g) {
  const std::int64_t len = static_cast<std::int64_t>(sequence_len) + 64;
  const std::int64_t abs_g = std::abs(g);
  const std::int64_t max_penalty = std::max<std::int64_t>(std::abs(n), abs_g);
  const std::int64_t lowest = -(len * abs_g + max_penalty);
  const std::int64_t highest = len * std::max(m, 0) + std::abs(m);
  return lowest > ScanNegativeInfinity<std::int16_t>() + 1000 &&
         highest < std::numeric_limits<std::int16_t>::max() - 1000;
}

// Compute score row H from its predecessor rows, returning the maximum score in the row
//
// Each cell is first set to its best diagonal or vertical move score over all predecessor rows:
//
//   H[j] = max over p of (P[j - 1] + profile[j], P[j] + g)
//
// Then the left moves are added, H[j] = max(H[j], H[j - 1] + g) in sequence, with H[0] = 0.
// Subtracting j * g from each cell turns this into a prefix maximum, which is found for each
// vector of cells in registers, then combined with the running maximum over all preceding
// cells. For kSW, clamping each cell to zero in sequence is the same as clamping after the left
// moves are added.
//
template <typename T, bool kIsSW>
std::int32_t ComputeRowScalar(
    const T* const* pred_rows,
    std::uint32_t num_preds,
    const T* __restrict profile,
    std::int32_t g,
    std::uint32_t width,
    T* __restrict H,
    std::uint32_t j,
    std::int32_t running,
    std::int32_t max_score) {
  for (; j < width; ++j) {
    std::int32_t t = ScanNegativeInfinity<std::int32_t>();
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      const T* P = pred_rows[p];
      t = std::max(t, std::max(P[j - 1] + profile[j], P[j] + g));
    }
    const std::int32_t jg = static_cast<std::int32_t>(j) * g;
    running = std::max(running, t - jg);
    std::int32_t h = running + jg;
    if (kIsSW) {
      h = std::max(h, 0);
    }
    // Only cells derived from -infinity predecessors can be below -infinity
    h = std::max<std::int32_t>(h, ScanNegativeInfinity<T>());
    H[j] = static_cast<T>(h);
    max_score = std::max(max_score, h);
  }
  return max_score;
}

template <typename T, bool kIsSW>
std::int32_t ComputeRow(
    const T* const* pred_rows,
    std::uint32_t num_preds,
    const T* __restrict profile,
    std::int32_t g,
    std::uint32_t width,
    T* __restrict H) {
  H[0] = 0;
  return ComputeRowScalar<T, kIsSW>(pred_rows, num_preds, profile, g, width, H, 1, 0,
      ScanNegativeInfinity<T>());
}

// Vectorized ComputeRow for vector type V, with the same results as the scalar version
//
// V provides kLanes scores of type Score, and these operations:
// Load, Store, Set, Add, Sub, Max, Ramp(g) = {g, 2g, ..., kLanes * g}, PrefixMax(x, fill) = the
// running maximum over lanes with fill shifted in, BroadcastLast, First and HorizontalMax.
//
#define SPOA_COMPUTE_ROW_VEC_BODY                                                     \
  using T = typename V::Score;                                                        \
  using VT = typename V::Type;                                                        \
  H[0] = 0;                                                                           \
  const VT fill = V::Set(static_cast<T>(ScanNegativeInfinity<T>()));                  \
  const VT zero = V::Set(0);                                                          \
  const VT gap = V::Set(static_cast<T>(g));                                           \
  const VT ramp = V::Ramp(g);                                                         \
  const VT block_gap = V::Set(static_cast<T>(static_cast<std::int32_t>(V::kLanes) * g)); \
  /* Score of the cell to the left of the current vector, in every lane */            \
  VT carry = zero;                                                                    \
  VT row_max = fill;                                                                  \
  std::uint32_t j = 1;                                                                \
  for (; j + V::kLanes <= width; j += V::kLanes) {                                    \
    const VT char_profile = V::Load(profile + j);                                     \
    const T* P = pred_rows[0];                                                        \
    VT x = V::Max(                                                                    \
        V::Add(V::Load(P + j - 1), char_profile),                                     \
        V::Add(V::Load(P + j), gap));                                                 \
    for (std::uint32_t p = 1; p < num_preds; ++p) {                                   \
      P = pred_rows[p];                                                               \
      x = V::Max(x, V::Max(                                                           \
          V::Add(V::Load(P + j - 1), char_profile),                                   \
          V::Add(V::Load(P + j), gap)));                                              \
    }                                                                                 \
    x = V::PrefixMax(V::Sub(x, ramp), fill);                                          \
    VT h = V::Add(V::Max(x, carry), ramp);                                            \
    if (kIsSW) {                                                                      \
      h = V::Max(h, zero);                                                            \
    }                                                                                 \
    V::Store(H + j, h);                                                               \
    row_max = V::Max(row_max, h);                                                     \
    carry = V::Add(V::Max(carry, V::BroadcastLast(x)), block_gap);                    \
  }                                                                                   \
  /* carry has the same value in every lane */                                        \
  return ComputeRowScalar<T, kIsSW>(pred_rows, num_preds, profile, g, width, H, j,     \
      V::First(carry) - static_cast<std::int32_t>(j - 1) * g, V::HorizontalMax(row_max));

// Range version of SPOA_COMPUTE_ROW_VEC_BODY: computes H[lo..end) given the running maximum
// (H[lo - 1] - (lo - 1) * g) of the cells to the left, without a row maximum
// Range version of SPOA_COMPUTE_ROW_VEC_BODY: computes H[lo..end) given the running maximum
// (H[lo - 1] - (lo - 1) * g) of the cells to the left, without a row maximum. Predecessor scores
// may be -infinity, so additions saturate (scores of valid cells never come close to saturating).
#define SPOA_COMPUTE_ROW_RANGE_VEC_BODY                                               \
  using T = typename V::Score;                                                        \
  using VT = typename V::Type;                                                        \
  const VT fill = V::Set(static_cast<T>(ScanNegativeInfinity<T>()));                  \
  const VT zero = V::Set(0);                                                          \
  const VT gap = V::Set(static_cast<T>(g));                                           \
  const VT ramp = V::Ramp(g);                                                         \
  const VT block_gap = V::Set(static_cast<T>(static_cast<std::int32_t>(V::kLanes) * g)); \
  /* Score of the cell to the left of the current vector, in every lane */            \
  VT carry = V::Set(static_cast<T>(std::max<std::int32_t>(                            \
      running + static_cast<std::int32_t>(lo - 1) * g, ScanNegativeInfinity<T>())));  \
  std::uint32_t j = lo;                                                               \
  for (; j + V::kLanes <= end; j += V::kLanes) {                                      \
    const VT char_profile = V::Load(profile + j);                                     \
    const T* P = pred_rows[0];                                                        \
    VT x = V::Max(                                                                    \
        V::AddSat(V::Load(P + j - 1), char_profile),                                  \
        V::AddSat(V::Load(P + j), gap));                                              \
    for (std::uint32_t p = 1; p < num_preds; ++p) {                                   \
      P = pred_rows[p];                                                               \
      x = V::Max(x, V::Max(                                                           \
          V::AddSat(V::Load(P + j - 1), char_profile),                                \
          V::AddSat(V::Load(P + j), gap)));                                           \
    }                                                                                 \
    x = V::PrefixMax(V::Sub(x, ramp), fill);                                          \
    VT h = V::AddSat(V::Max(x, carry), ramp);                                         \
    if (kIsSW) {                                                                      \
      h = V::Max(h, zero);                                                            \
    }                                                                                 \
    V::Store(H + j, h);                                                               \
    carry = V::AddSat(V::Max(carry, V::BroadcastLast(x)), block_gap);                 \
  }                                                                                   \
  ComputeRowScalar<T, kIsSW>(pred_rows, num_preds, profile, g, end, H, j,             \
      V::First(carry) - static_cast<std::int32_t>(j - 1) * g, 0);

#if defined(__aarch64__) || defined(__SSE2__)

#if defined(__aarch64__)

struct VecI32 {
  using Score = std::int32_t;
  using Type = int32x4_t;
  static constexpr std::uint32_t kLanes = 4;
  static Type Load(const Score* p) { return vld1q_s32(p); }
  static void Store(Score* p, Type x) { vst1q_s32(p, x); }
  static Type Set(Score x) { return vdupq_n_s32(x); }
  static Type Add(Type a, Type b) { return vaddq_s32(a, b); }
  static Type AddSat(Type a, Type b) { return vaddq_s32(a, b); }
  static Type Sub(Type a, Type b) { return vsubq_s32(a, b); }
  static Type Max(Type a, Type b) { return vmaxq_s32(a, b); }
  static Type Ramp(std::int32_t g) {
    const Score x[kLanes] = {g, 2 * g, 3 * g, 4 * g};
    return vld1q_s32(x);
  }
  static Type PrefixMax(Type x, Type fill) {
    x = Max(x, vextq_s32(fill, x, 3));
    return Max(x, vextq_s32(fill, x, 2));
  }
  static Type BroadcastLast(Type x) { return vdupq_laneq_s32(x, 3); }
  static std::int32_t First(Type x) { return vgetq_lane_s32(x, 0); }
  static std::int32_t HorizontalMax(Type x) { return vmaxvq_s32(x); }
};

struct VecI16 {
  using Score = std::int16_t;
  using Type = int16x8_t;
  static constexpr std::uint32_t kLanes = 8;
  static Type Load(const Score* p) { return vld1q_s16(p); }
  static void Store(Score* p, Type x) { vst1q_s16(p, x); }
  static Type Set(Score x) { return vdupq_n_s16(x); }
  static Type Add(Type a, Type b) { return vaddq_s16(a, b); }
  static Type AddSat(Type a, Type b) { return vqaddq_s16(a, b); }
  static Type Sub(Type a, Type b) { return vsubq_s16(a, b); }
  static Type Max(Type a, Type b) { return vmaxq_s16(a, b); }
  static Type Ramp(std::int32_t g) {
    Score x[kLanes];
    for (std::uint32_t i = 0; i < kLanes; ++i) {
      x[i] = static_cast<Score>((i + 1) * g);
    }
    return vld1q_s16(x);
  }
  static Type PrefixMax(Type x, Type fill) {
    x = Max(x, vextq_s16(fill, x, 7));
    x = Max(x, vextq_s16(fill, x, 6));
    return Max(x, vextq_s16(fill, x, 4));
  }
  static Type BroadcastLast(Type x) { return vdupq_laneq_s16(x, 7); }
  static std::int32_t First(Type x) { return vgetq_lane_s16(x, 0); }
  static std::int32_t HorizontalMax(Type x) { return vmaxvq_s16(x); }
};

#else

struct VecI32 {
  using Score = std::int32_t;
  using Type = __m128i;
  static constexpr std::uint32_t kLanes = 4;
  static Type Load(const Score* p) {
    return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
  }
  static void Store(Score* p, Type x) {
    _mm_storeu_si128(reinterpret_cast<__m128i*>(p), x);
  }
  static Type Set(Score x) { return _mm_set1_epi32(x); }
  static Type Add(Type a, Type b) { return _mm_add_epi32(a, b); }
  static Type AddSat(Type a, Type b) { return _mm_add_epi32(a, b); }
  static Type Sub(Type a, Type b) { return _mm_sub_epi32(a, b); }
  static Type Max(Type a, Type b) {
#if defined(__SSE4_1__)
    return _mm_max_epi32(a, b);
#else
    const __m128i a_greater = _mm_cmpgt_epi32(a, b);
    return _mm_or_si128(
        _mm_and_si128(a_greater, a),
        _mm_andnot_si128(a_greater, b));
#endif
  }
  static Type Ramp(std::int32_t g) {
    return _mm_setr_epi32(g, 2 * g, 3 * g, 4 * g);
  }
  // {fill, x[0], x[1], x[2]} and {fill, fill, x[0], x[1]}, where fill has the same value in
  // every lane
  static Type PrefixMax(Type x, Type fill) {
    x = Max(x, _mm_or_si128(_mm_slli_si128(x, 4), _mm_srli_si128(fill, 12)));
    return Max(x, _mm_or_si128(_mm_slli_si128(x, 8), _mm_srli_si128(fill, 8)));
  }
  static Type BroadcastLast(Type x) { return _mm_shuffle_epi32(x, 0xFF); }
  static std::int32_t First(Type x) { return _mm_cvtsi128_si32(x); }
  static std::int32_t HorizontalMax(Type x) {
    x = Max(x, _mm_shuffle_epi32(x, 0x4E));
    x = Max(x, _mm_shuffle_epi32(x, 0xB1));
    return First(x);
  }
};

struct VecI16 {
  using Score = std::int16_t;
  using Type = __m128i;
  static constexpr std::uint32_t kLanes = 8;
  static Type Load(const Score* p) {
    return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
  }
  static void Store(Score* p, Type x) {
    _mm_storeu_si128(reinterpret_cast<__m128i*>(p), x);
  }
  static Type Set(Score x) { return _mm_set1_epi16(x); }
  static Type Add(Type a, Type b) { return _mm_add_epi16(a, b); }
  static Type AddSat(Type a, Type b) { return _mm_adds_epi16(a, b); }
  static Type Sub(Type a, Type b) { return _mm_sub_epi16(a, b); }
  static Type Max(Type a, Type b) { return _mm_max_epi16(a, b); }
  static Type Ramp(std::int32_t g) {
    return _mm_setr_epi16(
        static_cast<Score>(g), static_cast<Score>(2 * g), static_cast<Score>(3 * g),
        static_cast<Score>(4 * g), static_cast<Score>(5 * g), static_cast<Score>(6 * g),
        static_cast<Score>(7 * g), static_cast<Score>(8 * g));
  }
  static Type PrefixMax(Type x, Type fill) {
    x = Max(x, _mm_or_si128(_mm_slli_si128(x, 2), _mm_srli_si128(fill, 14)));
    x = Max(x, _mm_or_si128(_mm_slli_si128(x, 4), _mm_srli_si128(fill, 12)));
    return Max(x, _mm_or_si128(_mm_slli_si128(x, 8), _mm_srli_si128(fill, 8)));
  }
  static Type BroadcastLast(Type x) {
    const __m128i t = _mm_shufflehi_epi16(x, 0xFF);
    return _mm_unpackhi_epi64(t, t);
  }
  static std::int32_t First(Type x) {
    return static_cast<Score>(_mm_cvtsi128_si32(x));
  }
  static std::int32_t HorizontalMax(Type x) {
    x = Max(x, _mm_shuffle_epi32(x, 0x4E));
    x = Max(x, _mm_shuffle_epi32(x, 0xB1));
    x = Max(x, _mm_shufflelo_epi16(x, 0xB1));
    return First(x);
  }
};

#endif

template <typename V, bool kIsSW>
std::int32_t ComputeRowVec(
    const typename V::Score* const* pred_rows,
    std::uint32_t num_preds,
    const typename V::Score* __restrict profile,
    std::int32_t g,
    std::uint32_t width,
    typename V::Score* __restrict H) {
  SPOA_COMPUTE_ROW_VEC_BODY
}

template <typename V, bool kIsSW>
void ComputeRowRangeVec(
    const typename V::Score* const* pred_rows,
    std::uint32_t num_preds,
    const typename V::Score* __restrict profile,
    std::int32_t g,
    std::uint32_t lo,
    std::uint32_t end,
    std::int32_t running,
    typename V::Score* __restrict H) {
  SPOA_COMPUTE_ROW_RANGE_VEC_BODY
}

// Vector type for each score type
template <typename T>
struct VecFor;
template <>
struct VecFor<std::int32_t> { using V = VecI32; };
template <>
struct VecFor<std::int16_t> { using V = VecI16; };

#if defined(SPOA_LOW_MEMORY_AVX2)

// Vector types with twice as many lanes, for CPUs found to support AVX2 at runtime
#define SPOA_AVX2 __attribute__((target("avx2")))

// Shift x up by one 128-bit half, with fill shifted into the lower half
SPOA_AVX2 inline __m256i ShiftUpHalfAvx2(__m256i x, __m256i fill) {
  return _mm256_permute2x128_si256(x, fill, 0x02);  // {fill_lo, x_lo}
}

struct Avx2I32 {
  using Score = std::int32_t;
  using Type = __m256i;
  static constexpr std::uint32_t kLanes = 8;
  SPOA_AVX2 static Type Load(const Score* p) {
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
  }
  SPOA_AVX2 static void Store(Score* p, Type x) {
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p), x);
  }
  SPOA_AVX2 static Type Set(Score x) { return _mm256_set1_epi32(x); }
  SPOA_AVX2 static Type Add(Type a, Type b) { return _mm256_add_epi32(a, b); }
  SPOA_AVX2 static Type AddSat(Type a, Type b) { return _mm256_add_epi32(a, b); }
  SPOA_AVX2 static Type Sub(Type a, Type b) { return _mm256_sub_epi32(a, b); }
  SPOA_AVX2 static Type Max(Type a, Type b) { return _mm256_max_epi32(a, b); }
  SPOA_AVX2 static Type Ramp(std::int32_t g) {
    return _mm256_setr_epi32(g, 2 * g, 3 * g, 4 * g, 5 * g, 6 * g, 7 * g, 8 * g);
  }
  SPOA_AVX2 static Type PrefixMax(Type x, Type fill) {
    x = Max(x, _mm256_alignr_epi8(x, ShiftUpHalfAvx2(x, fill), 12));
    x = Max(x, _mm256_alignr_epi8(x, ShiftUpHalfAvx2(x, fill), 8));
    return Max(x, ShiftUpHalfAvx2(x, fill));
  }
  SPOA_AVX2 static Type BroadcastLast(Type x) {
    return _mm256_permutevar8x32_epi32(x, _mm256_set1_epi32(7));
  }
  SPOA_AVX2 static std::int32_t First(Type x) {
    return _mm_cvtsi128_si32(_mm256_castsi256_si128(x));
  }
  SPOA_AVX2 static std::int32_t HorizontalMax(Type x) {
    __m128i m = _mm_max_epi32(_mm256_castsi256_si128(x), _mm256_extracti128_si256(x, 1));
    m = _mm_max_epi32(m, _mm_shuffle_epi32(m, 0x4E));
    m = _mm_max_epi32(m, _mm_shuffle_epi32(m, 0xB1));
    return _mm_cvtsi128_si32(m);
  }
};

struct Avx2I16 {
  using Score = std::int16_t;
  using Type = __m256i;
  static constexpr std::uint32_t kLanes = 16;
  SPOA_AVX2 static Type Load(const Score* p) {
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
  }
  SPOA_AVX2 static void Store(Score* p, Type x) {
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p), x);
  }
  SPOA_AVX2 static Type Set(Score x) { return _mm256_set1_epi16(x); }
  SPOA_AVX2 static Type Add(Type a, Type b) { return _mm256_add_epi16(a, b); }
  SPOA_AVX2 static Type AddSat(Type a, Type b) { return _mm256_adds_epi16(a, b); }
  SPOA_AVX2 static Type Sub(Type a, Type b) { return _mm256_sub_epi16(a, b); }
  SPOA_AVX2 static Type Max(Type a, Type b) { return _mm256_max_epi16(a, b); }
  SPOA_AVX2 static Type Ramp(std::int32_t g) {
    Score x[kLanes];
    for (std::uint32_t i = 0; i < kLanes; ++i) {
      x[i] = static_cast<Score>((i + 1) * g);
    }
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x));
  }
  SPOA_AVX2 static Type PrefixMax(Type x, Type fill) {
    x = Max(x, _mm256_alignr_epi8(x, ShiftUpHalfAvx2(x, fill), 14));
    x = Max(x, _mm256_alignr_epi8(x, ShiftUpHalfAvx2(x, fill), 12));
    x = Max(x, _mm256_alignr_epi8(x, ShiftUpHalfAvx2(x, fill), 8));
    return Max(x, ShiftUpHalfAvx2(x, fill));
  }
  SPOA_AVX2 static Type BroadcastLast(Type x) {
    const __m256i t = _mm256_shufflehi_epi16(_mm256_permute4x64_epi64(x, 0xFF), 0xFF);
    return _mm256_unpackhi_epi64(t, t);
  }
  SPOA_AVX2 static std::int32_t First(Type x) {
    return static_cast<Score>(_mm_cvtsi128_si32(_mm256_castsi256_si128(x)));
  }
  SPOA_AVX2 static std::int32_t HorizontalMax(Type x) {
    __m128i m = _mm_max_epi16(_mm256_castsi256_si128(x), _mm256_extracti128_si256(x, 1));
    m = _mm_max_epi16(m, _mm_shuffle_epi32(m, 0x4E));
    m = _mm_max_epi16(m, _mm_shuffle_epi32(m, 0xB1));
    m = _mm_max_epi16(m, _mm_shufflelo_epi16(m, 0xB1));
    return static_cast<Score>(_mm_cvtsi128_si32(m));
  }
};

template <typename V, bool kIsSW>
SPOA_AVX2 std::int32_t ComputeRowVecAvx2(
    const typename V::Score* const* pred_rows,
    std::uint32_t num_preds,
    const typename V::Score* __restrict profile,
    std::int32_t g,
    std::uint32_t width,
    typename V::Score* __restrict H) {
  SPOA_COMPUTE_ROW_VEC_BODY
}

template <typename V, bool kIsSW>
SPOA_AVX2 void ComputeRowRangeVecAvx2(
    const typename V::Score* const* pred_rows,
    std::uint32_t num_preds,
    const typename V::Score* __restrict profile,
    std::int32_t g,
    std::uint32_t lo,
    std::uint32_t end,
    std::int32_t running,
    typename V::Score* __restrict H) {
  SPOA_COMPUTE_ROW_RANGE_VEC_BODY
}

template <typename T>
struct Avx2VecFor;
template <>
struct Avx2VecFor<std::int32_t> { using V = Avx2I32; };
template <>
struct Avx2VecFor<std::int16_t> { using V = Avx2I16; };

#undef SPOA_AVX2

bool HasAvx2() {
  static const bool has_avx2 =
      __builtin_cpu_supports("avx2") && std::getenv("SPOA_DISABLE_AVX2") == nullptr;
  return has_avx2;
}

#endif

#endif

#undef SPOA_COMPUTE_ROW_VEC_BODY
#undef SPOA_COMPUTE_ROW_RANGE_VEC_BODY

// Compute a score row with the fastest available implementation
template <typename T, bool kIsSW>
std::int32_t ComputeRowFastest(
    const T* const* pred_rows,
    std::uint32_t num_preds,
    const T* __restrict profile,
    std::int32_t g,
    std::uint32_t width,
    T* __restrict H) {
#if defined(__aarch64__) || defined(__SSE2__)
#if defined(SPOA_LOW_MEMORY_AVX2)
  if (HasAvx2()) {
    return ComputeRowVecAvx2<typename Avx2VecFor<T>::V, kIsSW>(
        pred_rows, num_preds, profile, g, width, H);
  }
#endif
  return ComputeRowVec<typename VecFor<T>::V, kIsSW>(
      pred_rows, num_preds, profile, g, width, H);
#else
  return ComputeRow<T, kIsSW>(pred_rows, num_preds, profile, g, width, H);
#endif
}

// Compute score row cells [lo, end) with the fastest available implementation, see
// SPOA_COMPUTE_ROW_RANGE_VEC_BODY
template <typename T, bool kIsSW>
void ComputeRowRangeFastest(
    const T* const* pred_rows,
    std::uint32_t num_preds,
    const T* __restrict profile,
    std::int32_t g,
    std::uint32_t lo,
    std::uint32_t end,
    std::int32_t running,
    T* __restrict H) {
#if defined(__aarch64__) || defined(__SSE2__)
#if defined(SPOA_LOW_MEMORY_AVX2)
  if (HasAvx2()) {
    ComputeRowRangeVecAvx2<typename Avx2VecFor<T>::V, kIsSW>(
        pred_rows, num_preds, profile, g, lo, end, running, H);
    return;
  }
#endif
  ComputeRowRangeVec<typename VecFor<T>::V, kIsSW>(
      pred_rows, num_preds, profile, g, lo, end, running, H);
#else
  ComputeRowScalar<T, kIsSW>(pred_rows, num_preds, profile, g, end, H, lo, running, 0);
#endif
}

}  // namespace

bool SisdAlignmentEngine::LinearLowMemory(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score,
    Alignment* alignment) {
  if (g_ >= 0) {
    return false;
  }
  bool ok;
  // Use 16-bit scores where they can't overflow, as these are faster to compute
  if (IsInt16ScoreRange(sequence_len, m_, n_, g_)) {
    ok = LinearLowMemoryT<std::int16_t>(sequence, sequence_len, graph, score, alignment,
        &(pimpl_->sequence_profile16));
  } else {
    ok = LinearLowMemoryT<std::int32_t>(sequence, sequence_len, graph, score, alignment,
      &(pimpl_->sequence_profile));
  }
  return ok;
}

template <typename T>
bool SisdAlignmentEngine::LinearLowMemoryT(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score,
    Alignment* alignment,
    std::vector<T>* sequence_profile) {
  const auto& rank_to_node = graph.rank_to_node();
  const std::uint32_t num_nodes = rank_to_node.size();
  const std::uint32_t width = sequence_len + 1;
  const bool is_sw = type_ == AlignmentType::kSW;
  const std::int32_t g = g_;

  if (pimpl_->node_id_to_rank.size() < graph.nodes().size()) {
    pimpl_->node_id_to_rank.resize(graph.nodes().size(), 0);
  }
  for (std::uint32_t i = 0; i < num_nodes; ++i) {
    pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
  }

  if (sequence_profile->size() < graph.num_codes() * width) {
    sequence_profile->resize(graph.num_codes() * width, 0);
  }
  for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
    char c = graph.decoder(i);
    T* profile = &((*sequence_profile)[i * width]);
    profile[0] = 0;
    for (std::uint32_t j = 0; j < sequence_len; ++j) {
      profile[j + 1] = static_cast<T>(c == sequence[j] ? m_ : n_);
    }
  }
  auto row_profile = [&] (std::uint32_t i) -> const T* {
    return &((*sequence_profile)[rank_to_node[i - 1]->code * width]);
  };

  // Predecessor rows of each row i (node rank i - 1), in in-edge order, with
  // row 0 standing in for nodes without in-edges, as in Linear
  std::vector<std::uint32_t> pred_offsets(num_nodes + 2, 0);
  std::vector<std::uint32_t> preds;
  // Number of rows still to be computed which use each row as a predecessor
  std::vector<std::uint32_t> pending(num_nodes + 1, 0);
  // Rows needed in full by the backtrack
  std::vector<bool> is_stored_row(num_nodes + 1, false);
  auto block = [] (std::uint32_t i) { return (i - 1) / kRankBlockSize; };
  for (std::uint32_t i = 1; i <= num_nodes; ++i) {
    const auto& it = rank_to_node[i - 1];
    pred_offsets[i] = preds.size();
    if (it->inedges.empty()) {
      preds.emplace_back(0);
      ++pending[0];
    } else {
      for (const auto& jt : it->inedges) {
        std::uint32_t pred_i = pimpl_->node_id_to_rank[jt->tail->id] + 1;
        preds.emplace_back(pred_i);
        ++pending[pred_i];
        if (block(pred_i) != block(i)) {
          is_stored_row[pred_i] = true;
        }
      }
    }
  }
  pred_offsets[num_nodes + 1] = preds.size();

  const std::uint32_t num_stored_columns = (width - 1) / kColumnStride + 1;
  ScoreBuffer<T> stored_columns(
      static_cast<std::uint64_t>(num_nodes + 1) * num_stored_columns);
  std::vector<ScoreBuffer<T>> stored_rows(num_nodes + 1);

  // Score rows, recycled once all successors are computed
  std::vector<ScoreBuffer<T>> free_rows;
  std::vector<ScoreBuffer<T>> rows(num_nodes + 1);
  auto get_row = [&] () -> ScoreBuffer<T> {
    if (free_rows.empty()) {
      return ScoreBuffer<T>(width);
    }
    auto row = std::move(free_rows.back());
    free_rows.pop_back();
    return row;
  };
  auto release_row = [&] (std::uint32_t i) {
    if (is_stored_row[i]) {
      stored_rows[i] = std::move(rows[i]);
    } else {
      free_rows.emplace_back(std::move(rows[i]));
    }
  };

  // Row 0 scores, also used directly in the backtrack
  auto row_zero_score = [is_sw, g] (std::uint32_t j) -> std::int32_t {
    return is_sw ? 0 : static_cast<std::int32_t>(j) * g;
  };
  rows[0] = get_row();
  for (std::uint32_t j = 0; j < width; ++j) {
    rows[0][j] = static_cast<T>(row_zero_score(j));
  }

  std::int32_t max_score = is_sw ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;

  std::vector<const T*> pred_rows;
  for (std::uint32_t i = 1; i <= num_nodes; ++i) {
    const auto& it = rank_to_node[i - 1];
    const T* char_profile = row_profile(i);
    const std::uint32_t* row_preds = &preds[pred_offsets[i]];
    const std::uint32_t num_preds = pred_offsets[i + 1] - pred_offsets[i];

    rows[i] = get_row();
    T* H_row = rows[i].get();

    pred_rows.clear();
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      pred_rows.emplace_back(rows[row_preds[p]].get());
    }
    const std::int32_t row_max = is_sw ?
        ComputeRowFastest<T, true>(pred_rows.data(), num_preds, char_profile, g, width, H_row) :  // NOLINT
        ComputeRowFastest<T, false>(pred_rows.data(), num_preds, char_profile, g, width, H_row);  // NOLINT

    // end cell search
    if ((is_sw || it->outedges.empty()) && max_score < row_max) {
      max_score = row_max;
      max_i = i;
      for (std::uint32_t j = 1; j < width; ++j) {
        if (H_row[j] == row_max) {
          max_j = j;
          break;
        }
      }
    }

    T* row_stored_columns =
        &stored_columns[static_cast<std::uint64_t>(i) * num_stored_columns];
    for (std::uint32_t k = 0; k < num_stored_columns; ++k) {
      row_stored_columns[k] = H_row[k * kColumnStride];
    }

    // recycle score rows which are no longer needed
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      std::uint32_t pred_i = row_preds[p];
      if (--pending[pred_i] == 0) {
        release_row(pred_i);
      }
    }
    if (pending[i] == 0) {
      release_row(i);
    }
  }
  rows.clear();
  free_rows.clear();

  alignment->clear();
  if (max_i == 0 && max_j == 0) {
    return true;
  }
  *score = max_score;

  // Backtrack tile: rows [tile_first_row, tile_last_row] and columns
  // [tile_first_column, tile_last_column]
  std::uint32_t tile_first_row = 1, tile_last_row = 0;
  std::uint32_t tile_first_column = 1, tile_last_column = 0;
  const std::uint32_t tile_width = kColumnStride + 1;
  std::vector<std::int32_t> tile(kRankBlockSize * tile_width);
  auto tile_score = [&] (std::uint32_t i, std::uint32_t j) -> std::int32_t& {
    return tile[(i - tile_first_row) * tile_width + (j - tile_first_column)];
  };

  // Score of a cell in the tile, row zero, or a stored row
  auto get_score = [&] (std::uint32_t i, std::uint32_t j) -> std::int32_t {
    if (i == 0) {
      return row_zero_score(j);
    } else if (i >= tile_first_row) {
      return tile_score(i, j);
    } else {
      return stored_rows[i][j];
    }
  };

  // Recompute the tile containing the cells needed to backtrack from H[i][j]
  auto compute_tile = [&] (std::uint32_t i, std::uint32_t j) {
    tile_first_row = block(i) * kRankBlockSize + 1;
    tile_last_row = i;
    tile_first_column = ((j - 1) / kColumnStride) * kColumnStride;
    tile_last_column = j;
    for (std::uint32_t r = tile_first_row; r <= tile_last_row; ++r) {
      const T* char_profile = row_profile(r);
      const std::uint32_t* row_preds = &preds[pred_offsets[r]];
      const std::uint32_t num_preds = pred_offsets[r + 1] - pred_offsets[r];
      tile_score(r, tile_first_column) = stored_columns[
          static_cast<std::uint64_t>(r) * num_stored_columns +
          tile_first_column / kColumnStride];
      for (std::uint32_t c = tile_first_column + 1; c <= tile_last_column; ++c) {  // NOLINT
        std::int32_t h = kNegativeInfinity;
        for (std::uint32_t p = 0; p < num_preds; ++p) {
          h = std::max(h, std::max(
              get_score(row_preds[p], c - 1) + char_profile[c],
              get_score(row_preds[p], c) + g));
        }
        h = std::max(h, tile_score(r, c - 1) + g);
        if (is_sw) {
          h = std::max(h, 0);
        }
        tile_score(r, c) = h;
      }
    }
  };

  std::uint32_t i = max_i;
  std::uint32_t j = max_j;
  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;
  while (i != 0 && j != 0) {
    if (i < tile_first_row || i > tile_last_row ||
        j - 1 < tile_first_column || j > tile_last_column) {
      compute_tile(i, j);
    }
    const std::int32_t H_ij = tile_score(i, j);
    if (is_sw && H_ij == 0) {
      break;
    }

    const std::uint32_t* row_preds = &preds[pred_offsets[i]];
    const std::uint32_t num_preds = pred_offsets[i + 1] - pred_offsets[i];
    const std::int32_t match_cost = row_profile(i)[j];
    bool predecessor_found = false;
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      if (H_ij == get_score(row_preds[p], j - 1) + match_cost) {
        prev_i = row_preds[p];
        prev_j = j - 1;
        predecessor_found = true;
        break;
      }
    }
    if (!predecessor_found) {
      for (std::uint32_t p = 0; p < num_preds; ++p) {
        if (H_ij == get_score(row_preds[p], j) + g) {
          prev_i = row_preds[p];
          prev_j = j;
          predecessor_found = true;
          break;
        }
      }
    }
    if (!predecessor_found && H_ij == tile_score(i, j - 1) + g) {
      prev_i = i;
      prev_j = j - 1;
    }

    alignment->emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;
  }

  std::reverse(alignment->begin(), alignment->end());

  return true;
}

// Linear gap kSW or kOV alignment which skips cells that can't be on an alignment scoring at least
// min_score
//
// Each cell's score H plus an upper bound on the score that can still be added after it,
// m * min(remaining sequence length, longest remaining graph path), is compared to min_score, and
// cells below min_score are pruned: they are treated as -infinity when used as predecessors. Each
// row is only computed over the columns that can be reached from unpruned cells. Rows store only
// their unpruned cells, so the backtrack can use the stored scores directly.
//
// Pruned cells can only lower the scores of other cells, to the score of another valid path, and
// no cell on an alignment with score >= min_score is ever pruned, since the bound is an upper
// bound. So if the best score found is at least min_score, it is the optimal score, every cell on
// an optimal path has its exact score, and the end cell search and backtrack (which only compare
// scores of cells against the optimal path) give exactly the same alignment as Linear.
//
// Returns false if the best score is below min_score, or too many cells can't be pruned, in which
// case the alignment isn't valid.
//
namespace {

// Number of nodes on the longest path through the graph
std::uint32_t GraphLongestPath(const Graph& graph) {
  const auto& rank_to_node = graph.rank_to_node();
  std::vector<std::uint32_t> node_id_to_rank(graph.nodes().size(), 0);
  for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
    node_id_to_rank[rank_to_node[i]->id] = i;
  }
  std::vector<std::uint32_t> longest(rank_to_node.size(), 0);
  std::uint32_t max_longest = 0;
  for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
    std::uint32_t x = 0;
    for (const auto& jt : rank_to_node[i]->inedges) {
      x = std::max(x, longest[node_id_to_rank[jt->tail->id]]);
    }
    longest[i] = x + 1;
    max_longest = std::max(max_longest, longest[i]);
  }
  return max_longest;
}

}  // namespace

bool SisdAlignmentEngine::LinearPruned(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t min_score,
    std::int32_t* score,
    Alignment* alignment,
    bool* aborted) {
  // Use 16-bit scores where they can't overflow, as these are faster to compute
  if (IsInt16ScoreRange(sequence_len, m_, n_, g_)) {
    return LinearPrunedT<std::int16_t>(sequence, sequence_len, graph, min_score, score, alignment,
        aborted, &(pimpl_->sequence_profile16));
  }
  return LinearPrunedT<std::int32_t>(sequence, sequence_len, graph, min_score, score, alignment,
      aborted, &(pimpl_->sequence_profile));
}

template <typename T>
PrunedBuffers<T>& SelectPrunedBuffers(
    PrunedBuffers<std::int16_t>& b16, PrunedBuffers<std::int32_t>& b32);

template <>
PrunedBuffers<std::int16_t>& SelectPrunedBuffers(
    PrunedBuffers<std::int16_t>& b16, PrunedBuffers<std::int32_t>&) {
  return b16;
}

template <>
PrunedBuffers<std::int32_t>& SelectPrunedBuffers(
    PrunedBuffers<std::int16_t>&, PrunedBuffers<std::int32_t>& b32) {
  return b32;
}

template <typename T>
bool SisdAlignmentEngine::LinearPrunedT(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t min_score,
    std::int32_t* score,
    Alignment* alignment,
    bool* aborted,
    std::vector<T>* sequence_profile_ptr) {
  *aborted = false;
  const T kNegInf = ScanNegativeInfinity<T>();
  const auto& rank_to_node = graph.rank_to_node();
  const std::uint32_t num_nodes = rank_to_node.size();
  const std::uint32_t width = sequence_len + 1;
  const bool is_sw = type_ == AlignmentType::kSW;
  const std::int32_t g = g_;
  const std::int32_t m = m_;
  if (m <= 0 || g >= 0 || num_nodes == 0) {
    *aborted = true;
    return false;
  }

  if (pimpl_->node_id_to_rank.size() < graph.nodes().size()) {
    pimpl_->node_id_to_rank.resize(graph.nodes().size(), 0);
  }
  for (std::uint32_t i = 0; i < num_nodes; ++i) {
    pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
  }

  // Profile with padding for vector loads
  std::vector<T>& sequence_profile = *sequence_profile_ptr;
  const std::uint32_t padded_width = width + 64;
  if (sequence_profile.size() < graph.num_codes() * padded_width) {
    sequence_profile.resize(graph.num_codes() * padded_width, 0);
  }
  for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
    char c = graph.decoder(i);
    T* profile = &sequence_profile[i * padded_width];
    profile[0] = 0;
    for (std::uint32_t j = 0; j < sequence_len; ++j) {
      profile[j + 1] = static_cast<T>(c == sequence[j] ? m_ : n_);
    }
  }
  auto row_profile = [&] (std::uint32_t i) -> const T* {
    return &sequence_profile[rank_to_node[i - 1]->code * padded_width];
  };

  // Predecessor rows, as in LinearLowMemoryT
  std::vector<std::uint32_t> pred_offsets(num_nodes + 2, 0);
  std::vector<std::uint32_t> preds;
  std::uint32_t max_preds = 1;
  for (std::uint32_t i = 1; i <= num_nodes; ++i) {
    const auto& it = rank_to_node[i - 1];
    pred_offsets[i] = preds.size();
    if (it->inedges.empty()) {
      preds.emplace_back(0);
    } else {
      for (const auto& jt : it->inedges) {
        preds.emplace_back(pimpl_->node_id_to_rank[jt->tail->id] + 1);
      }
      max_preds = std::max<std::uint32_t>(max_preds, it->inedges.size());
    }
  }
  pred_offsets[num_nodes + 1] = preds.size();

  // Longest path length (in nodes) after each row's node, to the end of the graph
  std::vector<std::uint32_t> remaining_path(num_nodes + 1, 0);
  for (std::uint32_t i = num_nodes; i >= 1; --i) {
    std::uint32_t x = 0;
    for (const auto& jt : rank_to_node[i - 1]->outedges) {
      const std::uint32_t succ_i = pimpl_->node_id_to_rank[jt->head->id] + 1;
      x = std::max(x, 1 + remaining_path[succ_i]);
    }
    remaining_path[i] = x;
  }
  // Upper bound of the score which can be added after cell (i, j)
  auto remaining_bound = [&] (std::uint32_t i, std::uint32_t j) -> std::int64_t {
    return static_cast<std::int64_t>(m) *
        std::min<std::int64_t>(sequence_len - j, remaining_path[i]);
  };

  PrunedBuffers<T>& buf = SelectPrunedBuffers<T>(pimpl_->pruned16, pimpl_->pruned32);

  // Row 0 is computed in full and never pruned
  std::vector<T>& row0 = buf.row0;
  row0.assign(padded_width, kNegInf);
  for (std::uint32_t j = 0; j < width; ++j) {
    row0[j] = is_sw ? 0 : static_cast<T>(j) * g;
  }

  // Unpruned cells of each row i: columns [kept_lo[i], kept_hi[i]] stored from vals[start[i]]
  std::vector<std::uint32_t> kept_lo(num_nodes + 1, 1), kept_hi(num_nodes + 1, 0);
  std::vector<std::uint64_t> start(num_nodes + 1, 0);
  std::vector<T>& vals = buf.vals;
  vals.clear();
  // Pruned alignment costs a few times more per cell than the full alignment, so it's aborted once
  // too many cells are kept for it to be worthwhile
  const std::uint64_t max_cells = std::min<std::uint64_t>(
      std::max<std::uint64_t>(1u << 16, static_cast<std::uint64_t>(num_nodes) * width / 16),
      1u << 25);

  // Full width score rows, with -infinity outside unpruned cells and 0 in column 0, recycled once
  // all successors are computed
  std::vector<std::uint32_t> pending(num_nodes + 1, 0);
  for (std::uint32_t i = 1; i <= num_nodes; ++i) {
    for (std::uint32_t p = pred_offsets[i]; p < pred_offsets[i + 1]; ++p) {
      ++pending[preds[p]];
    }
  }
  const std::uint32_t kNoRow = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> row_index(num_nodes + 1, kNoRow);
  auto get_row = [&] () -> std::uint32_t {
    std::uint32_t r;
    if (buf.free_rows.empty()) {
      r = buf.rows.size();
      buf.rows.emplace_back();
    } else {
      r = buf.free_rows.back();
      buf.free_rows.pop_back();
    }
    auto& row = buf.rows[r];
    if (row.size() < padded_width) {
      row.resize(padded_width, kNegInf);
    }
    row[0] = 0;
    return r;
  };
  auto release_row = [&] (std::uint32_t i) {
    auto& row = buf.rows[row_index[i]];
    if (kept_lo[i] <= kept_hi[i]) {
      std::fill(&row[kept_lo[i]], &row[kept_hi[i]] + 1, kNegInf);
    }
    buf.free_rows.emplace_back(row_index[i]);
    row_index[i] = kNoRow;
  };
  // Release all rows still in use, leaving all row buffers set to -infinity
  auto release_all_rows = [&] () {
    for (std::uint32_t i = 1; i <= num_nodes; ++i) {
      if (row_index[i] != kNoRow) {
        release_row(i);
      }
    }
  };
  std::vector<const T*> pred_rows;

  // Minimum number of sequence bases which must be matched to reach min_score
  const std::int64_t min_matches = (static_cast<std::int64_t>(min_score) + m - 1) / m;

  std::int32_t max_score = is_sw ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;

  for (std::uint32_t i = 1; i <= num_nodes; ++i) {
    const auto& it = rank_to_node[i - 1];
    const T* char_profile = row_profile(i);
    const std::uint32_t* row_preds = &preds[pred_offsets[i]];
    const std::uint32_t num_preds = pred_offsets[i + 1] - pred_offsets[i];

    // Columns which may hold unpruned cells: those reachable from unpruned predecessor cells,
    // column 1 (reachable from column 0), and for kSW any column where a local alignment can
    // start and still reach min_score
    std::int64_t lo = std::numeric_limits<std::int64_t>::max();
    std::int64_t hi = -1;
    auto add_range = [&] (std::int64_t a, std::int64_t b) {
      a = std::max<std::int64_t>(a, 1);
      b = std::min<std::int64_t>(b, width - 1);
      if (a <= b) {
        lo = std::min(lo, a);
        hi = std::max(hi, b);
      }
    };
    pred_rows.clear();
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      const std::uint32_t pi = row_preds[p];
      if (pi == 0) {
        add_range(1, width - 1);
        pred_rows.emplace_back(row0.data());
      } else {
        if (kept_lo[pi] <= kept_hi[pi]) {
          add_range(kept_lo[pi], static_cast<std::int64_t>(kept_hi[pi]) + 1);
        }
        pred_rows.emplace_back(buf.rows[row_index[pi]].data());
      }
    }
    const std::int32_t L = static_cast<std::int32_t>(
        std::min<std::uint32_t>(remaining_path[i], sequence_len));
    if (static_cast<std::int64_t>(m) + remaining_bound(i, 1) >= min_score) {
      add_range(1, 1);
    }
    if (is_sw && static_cast<std::int64_t>(remaining_path[i]) >= min_matches) {
      add_range(1, static_cast<std::int64_t>(sequence_len) - min_matches);
    }

    row_index[i] = get_row();
    T* cur = buf.rows[row_index[i]].data();
    std::uint32_t row_lo = 1, row_hi = 0;
    if (lo <= hi) {
      // Compute the candidate range, then extend it to the right while cells stay unpruned
      const std::int32_t running = (lo == 1) ? 0 : kNegInf;
      if (is_sw) {
        ComputeRowRangeFastest<T, true>(pred_rows.data(), num_preds, char_profile, g, lo, hi + 1,
            running, cur);
      } else {
        ComputeRowRangeFastest<T, false>(pred_rows.data(), num_preds, char_profile, g, lo, hi + 1,
            running, cur);
      }
      std::uint32_t j = hi + 1;
      for (; j < width; ++j) {
        std::int32_t h = cur[j - 1] + g;
        for (std::uint32_t p = 0; p < num_preds; ++p) {
          const T* P = pred_rows[p];
          h = std::max(h, std::max(P[j - 1] + char_profile[j], P[j] + g));
        }
        if (is_sw) {
          h = std::max(h, 0);
        }
        if (h + remaining_bound(i, j) < min_score) {
          break;
        }
        cur[j] = h;
      }
      hi = j - 1;

      // Prune: cell j is kept if cur[j] + m * min(sequence_len - j, L) >= min_score
      const std::uint32_t ulo = lo, uhi = hi;
      const std::int32_t len = sequence_len;
      {
        // Cells up to column len - L have the same bound, and the bound decreases after that. The
        // loops are written so that they can be vectorized.
        std::uint32_t jj = ulo;
        const std::int64_t split = static_cast<std::int64_t>(len) - L;
        if (static_cast<std::int64_t>(jj) <= split) {
          const std::uint32_t end1 = (static_cast<std::int64_t>(uhi) <= split)
              ? uhi : static_cast<std::uint32_t>(split);
          const std::int32_t threshold = min_score - m * L;
          T* __restrict row = cur;
          for (; jj <= end1; ++jj) {
            const T c = row[jj];
            row[jj] = (static_cast<std::int32_t>(c) < threshold) ? kNegInf : c;
          }
        }
        const std::int32_t threshold_base = min_score - m * len;
        T* __restrict row = cur;
        for (; jj <= uhi; ++jj) {
          const T c = row[jj];
          row[jj] = (static_cast<std::int32_t>(c) <
              threshold_base + m * static_cast<std::int32_t>(jj)) ? kNegInf : c;
        }
      }
      std::uint32_t a = ulo, b = uhi;
      while (a <= b && cur[a] == kNegInf) {
        ++a;
      }
      while (b >= a && b > 0 && cur[b] == kNegInf) {
        --b;
      }
      if (a <= b) {
        row_lo = a;
        row_hi = b;
      }
      // Clear computed cells outside the unpruned range, so that the row is -infinity outside
      // [row_lo, row_hi]
      if (row_lo > row_hi) {
        std::fill(&cur[ulo], &cur[uhi] + 1, kNegInf);
      }
    }

    // Store unpruned cells, and search for the end cell
    if (row_lo <= row_hi) {
      kept_lo[i] = row_lo;
      kept_hi[i] = row_hi;
      start[i] = vals.size();
      vals.insert(vals.end(), &cur[row_lo], &cur[row_hi] + 1);
      // Abort early if cells are kept at a higher rate than the limit allows
      if (vals.size() > max_cells ||
          vals.size() > (1u << 16) + 3 * (max_cells / 2) * static_cast<std::uint64_t>(i) /
              num_nodes) {
        release_all_rows();
        *aborted = true;
        return false;
      }
      if (is_sw || it->outedges.empty()) {
        T row_max_t = kNegInf;
        for (std::uint32_t jj = row_lo; jj <= row_hi; ++jj) {
          row_max_t = std::max<T>(row_max_t, cur[jj]);
        }
        const std::int32_t row_max = row_max_t;
        if (max_score < row_max) {
          std::uint32_t jj = row_lo;
          while (cur[jj] != row_max) {
            ++jj;
          }
          max_score = row_max;
          max_i = i;
          max_j = jj;
        }
      }
    }

    // Recycle score rows which are no longer needed
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      const std::uint32_t pi = row_preds[p];
      if (pi != 0 && --pending[pi] == 0) {
        release_row(pi);
      }
    }
    if (pending[i] == 0) {
      release_row(i);
    }
  }

  alignment->clear();
  if ((max_i == 0 && max_j == 0) || max_score < min_score) {
    return false;
  }
  *score = max_score;

  auto get_score = [&] (std::uint32_t i, std::uint32_t j) -> std::int32_t {
    if (i == 0) {
      return row0[j];
    }
    if (j == 0) {
      return 0;
    }
    if (j < kept_lo[i] || j > kept_hi[i]) {
      return kNegInf;
    }
    return vals[start[i] + (j - kept_lo[i])];
  };

  // Backtrack, as in LinearLowMemoryT
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;
  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;
  while (i != 0 && j != 0) {
    const std::int32_t H_ij = get_score(i, j);
    if (is_sw && H_ij == 0) {
      break;
    }
    const std::uint32_t* row_preds = &preds[pred_offsets[i]];
    const std::uint32_t num_preds = pred_offsets[i + 1] - pred_offsets[i];
    const std::int32_t match_cost = row_profile(i)[j];
    bool predecessor_found = false;
    for (std::uint32_t p = 0; p < num_preds; ++p) {
      if (H_ij == get_score(row_preds[p], j - 1) + match_cost) {
        prev_i = row_preds[p];
        prev_j = j - 1;
        predecessor_found = true;
        break;
      }
    }
    if (!predecessor_found) {
      for (std::uint32_t p = 0; p < num_preds; ++p) {
        if (H_ij == get_score(row_preds[p], j) + g) {
          prev_i = row_preds[p];
          prev_j = j;
          predecessor_found = true;
          break;
        }
      }
    }
    if (!predecessor_found && H_ij == get_score(i, j - 1) + g) {
      prev_i = i;
      prev_j = j - 1;
    }

    alignment->emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;
  }

  std::reverse(alignment->begin(), alignment->end());
  return true;
}

Alignment SisdAlignmentEngine::AlignPruned(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t min_score,
    bool full_fallback,
    std::int32_t* score,
    bool* reached) {
  static const bool disable = std::getenv("SPOA_DISABLE_PRUNED") != nullptr;
  if (!disable && subtype_ == AlignmentSubtype::kLinear &&
      (type_ == AlignmentType::kSW || type_ == AlignmentType::kOV) &&
      sequence_len > 0 && !graph.nodes().empty() && m_ > 0 && n_ <= m_ && !full_fallback) {
    // Each aligned sequence base adds at most m to the score, so this is an upper bound
    const std::int64_t max_possible_score = static_cast<std::int64_t>(m_) *
        std::min<std::int64_t>(sequence_len, GraphLongestPath(graph));
    if (max_possible_score < min_score) {
      static const bool verify = std::getenv("SPOA_VERIFY_PRUNED") != nullptr;
      if (verify) {
        std::int32_t full_score = 0;
        Alignment full_alignment = Align(sequence, sequence_len, graph, &full_score);
        if (!full_alignment.empty() && full_score >= min_score) {
          std::fprintf(stderr,
              "[spoa::SisdAlignmentEngine::AlignPruned] error: score %d above bound %lld\n",
              full_score, static_cast<long long>(max_possible_score));
          std::abort();
        }
      }
      *reached = false;
      return Alignment();
    }
  }
  if (!disable && subtype_ == AlignmentSubtype::kLinear &&
      (type_ == AlignmentType::kSW || type_ == AlignmentType::kOV) &&
      sequence_len > 0 && !graph.nodes().empty()) {
    Alignment alignment;
    std::int32_t pruned_score = 0;
    bool aborted = false;
    const bool ok = LinearPruned(sequence, sequence_len, graph, min_score, &pruned_score,
        &alignment, &aborted);
    if (ok) {
      static const bool verify = std::getenv("SPOA_VERIFY_PRUNED") != nullptr;
      if (verify) {
        std::int32_t full_score = 0;
        Alignment full_alignment = Align(sequence, sequence_len, graph, &full_score);
        if (full_alignment != alignment || full_score != pruned_score) {
          std::fprintf(stderr,
              "[spoa::SisdAlignmentEngine::AlignPruned] error: pruned alignment mismatch "
              "(score %d vs %d, length %zu vs %zu)\n",
              pruned_score, full_score, alignment.size(), full_alignment.size());
          std::abort();
        }
      }
      if (score) {
        *score = pruned_score;
      }
      *reached = true;
      return alignment;
    }
    // Unless the pruned alignment was aborted, it found that the score is below min_score
    if (!aborted && !full_fallback) {
      static const bool verify = std::getenv("SPOA_VERIFY_PRUNED") != nullptr;
      if (verify) {
        std::int32_t full_score = 0;
        Alignment full_alignment = Align(sequence, sequence_len, graph, &full_score);
        if (!full_alignment.empty() && full_score >= min_score) {
          std::fprintf(stderr,
              "[spoa::SisdAlignmentEngine::AlignPruned] error: pruned alignment missed score %d "
              ">= %d\n", full_score, min_score);
          std::abort();
        }
      }
      *reached = false;
      return Alignment();
    }
  }
  std::int32_t full_score = 0;
  Alignment alignment = Align(sequence, sequence_len, graph, &full_score);
  if (score) {
    *score = full_score;
  }
  *reached = !alignment.empty() && full_score >= min_score;
  return alignment;
}

std::int64_t SisdAlignmentEngine::MaxPossibleScore(
    std::uint32_t sequence_len,
    const Graph& graph) {
  if (m_ <= 0 || n_ > m_ || graph.nodes().empty()) {
    return std::numeric_limits<std::int32_t>::max();
  }
  return static_cast<std::int64_t>(m_) *
      std::min<std::int64_t>(sequence_len, GraphLongestPath(graph));
}

Alignment SisdAlignmentEngine::AlignPrunedAuto(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    const double* expected_fracs,
    std::uint32_t num_expected_fracs,
    std::int32_t* score) {
  static const bool disable = std::getenv("SPOA_DISABLE_PRUNED") != nullptr;
  if (!disable && subtype_ == AlignmentSubtype::kLinear &&
      (type_ == AlignmentType::kSW || type_ == AlignmentType::kOV) &&
      sequence_len > 0 && !graph.nodes().empty() && m_ > 0 && n_ <= m_) {
    const std::int64_t max_possible_score = static_cast<std::int64_t>(m_) *
        std::min<std::int64_t>(sequence_len, GraphLongestPath(graph));
    for (std::uint32_t f = 0; f < num_expected_fracs; ++f) {
      const std::int32_t min_score = std::max<std::int32_t>(1,
          static_cast<std::int32_t>(expected_fracs[f] * max_possible_score));
      bool reached = false;
      Alignment alignment = AlignPruned(sequence, sequence_len, graph, min_score, false, score,
          &reached);
      // If pruning was abandoned the full alignment was computed, so it's the result even if its
      // score is below min_score
      if (reached || !alignment.empty()) {
        return alignment;
      }
    }
  }
  return Align(sequence, sequence_len, graph, score);
}

Alignment SisdAlignmentEngine::Affine(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    std::int32_t* F_row = &(pimpl_->F[i * matrix_width]);
    std::int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

    // update F and H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      F_row[j] = std::max(
          H_pred_row[j] + g_,
          F_pred_row[j] + e_);
      H_row[j] = H_pred_row[j - 1] + char_profile[j];
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        F_row[j] = std::max(
            F_row[j],
            std::max(
                H_pred_row[j] + g_,
                F_pred_row[j] + e_));
        H_row[j] = std::max(
            H_row[j],
            H_pred_row[j - 1] + char_profile[j]);
      }
    }

    // update E and H
    std::int32_t* E_row = &(pimpl_->E[i * matrix_width]);
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      E_row[j] = std::max(H_row[j - 1] + g_, E_row[j - 1] + e_);
      H_row[j] = std::max(H_row[j], std::max(F_row[j], E_row[j]));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          (it->outedges.empty() && j == matrix_width - 1)) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && (it->outedges.empty())) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;

  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
          pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                       H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                           H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && j != 0) {
      if ((extend_left = H_ij == pimpl_->E[i * matrix_width + j - 1] + e_) ||
                         H_ij == pimpl_->H[i * matrix_width + j - 1] + g_) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;

    if (extend_left) {
      while (true) {
        alignment.emplace_back(-1, j - 1);
        --j;
        if (pimpl_->E[i * matrix_width + j] + e_ !=
            pimpl_->E[i * matrix_width + j + 1]) {
          break;
        }
      }
    } else if (extend_up) {
      while (true) {
        bool stop = false;
        prev_i = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

          if ((stop = pimpl_->F[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + g_) ||  // NOLINT
                      pimpl_->F[i * matrix_width + j] == pimpl_->F[pred_i * matrix_width + j] + e_) {  // NOLINT
            prev_i = pred_i;
            break;
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;
        if (stop || i == 0) {
          break;
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

Alignment SisdAlignmentEngine::Convex(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    std::int32_t* F_row = &(pimpl_->F[i * matrix_width]);
    std::int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

    std::int32_t* O_row = &(pimpl_->O[i * matrix_width]);
    std::int32_t* O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

    // update F, O and H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      F_row[j] = std::max(H_pred_row[j] + g_, F_pred_row[j] + e_);
      O_row[j] = std::max(H_pred_row[j] + q_, O_pred_row[j] + c_);
      H_row[j] = H_pred_row[j - 1] + char_profile[j];
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);
      O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        F_row[j] = std::max(
            F_row[j],
            std::max(
                H_pred_row[j] + g_,
                F_pred_row[j] + e_));
        O_row[j] = std::max(
            O_row[j],
            std::max(
                H_pred_row[j] + q_,
                O_pred_row[j] + c_));
        H_row[j] = std::max(H_row[j], H_pred_row[j - 1] + char_profile[j]);
      }
    }

    // update E, Q and H
    std::int32_t* E_row = &(pimpl_->E[i * matrix_width]);
    std::int32_t* Q_row = &(pimpl_->Q[i * matrix_width]);
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      E_row[j] = std::max(H_row[j - 1] + g_, E_row[j - 1] + e_);
      Q_row[j] = std::max(H_row[j - 1] + q_, Q_row[j - 1] + c_);
      H_row[j] = std::max(
        H_row[j],
        std::max(
            std::max(F_row[j], E_row[j]),
            std::max(O_row[j], Q_row[j])));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          (it->outedges.empty() && j == matrix_width - 1)) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && it->outedges.empty()) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;

  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
        pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if ((extend_up |= H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                        H_ij == pimpl_->H[pred_i * matrix_width + j] + g_  ||
          (extend_up |= H_ij == pimpl_->O[pred_i * matrix_width + j] + c_) ||
                        H_ij == pimpl_->H[pred_i * matrix_width + j] + q_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if ((extend_up |= H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||  // NOLINT
                            H_ij == pimpl_->H[pred_i * matrix_width + j] + g_  ||  // NOLINT
              (extend_up |= H_ij == pimpl_->O[pred_i * matrix_width + j] + c_) ||  // NOLINT
                            H_ij == pimpl_->H[pred_i * matrix_width + j] + q_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && j != 0) {
      if ((extend_left |= H_ij == pimpl_->E[i * matrix_width + j - 1] + e_) ||
                          H_ij == pimpl_->H[i * matrix_width + j - 1] + g_  ||
          (extend_left |= H_ij == pimpl_->Q[i * matrix_width + j - 1] + c_) ||
                          H_ij == pimpl_->H[i * matrix_width + j - 1] + q_) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;

    if (extend_left) {
      while (true) {
        alignment.emplace_back(-1, j - 1);
        --j;
        if (pimpl_->E[i * matrix_width + j] + e_ != pimpl_->E[i * matrix_width + j + 1] &&  // NOLINT
            pimpl_->Q[i * matrix_width + j] + c_ != pimpl_->Q[i * matrix_width + j + 1]) {  // NOLINT
          break;
        }
      }
    } else if (extend_up) {
      while (true) {
        bool stop = true;
        prev_i = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

          if (pimpl_->F[i * matrix_width + j] == pimpl_->F[pred_i * matrix_width + j] + e_ ||  // NOLINT
              pimpl_->O[i * matrix_width + j] == pimpl_->O[pred_i * matrix_width + j] + c_) {  // NOLINT
            prev_i = pred_i;
            stop = false;
            break;
          }
        }
        if (stop == true) {
          for (const auto& it : rank_to_node[i - 1]->inedges) {
            std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

            if (pimpl_->F[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + g_ ||  // NOLINT
                pimpl_->O[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + q_) {  // NOLINT
              prev_i = pred_i;
              break;
            }
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;

        if (stop || i == 0) {
          break;
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

}  // namespace spoa
