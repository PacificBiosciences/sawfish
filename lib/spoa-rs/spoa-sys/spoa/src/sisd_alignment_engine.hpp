// Copyright (c) 2020 Robert Vaser

#ifndef SISD_ALIGNMENT_ENGINE_HPP_
#define SISD_ALIGNMENT_ENGINE_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "spoa/alignment_engine.hpp"

namespace spoa {

class SisdAlignmentEngine: public AlignmentEngine {
 public:
  SisdAlignmentEngine(const SisdAlignmentEngine&) = delete;
  SisdAlignmentEngine& operator=(const SisdAlignmentEngine&) = delete;

  SisdAlignmentEngine(SisdAlignmentEngine&&) = default;
  SisdAlignmentEngine& operator=(SisdAlignmentEngine&&) = default;

  ~SisdAlignmentEngine() = default;

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  void Prealloc(
      std::uint32_t max_sequence_len,
      std::uint8_t alphabet_size) override;

  Alignment Align(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) override;

  Alignment AlignPruned(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t min_score,
      bool full_fallback,
      std::int32_t* score,
      bool* reached) override;

  std::int64_t MaxPossibleScore(
      std::uint32_t sequence_len,
      const Graph& graph) override;

  Alignment AlignPrunedAuto(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      const double* expected_fracs,
      std::uint32_t num_expected_fracs,
      std::int32_t* score) override;

 private:
  SisdAlignmentEngine(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  Alignment Linear(
      std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) noexcept;

  // Same result as Linear for kSW and kOV alignment, without storing the full
  // score matrix. Returns false if the graph can't be handled, in which case
  // Linear should be used instead.
  bool LinearLowMemory(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score,
      Alignment* alignment);

  // LinearLowMemory with scores of type T, and sequence profile storage for T
  template <typename T>
  bool LinearLowMemoryT(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score,
      Alignment* alignment,
      std::vector<T>* sequence_profile);

  bool LinearPruned(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t min_score,
      std::int32_t* score,
      Alignment* alignment,
      bool* aborted);

  template <typename T>
  bool LinearPrunedT(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t min_score,
      std::int32_t* score,
      Alignment* alignment,
      bool* aborted,
      std::vector<T>* sequence_profile);

  Alignment AlignFullMatrix(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score);

  Alignment Affine(
      std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) noexcept;

  Alignment Convex(
      std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score) noexcept;

  void Realloc(
      std::uint64_t matrix_width,
      std::uint64_t matrix_height,
      std::uint8_t num_codes);

  void Initialize(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph) noexcept;

  struct Implementation;
  std::unique_ptr<Implementation> pimpl_;
};

}  // namespace spoa

#endif  // SISD_ALIGNMENT_ENGINE_HPP_
