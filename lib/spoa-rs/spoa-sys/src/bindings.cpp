#include <vector>

#include "spoa-sys/include/bindings.hpp"

namespace spoa {

std::unique_ptr<AlignmentEngine>
create_alignment_engine(AlignmentType typ, std::int8_t m, std::int8_t n,
                        std::int8_t g, std::int8_t e, std::int8_t q,
                        std::int8_t c) {
  return AlignmentEngine::Create(typ, m, n, g, e, q, c);
}

std::unique_ptr<Alignment> create_alignment() {
  return std::unique_ptr<Alignment>(new Alignment());
}

bool alignment_is_empty(const Alignment &alignment) { return alignment.empty(); }

std::unique_ptr<Alignment> align(AlignmentEngine &engine, const char *sequence,
                                 std::uint32_t sequence_len,
                                 const Graph &graph,
                                 int32_t& score) {
  Alignment alignment = engine.Align(sequence, sequence_len, graph, &score);
  return std::unique_ptr<Alignment>(new Alignment(std::move(alignment)));
}

std::unique_ptr<Alignment> align_pruned(AlignmentEngine &engine, const char *sequence,
                                        std::uint32_t sequence_len,
                                        const Graph &graph,
                                        int32_t min_score,
                                        bool full_fallback,
                                        int32_t& score,
                                        bool& reached) {
  Alignment alignment = engine.AlignPruned(sequence, sequence_len, graph, min_score,
                                           full_fallback, &score, &reached);
  return std::unique_ptr<Alignment>(new Alignment(std::move(alignment)));
}

std::unique_ptr<Alignment> align_pruned_auto(AlignmentEngine &engine, const char *sequence,
                                             std::uint32_t sequence_len,
                                             const Graph &graph,
                                             const double *expected_fracs,
                                             std::uint32_t num_expected_fracs,
                                             int32_t& score) {
  Alignment alignment = engine.AlignPrunedAuto(sequence, sequence_len, graph, expected_fracs,
                                               num_expected_fracs, &score);
  return std::unique_ptr<Alignment>(new Alignment(std::move(alignment)));
}

int64_t max_possible_score(AlignmentEngine &engine, std::uint32_t sequence_len,
                           const Graph &graph) {
  return engine.MaxPossibleScore(sequence_len, graph);
}

uint32_t get_alignment_overlap_size(const Alignment& alignment) {
    const int32_t start(get_alignment_clip_size(alignment));

    int32_t end(-1);
    for (auto i = alignment.rbegin(); i != alignment.rend(); ++i) {
        if (i->second >= 0) {
            end = i->second;
            break;
        }
    }

    if (start >= 0 && end >= 0) {
        return 1 + (end - start);
    } else {
        return 0;
    }
}

int32_t get_alignment_clip_size(const Alignment& alignment) {
    for (auto i = alignment.begin(); i != alignment.end(); ++i) {
        if (i->second >= 0) {
            return i->second;
        }
    }
    return -1;
}


std::unique_ptr<Graph> create_graph() {
  return std::unique_ptr<Graph>(new Graph());
}

void clear(Graph &graph) {
  graph.Clear();
}

void add_alignment(Graph &graph, const Alignment &alignment,
                   const char *sequence, std::uint32_t sequence_len,
                   const char *quality, std::uint32_t quality_len) {
  graph.AddAlignment(alignment, sequence, sequence_len, quality, quality_len);
}

void add_alignment_noqual(Graph &graph, const Alignment &alignment,
                   const char *sequence, std::uint32_t sequence_len) {
  graph.AddAlignment(alignment, sequence, sequence_len);
}

std::unique_ptr<std::string> generate_consensus(Graph &graph) {
  std::string consensus = graph.GenerateConsensus();
  return std::unique_ptr<std::string>(new std::string(std::move(consensus)));
}

std::unique_ptr<std::vector<std::string>>
generate_multiple_sequence_alignment(Graph &graph, bool include_consensus) {
  std::vector<std::string> msa =
      graph.GenerateMultipleSequenceAlignment(include_consensus);
  return std::unique_ptr<std::vector<std::string>>(
      new std::vector<std::string>(std::move(msa)));
}

// Note Graph& could be const here, but I can't figure out how to do this and also pin it on the rust side
uint32_t
get_sequence_count(const Graph &graph) {
  return graph.sequences().size();
}

} // namespace spoa
