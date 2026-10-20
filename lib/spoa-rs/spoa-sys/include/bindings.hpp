#pragma once

#include <memory>
#include <string>
#include <vector>

#include "spoa-sys/spoa/include/spoa/spoa.hpp"

namespace spoa {

std::unique_ptr<spoa::AlignmentEngine>
create_alignment_engine(spoa::AlignmentType type, std::int8_t m, std::int8_t n,
                        std::int8_t g, std::int8_t e, std::int8_t q,
                        std::int8_t c);

std::unique_ptr<spoa::Alignment> create_alignment();

bool alignment_is_empty(const spoa::Alignment &alignment);

std::unique_ptr<spoa::Alignment> align(spoa::AlignmentEngine &engine,
                                       const char *sequence,
                                       std::uint32_t sequence_len,
                                       const Graph &graph,
                                       int32_t& score);

std::unique_ptr<spoa::Alignment> align_pruned(spoa::AlignmentEngine &engine,
                                              const char *sequence,
                                              std::uint32_t sequence_len,
                                              const Graph &graph,
                                              int32_t min_score,
                                              bool full_fallback,
                                              int32_t& score,
                                              bool& reached);

std::unique_ptr<spoa::Alignment> align_pruned_auto(spoa::AlignmentEngine &engine,
                                                   const char *sequence,
                                                   std::uint32_t sequence_len,
                                                   const Graph &graph,
                                                   const double *expected_fracs,
                                                   std::uint32_t num_expected_fracs,
                                                   int32_t& score);

int64_t max_possible_score(spoa::AlignmentEngine &engine, std::uint32_t sequence_len,
                           const Graph &graph);

/// Return the number of read bases used in the alignment (subtracting any clipping on the left and right sides)
///
uint32_t get_alignment_overlap_size(const Alignment& alignment);

/// Return the read base index at the start of the alignment, or -1 if no aligned base can be found
///
int32_t get_alignment_clip_size(const Alignment& alignment);

std::unique_ptr<spoa::Graph> create_graph();

void clear(spoa::Graph &graph);

void add_alignment(spoa::Graph &graph, const spoa::Alignment &alignment,
                   const char *sequence, std::uint32_t sequence_len,
                   const char *quality, std::uint32_t quality_len);

void add_alignment_noqual(spoa::Graph &graph, const spoa::Alignment &alignment,
                   const char *sequence, std::uint32_t sequence_len);

std::unique_ptr<std::string> generate_consensus(spoa::Graph &graph);

std::unique_ptr<std::vector<std::string> >
generate_multiple_sequence_alignment(spoa::Graph &graph,
                                     bool include_consensus);

uint32_t
get_sequence_count(const spoa::Graph &graph);

} // namespace spoa
