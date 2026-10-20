// Copyright (c) 2020 Robert Vaser

#ifndef SPOA_ALIGNMENT_ENGINE_HPP_
#define SPOA_ALIGNMENT_ENGINE_HPP_

#pragma once

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace spoa {

enum class AlignmentType {
  kSW,  // Smith Waterman
  kNW,  // Needleman Wunsch
  kOV   // Overlap
};

enum class AlignmentSubtype {
  kLinear,  // g * i
  kAffine,  // g + (i - 1) * e
  kConvex   // min(g1 + (i - 1) * e1, g2 + (i - 1) * e2)
};

class Graph;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class AlignmentEngine {
 public:
  virtual ~AlignmentEngine() = default;

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,   // match
      std::int8_t n,   // mismatch
      std::int8_t g);  // gap

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,   // gap open
      std::int8_t e);  // gap extend

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,   // gap open of second affine
      std::int8_t c);  // gap extend of second affine

  virtual void Prealloc(
      std::uint32_t max_sequence_len,
      std::uint8_t alphabet_size) = 0;

  Alignment Align(
      const std::string& sequence,
      const Graph& graph,
      std::int32_t* score = nullptr);

  virtual Alignment Align(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score = nullptr) = 0;

  // Same result as Align, computed more efficiently if the alignment score is at least min_score.
  // reached is set to true if the alignment is non-empty with a score of at least min_score.
  // Without full_fallback, the alignment and score are only valid if reached is true.
  //
  // A non-empty alignment is always exactly the alignment Align would return (and so is its
  // score), whether or not reached is set: pruning only discards cells that can't be on a path
  // scoring at least min_score, and a non-empty result means such a path was found, so the end
  // cell and backtrack equal those of the full alignment. Callers may rely on this.
  //
  // Pruning requires a positive match score, a negative gap score and a mismatch score no larger
  // than the match score. Otherwise this is the same as Align.
  virtual Alignment AlignPruned(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t min_score,
      bool full_fallback,
      std::int32_t* score,
      bool* reached) {
    (void)full_fallback;
    std::int32_t full_score = 0;
    Alignment alignment = Align(sequence, sequence_len, graph, &full_score);
    if (score) {
      *score = full_score;
    }
    *reached = !alignment.empty() && full_score >= min_score;
    return alignment;
  }

  // Upper bound of the alignment score of a sequence of this length to graph
  virtual std::int64_t MaxPossibleScore(
      std::uint32_t sequence_len,
      const Graph& graph) {
    (void)sequence_len;
    (void)graph;
    return std::numeric_limits<std::int32_t>::max();
  }

  // Same result as Align, computed more efficiently if the alignment score is close to one of the
  // expected fractions of the maximum possible score (tried in the given order)
  virtual Alignment AlignPrunedAuto(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      const double* expected_fracs,
      std::uint32_t num_expected_fracs,
      std::int32_t* score) {
    (void)expected_fracs;
    (void)num_expected_fracs;
    return Align(sequence, sequence_len, graph, score);
  }

 protected:
  AlignmentEngine(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  std::int64_t WorstCaseAlignmentScore(
      std::int64_t sequence_len,
      std::int64_t graph_len) const;

  AlignmentType type_;
  AlignmentSubtype subtype_;
  std::int8_t m_;
  std::int8_t n_;
  std::int8_t g_;
  std::int8_t e_;
  std::int8_t q_;
  std::int8_t c_;
};

}  // namespace spoa

#endif  // SPOA_ALIGNMENT_ENGINE_HPP_
