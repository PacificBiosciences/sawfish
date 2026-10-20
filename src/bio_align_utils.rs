use bio::alignment::pairwise::{MatchFunc, MatchParams, Scoring, banded};
use bio::alignment::{Alignment, sparse};
use rust_htslib::bam::record::Cigar;

use crate::banded_score::BandedScoreAligner;

pub struct AlignmentWeights {
    pub match_: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

/// A simple alignment re-scoring utility
///
/// Assumes that cigar uses equal/diff and not match (e.g. the output of bio_alignment_to_simple_alignment)
///
#[allow(unused)]
pub fn rescore_alignment(cigar: &[Cigar], weights: &AlignmentWeights) -> i32 {
    let mut score = 0;

    use Cigar::*;
    for c in cigar.iter() {
        match c {
            Equal(len) => {
                score += weights.match_ * *len as i32;
            }
            Diff(len) => {
                score += weights.mismatch * *len as i32;
            }
            Match(_) => {
                panic!("Can't accept match cigar element");
            }
            Del(len) | Ins(len) => {
                score += weights.gap_open + weights.gap_extend * *len as i32;
            }
            _ => (),
        }
    }

    score
}

/// Hide details behind a standard-ish interface to make it easy to change between aligners
pub struct PairwiseAligner<F: MatchFunc> {
    k: usize,

    aligner: banded::Aligner<F>,

    /// Store the last alignment result
    alignment: Option<Alignment>,
}

/// Gap open, gap extend, match and mismatch scores of the merge aligner
const MERGE_ALIGNER_SCORES: (i32, i32, i32, i32) = (0, -2, 1, -3);

/// Return false only if no alignment of all of `pattern` to part of `text` with the merge aligner
/// can reach a score of `min_norm_score` times the pattern length, so the alignment can be skipped
///
/// Each edit in an alignment lowers its score from the all-match score by at least the smallest
/// of (match - mismatch) for a mismatch, (match - gap_extend) for a deleted pattern base and
/// -gap_extend for an inserted text base, with no gap open score. So an alignment with E edits
/// scores at most n * match - E * min_edit_loss for a pattern of length n. If the pattern is within
/// E edits of a substring of the text, at least n - k + 1 - k * E of the pattern's k-mers occur in
/// the text (the q-gram lemma), so counting those bounds E from below.
///
/// The merge aligner is banded, so its score is at most the best score of any alignment.
///
pub fn can_merge_alignment_reach_norm_score(
    text: &[u8],
    pattern: &[u8],
    min_norm_score: f64,
) -> bool {
    const K: usize = 16;
    let (gap_open, gap_extend, match_, mismatch) = MERGE_ALIGNER_SCORES;
    assert_eq!(gap_open, 0);
    let min_edit_loss = (match_ - mismatch)
        .min(match_ - gap_extend)
        .min(-gap_extend) as i64;
    assert!(min_edit_loss > 0);

    let n = pattern.len();
    if n < K {
        return true;
    }

    // 2-bit k-mer codes, or None for k-mers with a base other than ACGT
    fn kmers(seq: &[u8]) -> impl Iterator<Item = Option<u32>> + '_ {
        let mut code = 0u32;
        let mut valid_len = 0usize;
        seq.iter().enumerate().filter_map(move |(i, &b)| {
            let x = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };
            if x == 4 {
                valid_len = 0;
            } else {
                code = (code << 2) | x;
                valid_len += 1;
            }
            if i + 1 >= K {
                Some((valid_len >= K).then_some(code))
            } else {
                None
            }
        })
    }

    let text_kmers = kmers(text)
        .flatten()
        .collect::<std::collections::HashSet<_>>();

    // Pattern k-mers containing other bases are counted as found, which only weakens the bound
    let found = kmers(pattern)
        .filter(|x| x.is_none_or(|x| text_kmers.contains(&x)))
        .count();
    let missing = (n - K + 1 - found) as i64;
    let min_edits = (missing + K as i64 - 1) / K as i64;
    let max_score = n as i64 * match_ as i64 - min_edits * min_edit_loss;
    max_score as f64 / n as f64 >= min_norm_score
}

/// If true, run each merge alignment skipped by can_merge_alignment_reach_norm_score, and abort if
/// it reaches the merge score
pub fn verify_merge_alignment_filter() -> bool {
    static VERIFY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *VERIFY.get_or_init(|| std::env::var("SAWFISH_VERIFY_MERGE_FILTER").is_ok())
}

/// k-mer size and band width of the merge aligner
const MERGE_ALIGNER_K: usize = 20;
const MERGE_ALIGNER_W: usize = 10;

fn get_merge_aligner_scoring() -> Scoring<MatchParams> {
    let (gap_open, gap_extend, match_, mismatch) = MERGE_ALIGNER_SCORES;
    Scoring::from_scores(gap_open, gap_extend, match_, mismatch).yclip(0)
}

/// Score of `PairwiseAligner::new_merge_aligner().align(text, pattern)`, computed without the
/// alignment
///
/// rust-bio's banded aligner allocates a full traceback matrix, so this saves a lot of memory and
/// time for long sequences. See BandedScoreAligner.
pub fn merge_alignment_score(text: &[u8], pattern: &[u8]) -> i32 {
    assert!(!text.is_empty());
    assert!(!pattern.is_empty());
    let (k, w) = (MERGE_ALIGNER_K, MERGE_ALIGNER_W);
    let matches = sparse::find_kmer_matches(pattern, text, k);
    let mut aligner = BandedScoreAligner::with_scoring(get_merge_aligner_scoring(), k, w);
    let score = aligner.custom_with_matches(pattern, text, &matches);
    if verify_merge_alignment_filter() {
        let mut aligner = PairwiseAligner::new_merge_aligner();
        let expected = aligner.align(text, pattern);
        assert_eq!(
            score, expected,
            "Merge alignment score differs from rust-bio"
        );
    }
    score
}

impl PairwiseAligner<MatchParams> {
    /// Aligner used for determining whether to merge similar haplotypes across samples
    ///
    pub fn new_merge_aligner() -> Self {
        let k = MERGE_ALIGNER_K;
        let w = MERGE_ALIGNER_W;
        let scoring = get_merge_aligner_scoring();
        Self {
            k,
            aligner: banded::Aligner::with_scoring(scoring, k, w),
            alignment: None,
        }
    }

    /// Align pattern to text and return alignment score
    pub fn align(&mut self, text: &[u8], pattern: &[u8]) -> i32 {
        let matches = sparse::find_kmer_matches(pattern, text, self.k);
        self.align_with_matches(text, pattern, &matches)
    }

    /// Align pattern to text and return alignment score
    pub fn align_with_matches(
        &mut self,
        text: &[u8],
        pattern: &[u8],
        matches: &[(u32, u32)],
    ) -> i32 {
        assert!(!text.is_empty());
        assert!(!pattern.is_empty());
        self.alignment = Some(self.aligner.custom_with_matches(pattern, text, matches));
        self.alignment.as_ref().unwrap().score
    }
}

/// Handles wrapping long alignments into multiple rows
///
pub fn print_pairwise_alignment<T: MatchFunc>(
    aligner: &PairwiseAligner<T>,
    text: &[u8],
    pattern: &[u8],
) {
    const WIDTH: usize = 150;

    if let Some(alignment) = aligner.alignment.as_ref() {
        let ap = alignment.pretty(pattern, text, WIDTH);
        eprintln!("{ap}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Check that the merge alignment filter never rejects a pair that the merge aligner accepts
    #[test]
    fn test_merge_alignment_filter() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut rand = move |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let min_norm_score = 0.97;
        let mut skipped = 0;
        for iter in 0..2000 {
            let base_len = 20 + rand(3000);
            let base = (0..base_len)
                .map(|_| b"ACGTN"[rand(if iter % 7 == 0 { 5 } else { 4 })])
                .collect::<Vec<_>>();
            let start = rand(base_len / 2);
            let end = start + 1 + rand(base_len - start);
            let rate = if iter % 5 == 0 { rand(1000) } else { rand(40) };
            let mut query = Vec::new();
            for &b in base[start..end].iter() {
                if rand(1000) < rate {
                    match rand(3) {
                        0 => {}
                        1 => {
                            query.push(b);
                            query.push(b"ACGT"[rand(4)]);
                        }
                        _ => query.push(b"ACGT"[rand(4)]),
                    }
                } else {
                    query.push(b);
                }
            }
            if query.is_empty() {
                continue;
            }
            let mut aligner = PairwiseAligner::new_merge_aligner();
            let score = aligner.align(&base, &query);
            assert_eq!(merge_alignment_score(&base, &query), score, "iter {iter}");
            let norm_score = score as f64 / query.len() as f64;
            if !can_merge_alignment_reach_norm_score(&base, &query, min_norm_score) {
                skipped += 1;
                assert!(norm_score < min_norm_score, "iter {iter}");
            }
        }
        assert!(skipped > 0);
    }

    #[test]
    fn test_merge_aligner() {
        let mut aligner = PairwiseAligner::new_merge_aligner();

        // Text:    AACTGTATAA
        // Pattern: -ACT-TAT--
        let score = aligner.align(b"AACTGTATAA", b"ACTTAT");
        assert_eq!(score, 4);
    }

    #[test]
    fn test_merge_alignment_score_matches_rust_bio() {
        let mut state = 0x6a09_e667_f3bc_c908u64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let bases = b"ACGT";
        let mut aligner = PairwiseAligner::new_merge_aligner();
        for case in 0..60 {
            let len = 40 + rand(3000);
            let text = (0..len).map(|_| bases[rand(4)]).collect::<Vec<_>>();
            let pattern = match case % 5 {
                // Same sequence
                0 => text.clone(),
                // SNVs and small indels
                1 | 2 => {
                    let mut x = Vec::new();
                    let mut i = 0;
                    while i < text.len() {
                        match rand(100) {
                            0 => x.push(bases[rand(4)]),
                            1 => (0..1 + rand(4)).for_each(|_| x.push(bases[rand(4)])),
                            2 => i += 1 + rand(4),
                            _ => {}
                        }
                        if i < text.len() {
                            x.push(text[i]);
                        }
                        i += 1;
                    }
                    x
                }
                // A large indel and a different end
                3 => {
                    let cut = len / 3;
                    let mut x = text[..cut].to_vec();
                    x.extend_from_slice(&text[cut + len / 4..]);
                    x.extend((0..rand(50)).map(|_| bases[rand(4)]));
                    x
                }
                // Unrelated, or too short for any k-mer match
                _ if case % 10 == 4 => (0..len).map(|_| bases[rand(4)]).collect(),
                _ => text[..10 + rand(10)].to_vec(),
            };
            if pattern.is_empty() {
                continue;
            }
            let expected = aligner.align(&text, &pattern);
            assert_eq!(
                merge_alignment_score(&text, &pattern),
                expected,
                "case {case}"
            );
        }
    }
}
