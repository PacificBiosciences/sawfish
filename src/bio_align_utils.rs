use bio::alignment::pairwise::{banded, MatchFunc, MatchParams, Scoring};
use bio::alignment::{sparse, Alignment, AlignmentOperation};
use rust_htslib::bam::record::Cigar;

use crate::simple_alignment::SimpleAlignment;

pub struct AlignmentWeights {
    pub match_: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

/// Convert rust bio alignments into sawfish's internal format
///
pub fn bio_alignment_to_simple_alignment(rb_aln: &Alignment) -> SimpleAlignment {
    fn process_last_op(
        last_op: Option<&AlignmentOperation>,
        last_op_len: u32,
        cigar: &mut Vec<Cigar>,
    ) {
        if last_op.is_none() {
            return;
        }
        let last_op = last_op.unwrap();
        if let AlignmentOperation::Yclip(_) = last_op {
            assert_eq!(last_op_len, 1);
        } else {
            use AlignmentOperation::*;
            let cigar_op = match last_op {
                Del => Cigar::Del(last_op_len),
                Ins => Cigar::Ins(last_op_len),
                Match => Cigar::Equal(last_op_len),
                Subst => Cigar::Diff(last_op_len),
                Xclip(len) => {
                    assert_eq!(last_op_len, 1);
                    Cigar::SoftClip(*len as u32)
                }
                Yclip(_) => {
                    unreachable!("Yclip state should be filtered out");
                }
            };
            cigar.push(cigar_op);
        }
    }

    // x is pattern (read) and y is text (ref)
    let ref_offset = rb_aln.ystart as i64;
    let mut cigar = Vec::new();
    let mut last_op = None;
    let mut last_op_len = 0;
    for op in rb_aln.operations.iter() {
        if Some(op) == last_op {
            last_op_len += 1;
        } else {
            process_last_op(last_op, last_op_len, &mut cigar);
            last_op = Some(op);
            last_op_len = 1;
        }
    }
    process_last_op(last_op, last_op_len, &mut cigar);

    SimpleAlignment { ref_offset, cigar }
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

impl PairwiseAligner<MatchParams> {
    /// Aligner used for determining whether to merge similar haplotypes across samples
    ///
    pub fn new_merge_aligner() -> Self {
        let k = 20;
        let w = 10;
        let scoring = Scoring::from_scores(0, -2, 1, -3).yclip(0);
        Self {
            k,
            aligner: banded::Aligner::with_scoring(scoring, k, w),
            alignment: None,
        }
    }

    /// Supports global alignment of pattern, local alignment of text
    ///
    pub fn new_glocal_aligner(weights: &AlignmentWeights) -> Self {
        Self::new_glocal_kw_aligner(weights, 20, 10)
    }

    /// Supports global alignment of pattern, local alignment of text
    ///
    pub fn new_glocal_kw_aligner(weights: &AlignmentWeights, k: usize, w: usize) -> Self {
        let scoring = Scoring::from_scores(
            weights.gap_open,
            weights.gap_extend,
            weights.match_,
            weights.mismatch,
        );
        let scoring = scoring.yclip(0);
        Self {
            k,
            aligner: banded::Aligner::with_scoring(scoring, k, w),
            alignment: None,
        }
    }

    /// First (of two) aligners used for an 'overlap' alignment where the pattern suffix extends the text.
    ///
    pub fn new_pattern_suffix_aligner(scoring: Scoring<MatchParams>, k: usize, w: usize) -> Self {
        let scoring = scoring.yclip_prefix(0).xclip_suffix(0);
        Self {
            k,
            aligner: banded::Aligner::with_scoring(scoring, k, w),
            alignment: None,
        }
    }

    /// Second (of two) aligners used for an 'overlap' alignment where the pattern prefix preceeds the text.
    ///
    pub fn new_pattern_prefix_aligner(scoring: Scoring<MatchParams>, k: usize, w: usize) -> Self {
        let scoring = scoring.yclip_suffix(0).xclip_prefix(0);
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

    pub fn get_last_alignment(&self) -> Option<SimpleAlignment> {
        let alignment = self.alignment.as_ref()?;
        Some(bio_alignment_to_simple_alignment(alignment))
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

/// Make a special interface for the overlap aligner, which needs to be run as 2 pairwise aligners
///
pub struct OverlapPairwiseAligner<F: MatchFunc> {
    k: usize,
    is_suffix: bool,
    pattern_suffix_aligner: PairwiseAligner<F>,
    pattern_prefix_aligner: PairwiseAligner<F>,
}

impl OverlapPairwiseAligner<MatchParams> {
    /// Aligner used for building a large haplotype backbone
    ///
    /// This aligner allows for an 'overlap' alignment pattern
    ///
    pub fn new_aligner(weights: &AlignmentWeights) -> Self {
        let k = 20;
        let w = 10;
        let scoring = Scoring::from_scores(
            weights.gap_open,
            weights.gap_extend,
            weights.match_,
            weights.mismatch,
        );
        let pattern_suffix_aligner = PairwiseAligner::new_pattern_suffix_aligner(scoring, k, w);
        let pattern_prefix_aligner = PairwiseAligner::new_pattern_prefix_aligner(scoring, k, w);
        Self {
            k,
            is_suffix: true,
            pattern_suffix_aligner,
            pattern_prefix_aligner,
        }
    }

    /// Return alignment score and optional alignment of pattern to text
    pub fn align(&mut self, text: &[u8], pattern: &[u8]) -> i32 {
        let matches = sparse::find_kmer_matches(pattern, text, self.k);
        //let set = hash_kmers(seq2, k);
        //find_kmer_matches_seq2_hashed(seq1, &set, k)
        let score1 = self
            .pattern_suffix_aligner
            .align_with_matches(text, pattern, &matches);
        let score2 = self
            .pattern_prefix_aligner
            .align_with_matches(text, pattern, &matches);
        self.is_suffix = score1 >= score2;
        std::cmp::max(score1, score2)
    }

    pub fn get_last_alignment(&self) -> Option<SimpleAlignment> {
        if self.is_suffix {
            &self.pattern_suffix_aligner
        } else {
            &self.pattern_prefix_aligner
        }
        .get_last_alignment()
    }

    pub fn print_pairwise_alignment(&self, text: &[u8], pattern: &[u8]) {
        let aligner = if self.is_suffix {
            &self.pattern_suffix_aligner
        } else {
            &self.pattern_prefix_aligner
        };
        print_pairwise_alignment(aligner, text, pattern);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use unwrap::unwrap;

    fn get_test_weights() -> AlignmentWeights {
        AlignmentWeights {
            match_: 1,
            mismatch: -3,
            gap_open: -1,
            gap_extend: -1,
        }
    }

    fn new_backbone_aligner() -> OverlapPairwiseAligner<MatchParams> {
        OverlapPairwiseAligner::new_aligner(&get_test_weights())
    }

    #[test]
    fn test_merge_aligner() {
        let mut aligner = PairwiseAligner::new_merge_aligner();

        // Text:    AACTGTATAA
        // Pattern: -ACT-TAT--
        let score = aligner.align("AACTGTATAA".as_bytes(), "ACTTAT".as_bytes());
        assert_eq!(score, 4);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 1);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "3=1D3=");
    }

    #[test]
    fn test_glocal_aligner() {
        let mut aligner = PairwiseAligner::new_glocal_aligner(&get_test_weights());

        // Text:    AACTGTATAA
        // Pattern: --CTGAATA-
        let score = aligner.align("AACTGTATAA".as_bytes(), "CTGAATA".as_bytes());
        assert_eq!(score, 3);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 2);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "3=1X3=");

        // Text:    ---AACTGTATAA
        // Pattern: CCGAACT------
        let score = aligner.align("AACTGTATAA".as_bytes(), "CCGAACT".as_bytes());
        assert_eq!(score, 0); // Expected to align incorrectly

        // Text:    AACTGTATAA
        // Pattern: -ACTG-ATA-
        let score = aligner.align("AACTGTATAA".as_bytes(), "ACTGATA".as_bytes());
        assert_eq!(score, 5);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 1);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "4=1D3=");
    }

    #[test]
    fn test_suffix_aligner() {
        let weights = get_test_weights();
        let scoring = Scoring::from_scores(
            weights.gap_open,
            weights.gap_extend,
            weights.match_,
            weights.mismatch,
        );
        let k = 20;
        let w = 10;
        let mut aligner = PairwiseAligner::new_pattern_suffix_aligner(scoring, k, w);

        // Text:    AACTGTATAA---
        // Pattern: ------ATAACCG
        let score = aligner.align("AACTGTATAA".as_bytes(), "ATAACCG".as_bytes());
        assert_eq!(score, 4);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 6);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "4=3S");

        // Text:    ---AACTGTATAA
        // Pattern: CCGAACT------
        let score = aligner.align("AACTGTATAA".as_bytes(), "CCGAACT".as_bytes());
        assert_eq!(score, 0); // Expected to align incorrectly
    }

    #[test]
    fn test_prefix_aligner() {
        let weights = get_test_weights();
        let scoring = Scoring::from_scores(
            weights.gap_open,
            weights.gap_extend,
            weights.match_,
            weights.mismatch,
        );
        let k = 20;
        let w = 10;
        let mut aligner = PairwiseAligner::new_pattern_prefix_aligner(scoring, k, w);

        // Text:    ---AACTGTATAA
        // Pattern: CCGAACT------
        let score = aligner.align("AACTGTATAA".as_bytes(), "CCGAACT".as_bytes());
        assert_eq!(score, 4);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 0);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "3S4=");

        // Text:    AACTGTATAA---
        // Pattern: ------ATAACCG
        let score = aligner.align("AACTGTATAA".as_bytes(), "ATAACCG".as_bytes());
        assert_eq!(score, 1); // Expected to align incorrectly
    }

    #[test]
    fn test_overlap_aligner() {
        let mut aligner = new_backbone_aligner();

        // Text:    AACTGTATAA---
        // Pattern: ------ATAACCG
        let score = aligner.align("AACTGTATAA".as_bytes(), "ATAACCG".as_bytes());
        assert_eq!(score, 4);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 6);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "4=3S");

        // Text:    ---AACTGTATAA
        // Pattern: CCGAACT------
        let score = aligner.align("AACTGTATAA".as_bytes(), "CCGAACT".as_bytes());
        assert_eq!(score, 4);
        let alignment = unwrap!(aligner.get_last_alignment());
        assert_eq!(alignment.ref_offset, 0);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "3S4=");
    }
}
