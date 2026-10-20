pub mod affine2p;
pub mod aligner;
pub mod linear_dp;

pub mod wfa2 {
    //! Re-export wfa2-sys bindings
    pub use wfa2_sys::*;
}

/// Set a process-wide memory budget shared by all aligners, or zero to disable it (the default)
///
/// When an allocation would exceed the budget, the aligner waits for other aligners to free
/// memory, unless it holds the oldest claim on budget memory. This only changes when alignments
/// run, never their results.
pub fn set_memory_budget(num_bytes: u64) {
    unsafe extern "C" {
        fn mm_allocator_set_memory_budget(num_bytes: u64);
    }
    unsafe { mm_allocator_set_memory_budget(num_bytes) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use aligner::{
        AlignmentScope, AlignmentStatus, MemoryModel, WFAlignerGapAffine,
        WFAlignerGapAffine2Pieces, WFAlignerGapLinear,
    };

    /// Compress WFA2 cigar so that it's easier to read
    fn compress_cigar(cigar: &Vec<u8>) -> String {
        let mut out = String::new();

        let mut runlen = 0;
        let mut current_type = 0;

        for c in cigar {
            if current_type == *c {
                runlen += 1;
            } else {
                if runlen > 0 {
                    out += format!(
                        "{}{}",
                        runlen,
                        std::str::from_utf8(&[current_type]).unwrap()
                    )
                    .as_str();
                }
                runlen = 1;
                current_type = *c;
            }
        }
        if runlen > 0 {
            out += format!(
                "{}{}",
                runlen,
                std::str::from_utf8(&[current_type]).unwrap()
            )
            .as_str();
        }
        out
    }

    /// Reproduce basic test from library README
    #[test]
    fn test_end_to_end() {
        let mut aligner =
            WFAlignerGapAffine::new(6, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryLow);

        let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar(), b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM");
        let (a, b, c) = aligner.matching(pattern, text);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "TCTTTACTCGCGCGTT-GGAGAAATACAATAGT\n|||||||| ||||||| ||||||||||||||||\nTCTATACT-GCGCGTTTGGAGAAATAAAATAGT"
        );
    }

    /// Test align_ends_free method and new_with_match ctor
    #[test]
    fn test_ends_free() {
        let mut aligner = WFAlignerGapAffine::new_with_match(
            -1,
            3,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"CGCGTTTGGAGAA";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            text,
            pattern_size,
            pattern_size,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(compress_cigar(&aligner.cigar()), "9I13M10I");

        let pattern = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let text = b"CGCGTTTGGAGAA";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            text,
            pattern_size,
            pattern_size,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 13);
        assert_eq!(compress_cigar(&aligner.cigar()), "9D13M10D");
    }

    /// Change pattern to test test the left-right shift behavior of this library
    #[test]
    fn test_ends_free_shift() {
        let mut aligner = WFAlignerGapAffine::new_with_match(
            -1,
            3,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"TATATTTTTTTTGGAGAAATAAAATA";
        let text = b"TCTATATTTTTTTTTGGAGAAATAAAATAGT";
        let pattern_size = pattern.len() as i32;
        let text_size = text.len() as i32;
        let status = aligner.align_ends_free(
            pattern,
            text,
            pattern_size,
            pattern_size,
            text_size,
            text_size,
        );
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        //assert_eq!(aligner.score(), 18);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(
            std::str::from_utf8(aligner.cigar().as_slice()).unwrap(),
            "IIMMMMMMMMMMMMIMMMMMMMMMMMMMMII"
        );
    }

    /// Test double affine mode, and with 0 gap extension to see if a long gap
    /// is created as expected
    #[test]
    fn test_end_to_end_affine2() {
        let mut aligner = WFAlignerGapAffine2Pieces::new_with_match(
            -1,
            3,
            3,
            3,
            10,
            0,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);
        assert_eq!(compress_cigar(&aligner.cigar()), "6M21I5M");
    }

    /// Test double affine mode, and with 0 gap open
    #[test]
    fn test_end_to_end_affine2_zero_open() {
        let mut aligner = WFAlignerGapAffine2Pieces::new_with_match(
            -1,
            3,
            0,
            4,
            0,
            10,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"TCTATAATAGT";
        let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -73);
        assert_eq!(
            std::str::from_utf8(aligner.cigar().as_slice()).unwrap(),
            "MMMMMMIIIIIIIIIIIIMIIIIMMIIIIIMM"
        );
    }

    /// This test reproduces a bug found in WFA2-lib main branch at 94bcccd.
    ///
    /// A version directly in C is submitted here:
    /// https://github.com/smarco/WFA2-lib/issues/73
    ///
    #[test]
    fn test_ends_free_bug() {
        let mut aligner = WFAlignerGapLinear::new_with_match(
            -1,
            1,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );

        let pattern = b"A";
        let text = b"ACG";
        let status = aligner.align_ends_free(pattern, text, 0, 0, 0, 2);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 1);

        // bug version output
        //assert_eq!(aligner.score(), 2);

        // CIGAR output is configured for a reversed notion of pattern/text:
        assert_eq!(
            std::str::from_utf8(aligner.cigar().as_slice()).unwrap(),
            "MII"
        );
    }

    /// Test simple end to end affine gap example to verify that large gap is created as expected
    ///
    #[test]
    fn test_end_to_end_affine() {
        let mut aligner = WFAlignerGapAffine::new_with_match(
            -1,
            2,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -2);
        assert_eq!(
            std::str::from_utf8(aligner.cigar().as_slice()).unwrap(),
            "MMMIIIIIIMMM"
        );
    }

    /// Test case expected to have equal score
    #[test]
    fn test_linear() {
        let mut aligner = WFAlignerGapAffine::new_with_match(
            -1,
            2,
            0,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);

        let mut aligner = WFAlignerGapLinear::new_with_match(
            -1,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);
    }

    /// Test case expected to have equal score
    #[test]
    fn test_score_only() {
        let mut aligner = WFAlignerGapLinear::new_with_match(
            -1,
            2,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryLow,
        );

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);

        let mut aligner = WFAlignerGapLinear::new_with_match(
            -1,
            2,
            1,
            AlignmentScope::Score,
            MemoryModel::MemoryLow,
        );

        let pattern = b"ATAATA";
        let text = b"ATACATAAAATA";
        let status = aligner.align_end_to_end(pattern, text);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 0);
    }
}
