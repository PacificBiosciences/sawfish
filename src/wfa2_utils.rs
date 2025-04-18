//! Utilities to help interface with the WFA2-lib via the rust-wfa2 crate
//!

use rust_htslib::bam::record::Cigar;
use rust_wfa2::aligner::{
    AlignmentScope, AlignmentStatus, Heuristic, MemoryModel, WFAligner, WFAlignerGapAffine2Pieces,
    WFAlignerGapLinear,
};

use crate::simple_alignment::SimpleAlignment;
use crate::utils::pairwise_alignment_printer;

fn process_wfa2_cigar_string(wfa2_cigar: Vec<u8>) -> Option<SimpleAlignment> {
    // The WFA2 'cigar' is a string of characters from the set [MXID]
    // To turn this into a conventional CIGAR:
    // 'M' maps to '='
    // 'X' maps to 'X'
    // 'D' maps to an interior 'I' or if on the edge translate this into 'S'
    // 'I' maps to an interior 'D' or if on the edge translate this into ref offset
    struct Wfa2Processor {
        last_state: u8,
        state_count: u32,
        first_al_match: bool,
        al: SimpleAlignment,
    }

    impl Wfa2Processor {
        fn new() -> Self {
            Self {
                last_state: 0,
                state_count: 0,
                first_al_match: false,
                al: SimpleAlignment::new(),
            }
        }

        fn process_core(&mut self, final_state: bool) {
            if self.state_count == 0 {
                return;
            }
            use Cigar::*;
            let c = match self.last_state {
                b'M' => Some(Equal(self.state_count)),
                b'X' => Some(Diff(self.state_count)),
                b'D' => {
                    if !self.first_al_match || final_state {
                        Some(SoftClip(self.state_count))
                    } else {
                        Some(Ins(self.state_count))
                    }
                }
                b'I' => {
                    if !self.first_al_match {
                        self.al.ref_offset += self.state_count as i64;
                        None
                    } else if final_state {
                        None
                    } else {
                        Some(Del(self.state_count))
                    }
                }
                _ => {
                    panic!("Unexpected wfa2 cigar string content");
                }
            };
            if let Some(c) = c {
                self.al.cigar.push(c);
            }
        }

        fn process(&mut self, state: u8) {
            if state != self.last_state {
                self.process_core(state == 0);
                self.state_count = 0;
            }
            if !self.first_al_match && (state == b'M' || state == b'X') {
                self.first_al_match = true;
            }
            self.state_count += 1;
            self.last_state = state;
        }
    }

    let mut wp = Wfa2Processor::new();
    for &state in wfa2_cigar.iter().rev() {
        wp.process(state);
    }
    wp.process(0);

    if !wp.first_al_match {
        None
    } else {
        Some(wp.al)
    }
}

/// Translate wfa2 alignment results into a standard format used within sawfish
///
/// wfa2 can generate alignments without any matches (ie. just "D" and "I" ), so check for this and
/// return None if it occurs
///
fn process_wfa2_alignment(aligner: &WFAligner) -> Option<SimpleAlignment> {
    process_wfa2_cigar_string(aligner.cigar())
}

fn align_ends_free(aligner: &mut WFAligner, pattern: &[u8], text: &[u8]) {
    // Start with a conservative notion of 'ends-free', by allowing much of the length of
    // pattern or text to be clipped on the ends. This assumes that a good candidate will not
    // deviate far from a good alignment in the central region of the contig and reference
    // regions. Note that for some STRs there might be large offsets from the cluster 'center'
    //
    let clip_frac = 0.4;
    let p_clip = (pattern.len() as f64 * clip_frac) as i32;
    let t_clip = (text.len() as f64 * clip_frac) as i32;
    let clip = std::cmp::min(p_clip, t_clip);
    let status = aligner.align_ends_free(pattern, text, clip, clip, clip, clip);
    assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
}

fn align_ends_free_clip_frac(
    aligner: &mut WFAligner,
    pattern: &[u8],
    text: &[u8],
    text_clip_frac: f64,
    pattern_clip_frac: f64,
) {
    let p_clip = (pattern.len() as f64 * pattern_clip_frac) as i32;
    let t_clip = (text.len() as f64 * text_clip_frac) as i32;
    //let clip = std::cmp::min(p_clip, t_clip);
    let status = aligner.align_ends_free(pattern, text, p_clip, p_clip, t_clip, t_clip);
    assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
}

fn align_end_to_end(aligner: &mut WFAligner, pattern: &[u8], text: &[u8]) {
    let status = aligner.align_end_to_end(pattern, text);
    assert!(
        status == AlignmentStatus::StatusAlgCompleted
            || status == AlignmentStatus::StatusAlgPartial
    );
}

/// Hide details behind a standard-ish interface to make it easy to change aligners in future
pub struct PairwiseAligner {
    aligner: WFAligner,
    text: Vec<u8>,
}

impl PairwiseAligner {
    pub fn new_large_indel_aligner(is_long_contig: bool) -> Self {
        // For this case we need a convex alignment that allows very long indels.
        // We also lower the mismatch cost compared to assembly to allow for SNPs

        // Adjust memory mode depending on anticipated contig size
        let memory_model = if is_long_contig {
            MemoryModel::MemoryMed
        } else {
            MemoryModel::MemoryHigh
        };
        let mut aligner = WFAlignerGapAffine2Pieces::new_with_match(
            -1,
            3,
            1,
            2,
            60,
            0,
            AlignmentScope::Alignment,
            memory_model,
        );
        aligner.set_heuristic(Heuristic::None);

        Self {
            aligner,
            text: Vec::new(),
        }
    }

    /// Aligner used for consensus alignment where we only expect sequencing errors between reads
    /// sampled from the same haplotype
    ///
    /// Attempting to sync this with the scores used by spoa for assembly:
    ///
    pub fn new_consensus_aligner() -> Self {
        let mut aligner = WFAlignerGapLinear::new_with_match(
            -1,
            3,
            2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        aligner.set_heuristic(Heuristic::None);

        Self {
            aligner,
            text: Vec::new(),
        }
    }

    pub fn set_text(&mut self, text: &[u8]) {
        assert!(!text.is_empty());

        self.text = text.to_vec();
        self.text.reverse();
    }

    pub fn text(&self) -> &[u8] {
        &self.text
    }

    pub fn align(&mut self, pattern: &[u8]) -> (i32, Option<SimpleAlignment>) {
        assert!(!self.text.is_empty());
        assert!(!pattern.is_empty());

        let mut rp = pattern.to_vec();
        rp.reverse();
        align_ends_free(&mut self.aligner, &rp, &self.text);
        let score = self.aligner.score();

        (score, process_wfa2_alignment(&self.aligner))
    }

    #[allow(dead_code)]
    pub fn align_clip_frac(
        &mut self,
        pattern: &[u8],
        text_clip_frac: f64,
        pattern_clip_frac: f64,
    ) -> (i32, Option<SimpleAlignment>) {
        assert!(!self.text.is_empty());
        assert!(!pattern.is_empty());

        let mut rp = pattern.to_vec();
        rp.reverse();
        align_ends_free_clip_frac(
            &mut self.aligner,
            &rp,
            &self.text,
            text_clip_frac,
            pattern_clip_frac,
        );
        let score = self.aligner.score();
        (score, process_wfa2_alignment(&self.aligner))
    }

    #[allow(dead_code)]
    pub fn align_end_to_end(&mut self, pattern: &[u8]) -> (i32, Option<SimpleAlignment>) {
        let mut rp = pattern.to_vec();
        rp.reverse();
        align_end_to_end(&mut self.aligner, &rp, &self.text);
        let score = self.aligner.score();
        (score, process_wfa2_alignment(&self.aligner))
    }

    pub fn matching(&self, text: &[u8], pattern: &[u8]) -> (String, String, String) {
        self.aligner.matching(pattern, text)
    }
}

/// Handles wrapping long alignments into multiple rows
///
pub fn print_pairwise_alignment(aligner: &PairwiseAligner, pattern: &[u8]) {
    const WIDTH: usize = 150;

    let mut pattern = pattern.to_vec();
    pattern.reverse();
    let (l1, _, l3) = aligner.matching(aligner.text(), &pattern);

    let mut l1 = l1.as_bytes().to_vec();
    l1.reverse();
    let mut l3 = l3.as_bytes().to_vec();
    l3.reverse();

    pairwise_alignment_printer(WIDTH, &[&l3, &l1]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    #[test]
    fn test_process_wfa2_cigar_string() {
        let mut wfa2_cigar = b"DDDDMMMXXMMMMIIIMMMM".to_vec();
        wfa2_cigar.reverse();
        let al = process_wfa2_cigar_string(wfa2_cigar);

        assert!(al.is_some());

        let al = al.unwrap();

        assert_eq!(CigarString(al.cigar.clone()).to_string(), "4S3=2X4=3D4=");
        assert_eq!(al.ref_offset, 0);
    }

    #[test]
    fn test_consensus_aligner() {
        let mut aligner = PairwiseAligner::new_consensus_aligner();
        aligner.set_text(b"AACTGTATAA");
        let (score, alignment) = aligner.align(b"ACTTAT");

        // Expect:
        // AACTGTATAA
        // -ACT-TAT--
        assert_eq!(score, 4);
        assert!(alignment.is_some());
        let alignment = alignment.unwrap();
        assert_eq!(alignment.ref_offset, 1);
        assert_eq!(CigarString(alignment.cigar.clone()).to_string(), "3=1D3=");
    }
}
