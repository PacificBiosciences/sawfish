//! Utilities to help interface with the WFA2-lib via the rust-wfa2 crate
//!

use rust_htslib::bam::record::Cigar;
use rust_wfa2::affine2p::Affine2pAligner;
use rust_wfa2::aligner::{
    AlignmentScope, AlignmentStatus, Heuristic, MemoryModel, WFAligner, WFAlignerGapAffine2Pieces,
    WFAlignerGapLinear,
};
use rust_wfa2::linear_dp::LinearEndsFreeDp;

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
/// Ends-free clip lengths for an alignment of pattern to text
///
/// Start with a conservative notion of 'ends-free', by allowing much of the length of pattern or
/// text to be clipped on the ends. This assumes that a good candidate will not deviate far from a
/// good alignment in the central region of the contig and reference regions. Note that for some
/// STRs there might be large offsets from the cluster 'center'
///
fn get_ends_free_clip(pattern: &[u8], text: &[u8]) -> i32 {
    let clip_frac = 0.4;
    let p_clip = (pattern.len() as f64 * clip_frac) as i32;
    let t_clip = (text.len() as f64 * clip_frac) as i32;
    std::cmp::min(p_clip, t_clip)
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
    /// WFA2 aligner, used unless the fast aligner is set, or to verify the fast aligner
    aligner: Option<WFAligner>,
    /// Low-memory equivalent of the WFA2 aligner, if available for this alignment configuration
    fast: Option<Affine2pAligner>,
    /// Score and WFA2-style CIGAR of the last fast alignment
    fast_result: (i32, Vec<u8>),
    text: Vec<u8>,
    kind: u8,
}

/// If true, compare each fast alignment to WFA2, and abort on any difference
fn verify_fast_aligner() -> bool {
    static VERIFY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *VERIFY.get_or_init(|| std::env::var("SAWFISH_VERIFY_FAST_WFA").is_ok())
}

/// Trace memory at which the fast aligner takes a checkpoint, if set by SAWFISH_WFA_CHECKPOINT_BYTES
///
/// This is for testing only, since checkpoints never change the result.
fn fast_aligner_checkpoint_bytes() -> Option<u64> {
    static BYTES: std::sync::OnceLock<Option<u64>> = std::sync::OnceLock::new();
    *BYTES.get_or_init(|| {
        std::env::var("SAWFISH_WFA_CHECKPOINT_BYTES").ok().map(|x| {
            x.parse()
                .unwrap_or_else(|_| panic!("Invalid SAWFISH_WFA_CHECKPOINT_BYTES value '{x}'"))
        })
    })
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
        let (match_, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2) = (-1, 3, 1, 2, 60, 0);

        // The fast aligner gives the same alignments as WFA2 with less memory, see Affine2pAligner
        let piggyback = memory_model != MemoryModel::MemoryHigh;
        let mut fast = Affine2pAligner::new(
            match_, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2, piggyback,
        );
        if let Some(x) = fast_aligner_checkpoint_bytes() {
            fast.set_checkpoint_bytes(x);
        }

        let aligner = if verify_fast_aligner() {
            let mut aligner = WFAlignerGapAffine2Pieces::new_with_match(
                match_,
                mismatch,
                gap_open1,
                gap_ext1,
                gap_open2,
                gap_ext2,
                AlignmentScope::Alignment,
                memory_model,
            );
            aligner.set_heuristic(Heuristic::None);
            Some(aligner)
        } else {
            None
        };

        Self {
            aligner,
            fast: Some(fast),
            fast_result: (0, Vec::new()),
            text: Vec::new(),
            kind: if is_long_contig { b'L' } else { b'S' },
        }
    }

    /// Aligner used for consensus alignment where we only expect sequencing errors between reads
    /// sampled from the same haplotype
    ///
    /// Attempting to sync this with the scores used by spoa for assembly:
    ///
    pub fn new_consensus_aligner() -> Self {
        let aligner = new_consensus_wfa2_aligner();

        Self {
            aligner: Some(aligner),
            fast: None,
            fast_result: (0, Vec::new()),
            text: Vec::new(),
            kind: b'C',
        }
    }

    fn score(&self) -> i32 {
        if self.fast.is_some() {
            self.fast_result.0
        } else {
            self.aligner.as_ref().unwrap().score()
        }
    }

    /// WFA2-style CIGAR string of the last alignment
    fn cigar(&self) -> Vec<u8> {
        if self.fast.is_some() {
            self.fast_result.1.clone()
        } else {
            self.aligner.as_ref().unwrap().cigar()
        }
    }

    /// Ends-free alignment of the reversed pattern to the text, with the fast aligner if available
    fn align_ends_free_clips(&mut self, rp: &[u8], clips: [i32; 4]) {
        if let Some(fast) = self.fast.as_mut() {
            let x = fast.align_ends_free(rp, &self.text, clips[0], clips[1], clips[2], clips[3]);
            self.fast_result = (x.score, x.cigar);
            if let Some(aligner) = self.aligner.as_mut() {
                let status =
                    aligner.align_ends_free(rp, &self.text, clips[0], clips[1], clips[2], clips[3]);
                assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
                if aligner.score() != self.fast_result.0 || aligner.cigar() != self.fast_result.1 {
                    panic!(
                        "Fast aligner result differs from WFA2: kind {} clips {clips:?} pattern {} text {}",
                        self.kind as char,
                        String::from_utf8_lossy(rp),
                        String::from_utf8_lossy(&self.text)
                    );
                }
            }
        } else {
            let text = &self.text;
            let aligner = self.aligner.as_mut().unwrap();
            let status = aligner.align_ends_free(rp, text, clips[0], clips[1], clips[2], clips[3]);
            assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
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
        let clip = get_ends_free_clip(&rp, &self.text);
        let clips = [clip, clip, clip, clip];
        self.align_ends_free_clips(&rp, clips);

        (self.score(), process_wfa2_cigar_string(self.cigar()))
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
        let p_clip = (rp.len() as f64 * pattern_clip_frac) as i32;
        let t_clip = (self.text.len() as f64 * text_clip_frac) as i32;
        let clips = [p_clip, p_clip, t_clip, t_clip];
        self.align_ends_free_clips(&rp, clips);
        (self.score(), process_wfa2_cigar_string(self.cigar()))
    }

    #[allow(dead_code)]
    pub fn align_end_to_end(&mut self, pattern: &[u8]) -> (i32, Option<SimpleAlignment>) {
        let mut rp = pattern.to_vec();
        rp.reverse();
        assert!(
            self.fast.is_none(),
            "End-to-end alignment isn't supported by the fast aligner"
        );
        let text = &self.text;
        align_end_to_end(self.aligner.as_mut().unwrap(), &rp, text);
        (self.score(), process_wfa2_cigar_string(self.cigar()))
    }

    pub fn matching(&self, text: &[u8], pattern: &[u8]) -> (String, String, String) {
        let mut pattern_iter = pattern.iter();
        let mut text_iter = text.iter();
        let mut pattern_match = String::new();
        let mut middle_match = String::new();
        let mut text_match = String::new();
        for c in self.cigar().iter() {
            match c {
                b'M' | b'X' => {
                    pattern_match.push(*pattern_iter.next().unwrap() as char);
                    middle_match.push('|');
                    text_match.push(*text_iter.next().unwrap() as char);
                }
                b'D' => {
                    pattern_match.push(*pattern_iter.next().unwrap() as char);
                    middle_match.push(' ');
                    text_match.push('-');
                }
                b'I' => {
                    pattern_match.push('-');
                    middle_match.push(' ');
                    text_match.push(*text_iter.next().unwrap() as char);
                }
                _ => panic!("Unknown cigar operation: {c}"),
            }
        }
        (pattern_match, middle_match, text_match)
    }
}

/// Match, mismatch and indel scores of the consensus aligner
const CONSENSUS_SCORES: (i32, i32, i32) = (-1, 3, 2);

fn new_consensus_wfa2_aligner() -> WFAligner {
    let (match_, mismatch, indel) = CONSENSUS_SCORES;
    let mut aligner = WFAlignerGapLinear::new_with_match(
        match_,
        mismatch,
        indel,
        AlignmentScope::Alignment,
        MemoryModel::MemoryHigh,
    );
    aligner.set_heuristic(Heuristic::None);
    aligner
}

/// WFA2 steps after which ConsensusScoreAligner switches to the DP matrix
///
/// WFA2 work grows with the square of the score, while the DP matrix takes a fixed ~30 us for the
/// 500 bp flanks scored in genotyping. Read flanks that match an allele finish well within this
/// limit, those that don't match it can take WFA2 up to ~1 ms.
const CONSENSUS_SCORE_MAX_WFA_STEPS: i32 = 60;

/// Alignment score of `PairwiseAligner::new_consensus_aligner`, computed without the alignment
///
/// The score is the same as `PairwiseAligner::align` returns. Each alignment runs WFA2 with a step
/// limit, and alignments that reach the limit are scored by filling the DP matrix instead, see
/// LinearEndsFreeDp.
pub struct ConsensusScoreAligner {
    aligner: WFAligner,
    /// WFA2 step limit
    max_steps: i32,
    dp: LinearEndsFreeDp,
    /// Reversed text, as in PairwiseAligner
    text: Vec<u8>,
    /// Reversed pattern of the current alignment
    pattern: Vec<u8>,
    /// Aligner used to check each score if SAWFISH_VERIFY_FAST_WFA is set
    verify: Option<PairwiseAligner>,
}

impl ConsensusScoreAligner {
    pub fn new() -> Self {
        let mut aligner = new_consensus_wfa2_aligner();
        let max_steps = CONSENSUS_SCORE_MAX_WFA_STEPS;
        aligner.set_max_alignment_steps(max_steps);
        Self {
            aligner,
            max_steps,
            dp: LinearEndsFreeDp::new(),
            text: Vec::new(),
            pattern: Vec::new(),
            verify: verify_fast_aligner().then(PairwiseAligner::new_consensus_aligner),
        }
    }

    pub fn set_text(&mut self, text: &[u8]) {
        assert!(!text.is_empty());
        self.text.clear();
        self.text.extend(text.iter().rev());
        if let Some(verify) = self.verify.as_mut() {
            verify.set_text(text);
        }
    }

    /// Set the reversed pattern for an alignment, and return its ends-free clip size
    fn set_pattern(&mut self, pattern: &[u8]) -> i32 {
        assert!(!self.text.is_empty());
        assert!(!pattern.is_empty());
        self.pattern.clear();
        self.pattern.extend(pattern.iter().rev());
        get_ends_free_clip(&self.pattern, &self.text)
    }

    /// Score of the alignment of the current pattern by the DP matrix
    fn get_dp_score(&mut self, clip: i32) -> i32 {
        let (match_, mismatch, indel) = CONSENSUS_SCORES;
        self.dp.score(
            match_,
            mismatch,
            indel,
            &self.pattern,
            &self.text,
            clip,
            clip,
            clip,
            clip,
        )
    }

    /// If SAWFISH_VERIFY_FAST_WFA is set, check a score of pattern against WFA2, or check that WFA2's
    /// score is no higher than an upper bound
    fn verify_score(&mut self, pattern: &[u8], score: Option<i32>, upper_bound: Option<i32>) {
        let Some(verify) = self.verify.as_mut() else {
            return;
        };
        let expected = verify.align(pattern).0;
        let clip = self.set_pattern(pattern);
        let dp = self.get_dp_score(clip);
        if score.is_some_and(|x| x != expected)
            || upper_bound.is_some_and(|x| expected > x)
            || dp != expected
        {
            panic!(
                "Consensus score {score:?} (upper bound {upper_bound:?}, DP {dp}) differs from WFA2 {expected}: pattern {} reversed text {}",
                String::from_utf8_lossy(pattern),
                String::from_utf8_lossy(&self.text)
            );
        }
    }

    /// Alignment score, or an upper bound on it if WFA2 reaches its step limit
    ///
    /// When WFA2 stops at the step limit, no end cell can be reached with a lower penalty, and the
    /// score of an end cell (v,h) with penalty s is (v + h - s) / 2.
    ///
    pub fn align_score_bounded(&mut self, pattern: &[u8]) -> BoundedScore {
        let clip = self.set_pattern(pattern);
        let status =
            self.aligner
                .align_ends_free(&self.pattern, &self.text, clip, clip, clip, clip);
        let result = match status {
            AlignmentStatus::StatusAlgCompleted => BoundedScore::Exact(self.aligner.score()),
            AlignmentStatus::StatusMaxStepsReached => {
                let max_end_sum = (self.pattern.len() + self.text.len()) as i32;
                BoundedScore::AtMost((max_end_sum - self.max_steps).div_euclid(2))
            }
            _ => panic!("Unexpected consensus alignment status: {status:?}"),
        };
        match result {
            BoundedScore::Exact(x) => self.verify_score(pattern, Some(x), None),
            BoundedScore::AtMost(x) => self.verify_score(pattern, None, Some(x)),
        }
        result
    }

    /// Alignment score from the DP matrix, without first trying WFA2
    pub fn dp_align_score(&mut self, pattern: &[u8]) -> i32 {
        let clip = self.set_pattern(pattern);
        let score = self.get_dp_score(clip);
        self.verify_score(pattern, Some(score), None);
        score
    }

    pub fn align_score(&mut self, pattern: &[u8]) -> i32 {
        match self.align_score_bounded(pattern) {
            BoundedScore::Exact(x) => x,
            BoundedScore::AtMost(_) => self.dp_align_score(pattern),
        }
    }
}

/// Alignment score, or an upper bound on it
#[derive(Clone, Copy, Debug)]
pub enum BoundedScore {
    Exact(i32),
    AtMost(i32),
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

    /// Check that the low-memory large indel aligner gives the same result as WFA2, for both
    /// memory models, on random sequence pairs with SNVs, small indels and large indels
    #[test]
    fn test_large_indel_aligner_matches_wfa2() {
        // Simple deterministic pseudo-random generator
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let bases = b"ACGT";
        for case in 0..60 {
            let ref_len = 50 + rand(1500);
            let text = (0..ref_len).map(|_| bases[rand(4)]).collect::<Vec<_>>();

            // Derive the pattern from a section of the text with edits
            let start = rand(ref_len / 3 + 1);
            let end = ref_len - rand(ref_len / 3 + 1);
            let mut pattern = Vec::new();
            let mut i = start;
            while i < end.max(start + 1).min(ref_len) {
                match rand(400) {
                    0..=3 => pattern.push(bases[rand(4)]),
                    4..=5 => {
                        for _ in 0..1 + rand(8) {
                            pattern.push(bases[rand(4)]);
                        }
                    }
                    6..=7 => i += 1 + rand(8),
                    8 => i += 50 + rand(600),
                    9 => {
                        for _ in 0..50 + rand(600) {
                            pattern.push(bases[rand(4)]);
                        }
                    }
                    _ => pattern.push(text[i]),
                }
                i += 1;
            }
            if pattern.is_empty() {
                pattern.push(b'A');
            }
            let clip = |len: usize| (len as f64 * 0.4) as i32;
            let c = std::cmp::min(clip(pattern.len()), clip(text.len()));

            for memory_model in [MemoryModel::MemoryHigh, MemoryModel::MemoryMed] {
                let mut wfa2 = WFAlignerGapAffine2Pieces::new_with_match(
                    -1,
                    3,
                    1,
                    2,
                    60,
                    0,
                    AlignmentScope::Alignment,
                    memory_model,
                );
                wfa2.set_heuristic(Heuristic::None);
                let status = wfa2.align_ends_free(&pattern, &text, c, c, c, c);
                assert_eq!(status, AlignmentStatus::StatusAlgCompleted);

                // Check with the default checkpoint size (no checkpoints at this size), and
                // with checkpoints after every score or every few scores
                let piggyback = memory_model != MemoryModel::MemoryHigh;
                let mut fast = Affine2pAligner::new(-1, 3, 1, 2, 60, 0, piggyback);
                for checkpoint_bytes in [None, Some(1), Some(5000)] {
                    if let Some(x) = checkpoint_bytes {
                        fast.set_checkpoint_bytes(x);
                    }
                    let result = fast.align_ends_free(&pattern, &text, c, c, c, c);
                    assert_eq!(result.score, wfa2.score(), "case {case}");
                    assert!(result.cigar == wfa2.cigar(), "case {case}");
                }
            }
        }
    }

    /// Check that checkpoints don't change the result of an alignment with i32 offsets
    #[test]
    fn test_large_indel_aligner_checkpoints() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let bases = b"ACGT";
        let text = (0..26_000).map(|_| bases[rand(4)]).collect::<Vec<_>>();
        let mut pattern = text[..12_000].to_vec();
        pattern.extend_from_slice(&text[14_000..]);
        for x in pattern.iter_mut() {
            if rand(100) == 0 {
                *x = bases[rand(4)];
            }
        }
        let c = (pattern.len().min(text.len()) as f64 * 0.4) as i32;
        for piggyback in [false, true] {
            let mut fast = Affine2pAligner::new(-1, 3, 1, 2, 60, 0, piggyback);
            let expected = fast.align_ends_free(&pattern, &text, c, c, c, c);
            fast.set_checkpoint_bytes(100_000);
            let result = fast.align_ends_free(&pattern, &text, c, c, c, c);
            assert_eq!(result.score, expected.score);
            assert!(result.cigar == expected.cigar);
        }
    }

    fn xorshift(seed: u64) -> impl FnMut(usize) -> usize {
        let mut state = seed;
        move |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        }
    }

    /// Assert that Affine2pAligner gives WFA2's score and cigar for both memory models, with
    /// checkpoints after every score and with the default checkpoint size
    ///
    /// penalties are (match, mismatch, gap open 1, gap extend 1, gap open 2, gap extend 2), and
    /// ends_free is (pattern begin, pattern end, text begin, text end). Returns the score and cigar.
    ///
    fn assert_aligner_matches_wfa2(
        pattern: &[u8],
        text: &[u8],
        ends_free: (i32, i32, i32, i32),
        penalties: (i32, i32, i32, i32, i32, i32),
        memory_models: &[MemoryModel],
        label: &str,
    ) -> (i32, Vec<u8>) {
        let (m, x, o1, e1, o2, e2) = penalties;
        let (pb, pe, tb, te) = ends_free;
        let mut result_score_cigar = (0, Vec::new());
        for &memory_model in memory_models {
            let mut wfa2 = WFAlignerGapAffine2Pieces::new_with_match(
                m,
                x,
                o1,
                e1,
                o2,
                e2,
                AlignmentScope::Alignment,
                memory_model,
            );
            wfa2.set_heuristic(Heuristic::None);
            let status = wfa2.align_ends_free(pattern, text, pb, pe, tb, te);
            assert_eq!(status, AlignmentStatus::StatusAlgCompleted, "{label}");

            let piggyback = memory_model != MemoryModel::MemoryHigh;
            let mut fast = Affine2pAligner::new(m, x, o1, e1, o2, e2, piggyback);
            for checkpoint_bytes in [None, Some(1)] {
                if let Some(x) = checkpoint_bytes {
                    fast.set_checkpoint_bytes(x);
                }
                let result = fast.align_ends_free(pattern, text, pb, pe, tb, te);
                assert_eq!(result.score, wfa2.score(), "{label}");
                assert!(result.cigar == wfa2.cigar(), "{label}");
                result_score_cigar = (result.score, result.cigar);
            }
        }
        result_score_cigar
    }

    /// The score the aligner works with (the score of WFA2's adjusted penalties) for an alignment
    /// scored with match score -1: sequence lengths aligned minus twice the alignment score
    fn internal_score(score: i32, cigar: &[u8]) -> i32 {
        let count = |ops: &[u8]| cigar.iter().filter(|x| ops.contains(x)).count() as i32;
        // D consumes the pattern and I the text
        count(b"MXD") + count(b"MXI") - 2 * score
    }

    /// Copy of text with SNVs and small indels, at a rate of per_thousand per 1000 bases
    fn mutate_seq(
        text: &[u8],
        per_thousand: usize,
        rand: &mut impl FnMut(usize) -> usize,
    ) -> Vec<u8> {
        let bases = b"ACGT";
        let mut x = Vec::with_capacity(text.len());
        let mut i = 0;
        while i < text.len() {
            if rand(1000) < per_thousand {
                match rand(4) {
                    0 | 1 => x.push(bases[rand(4)]),
                    2 => (0..1 + rand(5)).for_each(|_| x.push(bases[rand(4)])),
                    _ => i += rand(5),
                }
            }
            if i < text.len() {
                x.push(text[i]);
            }
            i += 1;
        }
        x
    }

    /// Sequences too long for 16-bit offsets, with long exact matches between edits and large indels
    #[test]
    fn test_large_indel_aligner_i32_lengths_matches_wfa2() {
        let mut rand = xorshift(0x1234_5678_9abc_def1);
        let bases = b"ACGT";
        for case in 0..3 {
            let text = (0..26_000 + rand(2000))
                .map(|_| bases[rand(4)])
                .collect::<Vec<_>>();
            // Few edits, so that matches between them are far longer than 255 bases, plus a
            // deletion and an insertion of 300 bases (WFA2 memory grows with the indel size)
            let mut pattern = mutate_seq(&text[..11_000], 1, &mut rand);
            pattern.extend((0..300).map(|_| bases[rand(4)]));
            pattern.extend(mutate_seq(&text[11_300..], 1, &mut rand));
            // The aligner uses 16-bit offsets only if the lengths sum to at most 48 kb
            assert!(pattern.len() + text.len() > 48_000);
            // Ends are free only over a short distance, because WFA2 uses tens of GB when they
            // are free over sequences this long
            let c = 500;
            assert_aligner_matches_wfa2(
                &pattern,
                &text,
                (c, c, c, c),
                (-1, 3, 1, 2, 60, 0),
                &[MemoryModel::MemoryHigh, MemoryModel::MemoryMed],
                &format!("i32 lengths case {case}"),
            );
        }
    }

    /// Score of an alignment from its cigar, for penalties (match, mismatch, gap open 1, gap extend
    /// 1, gap open 2, gap extend 2), where a match scores -match and a gap of length n costs the
    /// lower of the two gap models
    fn score_from_cigar(cigar: &[u8], penalties: (i32, i32, i32, i32, i32, i32)) -> i32 {
        let (m, x, o1, e1, o2, e2) = penalties;
        let mut score = 0;
        let mut i = 0;
        while i < cigar.len() {
            let op = cigar[i];
            let run = cigar[i..].iter().take_while(|&&y| y == op).count() as i32;
            score += match op {
                b'M' => -m * run,
                b'X' => -x * run,
                _ => -std::cmp::min(o1 + e1 * run, o2 + e2 * run),
            };
            i += run as usize;
        }
        score
    }

    /// Sequences short enough for 16-bit offsets, but with a score above the 16-bit limit, which
    /// restarts the alignment with 32-bit offsets
    ///
    /// WFA2 needs tens of GB for alignments like this, so the restarted result is checked without
    /// it. The score must equal the score of its own cigar, and both backtrace models and
    /// checkpoints must agree.
    ///
    /// The limit is 32600 - text length, about 10,000 for the lengths here, and about 5% edits
    /// gives a score the aligner works with of about 11,000.
    #[test]
    fn test_large_indel_aligner_restart_to_i32() {
        let mut rand = xorshift(0x0fed_cba9_8765_4321);
        let bases = b"ACGT";
        let penalties = (-1, 3, 1, 2, 60, 0);
        let (m, x, o1, e1, o2, e2) = penalties;
        for case in 0..2 {
            let text = (0..22_000 + 200 * case)
                .map(|_| bases[rand(4)])
                .collect::<Vec<_>>();
            let pattern = mutate_seq(&text, 50, &mut rand);
            assert!(pattern.len() + text.len() <= 48_000);

            let mut results = Vec::new();
            for piggyback in [false, true] {
                let mut aligner = Affine2pAligner::new(m, x, o1, e1, o2, e2, piggyback);
                // A checkpoint after every score would keep gigabytes of them
                for checkpoint_bytes in [None, Some(1_000_000)] {
                    if let Some(bytes) = checkpoint_bytes {
                        aligner.set_checkpoint_bytes(bytes);
                    }
                    results.push(aligner.align_ends_free(&pattern, &text, 500, 500, 500, 500));
                }
            }
            let first = &results[0];
            assert!(
                first.score == score_from_cigar(&first.cigar, penalties),
                "case {case}"
            );
            assert!(internal_score(first.score, &first.cigar) > 32_600 - text.len() as i32);
            for result in results.iter() {
                assert_eq!(result.score, first.score, "case {case}");
                assert!(result.cigar == first.cigar, "case {case}");
            }
        }
    }

    /// Penalties other than sawfish's, and empty and tiny inputs
    #[test]
    fn test_large_indel_aligner_penalties_and_tiny_inputs_match_wfa2() {
        let mut rand = xorshift(0x7777_aaaa_5555_3333);
        let bases = b"ACGT";
        let penalties = [
            (-1, 3, 1, 2, 60, 0),
            (-2, 5, 4, 3, 40, 1),
            (-1, 4, 6, 2, 24, 1),
            (-3, 2, 2, 2, 20, 2),
        ];
        for (pi, pen) in penalties.iter().enumerate() {
            for case in 0..20 {
                let len = 1 + rand(300);
                let text = (0..len).map(|_| bases[rand(4)]).collect::<Vec<_>>();
                let mut pattern = mutate_seq(&text, 100, &mut rand);
                if pattern.is_empty() {
                    pattern.push(b'A');
                }
                let c = if case % 2 == 0 {
                    0
                } else {
                    (pattern.len().min(text.len()) / 3) as i32
                };
                assert_aligner_matches_wfa2(
                    &pattern,
                    &text,
                    (c, c, c, c),
                    *pen,
                    &[MemoryModel::MemoryHigh, MemoryModel::MemoryMed],
                    &format!("penalties {pi} case {case}"),
                );
            }
            // One base against one base, and one base against many
            for (p, t) in [
                (&b"A"[..], &b"A"[..]),
                (b"A", b"C"),
                (b"A", b"ACGTACGTAC"),
                (b"ACGTACGTAC", b"T"),
            ] {
                assert_aligner_matches_wfa2(
                    p,
                    t,
                    (0, 0, 0, 0),
                    *pen,
                    &[MemoryModel::MemoryHigh, MemoryModel::MemoryMed],
                    &format!("tiny penalties {pi}"),
                );
            }
        }
    }

    /// Check ConsensusScoreAligner against WFA2 and its own DP matrix, for patterns that match
    /// the text, are mutated, or are unrelated, which makes WFA2 reach its step limit
    #[test]
    fn test_consensus_score_aligner_matches_wfa2() {
        let mut rand = xorshift(0x3141_5926_5358_9793);
        let bases = b"ACGT";
        let mut scorer = ConsensusScoreAligner::new();
        let mut reference = PairwiseAligner::new_consensus_aligner();
        let mut at_most_count = 0;
        for case in 0..60 {
            let text_len = 60 + rand(500);
            let text = (0..text_len).map(|_| bases[rand(4)]).collect::<Vec<_>>();
            scorer.set_text(&text);
            reference.set_text(&text);
            for kind in 0..4 {
                let pattern = match kind {
                    0 => {
                        let start = rand(text_len / 2);
                        text[start..start + 20 + rand(text_len / 2 - 20)].to_vec()
                    }
                    1 => mutate_seq(&text, 20 + rand(60), &mut rand),
                    2 => (0..20 + rand(text_len))
                        .map(|_| bases[rand(4)])
                        .collect::<Vec<_>>(),
                    _ => mutate_seq(&text[..text_len / 2 + 1], 150, &mut rand),
                };
                let expected = reference.align(&pattern).0;
                let label = format!("case {case} kind {kind}");
                match scorer.align_score_bounded(&pattern) {
                    BoundedScore::Exact(x) => assert_eq!(x, expected, "{label}"),
                    BoundedScore::AtMost(x) => {
                        at_most_count += 1;
                        assert!(expected <= x, "{label}: bound {x} below score {expected}");
                    }
                }
                assert_eq!(scorer.align_score(&pattern), expected, "{label}");
                assert_eq!(scorer.dp_align_score(&pattern), expected, "{label}");
            }
        }
        // The cases must include scores found only by the DP matrix
        assert!(at_most_count > 10, "{at_most_count}");
    }
}
