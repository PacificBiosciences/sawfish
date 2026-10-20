use serde::{Deserialize, Serialize};

use crate::wfa2;
use std::ptr;
use std::slice;

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryModel {
    MemoryHigh,
    MemoryMed,
    MemoryLow,
    MemoryUltraLow,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentScope {
    Score,
    Alignment,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Heuristic {
    None,
    BandedStatic(i32, i32),
    BandedAdaptive(i32, i32, i32),
    WFadaptive(i32, i32, i32),
    XDrop(i32, i32),
    ZDrop(i32, i32),
}

#[derive(Debug, PartialEq, Eq)]
pub enum AlignmentStatus {
    // OK Status (>=0)
    StatusAlgCompleted = wfa2::WF_STATUS_ALG_COMPLETED as isize,
    StatusAlgPartial = wfa2::WF_STATUS_ALG_PARTIAL as isize,
    // FAILED Status (<0)
    StatusMaxStepsReached = wfa2::WF_STATUS_MAX_STEPS_REACHED as isize,
    StatusOOM = wfa2::WF_STATUS_OOM as isize,
}

impl From<i32> for AlignmentStatus {
    fn from(value: i32) -> Self {
        match value {
            x if x == wfa2::WF_STATUS_ALG_COMPLETED as i32 => AlignmentStatus::StatusAlgCompleted,
            x if x == wfa2::WF_STATUS_ALG_PARTIAL as i32 => AlignmentStatus::StatusAlgPartial,
            wfa2::WF_STATUS_MAX_STEPS_REACHED => AlignmentStatus::StatusMaxStepsReached,
            wfa2::WF_STATUS_OOM => AlignmentStatus::StatusOOM,
            _ => panic!("Unknown alignment status: {}", value),
        }
    }
}

struct WFAttributes {
    inner: wfa2::wavefront_aligner_attr_t,
}

impl WFAttributes {
    fn default() -> Self {
        Self {
            inner: unsafe { wfa2::wavefront_aligner_attr_default },
        }
    }

    fn memory_model(mut self, memory_model: MemoryModel) -> Self {
        let memory_mode = match memory_model {
            MemoryModel::MemoryHigh => wfa2::wavefront_memory_t_wavefront_memory_high,
            MemoryModel::MemoryMed => wfa2::wavefront_memory_t_wavefront_memory_med,
            MemoryModel::MemoryLow => wfa2::wavefront_memory_t_wavefront_memory_low,
            MemoryModel::MemoryUltraLow => wfa2::wavefront_memory_t_wavefront_memory_ultralow,
        };
        self.inner.memory_mode = memory_mode;
        self
    }

    fn alignment_scope(mut self, alignment_scope: AlignmentScope) -> Self {
        let alignment_scope = match alignment_scope {
            AlignmentScope::Score => wfa2::alignment_scope_t_compute_score,
            AlignmentScope::Alignment => wfa2::alignment_scope_t_compute_alignment,
        };
        self.inner.alignment_scope = alignment_scope;
        self
    }

    fn linear_penalties(mut self, match_: i32, mismatch: i32, indel: i32) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_linear;
        self.inner.linear_penalties.match_ = match_;
        self.inner.linear_penalties.mismatch = mismatch;
        self.inner.linear_penalties.indel = indel;
        self
    }

    fn affine_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine;
        self.inner.affine_penalties.match_ = match_;
        self.inner.affine_penalties.mismatch = mismatch;
        self.inner.affine_penalties.gap_opening = gap_opening;
        self.inner.affine_penalties.gap_extension = gap_extension;
        self
    }

    fn affine2p_penalties(
        mut self,
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
    ) -> Self {
        self.inner.distance_metric = wfa2::distance_metric_t_gap_affine_2p;
        self.inner.affine2p_penalties.match_ = match_;
        self.inner.affine2p_penalties.mismatch = mismatch;
        self.inner.affine2p_penalties.gap_opening1 = gap_opening1;
        self.inner.affine2p_penalties.gap_extension1 = gap_extension1;
        self.inner.affine2p_penalties.gap_opening2 = gap_opening2;
        self.inner.affine2p_penalties.gap_extension2 = gap_extension2;
        self
    }
}

pub struct WFAligner {
    attributes: WFAttributes,
    inner: *mut wfa2::wavefront_aligner_t,
}

impl WFAligner {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> Self {
        let attributes = WFAttributes::default()
            .memory_model(memory_model)
            .alignment_scope(alignment_scope);
        Self {
            attributes,
            inner: ptr::null_mut(),
        }
    }
}

impl Drop for WFAligner {
    fn drop(&mut self) {
        unsafe {
            if !self.inner.is_null() {
                wfa2::wavefront_aligner_delete(self.inner);
            }
        }
    }
}

impl WFAligner {
    pub fn align_end_to_end(&mut self, pattern: &[u8], text: &[u8]) -> AlignmentStatus {
        let status = unsafe {
            // Configure
            wfa2::wavefront_aligner_set_alignment_end_to_end(self.inner);
            // Align
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const i8,
                pattern.len() as i32,
                text.as_ptr() as *const i8,
                text.len() as i32,
            )
        };
        AlignmentStatus::from(status)
    }

    pub fn align_ends_free(
        &mut self,
        pattern: &[u8],
        text: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) -> AlignmentStatus {
        let status = unsafe {
            // Configure
            wfa2::wavefront_aligner_set_alignment_free_ends(
                self.inner,
                pattern_begin_free,
                pattern_end_free,
                text_begin_free,
                text_end_free,
            );
            // Align
            wfa2::wavefront_align(
                self.inner,
                pattern.as_ptr() as *const i8,
                pattern.len() as i32,
                text.as_ptr() as *const i8,
                text.len() as i32,
            )
        };
        AlignmentStatus::from(status)
    }

    pub fn score(&self) -> i32 {
        unsafe { (*(*self.inner).cigar).score }
    }

    pub fn cigar(&self) -> Vec<u8> {
        let cigar_str = unsafe {
            let begin_offset = (*(*self.inner).cigar).begin_offset;
            let cigar_operations = (*(*self.inner).cigar)
                .operations
                .offset(begin_offset as isize) as *const u8;
            let cigar_length = ((*(*self.inner).cigar).end_offset - begin_offset) as usize;
            slice::from_raw_parts(cigar_operations, cigar_length)
        };
        cigar_str.to_vec()
    }

    pub fn matching(&self, pattern: &[u8], text: &[u8]) -> (String, String, String) {
        let cigar = self.cigar();
        let mut pattern_iter = pattern.iter();
        let mut text_iter = text.iter();
        let mut pattern_match = String::new();
        let mut middle_match = String::new();
        let mut text_match = String::new();
        for c in cigar.iter() {
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
                _ => panic!("Unknown cigar operation: {}", c),
            }
        }
        (pattern_match, middle_match, text_match)
    }

    /// Stop alignments that reach this score in WFA steps, with status StatusMaxStepsReached
    pub fn set_max_alignment_steps(&mut self, max_alignment_steps: i32) {
        unsafe { wfa2::wavefront_aligner_set_max_alignment_steps(self.inner, max_alignment_steps) }
    }

    pub fn set_heuristic(&mut self, heuristic: Heuristic) {
        unsafe {
            match heuristic {
                Heuristic::None => wfa2::wavefront_aligner_set_heuristic_none(self.inner),
                Heuristic::BandedStatic(band_min_k, band_max_k) => {
                    wfa2::wavefront_aligner_set_heuristic_banded_static(
                        self.inner, band_min_k, band_max_k,
                    )
                }
                Heuristic::BandedAdaptive(band_min_k, band_max_k, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_banded_adaptive(
                        self.inner,
                        band_min_k,
                        band_max_k,
                        steps_between_cutoffs,
                    )
                }
                Heuristic::WFadaptive(
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ) => wfa2::wavefront_aligner_set_heuristic_wfadaptive(
                    self.inner,
                    min_wavefront_length,
                    max_distance_threshold,
                    steps_between_cutoffs,
                ),
                Heuristic::XDrop(xdrop, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_xdrop(
                        self.inner,
                        xdrop,
                        steps_between_cutoffs,
                    )
                }
                Heuristic::ZDrop(zdrop, steps_between_cutoffs) => {
                    wfa2::wavefront_aligner_set_heuristic_zdrop(
                        self.inner,
                        zdrop,
                        steps_between_cutoffs,
                    )
                }
            }
        }
    }
}

/// Indel Aligner (a.k.a Longest Common Subsequence - LCS)
pub struct WFAlignerIndel;

impl WFAlignerIndel {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes.inner.distance_metric = wfa2::distance_metric_t_indel;
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }
}

/// Edit Aligner (a.k.a Levenshtein)
pub struct WFAlignerEdit;

impl WFAlignerEdit {
    pub fn new(alignment_scope: AlignmentScope, memory_model: MemoryModel) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes.inner.distance_metric = wfa2::distance_metric_t_edit;
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }
}

/// Gap-Linear Aligner (a.k.a Needleman-Wunsch)
pub struct WFAlignerGapLinear;

impl WFAlignerGapLinear {
    pub fn new_with_match(
        match_: i32,
        mismatch: i32,
        indel: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .linear_penalties(match_, mismatch, indel)
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn new(
        mismatch: i32,
        indel: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::new_with_match(0, mismatch, indel, alignment_scope, memory_model)
    }
}

/// Gap-Affine Aligner (a.k.a Smith-Waterman-Gotoh)
pub struct WFAlignerGapAffine;

impl WFAlignerGapAffine {
    pub fn new_with_match(
        match_: i32,
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .affine_penalties(match_, mismatch, gap_opening, gap_extension)
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn new(
        mismatch: i32,
        gap_opening: i32,
        gap_extension: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::new_with_match(
            0,
            mismatch,
            gap_opening,
            gap_extension,
            alignment_scope,
            memory_model,
        )
    }
}

/// Gap-Affine Dual-Cost Aligner (a.k.a. concave 2-pieces)
pub struct WFAlignerGapAffine2Pieces;

impl WFAlignerGapAffine2Pieces {
    pub fn new_with_match(
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        let mut aligner = WFAligner::new(alignment_scope, memory_model);
        aligner.attributes = WFAttributes::default()
            .affine2p_penalties(
                match_,
                mismatch,
                gap_opening1,
                gap_extension1,
                gap_opening2,
                gap_extension2,
            )
            .alignment_scope(alignment_scope)
            .memory_model(memory_model);
        unsafe {
            aligner.inner = wfa2::wavefront_aligner_new(&mut aligner.attributes.inner);
        }
        aligner
    }

    pub fn new(
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
        alignment_scope: AlignmentScope,
        memory_model: MemoryModel,
    ) -> WFAligner {
        Self::new_with_match(
            0,
            mismatch,
            gap_opening1,
            gap_extension1,
            gap_opening2,
            gap_extension2,
            alignment_scope,
            memory_model,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PATTERN: &[u8] = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
    const TEXT: &[u8] = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";

    fn aligner_indel() -> WFAligner {
        WFAlignerIndel::new(AlignmentScope::Alignment, MemoryModel::MemoryHigh)
    }

    fn aligner_edit() -> WFAligner {
        WFAlignerEdit::new(AlignmentScope::Alignment, MemoryModel::MemoryMed)
    }

    fn aligner_gap_linear() -> WFAligner {
        WFAlignerGapLinear::new(6, 2, AlignmentScope::Alignment, MemoryModel::MemoryHigh)
    }

    fn aligner_gap_affine() -> WFAligner {
        WFAlignerGapAffine::new(6, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryLow)
    }

    fn aligner_gap_affine_2pieces() -> WFAligner {
        WFAlignerGapAffine2Pieces::new(
            6,
            4,
            2,
            12,
            1,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        )
    }

    #[test]
    fn test_aligner_indel() {
        let mut aligner = aligner_indel();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 6);
    }

    #[test]
    fn test_aligner_edit() {
        let mut aligner = aligner_edit();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), 4);
    }

    #[test]
    fn test_aligner_gap_linear() {
        let mut aligner = aligner_gap_linear();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -12);
    }

    #[test]
    fn test_aligner_gap_affine() {
        let mut aligner = aligner_gap_affine();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
        assert_eq!(aligner.cigar(), b"MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM");
        let (a, b, c) = aligner.matching(PATTERN, TEXT);
        assert_eq!(
            format!("{}\n{}\n{}", a, b, c),
            "TCTTTACTCGCGCGTT-GGAGAAATACAATAGT\n|||||||| ||||||| ||||||||||||||||\nTCTATACT-GCGCGTTTGGAGAAATAAAATAGT"
        );
    }

    #[test]
    fn test_aligner_gap_affine_2pieces() {
        let mut aligner = aligner_gap_affine_2pieces();
        let status = aligner.align_end_to_end(PATTERN, TEXT);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        assert_eq!(aligner.score(), -24);
    }

    #[test]
    fn test_set_heuristic() {
        let mut aligner = aligner_gap_affine();
        aligner.set_heuristic(Heuristic::BandedStatic(1, 2));
        aligner.set_heuristic(Heuristic::BandedAdaptive(1, 2, 3));
        aligner.set_heuristic(Heuristic::WFadaptive(1, 2, 3));
        aligner.set_heuristic(Heuristic::XDrop(1, 2));
        aligner.set_heuristic(Heuristic::ZDrop(1, 2));
        aligner.set_heuristic(Heuristic::None);
    }
}
