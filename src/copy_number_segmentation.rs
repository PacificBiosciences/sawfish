use std::collections::BTreeMap;

use camino::Utf8Path;
use itertools::Itertools;
use log::info;
use rust_vc_utils::ChromList;
use serde::{Deserialize, Serialize};
use statrs::distribution::Discrete;
use strum::EnumCount;
use unwrap::unwrap;

use crate::bam_scanner::SampleAlignmentScanResult;
use crate::cli;
use crate::depth_bins::{DepthBin, GenomeDepthBins, get_bin_index, get_complete_bin_count};
use crate::discover::{COPYNUM_SEGMENT_BEDGRAPH_FILENAME, COPYNUM_SEGMENT_MESSAGEPACK_FILENAME};
use crate::expected_ploidy::get_expected_copy_number_info_for_regions;
use crate::gc_correction::{GCBiasCorrectionData, GenomeGCLevels, SampleGCBiasCorrectionData};
use crate::genome_regions::GenomeRegionsByChromIndex;
use crate::genome_segment::GenomeSegment;
use crate::int_range::IntRange;

#[derive(Debug)]
pub struct HaploidCoverage {
    pub depth: f64,
    pub gc_corrected_depth: f64,
}

/// Get whole genome haploid coverage estimates
///
/// Averages the coverage over non-excluded non-zero depth bins in chromosomes matching
/// coverage_est_regex.
///
/// GC-correction factor is applied such that haploid depth estimate should
/// reflect the depth of the least biased GC window.
///
/// Also produces the corresponding midpoint_depth estimates, but uncorrected and gc-corrected
///
/// Returns full HaploidCoverage struct with all haploid depth estimates
///
pub fn get_haploid_genome_coverage(
    coverage_est_regex: &str,
    chrom_list: &ChromList,
    genome_depth_bins: &GenomeDepthBins,
    genome_gc_levels: &GenomeGCLevels,
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
    copy_number_segments: Option<&SampleCopyNumberSegments>,
) -> HaploidCoverage {
    use regex::Regex;
    let chrom_include_regex = Regex::new(coverage_est_regex).unwrap();

    // This application requires the reciprocal of the depth reduction factor used for the
    // emission prob correction, so go ahead and setup the reciprocal value array here:
    let gc_depth_correction = sample_gc_bias_data
        .gc_depth_reduction
        .iter()
        .map(|x| 1.0 / x)
        .collect::<Vec<_>>();

    let mut total = 0.0;
    let mut gc_corrected_total = 0.0;

    // True if any chromosome names match the chromosome regex
    //
    // This flag is used to improve error message if user specifies an invalid chrom regex.
    //
    let mut chrom_match = false;

    // The total of the predicted copy number at each bin, meaning for a bin covering a 2-copy region in the
    // sample we would add 2.
    let mut count = 0.0;

    for (chrom_index, chrom_entry) in chrom_list.data.iter().enumerate() {
        if !chrom_include_regex.is_match(chrom_entry.label.as_str()) {
            continue;
        }
        chrom_match = true;

        let chrom_depth_bins = &genome_depth_bins.depth_bins[chrom_index];
        let chrom_gc_levels = &genome_gc_levels[chrom_index];
        let bin_count = chrom_depth_bins.len();
        let chrom_cn_segments = {
            match copy_number_segments {
                Some(x) => &x.cn_segments[chrom_index],
                None => {
                    let seg = CopyNumberSegment {
                        begin_bin: 0,
                        end_bin: bin_count,
                        copy_number_info: ExtendedCopyNumberState {
                            state: CopyNumberState::Two,
                            high_state_copy_number: None,
                        },
                    };
                    &vec![seg]
                }
            }
        };
        for cn_segment in chrom_cn_segments
            .iter()
            .filter(|x| x.copy_number_info.state != CopyNumberState::Unknown)
        {
            let segment_copy_number = cn_segment.copy_number_info.state as usize as f64;
            for bin_index in cn_segment.begin_bin..cn_segment.end_bin {
                let depth = &chrom_depth_bins[bin_index];
                let gc_level = chrom_gc_levels[bin_index];
                if let DepthBin::Depth(depth) = depth {
                    if *depth > 0.0 {
                        let gc_correction = gc_depth_correction[gc_level];
                        total += depth;
                        gc_corrected_total += depth * gc_correction;

                        count += segment_copy_number;
                    }
                };
            }
        }
    }

    assert!(
        chrom_match,
        "Diploid chromosome regex '{}' does not match any sample chromosome names.",
        coverage_est_regex
    );
    assert!(
        count > 0.0,
        "No usable bins found for haploid depth estimation."
    );

    let depth = total / count;
    let gc_corrected_depth = gc_corrected_total / count;

    HaploidCoverage {
        depth,
        gc_corrected_depth,
    }
}

#[allow(clippy::needless_range_loop)]
/// Compute and return the HMM state transition matrix as log prob values
///
/// # Arguments
/// * `go_prob` - probability of changing from one state to any other state
///
/// Matrix returned with lookup format `matrix[fromState][toState]`
///
fn get_copynum_transition_probs(transition_prob: f64) -> Vec<Vec<f64>> {
    let state_count = CopyNumberState::COUNT;

    let stay_prob = 1.0 - (transition_prob * ((state_count - 1) as f64));
    let transition_lnprob = transition_prob.ln();
    let stay_lnprob = stay_prob.ln();

    //clippy doesn't like this state indexing, but changing it would be counter-intuitive IMO
    let mut tm = vec![vec![0.0; state_count]; state_count];
    for from_index in 0..state_count {
        for to_index in 0..state_count {
            tm[from_index][to_index] = {
                if from_index == to_index {
                    stay_lnprob
                } else {
                    transition_lnprob
                }
            };
        }
    }
    tm
}

#[derive(PartialEq)]
enum CopyNumTransitionType {
    Loss,
    Gain,
}

/// Get table of transition probabilities between states for bins which have been annotated with a breakpoint
///
/// Breakpoints can either suggest that we expect a CN gain or a CN loss. The `cn_trans_type` argument here
/// is used to indicate which table to compute and return.
///
#[allow(clippy::needless_range_loop)]
fn get_breakpoint_copynum_transition_probs(
    transition_prob: f64,
    breakpoint_transition_prob: f64,
    cn_trans_type: CopyNumTransitionType,
) -> Vec<Vec<f64>> {
    let state_count = CopyNumberState::COUNT;

    let mut tm = vec![vec![0.0; state_count]; state_count];
    for from_index in 0..state_count {
        for to_index in 0..state_count {
            if from_index == to_index {
                // fill in the diagonal later
                continue;
            }
            // Determine if from->to represents the corresponding CN change:
            let is_target_cn_type = if from_index == CopyNumberState::Unknown as usize {
                false
            } else {
                match cn_trans_type {
                    CopyNumTransitionType::Loss => from_index > to_index,
                    CopyNumTransitionType::Gain => to_index > from_index,
                }
            };

            tm[from_index][to_index] = if is_target_cn_type {
                breakpoint_transition_prob
            } else {
                transition_prob
            };
        }
    }

    // Fill in diagonal
    for from_index in 0..state_count {
        let mut sum = 0.0;
        for to_index in 0..state_count {
            sum += tm[from_index][to_index];
        }
        tm[from_index][from_index] = 1.0 - sum;
    }

    // Finally transform the whole to log probs
    for from_index in 0..state_count {
        for to_index in 0..state_count {
            tm[from_index][to_index] = tm[from_index][to_index].ln();
        }
    }
    tm
}

/// This struct manages all transition scores that the segmentation routine needs to access
///
/// All transition matrices are x[from state][to state]
///
pub struct TransitionProbInfo {
    pub standard_transition_lnprob_matrix: Vec<Vec<f64>>,
    pub breakend_copynum_loss_transition_lnprob_matrix: Vec<Vec<f64>>,
    pub breakend_copynum_gain_transition_lnprob_matrix: Vec<Vec<f64>>,
}

impl TransitionProbInfo {
    /// * `transition_prob` - Standard probability of transitioning between copy number states
    ///
    /// * `breakpoint_transition_factor` - At breakpoint bins, the transition prob is raised to this power
    ///
    pub fn new(transition_prob: f64, breakpoint_transition_factor: f64) -> Self {
        assert!((0.0..=1.0).contains(&transition_prob));
        assert!((0.0..=1.0).contains(&breakpoint_transition_factor));

        let breakpoint_transition_prob = f64::powf(transition_prob, breakpoint_transition_factor);

        let std_transition_lnprob_matrix = get_copynum_transition_probs(transition_prob);
        let breakend_copynum_loss_transition_lnprob_matrix =
            get_breakpoint_copynum_transition_probs(
                transition_prob,
                breakpoint_transition_prob,
                CopyNumTransitionType::Loss,
            );
        let breakend_copynum_gain_transition_lnprob_matrix =
            get_breakpoint_copynum_transition_probs(
                transition_prob,
                breakpoint_transition_prob,
                CopyNumTransitionType::Gain,
            );

        Self {
            standard_transition_lnprob_matrix: std_transition_lnprob_matrix,
            breakend_copynum_loss_transition_lnprob_matrix,
            breakend_copynum_gain_transition_lnprob_matrix,
        }
    }
}

/// A state for every copy number bin describing the transition probability to use when
/// transitioning from the previous bin to the indexed bin.
///
#[derive(Clone, PartialEq)]
pub enum CopyNumBinTransitionType {
    Standard,
    BreakendGain,
    BreakendLoss,

    // This state is set if there's a conflicting gain/loss breakend on the same bin
    Complex,
}

/// Vector over all bins in one chromosome
pub type ChromTransitionTypes = Vec<CopyNumBinTransitionType>;

/// SV breakend-based CN bin transition modifications for one sample
///
/// The bin-index where a transition type is annotated should be the bin after the boosted copy number transition
///
#[derive(Clone)]
pub struct SampleTransitionTypes {
    bin_size: u32,

    /// Vector over all chromosomes
    pub sample_transition_types: Vec<ChromTransitionTypes>,
}

impl SampleTransitionTypes {
    pub fn new(chrom_list: &ChromList, bin_size: u32) -> Self {
        let mut sample_transition_types = Vec::new();
        for chrom_info in chrom_list.data.iter() {
            let bin_count = get_complete_bin_count(chrom_info.length, bin_size);
            let chrom_transition_types = vec![CopyNumBinTransitionType::Standard; bin_count];
            sample_transition_types.push(chrom_transition_types);
        }

        Self {
            bin_size,
            sample_transition_types,
        }
    }

    fn get_bin_type(
        &mut self,
        chrom_index: usize,
        pos: i64,
    ) -> Option<&mut CopyNumBinTransitionType> {
        let bin_index = get_bin_index(pos as u64, self.bin_size);

        let ctt = &mut self.sample_transition_types[chrom_index];
        if bin_index >= ctt.len() {
            None
        } else {
            Some(&mut ctt[bin_index])
        }
    }

    pub fn set_gain_pos(&mut self, chrom_index: usize, pos: i64) {
        // Adjust pos-to-bin mapping so that transition score adjustment results in a transitioning bin boundary
        // as close to the input pos value as possible.
        let pos = pos + self.bin_size as i64 / 2;
        if let Some(btt) = self.get_bin_type(chrom_index, pos) {
            use CopyNumBinTransitionType::*;
            if *btt == Standard {
                *btt = BreakendGain;
            } else if *btt == BreakendLoss {
                *btt = Complex;
            }
        }
    }

    pub fn set_loss_pos(&mut self, chrom_index: usize, pos: i64) {
        // Adjust pos-to-bin mapping so that transition score adjustment results in a transitioning bin boundary
        // as close to the input pos value as possible.
        let pos = pos + self.bin_size as i64 / 2;
        if let Some(btt) = self.get_bin_type(chrom_index, pos) {
            use CopyNumBinTransitionType::*;
            if *btt == Standard {
                *btt = BreakendLoss;
            } else if *btt == BreakendGain {
                *btt = Complex;
            }
        }
    }
}

/// * bin_index - refers to the bin following the copy number change
pub fn get_transition_lnprob(
    transition_info: &TransitionProbInfo,
    sample_transition_types: Option<&SampleTransitionTypes>,
    chrom_index: usize,
    bin_index: usize,
    from_state: usize,
    to_state: usize,
) -> f64 {
    let trans_matrix = if let Some(sample_trans_types) = sample_transition_types {
        use CopyNumBinTransitionType::*;
        match sample_trans_types.sample_transition_types[chrom_index][bin_index] {
            Standard | Complex => &transition_info.standard_transition_lnprob_matrix,
            BreakendGain => &transition_info.breakend_copynum_gain_transition_lnprob_matrix,
            BreakendLoss => &transition_info.breakend_copynum_loss_transition_lnprob_matrix,
        }
    } else {
        &transition_info.standard_transition_lnprob_matrix
    };

    trans_matrix[from_state][to_state]
}

fn get_expected_depth(
    gc_corrected_haploid_coverage: f64,
    copy_num_state: CopyNumberState,
    depth: f64,
    gc_depth_reduction: f64,
) -> f64 {
    let default_expected_depth =
        || gc_corrected_haploid_coverage * gc_depth_reduction * copy_num_state as usize as f64;

    // Special handlers to modify expectedDepth for certain states:
    match copy_num_state {
        CopyNumberState::High => {
            // This state accounts for anything at copy number HIGH or higher.
            // A hacky way to handle this state is to set expected depth to match observed depth,
            // so long as this does not lower expected depth below that for copy number HIGH.
            default_expected_depth().max(depth.round())
        }
        CopyNumberState::Zero => {
            // Don't use zero for expected depth, instead give the model a way to explain occasional
            // garbage alignments inside zero copy regions.
            let min_copy_num = 0.01;
            gc_corrected_haploid_coverage * min_copy_num
        }
        _ => default_expected_depth(),
    }
}

/// Return the log emission prob of the depth bin observation given the copy number state
///
/// # Arguments
///
/// * `gc_bias_reduction` - Depth is multiplied by this factor to reflect reduction from theoretical full depth associated
///   with the locale gc-bias level.
///
pub fn get_depth_emission_lnprob(
    gc_corrected_haploid_coverage: f64,
    copy_num_state: CopyNumberState,
    depth_bin_value: &DepthBin,
    gc_depth_reduction: f64,
) -> f64 {
    use statrs::distribution::Poisson;

    // For excluded regions, set unknown state emission probability to 1.0, and all other states to
    // some factor below that. The factor will determine how long of a exclusion gap can be spanned
    // without interruption.
    //
    let excluded_region_unknown_state_prob = 0f64;
    let excluded_region_known_state_prob = (0.5f64).ln();

    let depth = match depth_bin_value {
        DepthBin::Depth(depth) => {
            match copy_num_state {
                CopyNumberState::Unknown => {
                    // Disallow Unknown state unless the depth bin is excluded
                    return f64::NEG_INFINITY;
                }
                _ => *depth,
            }
        }
        DepthBin::Excluded => {
            return match copy_num_state {
                CopyNumberState::Unknown => excluded_region_unknown_state_prob,
                _ => excluded_region_known_state_prob,
            };
        }
    };

    assert!(gc_corrected_haploid_coverage >= 0.0);
    assert!(depth >= 0.0);

    // Set a max copy number to retain numerical stability in depth spikes:
    let max_copy_number = 20;
    let max_depth = gc_corrected_haploid_coverage * max_copy_number as f64;
    let depth = depth.min(max_depth);
    let expected_depth = get_expected_depth(
        gc_corrected_haploid_coverage,
        copy_num_state,
        depth,
        gc_depth_reduction,
    );

    let pd = Poisson::new(expected_depth).unwrap();

    let low_depth = depth as u64;
    let low_prob = pd.pmf(low_depth);
    let delta = depth - low_depth as f64;

    match delta {
        x if x < 0.001 => low_prob,
        _ => {
            let high_depth = low_depth + 1;
            let high_prob = pd.pmf(high_depth);
            delta * high_prob + (1.0 - delta) * low_prob
        }
    }
    .ln()
}

fn max_index<T: std::cmp::PartialOrd>(x: &[T]) -> usize {
    assert!(!x.is_empty());
    let mut mi = 0;
    for i in 1..x.len() {
        if x[i] > x[mi] {
            mi = i;
        }
    }
    mi
}

/// Backtrace to get viterbi parse
fn get_backtrace(last_row: &[f64], back_pointer: &[Vec<u8>]) -> Vec<u8> {
    let mut max_state = max_index(last_row);

    let obs_count = back_pointer.len();
    let mut max_path: Vec<u8> = vec![0; obs_count];
    for obs_index in (0..obs_count).rev() {
        max_path[obs_index] = max_state as u8;
        max_state = back_pointer[obs_index][max_state] as usize;
    }
    max_path
}

/// Get the init probs over all copy number states for the given sample and chromosome
fn get_sample_chrom_init_probs<'a>(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    init_probs: &'a [Vec<f64>],
    depth_bin_size: u32,
    chrom_index: usize,
) -> &'a Vec<f64> {
    let first_segment = GenomeSegment {
        chrom_index,
        range: IntRange::from_pair(0, depth_bin_size as i64),
    };
    let expected_cn_info =
        get_expected_copy_number_info_for_regions(expected_copy_number_regions, &[first_segment]);

    &init_probs[expected_cn_info.expected_copy_number as usize]
}

/// A Viterbi parse for copy number states of a single chromosome over all input samples
///
/// All prob values are tracked in log space. This parse is extensively customized for the specific
/// copy number inference case.
///
/// Returns the viterbi bin parse converted into contiguous copy number segments for all input samples
///
#[allow(clippy::needless_range_loop)]
fn viterbi_copy_number_parse(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    sample_seg_input: &SampleCopyNumberSegmentationInput,
    init_probs: &[Vec<f64>],
    depth_bin_size: u32,
    chrom_index: usize,
) -> Vec<CopyNumberSegment> {
    let chrom_gc_levels = &genome_gc_levels[chrom_index];

    // # Arguments
    //
    // * `gc_corrected_haploid_coverage` - Estimated depth per chromosome in this sample after GC-bias corrections are applied
    //
    // * `init` - Log prob prior on copy number states [S]
    //
    // * `sample_transition_types` - Describes different transition states for each bin, which can be used to look up different
    //   transition probabilities from `transtition_info` for each bin.
    //
    let observations = &sample_seg_input.genome_depth_bins.depth_bins[chrom_index];
    let gc_corrected_haploid_coverage = sample_seg_input.gc_corrected_haploid_coverage;
    let sample_gc_bias_data = sample_seg_input.sample_gc_bias_data;
    let sample_transition_types = sample_seg_input.sample_transition_types;
    let bin_dependency_correction_factor =
        sample_seg_input.settings.bin_dependency_correction_factor;
    let init = get_sample_chrom_init_probs(
        sample_seg_input.expected_copy_number_regions,
        init_probs,
        depth_bin_size,
        chrom_index,
    );

    assert!(gc_corrected_haploid_coverage >= 0.0);

    let state_count = CopyNumberState::COUNT;
    assert_eq!(init.len(), state_count);

    let obs_count = observations.len();
    if obs_count == 0 {
        return Vec::new();
    }

    // Instead of having a full SxO DP matrix, just ping-pong on two rows
    let mut max_pr_row1 = vec![0.0; state_count];
    let mut max_pr_row2 = vec![0.0; state_count];

    let mut back_pointer = vec![vec![0u8; state_count]; obs_count];

    for state_index in 0..state_count {
        let emit_lnprob = get_depth_emission_lnprob(
            gc_corrected_haploid_coverage,
            CopyNumberState::from_repr(state_index).unwrap(),
            &observations[0],
            sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[0]],
        );
        max_pr_row1[state_index] = init[state_index] + emit_lnprob;
    }

    let this_row = &mut max_pr_row1;
    let last_row = &mut max_pr_row2;
    for obs_index in 1..obs_count {
        std::mem::swap(this_row, last_row);
        for (to_state_index, row_value) in this_row.iter_mut().enumerate() {
            let emit_lnprob = get_depth_emission_lnprob(
                gc_corrected_haploid_coverage,
                CopyNumberState::from_repr(to_state_index).unwrap(),
                &observations[obs_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[obs_index]],
            );

            let emit_lnprob = emit_lnprob * bin_dependency_correction_factor;

            let mut max_index = 0;
            let mut max_lnprob = 0.0;
            for from_state_index in 0..state_count {
                let transition_ln_prob = get_transition_lnprob(
                    transition_info,
                    sample_transition_types,
                    chrom_index,
                    obs_index,
                    from_state_index,
                    to_state_index,
                );

                let lnprob = last_row[from_state_index] + transition_ln_prob + emit_lnprob;

                if (from_state_index == 0) || (lnprob > max_lnprob) {
                    max_index = from_state_index;
                    max_lnprob = lnprob;
                }
            }

            *row_value = max_lnprob;
            back_pointer[obs_index][to_state_index] = max_index as u8;
        }
    }

    let max_cn_path = get_backtrace(this_row, &back_pointer);

    get_chrom_copy_number_segments(&max_cn_path)
}

/// Translate the most likely copy number state path into copy number segments for one chromosome
///
fn get_chrom_copy_number_segments(max_cn_path: &[u8]) -> Vec<CopyNumberSegment> {
    let mut chrom_segments: Vec<CopyNumberSegment> = Vec::new();

    if max_cn_path.is_empty() {
        return chrom_segments;
    }

    let mut begin_bin: usize = 0;
    let chrom_segments_ref = &mut chrom_segments;

    let mut add_segment = |last_cn_state, end_bin| {
        assert!(end_bin > begin_bin);
        chrom_segments_ref.push(CopyNumberSegment {
            begin_bin,
            end_bin,
            copy_number_info: ExtendedCopyNumberState {
                state: last_cn_state,
                high_state_copy_number: None,
            },
        });
        begin_bin = end_bin;
    };

    let mut last_cn_state = CopyNumberState::Zero;
    for (bin_index, &bin_value) in max_cn_path.iter().enumerate() {
        let this_cn_state = CopyNumberState::from_repr(bin_value as usize).unwrap();
        if (bin_index > 0) && (this_cn_state != last_cn_state) {
            add_segment(last_cn_state, bin_index);
        }
        last_cn_state = this_cn_state;
    }

    let bin_count = max_cn_path.len();
    add_segment(last_cn_state, bin_count);

    chrom_segments
}

/// Get init probs for copy number segmentation
///
/// The first index of init probs is the expected initial state index, this can be used to select the correct init prob for
/// each chromosome.
///
/// The inner vector is the actual init_probs, selected based on the expected cn for a spectic sample+chromosome
///
fn get_init_probs() -> Vec<Vec<f64>> {
    let state_count = CopyNumberState::COUNT;
    let off_expected_init_state_prob = 0.001;
    let expected_init_state_prob = 1.0 - (state_count as f64 * off_expected_init_state_prob);
    let off_expected_init_state_lnprob = off_expected_init_state_prob.ln();
    let expected_init_state_lnprob = expected_init_state_prob.ln();

    let mut init_probs = Vec::new();
    for init_state_index in 0..state_count {
        let mut state_init_probs = vec![off_expected_init_state_lnprob; state_count];
        state_init_probs[init_state_index] = expected_init_state_lnprob;
        init_probs.push(state_init_probs);
    }
    init_probs
}

#[derive(Clone)]
struct CopyNumberBoundary {
    sample_index: usize,

    /// First bin affected by (ie. after) the copy-number change (zero-indexed)
    bin_index: usize,

    /// Index of the segment immediately before the bin change
    before_segment_index: usize,

    before_state: CopyNumberState,
    after_state: CopyNumberState,
}

/// Prune all cnbs that are too far away form any other cnbs to be adjusted
///
/// This assumes the cnb map is all one type (loss or gain)
///
fn prune_copy_number_boundaries(max_bin_shift: usize, cnbs: &mut CopyNumberBoundarySet) {
    let keys = cnbs.keys().cloned().collect::<Vec<_>>();

    let mut is_prev_key_far = true;
    for (key, next_key) in keys.into_iter().chain([0]).tuple_windows() {
        let is_next_key_far = next_key == 0 || (next_key - key) > max_bin_shift;
        if is_prev_key_far && is_next_key_far {
            cnbs.remove(&key);
        }
        is_prev_key_far = is_next_key_far;
    }
}

/// Return the segmentation score delta shifting from_cnb_set to to_bin, divided by the number of shifted samples
///
fn cnb_score_shift_per_sample(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    from_bin: usize,
    from_cnb_set: &[CopyNumberBoundary],
    to_bin: usize,
) -> f64 {
    // Hard code parameters for now
    let shifted_cnb_transition_factor = 0.5;

    if from_bin == to_bin {
        return 0.0;
    }

    let (start_bin, end_bin, after_state_is_baseline) = if from_bin < to_bin {
        (from_bin, to_bin, true)
    } else {
        (to_bin, from_bin, false)
    };

    let mut score_delta = 0.0;

    let chrom_gc_levels = &genome_gc_levels[chrom_index];

    // First find the difference in emit scores the difference between current and revised copy number
    for cnb in from_cnb_set {
        let sample_seg_input = &seg_input[cnb.sample_index];
        let chrom_depth_bins = &sample_seg_input.genome_depth_bins.depth_bins[chrom_index];
        let sample_gc_bias_data = sample_seg_input.sample_gc_bias_data;

        let mut sample_emit_score_delta = 0.0;
        for bin_index in start_bin..end_bin {
            let after_state_emit_lnprob = get_depth_emission_lnprob(
                sample_seg_input.gc_corrected_haploid_coverage,
                cnb.after_state,
                &chrom_depth_bins[bin_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]],
            );

            let before_state_emit_lnprob = get_depth_emission_lnprob(
                sample_seg_input.gc_corrected_haploid_coverage,
                cnb.before_state,
                &chrom_depth_bins[bin_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]],
            );

            let (baseline_emit_lnprob, shifted_emit_lnprob) = if after_state_is_baseline {
                (after_state_emit_lnprob, before_state_emit_lnprob)
            } else {
                (before_state_emit_lnprob, after_state_emit_lnprob)
            };

            sample_emit_score_delta += shifted_emit_lnprob - baseline_emit_lnprob;
        }
        // At this point we have the delta for emit values from all samples, so the bin_dependency correction factor applies to the entire score:
        let bin_dependency_correction_factor =
            sample_seg_input.settings.bin_dependency_correction_factor;
        sample_emit_score_delta *= bin_dependency_correction_factor;

        score_delta += sample_emit_score_delta;
    }

    // Next, find the difference in transition scores, accounting for the shared transition site bonus
    for cnb in from_cnb_set {
        let sample_seg_input = &seg_input[cnb.sample_index];
        let sample_transition_types = sample_seg_input.sample_transition_types;

        let baseline_transition_score = get_transition_lnprob(
            transition_info,
            sample_transition_types,
            chrom_index,
            from_bin,
            cnb.before_state as usize,
            cnb.after_state as usize,
        );
        let shifted_transition_score = get_transition_lnprob(
            transition_info,
            sample_transition_types,
            chrom_index,
            to_bin,
            cnb.before_state as usize,
            cnb.after_state as usize,
        );

        let sample_transition_score_delta =
            (shifted_transition_score * shifted_cnb_transition_factor) - baseline_transition_score;
        score_delta += sample_transition_score_delta;
    }

    // return the score change per shifted sample
    score_delta / from_cnb_set.len() as f64
}

/// Given a set of CNBs within the close distance threshold, test possible consolidating shifts and
/// report out the recommended changes as a set of from->to bin index transitions.
///
#[allow(clippy::needless_range_loop)]
fn refine_close_cnb_set(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    cnb_set: &[(usize, &Vec<CopyNumberBoundary>)],
    debug: bool,
) -> Vec<(usize, usize)> {
    if debug {
        eprintln!(
            "refine_close_cnb_set: start with cnb_set_len {}",
            cnb_set.len()
        );
        for (key, val) in cnb_set.iter() {
            eprintln!("cnb key {key} cnb_count {}", val.len());
        }
    }

    // Ideally we could assert on this being len 2 or more, as a practical solution just allow this input but bail out
    if cnb_set.len() < 2 {
        return Vec::new();
    }

    // Get full set of shift scores that are possible between cnb_set members
    #[derive(PartialEq, PartialOrd)]
    struct ShiftScore {
        score: f64,
        from_key_index: usize,
        to_key_index: usize,
    }
    let mut shift_scores = Vec::new();

    let cnb_location_count = cnb_set.len();
    for cnb_location_head_index in 0..cnb_location_count {
        let (key1, cnb1) = cnb_set[cnb_location_head_index];

        for cnb_location_compare_index in (cnb_location_head_index + 1)..cnb_location_count {
            let (key2, cnb2) = cnb_set[cnb_location_compare_index];

            let score_1to2 = cnb_score_shift_per_sample(
                genome_gc_levels,
                transition_info,
                seg_input,
                chrom_index,
                key1,
                cnb1,
                key2,
            );

            shift_scores.push(ShiftScore {
                score: score_1to2,
                from_key_index: cnb_location_head_index,
                to_key_index: cnb_location_compare_index,
            });

            let score_2to1 = cnb_score_shift_per_sample(
                genome_gc_levels,
                transition_info,
                seg_input,
                chrom_index,
                key2,
                cnb2,
                key1,
            );

            shift_scores.push(ShiftScore {
                score: score_2to1,
                from_key_index: cnb_location_compare_index,
                to_key_index: cnb_location_head_index,
            });
        }
    }

    // Convert shift_scores to cnb_shifts:
    shift_scores.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let mut cnb_shifts = Vec::new();
    let mut shift_from_locations = vec![false; cnb_location_count];
    let mut shift_to_locations = vec![false; cnb_location_count];

    for shift_score in shift_scores {
        if shift_score.score <= 0.0 {
            break;
        }

        if shift_from_locations[shift_score.from_key_index] {
            // We've shifted out of this location already so can't shift out again
            continue;
        }

        if shift_from_locations[shift_score.to_key_index] {
            // We've shifted out of this location already so don't shift into it now
            continue;
        }

        if shift_to_locations[shift_score.from_key_index] {
            // We've shifted into this location already so don't shift out of it now
            continue;
        }

        let from_key = cnb_set[shift_score.from_key_index].0;
        let to_key = cnb_set[shift_score.to_key_index].0;
        cnb_shifts.push((from_key, to_key));
        shift_from_locations[shift_score.from_key_index] = true;
        shift_to_locations[shift_score.to_key_index] = true;
    }

    cnb_shifts
}

/// Determine and report small syncronization shifts in the copy number boundaries over all samples
///
fn get_copy_number_boundary_shifts(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    max_bin_shift: usize,
    cnbs: &CopyNumberBoundarySet,
    debug: bool,
) -> Vec<(usize, usize)> {
    let mut cnb_shifts = Vec::new();

    // Build CNBs into close sets, then process each set
    let mut close_cnb_set = Vec::new();
    let keys = cnbs.keys().cloned().collect::<Vec<_>>();
    for (key, next_key) in keys.into_iter().chain([0]).tuple_windows() {
        close_cnb_set.push((key, cnbs.get(&key).unwrap()));
        if next_key == 0 || (next_key - key) > max_bin_shift {
            let mut close_cnb_set_shifts = refine_close_cnb_set(
                genome_gc_levels,
                transition_info,
                seg_input,
                chrom_index,
                &close_cnb_set,
                debug,
            );

            cnb_shifts.append(&mut close_cnb_set_shifts);
            close_cnb_set.clear();
        }
    }

    cnb_shifts
}

fn update_chrom_cn_segment_shifts(
    cnbs: &CopyNumberBoundarySet,
    cnb_shifts: &[(usize, usize)],
    chrom_cn_segments: &mut [Vec<CopyNumberSegment>],
    debug: bool,
) {
    for (from_bin_index, to_bin_index) in cnb_shifts.iter() {
        let from_bin_cnbs = cnbs.get(from_bin_index).unwrap();
        if debug {
            eprintln!(
                "update_chrom_cn_segment_shifts: formi/toi/cnbs_len: {}/{}/{}",
                from_bin_index,
                to_bin_index,
                from_bin_cnbs.len()
            );
        }

        for from_bin_cnb in from_bin_cnbs.iter() {
            let sample_chrom_cn_segments = &mut chrom_cn_segments[from_bin_cnb.sample_index];
            let before_cn_segment = &sample_chrom_cn_segments[from_bin_cnb.before_segment_index];
            let after_cn_segment = &sample_chrom_cn_segments[from_bin_cnb.before_segment_index + 1];

            if from_bin_index < to_bin_index {
                let bin_shift = to_bin_index - from_bin_index;
                // Check that after_cn_segment has room for the shift
                if after_cn_segment.begin_bin + bin_shift < after_cn_segment.end_bin {
                    sample_chrom_cn_segments[from_bin_cnb.before_segment_index].end_bin +=
                        bin_shift;
                    sample_chrom_cn_segments[from_bin_cnb.before_segment_index + 1].begin_bin +=
                        bin_shift;
                }
            } else {
                let bin_shift = from_bin_index - to_bin_index;
                // Check that before_cn_segment has room for the shift
                if before_cn_segment.end_bin > bin_shift
                    && before_cn_segment.end_bin - bin_shift > before_cn_segment.begin_bin
                {
                    sample_chrom_cn_segments[from_bin_cnb.before_segment_index].end_bin -=
                        bin_shift;
                    sample_chrom_cn_segments[from_bin_cnb.before_segment_index + 1].begin_bin -=
                        bin_shift;
                }
            }
        }
    }
}

/// map key should be CopyNumberBoundary::bin_index
type CopyNumberBoundarySet = BTreeMap<usize, Vec<CopyNumberBoundary>>;

/// 'losses' mean the boundary goes from copy number N to N-1, while transitioning from bin_index to bin_index+1
struct CopyNumberBoundaryInfo {
    losses: CopyNumberBoundarySet,
    gains: CopyNumberBoundarySet,
}

fn get_copy_number_boundary_info(
    chrom_cn_segments: &mut [Vec<CopyNumberSegment>],
) -> CopyNumberBoundaryInfo {
    let mut losses = BTreeMap::new();
    let mut gains = BTreeMap::new();

    for (sample_index, sample_chrom_cn_segments) in chrom_cn_segments.iter().enumerate() {
        let segment_count = sample_chrom_cn_segments.len();
        for segment_index in 1..segment_count {
            let before_segment_index = segment_index - 1;
            let after_segment_index = segment_index;
            let before_seg = &sample_chrom_cn_segments[before_segment_index];
            let after_seg = &sample_chrom_cn_segments[after_segment_index];

            let before_state = before_seg.copy_number_info.state;
            let after_state = after_seg.copy_number_info.state;
            if before_state == CopyNumberState::Unknown || after_state == CopyNumberState::Unknown {
                continue;
            }

            // This should never happen
            if before_state == after_state {
                continue;
            }

            let cnb = CopyNumberBoundary {
                sample_index,
                bin_index: after_seg.begin_bin,
                before_segment_index,
                before_state,
                after_state,
            };

            let cnbs = if before_state as usize > after_state as usize {
                &mut losses
            } else {
                &mut gains
            };

            let x = cnbs.entry(cnb.bin_index).or_insert_with(Vec::new);
            x.push(cnb);
        }
    }

    CopyNumberBoundaryInfo { losses, gains }
}

/// Refine copy number segmentation boundaries to sync them across multiple samples
///
/// The goal here is to give a very minor bonus for syncing the bin at which same copy
/// number shift direction occurs in multiple samples. The small bonus encourages
/// shifting boundary differences caused by sampling noise only, with the goal that true
/// differences can be preserved.
///
fn refine_multi_sample_chrom_cn_segments(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    depth_bin_size: u32,
    chrom_cn_segments: &mut [Vec<CopyNumberSegment>],
) {
    let debug = false;
    if debug {
        eprintln!(
            "rms: starting refine_multi_sample_chrom_cn_segments for chrom_index {chrom_index}"
        );
    }

    let sample_count = chrom_cn_segments.len();

    // The refine operation is only relevant to multi-sample analysis, so just pass-through
    // single-sample
    if sample_count < 2 {
        return;
    }

    // Hard code parameters:
    //

    // Max distance in a segment boundary can be shifted in a single sync operation:
    let max_shift = 5000;

    // Translate max shift into a maximum bin shift:
    let max_bin_shift = (max_shift / depth_bin_size) as usize;

    // Test for another early exist condition
    if max_bin_shift < 1 {
        return;
    }

    let mut cnb_info = get_copy_number_boundary_info(chrom_cn_segments);

    if debug {
        eprintln!("rms: before filtration");
        eprintln!("rms: cnb_loss_bin_count: {}", cnb_info.losses.len());
        eprintln!("rms: cnb_gain_bin_count: {}", cnb_info.gains.len());
    }

    prune_copy_number_boundaries(max_bin_shift, &mut cnb_info.losses);
    prune_copy_number_boundaries(max_bin_shift, &mut cnb_info.gains);

    if debug {
        eprintln!("rms: after filtration");
        eprintln!("rms: cnb_loss_bin_count: {}", cnb_info.losses.len());
        for x in cnb_info.losses.keys() {
            eprintln!("rms: cnb_loss_bin: {:?}", x);
        }
        eprintln!("rms: cnb_gain_bin_count: {}", cnb_info.gains.len());
        for x in cnb_info.gains.keys() {
            eprintln!("rms: cnb_gain_bin: {:?}", x);
        }
    }

    let cnb_loss_shifts = get_copy_number_boundary_shifts(
        genome_gc_levels,
        transition_info,
        seg_input,
        chrom_index,
        max_bin_shift,
        &cnb_info.losses,
        debug,
    );

    let cnb_gain_shifts = get_copy_number_boundary_shifts(
        genome_gc_levels,
        transition_info,
        seg_input,
        chrom_index,
        max_bin_shift,
        &cnb_info.gains,
        debug,
    );

    if debug {
        eprintln!("rms: cnb_loss_shifts: {}", cnb_loss_shifts.len());
        for x in cnb_loss_shifts.iter() {
            eprintln!("rms: cnb_loss_shift: {:?}", x);
        }
        eprintln!("rms: cnb_gain_shifts: {}", cnb_gain_shifts.len());
        for x in cnb_gain_shifts.iter() {
            eprintln!("rms: cnb_gain_shift: {:?}", x);
        }
    }

    // Final step is to adjust the original segments to reflect the cnb_shifts:
    update_chrom_cn_segment_shifts(&cnb_info.losses, &cnb_loss_shifts, chrom_cn_segments, debug);
    update_chrom_cn_segment_shifts(&cnb_info.gains, &cnb_gain_shifts, chrom_cn_segments, debug);
}

/// All per-sample inputs needed for copy-number segmentation
pub struct SampleCopyNumberSegmentationInput<'a> {
    pub settings: &'a cli::DiscoverSettings,
    pub gc_corrected_haploid_coverage: f64,
    pub genome_depth_bins: &'a GenomeDepthBins,
    pub sample_gc_bias_data: &'a SampleGCBiasCorrectionData,
    pub sample_transition_types: Option<&'a SampleTransitionTypes>,
    pub expected_copy_number_regions: Option<&'a GenomeRegionsByChromIndex>,
}

/// Run a joint copy-number segmentation over all samples to predict copy number state bins, and further reduce these into segments of continuous copy number
/// for each sample
///
/// This is the simple version of copy number segmentation which runs once on a single haploid coverage estimate
///
pub fn get_single_pass_sample_copy_number_segments(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
) -> Vec<SampleCopyNumberSegments> {
    assert!(!seg_input.is_empty());

    let init_probs = get_init_probs();
    let sample_count = seg_input.len();

    // Use settings fromt he first sample, leaning on the requirement that bin size is the same over all samples
    assert!(sample_count > 0);
    let depth_bin_size = seg_input[0].settings.depth_bin_size;

    //Run a parse for each chromosome:
    let chrom_count = genome_gc_levels.len();
    let mut cn_segments = (0..sample_count).map(|_| Vec::new()).collect::<Vec<_>>();

    for chrom_index in 0..chrom_count {
        let mut chrom_cn_segments = Vec::new();
        for sample_seg_input in seg_input.iter() {
            let sample_chrom_cn_segments = viterbi_copy_number_parse(
                genome_gc_levels,
                transition_info,
                sample_seg_input,
                &init_probs,
                depth_bin_size,
                chrom_index,
            );
            chrom_cn_segments.push(sample_chrom_cn_segments);
        }

        // Initial viterbi parse occurs independently for each sample over the chromosome.
        //
        // The next step adjusts the segmentation based on joint-sample information over the chromosome
        //
        refine_multi_sample_chrom_cn_segments(
            genome_gc_levels,
            transition_info,
            seg_input,
            chrom_index,
            depth_bin_size,
            &mut chrom_cn_segments,
        );

        for (sample_index, sample_chrom_cn_segments) in chrom_cn_segments.into_iter().enumerate() {
            cn_segments[sample_index].push(sample_chrom_cn_segments);
        }
    }

    let mut sample_results = Vec::new();
    for sample_cn_segments in cn_segments.into_iter() {
        let sample_result = SampleCopyNumberSegments {
            bin_size: depth_bin_size,
            cn_segments: sample_cn_segments,
        };
        sample_results.push(sample_result);
    }

    sample_results
}

/// Write out a bedgraph track for copy number segments (for use in IGV)
///
pub fn write_copy_number_segment_file(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    sample_cn_segments: &SampleCopyNumberSegments,
) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let filename = output_dir.join(COPYNUM_SEGMENT_BEDGRAPH_FILENAME);
    info!("Writing bedgraph copy number segment track to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create bedgraph copy number segment track file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    let chrom_count = chrom_list.data.len();
    let bin_size = sample_cn_segments.bin_size;

    for chrom_index in 0..chrom_count {
        let chrom_label = &chrom_list.data[chrom_index].label;
        let chrom_cn_segments = &sample_cn_segments.cn_segments[chrom_index];
        for s in chrom_cn_segments.iter() {
            if s.copy_number_info.state == CopyNumberState::Unknown {
                continue;
            }
            writeln!(
                f,
                "{}\t{}\t{}\t{}",
                chrom_label,
                s.begin_pos(bin_size),
                s.end_pos(bin_size),
                s.copy_number_info.as_u32_depth(),
            )
            .unwrap();
        }
    }
}

/// Serialize copy number segments to file so that sawfish can load this back in for joint-genotyping
pub fn serialize_copy_number_segments(
    discover_dir: &Utf8Path,
    sample_cn_segments: &SampleCopyNumberSegments,
) {
    let mut buf = Vec::new();
    sample_cn_segments
        .serialize(&mut rmp_serde::Serializer::new(&mut buf))
        .unwrap();

    let filename = discover_dir.join(COPYNUM_SEGMENT_MESSAGEPACK_FILENAME);

    info!("Writing copy number segments to binary file: '{filename}'");

    unwrap!(
        std::fs::write(&filename, buf.as_slice()),
        "Unable to open and write copy number segments to binary file: '{filename}'"
    );
}

/// Deserialize copy number segments from file
pub fn deserialize_copy_number_segments(discover_dir: &Utf8Path) -> SampleCopyNumberSegments {
    let filename = discover_dir.join(COPYNUM_SEGMENT_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read copy number segments from binary file: '{filename}'"
    );
    rmp_serde::from_slice(&buf).unwrap()
}

#[derive(
    Copy, Clone, Debug, Deserialize, Eq, PartialEq, Serialize, strum::EnumCount, strum::FromRepr,
)]
#[repr(usize)]
pub enum CopyNumberState {
    Zero,
    One,
    Two,
    Three,
    Four,
    /// All copy numbers greater than 4
    High,
    /// This state is used for excluded regions, if an excluded region is long enough, the weights
    /// are configured such that you'll eventually transition into it.
    Unknown,
}

/// Copy number information extended to include an optional depth estimate for the "High" state.
///
/// The idea here is that segmentation still uses a single High state to represent all copy numbers
/// higher than some threshold, but after segmentation, we can still go back and get a better copy
/// number estimate for reporting.
///
#[derive(Debug, Deserialize, Serialize)]
pub struct ExtendedCopyNumberState {
    /// Copy number state used for segmentation
    pub state: CopyNumberState,
    /// Supplemental copy number estimate for the "High" copy number state
    pub high_state_copy_number: Option<f64>,
}

impl ExtendedCopyNumberState {
    /// Convert copy number state information into an integer depth value, accounting for
    /// the optional high state depth estimate.
    ///
    pub fn as_u32_depth(&self) -> u32 {
        use CopyNumberState::*;
        assert!(self.state != Unknown);
        match (self.state, self.high_state_copy_number) {
            (High, Some(x)) => x.round() as u32,
            (s, _) => s as u32,
        }
    }
}

/// Copy Number state over a given bin interval
///
/// Interval defined by begin_bin,end_bin is zero-indexed, half-closed
///
#[derive(Debug, Deserialize, Serialize)]
pub struct CopyNumberSegment {
    pub begin_bin: usize,
    pub end_bin: usize,
    pub copy_number_info: ExtendedCopyNumberState,
}

impl CopyNumberSegment {
    pub fn begin_pos(&self, bin_size: u32) -> i64 {
        (self.begin_bin * bin_size as usize) as i64
    }

    pub fn end_pos(&self, bin_size: u32) -> i64 {
        (self.end_bin * bin_size as usize) as i64
    }

    pub fn to_range(&self, bin_size: u32) -> IntRange {
        IntRange::from_pair(self.begin_pos(bin_size), self.end_pos(bin_size))
    }
}

/// Genome copy-number segments for one sample
#[derive(Deserialize, Serialize)]
pub struct SampleCopyNumberSegments {
    pub bin_size: u32,
    /// Segments indexed on chromosome, then listed in order for each chromosome
    pub cn_segments: Vec<Vec<CopyNumberSegment>>,
}

/// Parse the depth track of one sample into copy number state segments
///
/// This method iterates the segmentation procedure, with the goal of refining the haploid depth estimate
/// in each iteration. The iterative approach has been found to be helpful in samples with at least
/// one very large CNV.
///
pub fn get_sample_copy_number_segments(
    settings: &cli::DiscoverSettings,
    chrom_list: &ChromList,
    sample_scan_result: &SampleAlignmentScanResult,
    gc_bias_data: &GCBiasCorrectionData,
) -> SampleCopyNumberSegments {
    info!("Segmenting copy number");

    let transition_info = TransitionProbInfo::new(settings.transition_prob, 1.0);

    // Iterate until the haploid coverage estimate converges or we reach max iteration count
    let max_iter = 8;
    let mut copy_number_segments = None;
    let mut last_coverage = 0.0;
    for iter_index in 0..max_iter {
        let haploid_coverage = get_haploid_genome_coverage(
            settings.coverage_est_regex.as_str(),
            chrom_list,
            &sample_scan_result.genome_depth_bins,
            &gc_bias_data.genome_gc_levels,
            &gc_bias_data.sample_gc_bias_data,
            copy_number_segments.as_ref(),
        );

        info!(
            "Haploid coverage estimate iteration {}. Uncorrected: {:.3} GC-Corrected: {:.3}",
            (iter_index + 1),
            haploid_coverage.depth,
            haploid_coverage.gc_corrected_depth
        );

        let sample_cns_input = SampleCopyNumberSegmentationInput {
            settings,
            gc_corrected_haploid_coverage: haploid_coverage.gc_corrected_depth,
            genome_depth_bins: &sample_scan_result.genome_depth_bins,
            sample_gc_bias_data: &gc_bias_data.sample_gc_bias_data,
            sample_transition_types: None,
            expected_copy_number_regions: None,
        };

        let mut cns = get_single_pass_sample_copy_number_segments(
            &gc_bias_data.genome_gc_levels,
            &transition_info,
            &[sample_cns_input],
        );
        copy_number_segments = Some(cns.pop().unwrap());
        if iter_index > 0 && (last_coverage - haploid_coverage.gc_corrected_depth).abs() < 0.001 {
            break;
        }
        last_coverage = haploid_coverage.gc_corrected_depth;
    }
    copy_number_segments.unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_emit() {
        // Value should match `log(dpois(5,10))` in R:
        let emit = get_depth_emission_lnprob(5.0, CopyNumberState::Two, &DepthBin::Depth(5.0), 1.0);
        approx::assert_ulps_eq!(emit, -3.2745662778118154, max_ulps = 4);

        let emit =
            get_depth_emission_lnprob(5.0, CopyNumberState::Unknown, &DepthBin::Excluded, 1.0);
        approx::assert_ulps_eq!(emit, 0.0, max_ulps = 4);
    }

    #[test]
    fn test_prune_copy_number_boundaries() {
        let mut cnbs = [3, 6, 8, 10, 13]
            .iter()
            .map(|&x| {
                (
                    x,
                    vec![CopyNumberBoundary {
                        sample_index: 0,
                        bin_index: x,
                        before_segment_index: 0,
                        before_state: CopyNumberState::Two,
                        after_state: CopyNumberState::Two,
                    }],
                )
            })
            .collect();

        prune_copy_number_boundaries(2, &mut cnbs);

        assert!(cnbs.keys().cloned().eq([6, 8, 10]));
    }
}
