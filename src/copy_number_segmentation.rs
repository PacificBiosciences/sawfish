use log::info;
use rust_vc_utils::ChromList;
use statrs::distribution::Discrete;

use crate::bam_scanner::SampleAlignmentScanResult;
use crate::cli;
use crate::depth_bins::{ChromDepthBins, DepthBin, GenomeDepthBins};
use crate::gc_correction::{GCBiasCorrectionData, SampleGCBiasCorrectionData};
use crate::CNState::Unknown;

/// Get whole genome haploid coverage estimates
///
/// Averages the coverage over non-excluded non-zero depth bins in chromosomes matching
/// coverage_est_regex.
///
/// GC-correction factor is applied such that haploid depth estimate should
/// reflect the depth of the least biased GC window.
///
/// Returns a 2-tuple of uncorrected and GC-corrected haploid coverage
fn get_haploid_genome_coverage(
    coverage_est_regex: &str,
    chrom_list: &ChromList,
    genome_depth_bins: &GenomeDepthBins,
    gc_bias_data: &GCBiasCorrectionData,
    input_cn_states: &Option<SampleCNStateBins>,
) -> (f64, f64) {
    use regex::Regex;
    let chrom_include_regex = Regex::new(coverage_est_regex).unwrap();

    // This application requires the reciprocal of the depth reduction factor used for the
    // emission prob correction, so go ahead and setup the reciprocal value array here:
    let gc_depth_correction = gc_bias_data
        .sample_gc_bias_data
        .gc_depth_reduction
        .iter()
        .map(|x| 1.0 / x)
        .collect::<Vec<_>>();

    let mut total = 0.0;
    let mut gc_corrected_total = 0.0;

    // This flag is used to improve error message if user specifies an invalid chrom regex:
    let mut chrom_match = false;

    // The number of observed chromosome bins, meaning for a bin covering a 2-copy region in the
    // sample we would add 2.
    let mut count = 0.0;

    for (chrom_index, chrom_entry) in chrom_list.data.iter().enumerate() {
        if !chrom_include_regex.is_match(chrom_entry.label.as_str()) {
            continue;
        }
        chrom_match = true;

        let chrom_depth_bins = &genome_depth_bins.chroms[chrom_index];
        let chrom_gc_levels = &gc_bias_data.genome_gc_levels[chrom_index];
        let chrom_cn_states = input_cn_states
            .as_ref()
            .map(|sample_cn_segments| &sample_cn_segments.data[chrom_index]);
        let bin_count = chrom_depth_bins.len();
        for bin_index in 0..bin_count {
            let depth = &chrom_depth_bins[bin_index];
            let gc_level = chrom_gc_levels[bin_index];
            if let DepthBin::Depth(x) = depth {
                // Get expected copy number (defaults to 2 unless copy number is available)
                let bin_copy_number = match chrom_cn_states {
                    Some(chrom_cn_states) => {
                        if chrom_cn_states[bin_index] == CNState::Unknown as u8 {
                            continue;
                        }
                        chrom_cn_states[bin_index] as f64
                    }
                    None => 2.0,
                };
                if *x > 0.0 {
                    total += x;
                    gc_corrected_total += x * gc_depth_correction[gc_level];
                    count += bin_copy_number;
                }
            };
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

    (total / count, gc_corrected_total / count)
}

#[allow(clippy::needless_range_loop)]
/// Compute and return the HMM state transition matrix as log prob values
///
/// Matrix returned with lookup format `matrix[fromState][toState]`
///
fn get_chrom_transition_probs(go_prob: f64) -> Vec<Vec<f64>> {
    // Start out with a really simple transition model.
    //
    // Note that transition probs, in addition to reflecting the biological phenomena, have to
    // account for the very high bin-to-bin dependency in the depth value from long reads.
    //

    let state_count = CNState::Unknown as usize + 1;

    // Note this is a critical parameter, and it is dependent on bin_size (and read size)
    //
    // At bin_size = 5000, go_prob of 1e-20 seemed to work well based on ad hoc optimization
    // At bin_size = 1000, go_prob of 1e-100 seemed to work well
    //
    // TODO: Partition this value into one core component that is independent of bin-to-bin
    // correlation, and another component which accounts for the correlation. We might be able
    // to model this dependency so that a change in bin-size or read-size is automatically
    // adjusted for.
    //
    //let go_prob = 1e-50;
    let stay_prob = 1.0 - (go_prob * ((state_count - 1) as f64));
    let go_ln_prob = go_prob.ln();
    let stay_ln_prob = stay_prob.ln();

    //clippy doesn't like this state indexing, but changing it would be counter-intuitive IMO
    let mut tm = vec![vec![0.0; state_count]; state_count];
    for from_index in 0..state_count {
        for to_index in 0..state_count {
            tm[from_index][to_index] = {
                if from_index == to_index {
                    stay_ln_prob
                } else {
                    go_ln_prob
                }
            };
        }
    }
    tm
}

fn get_expected_depth(
    haploid_coverage: f64,
    copy_num_state: CNState,
    depth: f64,
    gc_depth_reduction: f64,
) -> f64 {
    let default_expected_depth =
        || haploid_coverage * gc_depth_reduction * copy_num_state as usize as f64;

    // Special handlers to modify expectedDepth for certain states:
    match copy_num_state {
        CNState::High => {
            // This state accounts for anything at copy number HIGH or higher.
            // A hacky way to handle this state is to set expected depth to match observed depth,
            // so long as this does not lower expected depth below that for copy number HIGH.
            default_expected_depth().max(depth.round())
        }
        CNState::Zero => {
            // Don't use zero for expected depth, instead give the model a way to explain occasional
            // garbage alignments inside zero copy regions.
            let min_copy_num = 0.01;
            haploid_coverage * min_copy_num
        }
        _ => default_expected_depth(),
    }
}

/// Return the log emission prob of the depth bin observation given the copy number state
///
fn get_emission(
    haploid_coverage: f64,
    copy_num_state: CNState,
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
                Unknown => {
                    // Disallow Unknown state unless the depth bin is excluded
                    return f64::NEG_INFINITY;
                }
                _ => *depth,
            }
        }
        DepthBin::Excluded => {
            return match copy_num_state {
                Unknown => excluded_region_unknown_state_prob,
                _ => excluded_region_known_state_prob,
            };
        }
    };

    assert!(haploid_coverage >= 0.0);
    assert!(depth >= 0.0);

    // Set a max copy number to retain numerical stability in depth spikes:
    let max_copy_number = 20;
    let depth = depth.min(haploid_coverage * max_copy_number as f64);
    let expected_depth =
        get_expected_depth(haploid_coverage, copy_num_state, depth, gc_depth_reduction);

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

/// A Viterbi parse for copy number states
///
/// All prob values are tracked in log space. This parse is extensively customized for the specific
/// CN inference case.
///
/// # Arguments
///
/// * `init` - Log prob prior over state space
///
/// * `transition` - Log prob of state transition
///
fn viterbi_cn_parse(
    hapliod_coverage: f64,
    init: &[f64],
    transition: &[Vec<f64>],
    observations: &ChromDepthBins,
    chrom_gc_levels: &[usize],
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
) -> (Vec<u8>, Vec<f32>) {
    assert!(hapliod_coverage >= 0.0);

    let state_count = CNState::Unknown as usize + 1;
    assert_eq!(init.len(), state_count);

    let obs_count = observations.len();
    if obs_count == 0 {
        return (Vec::new(), Vec::new());
    }

    // Instead of having a full SxO DP matrix, just ping-pong on two rows
    let mut max_pr_row1 = vec![0.0; state_count];
    let mut max_pr_row2 = vec![0.0; state_count];

    let mut back_pointer = vec![vec![0u8; state_count]; obs_count];
    // note: this stores all probabilities, removing the ping-pong benefits; not sure we care though, still relatively small
    // we could consolidate/remove the max_pr_row# variables and just use this one
    let mut max_probs: Vec<Vec<f64>> = vec![vec![0.0; state_count]; obs_count];

    for state_index in 0..state_count {
        max_pr_row1[state_index] = init[state_index]
            + get_emission(
                hapliod_coverage,
                CNState::from_repr(state_index).unwrap(),
                &observations[0],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[0]],
            );
    }

    let this_row = &mut max_pr_row1;
    let last_row = &mut max_pr_row2;
    for obs_index in 0..obs_count {
        std::mem::swap(this_row, last_row);
        for (to_state_index, row_value) in this_row.iter_mut().enumerate() {
            let emit = get_emission(
                hapliod_coverage,
                CNState::from_repr(to_state_index).unwrap(),
                &observations[obs_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[obs_index]],
            );

            let mut max_index = 0;
            let mut max_log_prob = 0.0;
            for from_state_index in 0..state_count {
                let log_prob = last_row[from_state_index]
                    + transition[from_state_index][to_state_index]
                    + emit;

                if (from_state_index == 0) || (log_prob > max_log_prob) {
                    max_index = from_state_index;
                    max_log_prob = log_prob;
                }
            }

            *row_value = max_log_prob;
            back_pointer[obs_index][to_state_index] = max_index as u8;
            max_probs[obs_index][to_state_index] = max_log_prob;
        }
    }

    // Backtrace
    let mut max_state = 0;
    for (state, val) in this_row.iter().enumerate() {
        if (state == 0) || (*val > this_row[max_state]) {
            max_state = state;
        }
    }

    let mut max_path: Vec<u8> = vec![0; obs_count];
    let mut relative_prob: Vec<f32> = vec![0.0; obs_count];
    for obs_index in (0..obs_count).rev() {
        //store minimum relative probabilities deltas
        max_path[obs_index] = max_state as u8;
        let ms_prob = max_probs[obs_index][max_state];
        let mut min_prob_delta: f32 = f32::MAX;
        for (i, mrp) in max_probs[obs_index].iter().enumerate() {
            let rel_mrp = ms_prob - mrp;
            if i != max_state && rel_mrp < min_prob_delta as f64 {
                min_prob_delta = rel_mrp as f32;
            }
        }
        //this *can* be negative because at a given locus, another point may be more likely but not part of viterbi path
        relative_prob[obs_index] = min_prob_delta;

        //get next in backtrace
        max_state = back_pointer[obs_index][max_state] as usize;
    }
    (max_path, relative_prob)
}

/// Translate the most likely bin CN state path into CN segments
///
fn get_copy_segments(bin_size: u32, max_path: &[u8], quals: &[f32]) -> Vec<CNSegment> {
    assert!(max_path.len() == quals.len());
    let mut chrom_segments: Vec<CNSegment> = Vec::new();

    let bin_count = max_path.len();

    let mut begin_bin: usize = 0;
    let mut last_cn = CNState::Zero;

    let chrom_segments_ref = &mut chrom_segments;

    let mut add_segment = |last_cn, end_bin| {
        let begin = (begin_bin * bin_size as usize) as i64;
        let end = (end_bin * bin_size as usize) as i64;

        //mean quality where all negatives are converted to 0.0 first
        //TODO: consider other options for QUAL scoring
        let qual = quals[begin_bin..end_bin]
            .iter()
            .map(|&x| if x > 0.0 { x } else { 0.0 })
            .sum::<f32>()
            / (end_bin - begin_bin) as f32;
        chrom_segments_ref.push(CNSegment {
            begin,
            end,
            state: last_cn,
            qual,
        });
        begin_bin = end_bin;
    };

    for (bin_index, &bin_value) in max_path.iter().enumerate() {
        let this_cn = CNState::from_repr(bin_value as usize).unwrap();
        if (bin_index > 0) && (this_cn != last_cn) {
            add_segment(last_cn, bin_index);
        }
        last_cn = this_cn;
    }
    add_segment(last_cn, bin_count);

    chrom_segments
}

/// Genome copy-number state bins for one sample
///
/// This is the direct output of the segmentation method before bins have been consolidated.
struct SampleCNStateBins {
    /// Bins indexed on chromosome, then arranged to match depth bins on each chromosome
    data: Vec<Vec<u8>>,
    quals: Vec<Vec<f32>>,
}

/// Classify depth track bins into copy number state bins
///
/// This method can be iterated by specifying CN state results from a previous round as an optional
/// input parameter.
///
fn get_sample_copy_number_states(
    settings: &cli::DiscoverSettings,
    gc_corrected_haploid_coverage: f64,
    genome_depth_bins: &GenomeDepthBins,
    gc_bias_data: &GCBiasCorrectionData,
) -> (SampleCNStateBins, f64) {
    let transition_probs = get_chrom_transition_probs(settings.transition_prob);

    // Setup initial state log probs
    // TODO: change for haploid chromosomes:
    let diploid_init_probs = [0.001, 0.001, 0.994, 0.001, 0.001, 0.001, 0.001]
        .iter()
        .map(|x: &f64| x.ln())
        .collect::<Vec<_>>();

    //Run a parse for each chromosome:
    let mut sample_cn_states = SampleCNStateBins {
        data: Vec::new(),
        quals: Vec::new(),
    };

    for (chrom_index, chrom_depth_bins) in genome_depth_bins.chroms.iter().enumerate() {
        let chrom_gc_levels = &gc_bias_data.genome_gc_levels[chrom_index];
        let (max_path, relative_probs) = viterbi_cn_parse(
            gc_corrected_haploid_coverage,
            &diploid_init_probs,
            &transition_probs,
            chrom_depth_bins,
            chrom_gc_levels,
            &gc_bias_data.sample_gc_bias_data,
        );

        let bin_count = chrom_depth_bins.len();
        assert_eq!(bin_count, max_path.len());

        sample_cn_states.data.push(max_path);
        sample_cn_states.quals.push(relative_probs);
    }
    (sample_cn_states, gc_corrected_haploid_coverage)
}

#[derive(Copy, Clone, PartialEq, Eq, strum::FromRepr)]
#[repr(usize)]
pub enum CNState {
    Zero,
    One,
    Two,
    Three,
    Four,
    /// All copy numbers greater than 4
    High,
    /// This state is used for excluded regions, if an excluded region is long enough, the weights
    /// are configured such that you'll eventually transfer into it.
    /// transition into it.
    Unknown,
}

/// Copy Number state over a given interval
///
/// Interval defined by begin,end is zero-indexed, half-closed
pub struct CNSegment {
    pub begin: i64,
    pub end: i64,
    pub state: CNState,
    pub qual: f32,
}

/// Genome copy-number segments for one sample
pub struct SampleCopyNumberSegments {
    pub chrom_list: ChromList,
    pub sample_name: String,
    /// Segments indexed on chromosome, then listed in order for each chromosome
    pub data: Vec<Vec<CNSegment>>,
}

/// Segment the depth track of one sample
///
/// This method can iterate the segmentation multiple times, with the goal of refining the true
/// haploid depth estimate in each round. This should be helpful in samples with at least one very
/// large CNV.
///
pub fn segment_sample_copy_number(
    settings: &cli::DiscoverSettings,
    chrom_list: &ChromList,
    sample_scan_result: &SampleAlignmentScanResult,
    gc_bias_data: &GCBiasCorrectionData,
) -> SampleCopyNumberSegments {
    // Iterate until the haploid coverage estimate converges (or we reach max iteration count)
    let max_iter = 8;
    let mut cn_states = None;
    let mut last_coverage = 0.0;
    for iter_index in 0..max_iter {
        let (haploid_coverage, gc_corrected_haploid_coverage) = get_haploid_genome_coverage(
            settings.coverage_est_regex.as_str(),
            chrom_list,
            &sample_scan_result.genome_depth_bins,
            gc_bias_data,
            &cn_states,
        );

        info!(
            "Haploid coverage estimates for sample '{}', iteration {}. Uncorrected: {:.3} GC-Corrected: {:.3}",
            sample_scan_result.sample_name, (iter_index+1), haploid_coverage, gc_corrected_haploid_coverage
        );

        let (cn_states_bins, coverage) = get_sample_copy_number_states(
            settings,
            gc_corrected_haploid_coverage,
            &sample_scan_result.genome_depth_bins,
            gc_bias_data,
        );
        cn_states = Some(cn_states_bins);
        if iter_index > 0 && (last_coverage - coverage).abs() < 0.001 {
            break;
        }
        last_coverage = coverage;
    }
    let cn_states = cn_states.unwrap();

    let mut sample_cn_segments = SampleCopyNumberSegments {
        chrom_list: chrom_list.clone(),
        sample_name: sample_scan_result.sample_name.clone(),
        data: Vec::new(),
    };
    for (i, chrom_cn_states) in cn_states.data.iter().enumerate() {
        let quals = &cn_states.quals[i];
        sample_cn_segments.data.push(get_copy_segments(
            settings.depth_bin_size,
            chrom_cn_states,
            quals,
        ));
    }
    sample_cn_segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_emit() {
        // Value should match `log(dpois(5,10))` in R:
        let emit = get_emission(5.0, CNState::Two, &DepthBin::Depth(5.0), 1.0);
        approx::assert_ulps_eq!(emit, -3.2745662778118154, max_ulps = 4);

        let emit = get_emission(5.0, CNState::Unknown, &DepthBin::Excluded, 1.0);
        approx::assert_ulps_eq!(emit, 0.0, max_ulps = 4);
    }
}
