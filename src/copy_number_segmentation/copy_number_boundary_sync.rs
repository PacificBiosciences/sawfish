use std::collections::BTreeMap;

use itertools::Itertools;

use super::{
    CopyNumberSegment, CopyNumberState, SampleCopyNumberSegmentationInput, TransitionProbInfo,
    get_depth_emission_lnprob, get_transition_lnprob,
};
use crate::gc_correction::GenomeGCLevels;

/// Information to track a single copy-number shift boundary in one sample
#[derive(Clone)]
struct CopyNumberBoundary {
    /// This is the index among samples with CNV enabled, which is different than regular sample_index
    cnv_sample_index: usize,

    /// First bin affected by (ie. after) the copy-number change (zero-indexed)
    bin_index: usize,

    /// Index of the copy number segment immediately before the bin change
    before_segment_index: usize,

    /// Copy number before the boundary
    before_state: CopyNumberState,

    /// Copy number after the boundary
    after_state: CopyNumberState,
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

    assert!(!from_cnb_set.is_empty());

    let (start_bin, end_bin, after_state_is_baseline) = if from_bin < to_bin {
        (from_bin, to_bin, true)
    } else {
        (to_bin, from_bin, false)
    };

    let mut score_delta = 0.0;

    let chrom_gc_levels = &genome_gc_levels[chrom_index];

    // First find the difference in emit scores between current and revised copy number
    for copy_number_boundary in from_cnb_set {
        let sample_seg_input = &seg_input[copy_number_boundary.cnv_sample_index];

        let chrom_depth_bins = &sample_seg_input.genome_depth_bins.depth_bins[chrom_index];
        let sample_gc_bias_data = sample_seg_input.sample_gc_bias_data;
        let gc_corrected_haploid_coverage = sample_seg_input.gc_corrected_haploid_coverage;

        let mut sample_emit_score_delta = 0.0;
        for bin_index in start_bin..end_bin {
            let after_state_emit_lnprob = get_depth_emission_lnprob(
                gc_corrected_haploid_coverage,
                copy_number_boundary.after_state,
                &chrom_depth_bins[bin_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]],
            );

            let before_state_emit_lnprob = get_depth_emission_lnprob(
                gc_corrected_haploid_coverage,
                copy_number_boundary.before_state,
                &chrom_depth_bins[bin_index],
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]],
            );

            if after_state_emit_lnprob == f64::NEG_INFINITY
                || before_state_emit_lnprob == f64::NEG_INFINITY
            {
                return f64::NEG_INFINITY;
            }

            let (baseline_emit_lnprob, shifted_emit_lnprob) = if after_state_is_baseline {
                (after_state_emit_lnprob, before_state_emit_lnprob)
            } else {
                (before_state_emit_lnprob, after_state_emit_lnprob)
            };

            sample_emit_score_delta += shifted_emit_lnprob - baseline_emit_lnprob;
        }

        // At this point we have the delta for emit values from all samples, so the bin_dependency correction factor
        // applies to the entire score:
        let bin_dependency_correction_factor = sample_seg_input
            .discover_settings
            .bin_dependency_correction_factor;
        sample_emit_score_delta *= bin_dependency_correction_factor;

        score_delta += sample_emit_score_delta;
    }

    // Next, find the difference in transition scores, accounting for the shared transition site bonus
    for copy_number_boundary in from_cnb_set {
        let sample_seg_input = &seg_input[copy_number_boundary.cnv_sample_index];

        let sample_transition_types = sample_seg_input.sample_transition_types;

        let baseline_transition_score = get_transition_lnprob(
            transition_info,
            sample_transition_types,
            chrom_index,
            from_bin,
            copy_number_boundary.before_state as usize,
            copy_number_boundary.after_state as usize,
        );
        let shifted_transition_score = get_transition_lnprob(
            transition_info,
            sample_transition_types,
            chrom_index,
            to_bin,
            copy_number_boundary.before_state as usize,
            copy_number_boundary.after_state as usize,
        );

        let sample_transition_score_delta =
            (shifted_transition_score * shifted_cnb_transition_factor) - baseline_transition_score;
        score_delta += sample_transition_score_delta;
    }

    // score_delta can technically be -inf, but nan and inf should never occur:
    assert!(
        !score_delta.is_nan() && score_delta != f64::INFINITY,
        "Invalid copy number boundary shift assessment score: {score_delta} at chrom_index: {chrom_index} from_bin: {from_bin} to_bin: {to_bin}",
    );

    // return the score change per shifted sample
    score_delta / from_cnb_set.len() as f64
}

// Get full set of shift scores that are possible between cnb_set members
#[derive(Debug, PartialEq, PartialOrd)]
struct ShiftScore {
    score: f64,
    from_cluster_member_index: usize,
    to_cluster_member_index: usize,
}

#[allow(clippy::too_many_arguments)]
fn update_cnb_shift_scores(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    from_bin: usize,
    from_cnb_set: &[CopyNumberBoundary],
    to_bin: usize,
    from_cluster_member_index: usize,
    to_cluster_member_index: usize,
    shift_scores: &mut Vec<ShiftScore>,
) {
    let score = cnb_score_shift_per_sample(
        genome_gc_levels,
        transition_info,
        seg_input,
        chrom_index,
        from_bin,
        from_cnb_set,
        to_bin,
    );

    shift_scores.push(ShiftScore {
        score,
        from_cluster_member_index,
        to_cluster_member_index,
    });
}

/// Given a cluster of CNBs within the max_bin_shift threshold, test possible consolidating shifts and report out the
/// recommended changes as a set of from->to bin index transitions.
///
/// Returns a series of (from_index,to_index) tuples. Each index in the tuple is a CopyNumberBoundary.bin_index,
/// describing the bin_index boundaries that should be shifted and to where.
///
#[allow(clippy::needless_range_loop)]
fn get_cnb_shifts_for_cnb_cluster(
    genome_gc_levels: &GenomeGCLevels,
    transition_info: &TransitionProbInfo,
    seg_input: &[SampleCopyNumberSegmentationInput],
    chrom_index: usize,
    cnb_cluster: &[(usize, &Vec<CopyNumberBoundary>)],
    debug: bool,
) -> Vec<(usize, usize)> {
    // Singleton clusters will never have shifts and shouldn't be processed through this method.
    assert!(cnb_cluster.len() >= 2, "Unexpected singleton cluster");

    if debug {
        eprintln!(
            "get_cnb_shifts_for_cnb_cluster: start cnb_cluster member count: {}",
            cnb_cluster.len()
        );
        for (key, val) in cnb_cluster.iter() {
            eprintln!(
                "cluster member key/bin_index: {key} cnb_count: {}",
                val.len()
            );
        }
    }

    let mut cnb_shift_scores = Vec::new();

    let cluster_member_count = cnb_cluster.len();
    for cluster_member_index1 in 0..cluster_member_count {
        let (key1, cnb1) = cnb_cluster[cluster_member_index1];

        for cluster_member_index2 in (cluster_member_index1 + 1)..cluster_member_count {
            let (key2, cnb2) = cnb_cluster[cluster_member_index2];

            // Score key1 cnbs shifting to key2 location
            update_cnb_shift_scores(
                genome_gc_levels,
                transition_info,
                seg_input,
                chrom_index,
                key1,
                cnb1,
                key2,
                cluster_member_index1,
                cluster_member_index2,
                &mut cnb_shift_scores,
            );

            // Score key2 cnbs shifting to key1 location
            update_cnb_shift_scores(
                genome_gc_levels,
                transition_info,
                seg_input,
                chrom_index,
                key2,
                cnb2,
                key1,
                cluster_member_index2,
                cluster_member_index1,
                &mut cnb_shift_scores,
            );
        }
    }

    // Convert shift_scores to cnb_shifts:
    cnb_shift_scores.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let mut cnb_shifts = Vec::new();
    let mut shift_from_cluster_members = vec![false; cluster_member_count];
    let mut shift_to_cluster_members = vec![false; cluster_member_count];

    for cnb_shift_score in cnb_shift_scores {
        if cnb_shift_score.score <= 0.0 {
            break;
        }

        if shift_from_cluster_members[cnb_shift_score.from_cluster_member_index] {
            // We've shifted out of this location already so can't shift out again
            continue;
        }

        if shift_from_cluster_members[cnb_shift_score.to_cluster_member_index] {
            // We've shifted out of this location already so don't shift into it now
            continue;
        }

        if shift_to_cluster_members[cnb_shift_score.from_cluster_member_index] {
            // We've shifted into this location already so don't shift out of it now
            continue;
        }

        let from_key = cnb_cluster[cnb_shift_score.from_cluster_member_index].0;
        let to_key = cnb_cluster[cnb_shift_score.to_cluster_member_index].0;
        cnb_shifts.push((from_key, to_key));
        shift_from_cluster_members[cnb_shift_score.from_cluster_member_index] = true;
        shift_to_cluster_members[cnb_shift_score.to_cluster_member_index] = true;
    }

    cnb_shifts
}

/// Cluster copy number boundaries (CNBs), and then refine each cluster
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
    let mut all_cnb_cluster_shifts = Vec::new();

    // 1. Cluster CNBs to include all members with max_bin_shift (transitive)
    // 2. Obtain the recommended CNB shifts within each non-singleton cluster
    // 3. Return all recommended CNB shifts
    let mut cnb_cluster = Vec::new();
    let keys = cnbs.keys().cloned().collect::<Vec<_>>();
    for (key, next_key) in keys.into_iter().chain([0]).tuple_windows() {
        if next_key == 0 || (next_key - key) > max_bin_shift {
            // Filter out singleton clusters:
            if !cnb_cluster.is_empty() {
                cnb_cluster.push((key, cnbs.get(&key).unwrap()));
                let mut cnb_cluster_shifts = get_cnb_shifts_for_cnb_cluster(
                    genome_gc_levels,
                    transition_info,
                    seg_input,
                    chrom_index,
                    &cnb_cluster,
                    debug,
                );

                all_cnb_cluster_shifts.append(&mut cnb_cluster_shifts);
                cnb_cluster.clear();
            }
        } else {
            cnb_cluster.push((key, cnbs.get(&key).unwrap()));
        }
    }

    all_cnb_cluster_shifts
}

fn shift_copy_number_boundaries(
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
            let sample_chrom_cn_segments = &mut chrom_cn_segments[from_bin_cnb.cnv_sample_index];
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

/// Map key should exactly match CopyNumberBoundary::bin_index
type CopyNumberBoundarySet = BTreeMap<usize, Vec<CopyNumberBoundary>>;

/// All predicted copy number boundaries for one chromosome over multiple samples
struct CopyNumberBoundaryInfo {
    /// 'losses' mean the boundary goes from copy number N to N-1, while transitioning from bin_index to bin_index+1
    losses: CopyNumberBoundarySet,

    /// 'gains' mean the boundary goes from copy number N to N+1, while transitioning from bin_index to bin_index+1
    gains: CopyNumberBoundarySet,
}

/// Enumerate all predicted copy number boundaries for one chromosome over multiple samples
///
fn get_copy_number_boundary_info(
    chrom_cn_segments: &[Vec<CopyNumberSegment>],
) -> CopyNumberBoundaryInfo {
    let mut losses = BTreeMap::new();
    let mut gains = BTreeMap::new();

    for (cnv_sample_index, sample_chrom_cn_segments) in chrom_cn_segments.iter().enumerate() {
        let segment_count = sample_chrom_cn_segments.len();
        for segment_index in 1..segment_count {
            let before_segment_index = segment_index - 1;
            let after_segment_index = segment_index;
            let before_seg = &sample_chrom_cn_segments[before_segment_index];
            let after_seg = &sample_chrom_cn_segments[after_segment_index];

            let before_state = before_seg.copy_number_info.state;
            let after_state = after_seg.copy_number_info.state;

            if before_state == after_state {
                panic!(
                    "Adjacent copy number segments have the same copy number state {before_state:?}"
                );
            }

            if before_state == CopyNumberState::Unknown || after_state == CopyNumberState::Unknown {
                continue;
            }

            let cnb = CopyNumberBoundary {
                cnv_sample_index,
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

            cnbs.entry(cnb.bin_index).or_insert_with(Vec::new).push(cnb);
        }
    }

    CopyNumberBoundaryInfo { losses, gains }
}

/// Make small copy-number boundary shifts on one chromosome, to improve boundary synchronization across all samples
///
/// Over all samples, find instances where the copy number is shifting in the same direction over multiple samples
/// within a small window. Given these instances, test if copy boundaries can be moved by giving a small syncronizaition
/// bonus to either sample (of a pair) when one boundary is shifted to match the other. The small bonus encourages
/// shifting boundary differences caused by sampling noise only, with the goal that true differences can be preserved.
///
/// # Arguments
/// * `chrom_cn_segments` - An array over all samples from seg_input, where for each sample the entry contains all
///   predicted copy-number segments for the target chromosome
///
pub fn sync_chrom_copy_number_boundaries(
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
            "rms: refine_multi_sample_chrom_cn_segments starting on chrom_index: {chrom_index}"
        );
    }

    let sample_count = seg_input.len();
    assert_eq!(sample_count, chrom_cn_segments.len());

    // This refine operation is only relevant to multi-sample cnv analysis, so skip if there are fewer than 2 cnv-enabled samples:
    if sample_count < 2 {
        return;
    }

    // Hard code parameters:
    //

    // Max distance that a segment boundary can be shifted in a single sync operation:
    let max_base_shift = 5000;

    // Translate max shift from bases into depth bins
    let max_bin_shift = (max_base_shift / depth_bin_size) as usize;

    assert!(max_bin_shift >= 1, "Unexpected depth bin size");

    let cnb_info = get_copy_number_boundary_info(chrom_cn_segments);

    if debug {
        eprintln!("rms: cnb_loss_bin_count: {}", cnb_info.losses.len());
        eprintln!("rms: cnb_gain_bin_count: {}", cnb_info.gains.len());
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
            eprintln!("rms: cnb_loss_shift: {x:?}");
        }
        eprintln!("rms: cnb_gain_shifts: {}", cnb_gain_shifts.len());
        for x in cnb_gain_shifts.iter() {
            eprintln!("rms: cnb_gain_shift: {x:?}");
        }
    }

    // Adjust the original segments to reflect the cnb_shifts:
    shift_copy_number_boundaries(&cnb_info.losses, &cnb_loss_shifts, chrom_cn_segments, debug);
    shift_copy_number_boundaries(&cnb_info.gains, &cnb_gain_shifts, chrom_cn_segments, debug);
}
