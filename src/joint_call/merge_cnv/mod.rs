mod transition_types;

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use camino::Utf8Path;
use log::info;

use self::transition_types::get_copynum_transition_types_from_breakends;
use super::write_copy_number_info::write_copy_number_info;
use crate::breakpoint::BreakendDirection;
use crate::cli::JointCallSettings;
use crate::copy_number_segmentation::{
    CopyNumberSegment, CopyNumberState, ExtendedCopyNumberState, HaploidCoverage,
    SampleCopyNumberSegmentationInput, SampleCopyNumberSegments, SampleTransitionTypes,
    get_depth_emission_lnprob, get_haploid_genome_coverage,
    get_single_pass_sample_copy_number_segments, write_copy_number_segment_file,
};
use crate::depth_bins::DepthBin;
use crate::expected_ploidy::get_majority_expected_copy_number_for_region;
use crate::gc_correction::SampleGCBiasCorrectionData;
use crate::genome_segment::GenomeSegment;
use crate::int_range::{IntRange, get_recip_overlap};
use crate::joint_call::{SampleJointCallData, SharedJointCallData};
use crate::prob_utils::{ln_error_prob_to_qphred, normalize_ln_distro};
use crate::refine_sv::{Genotype, RefinedSV};
use crate::refined_cnv::{CNVSampleScoreInfo, CNVScoreInfo, RefinedCNV, SharedSampleCopyNumInfo};
use crate::sv_group::SVGroup;
use crate::sv_id::SVUniqueIdData;
use crate::utils::remove_sorted_indices;

/// Get the likelihood for each copy number level in a single sample
///
/// Return a vector of log likelihoods for each copy number state in the model
///
fn get_sample_cn_state_lnlhoods(
    cn_segment: &CopyNumberSegment,
    chrom_depth_bins: &[DepthBin],
    chrom_gc_levels: &[usize],
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
    gc_corrected_haploid_coverage: f64,
    bin_dependency_correction_factor: f64,
) -> Vec<f64> {
    let mut state_lnlhoods = Vec::new();

    let state_count = CopyNumberState::Unknown as usize;
    for state_index in 0..state_count {
        let copy_num_state = CopyNumberState::from_repr(state_index).unwrap();

        let mut state_lnlhood = 0.0;
        for bin_index in cn_segment.begin_bin..cn_segment.end_bin {
            let depth_bin_value = &chrom_depth_bins[bin_index];
            let gc_depth_reduction =
                sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]];
            let lnp = get_depth_emission_lnprob(
                gc_corrected_haploid_coverage,
                copy_num_state,
                depth_bin_value,
                gc_depth_reduction,
            );
            state_lnlhood += lnp * bin_dependency_correction_factor;
        }
        state_lnlhoods.push(state_lnlhood);
    }
    state_lnlhoods
}

/// Generate all copy number quality values for one CNV in one sample
///
/// All quality scores are approximated from pre-computed copy number state likelihoods, themselves mostly computed from
/// the emit probs for each cn state over the CNV segment.
///
fn get_cnv_sample_qualities_from_cn_state_lnlhoods(
    max_qscore: u32,
    copy_number_info: &ExtendedCopyNumberState,
    state_lnlhoods: &[f64],
    sample_score: &mut CNVSampleScoreInfo,
) {
    let copy_number = copy_number_info.as_u32_depth();

    // Compute CNQ:
    let copy_number_qscore = {
        // First generate the equiv of PL values for copy number (this doesn't exist in the VCF spec)
        //
        let mut cn_lhood_qscore = Vec::new();
        let cn_state_count = state_lnlhoods.len();
        for cn_index in 0..cn_state_count {
            let cnq = std::cmp::min(
                ln_error_prob_to_qphred(
                    state_lnlhoods[cn_index] - state_lnlhoods[copy_number_info.state as usize],
                ),
                max_qscore as i32,
            );

            // Max to zero since we aren't necessarily guaranteed that cn_state is the best lhood.
            // The max assumption has been broken at times by some non-probabilistic segmentation
            // schemes.
            //
            let cnq = std::cmp::max(cnq, 0);
            cn_lhood_qscore.push(cnq);
        }

        cn_lhood_qscore
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != copy_number_info.state as usize)
            .min_by(|(_, a), (_, b)| a.cmp(b))
            .map(|(_, &a)| a)
            .unwrap()
    };

    // Compute cnv_pprob, which is used to return the value we need for qual:
    //
    // TODO: Using naive prior for now but consider adding a more meaningful prior distro here
    //
    let mut cn_pprob = state_lnlhoods.to_vec();
    let _max_cn_pprob_index = normalize_ln_distro(&mut cn_pprob).unwrap();

    let expected_copy_number_prob = cn_pprob[sample_score.expected_copy_number as usize];

    sample_score.shared.copy_info = Some(SharedSampleCopyNumInfo {
        copy_number,
        copy_number_qscore,
        expected_copy_number_prob,
    });
}

/// Convert copy number segment to a refined CNV
///
/// A refined CNV is only produced if the segment is a copy number variant in at least one
/// sample
///
/// # Arguments
/// * `source_sample_index` - The sample this copy number segment was segmented from
///
#[allow(clippy::too_many_arguments)]
fn convert_cn_segment_to_refined_cnv(
    all_sample_data: &[SampleJointCallData],
    max_qscore: u32,
    source_sample_index: usize,
    chrom_index: usize,
    cn_segment: &CopyNumberSegment,
    bin_size: u32,
    chrom_gc_levels: &[usize],
    sample_seg_input: &SampleCopyNumberSegmentationInput,
) -> Option<RefinedCNV> {
    let segment = GenomeSegment {
        chrom_index,
        range: cn_segment.to_range(bin_size),
    };

    let mut score = CNVScoreInfo::default();

    for sample_data in all_sample_data.iter() {
        let expected_copy_number = get_majority_expected_copy_number_for_region(
            sample_data.expected_copy_number_regions.as_ref(),
            &segment,
        ) as u32;
        let sample_score = CNVSampleScoreInfo {
            expected_copy_number,
            ..Default::default()
        };
        score.samples.push(sample_score);
    }

    // Check whether this is a variant in the source sample. Return None if not.
    let source_sample_expected_copy_number =
        score.samples[source_sample_index].expected_copy_number;
    if cn_segment.copy_number_info.state as u32 == source_sample_expected_copy_number
        || cn_segment.copy_number_info.state == CopyNumberState::Unknown
    {
        return None;
    }

    // Add sample quality scores into the each CNV's source sample
    //
    // Scoring for other samples will be added later after all variant CNV records are consolidated downstream
    // of this method.
    //
    let sample_index = source_sample_index;
    let chrom_depth_bins = &sample_seg_input.genome_depth_bins.depth_bins[chrom_index];
    let sample_gc_bias_data = sample_seg_input.sample_gc_bias_data;
    let gc_corrected_haploid_coverage = sample_seg_input.gc_corrected_haploid_coverage;

    let state_lnlhoods = get_sample_cn_state_lnlhoods(
        cn_segment,
        chrom_depth_bins,
        chrom_gc_levels,
        sample_gc_bias_data,
        gc_corrected_haploid_coverage,
        sample_seg_input
            .discover_settings
            .bin_dependency_correction_factor,
    );

    get_cnv_sample_qualities_from_cn_state_lnlhoods(
        max_qscore,
        &cn_segment.copy_number_info,
        &state_lnlhoods,
        &mut score.samples[sample_index],
    );

    Some(RefinedCNV {
        segment,
        source_sample_index,
        score,
    })
}

/// Estimate haploid coverage for all samples
///
fn estimate_haploid_coverage_by_sample(
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    debug: bool,
) -> Vec<Option<HaploidCoverage>> {
    let mut x = Vec::new();
    for (sample_index, sample_data) in all_sample_data.iter().enumerate() {
        let haploid_coverage =
            sample_data
                .sample_gc_bias_data
                .as_ref()
                .map(|sample_gc_bias_data| {
                    let coverage_est_regex =
                        sample_data.discover_settings.coverage_est_regex.as_str();
                    let copy_number_segments = Some(&sample_data.copy_number_segments);
                    get_haploid_genome_coverage(
                        coverage_est_regex,
                        &shared_data.chrom_list,
                        &sample_data.genome_depth_bins,
                        &shared_data.genome_gc_levels,
                        sample_gc_bias_data,
                        copy_number_segments,
                    )
                });

        if debug {
            eprintln!("CNV haploid_coverage for sample {sample_index}:\n{haploid_coverage:?}");
        }
        x.push(haploid_coverage);
    }
    x
}

/// Input contract is that rcnv_merge_set contains a set of CNVs with exact segment match
///
fn consolidate_multi_sample_cnv(
    rcnv_merge_set: &mut Vec<RefinedCNV>,
    multisample_refined_cnvs: &mut Vec<RefinedCNV>,
) {
    if rcnv_merge_set.len() <= 1 {
        multisample_refined_cnvs.append(rcnv_merge_set);
    } else {
        let mut merged_rcnv = rcnv_merge_set.pop().unwrap();
        for rcnv in rcnv_merge_set.drain(..) {
            assert_eq!(rcnv.segment, merged_rcnv.segment);

            for (sample_index, sample_score) in rcnv.score.samples.into_iter().enumerate() {
                if sample_score.shared.copy_info.is_some() {
                    assert!(
                        merged_rcnv.score.samples[sample_index]
                            .shared
                            .copy_info
                            .is_none()
                    );
                    merged_rcnv.score.samples[sample_index] = sample_score;
                }
            }
        }

        multisample_refined_cnvs.push(merged_rcnv);
    }
}

/// Convert all copy number segmentation results into the "RefinedCNV" type which is ready for VCF output or merging with SV calls
///
/// Note that `all_sample_data` contains segmentation results from the discover step (in all_sample_data.copy_number_segments), but
/// that this function accepts segmentation results as a separate argument, allowing for newer segmentation results updated with
/// joint-sample and/or breakpoint information.
///
/// This method also computes the new copy number likelihoods and qualities for each CNV for eventual use in the VCF output.
///
/// # Arguments
/// * `max_qscore` - Maximum value to record for the copy number q-score
/// * `copy_number_segments` - segmentation results over all samples to be converted to RefinedCNV format
///
fn get_refined_cnvs(
    max_qscore: u32,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    seg_input: &[SampleCopyNumberSegmentationInput],
    cn_segments: &[SampleCopyNumberSegments],
    debug: bool,
) -> Vec<RefinedCNV> {
    let mut refined_cnvs = Vec::new();
    assert_eq!(seg_input.len(), cn_segments.len());
    for (sample_seg_input, sample_cn_segments) in seg_input.iter().zip(cn_segments.iter()) {
        let sample_index = sample_seg_input.sample_index;
        for (chrom_index, chrom_cn_segments) in sample_cn_segments.cn_segments.iter().enumerate() {
            let chrom_gc_levels = &shared_data.genome_gc_levels[chrom_index];

            for cn_segment in chrom_cn_segments.iter() {
                let rc_result = convert_cn_segment_to_refined_cnv(
                    all_sample_data,
                    max_qscore,
                    sample_index,
                    chrom_index,
                    cn_segment,
                    sample_cn_segments.bin_size,
                    chrom_gc_levels,
                    sample_seg_input,
                );

                if let Some(refined_cnv) = rc_result {
                    refined_cnvs.push(refined_cnv);
                }
            }
        }
    }

    let sample_count = all_sample_data.len();
    let mut multisample_refined_cnvs = if sample_count <= 1 {
        refined_cnvs
    } else {
        // First refined_cnv list is created with sample-independent CNV records, now merge exact segment matches into
        // multi-sample records and update QUAL scores
        refined_cnvs.sort_by_key(|x| x.segment.clone());

        let mut multisample_refined_cnvs = Vec::new();

        let mut rcnv_merge_set: Vec<RefinedCNV> = Vec::new();

        for rcnv in refined_cnvs.into_iter() {
            if !rcnv_merge_set.is_empty() && rcnv.segment != rcnv_merge_set[0].segment {
                consolidate_multi_sample_cnv(&mut rcnv_merge_set, &mut multisample_refined_cnvs);
            }
            rcnv_merge_set.push(rcnv);
        }
        consolidate_multi_sample_cnv(&mut rcnv_merge_set, &mut multisample_refined_cnvs);
        multisample_refined_cnvs
    };

    // Final step: Update the QUAL score for all CNVs
    //
    // This is delayed until after all sample consolidation
    //
    for rcnv in multisample_refined_cnvs.iter_mut() {
        rcnv.score.update_cnv_alt_score(max_qscore);
    }

    if debug {
        eprintln!("Input cnv count {}", multisample_refined_cnvs.len());
    }

    multisample_refined_cnvs
}

/// Return true if SV meets all criteria for the CNV matching routine
///
fn sv_eligible_for_cnv_matching(rsv: &RefinedSV, min_sv_size: usize, min_qual: f32) -> bool {
    if rsv.filter_sv() {
        // These cases will always be filtered from the VCF output entirely
        false
    } else if !rsv.bp.is_indel_or_dup() {
        // Filter out all types except indel/DUP:
        false
    } else if rsv.bp.deldup_size().unwrap() < min_sv_size {
        // Filter SVs below minimum size
        false
    } else {
        // These cases are in the VCF but non-PASS
        //
        // TODO: refactor to keep filter changes in sync between the VCF and here
        let is_max_scoring_depth = rsv.score.samples.iter().any(|x| x.is_max_scoring_depth);
        if is_max_scoring_depth {
            false
        } else if let Some(sv_alt_score) = rsv.score.sv_alt_score {
            sv_alt_score >= min_qual
        } else {
            false
        }
    }
}

/// Return true if the SV GT is compatible with the CNV CN in this sample, or false if they are incompatible.
///
/// Return None if either sample has missing information, or if both samples are non-variant.
///
fn is_gt_cn_match(sample_index: usize, rcnv: &RefinedCNV, rsv: &RefinedSV) -> Option<bool> {
    let rcnv_sample_score = &rcnv.score.samples[sample_index];
    let copy_number = match &rcnv_sample_score.shared.copy_info {
        Some(x) => x.copy_number,
        None => {
            return None;
        }
    };

    let be1 = &rsv.bp.breakend1;
    let expected_copy_number = rcnv_sample_score.expected_copy_number;
    let rsv_sample_score = &rsv.score.samples[sample_index];

    if be1.dir == BreakendDirection::LeftAnchor {
        // SV DEL:
        //
        // Return true if this is a CN loss matching the SV GT value
        match rsv_sample_score.shared.gt {
            Some(Genotype::Het) => {
                let is_gt_cn_match = copy_number == expected_copy_number.saturating_sub(1);
                Some(is_gt_cn_match)
            }
            Some(Genotype::Hom) => {
                let is_gt_cn_match = copy_number == 0;
                Some(is_gt_cn_match)
            }
            Some(Genotype::Ref) => {
                if copy_number == expected_copy_number {
                    None
                } else {
                    Some(false)
                }
            }
            _ => None,
        }
    } else {
        // SV DUP:
        //
        // Translating a GT to a CN gain match is a little too unreliable, so leave this for the DEL section only,
        // we just require any level of gain for het or hom
        match rsv_sample_score.shared.gt {
            Some(Genotype::Het) | Some(Genotype::Hom) => {
                let is_gt_cn_match = copy_number > expected_copy_number;
                Some(is_gt_cn_match)
            }
            Some(Genotype::Ref) => {
                if copy_number == expected_copy_number {
                    None
                } else {
                    Some(false)
                }
            }
            _ => None,
        }
    }
}

/// CNVToSVMatchTracker key is the cnv index, value is the list of matched sv indexes
///
type CNVToSVMatchTracker = BTreeMap<usize, BTreeSet<usize>>;

/// Describes an SV's candidate CNV match in one sample
#[derive(Clone, Debug, PartialEq)]
struct SVToCNVSampleMatchInfo {
    /// Index of the candidate CNV
    refined_cnv_index: usize,

    /// Total distance separating the left and right SV & CNV span edges
    total_breakend_dist: i64,

    /// True if the SV GT is compatible with the CNV CN shift in this sample
    gt_cn_match: bool,
}
struct SVToCNVMatchInfo {
    samples: Vec<Option<SVToCNVSampleMatchInfo>>,
}

/// key is the sv index
type SVToCNVMatchTracker = BTreeMap<usize, SVToCNVMatchInfo>;

/// Review SV2CNV entries and eliminate any SVs with less than min_gt_cn_match_fraction
///
/// Make corresponding updates to CNV2SV structure during SV entry elimination
///
fn remove_poorly_matching_svs(
    sv_to_cnv_matches: &mut SVToCNVMatchTracker,
    cnv_to_sv_matches: &mut CNVToSVMatchTracker,
    min_gt_cn_match_fraction: f64,
) {
    sv_to_cnv_matches.retain(|sv_group_index, sv2cnv_match_info| {
        let gt_cn_match_fraction = {
            let mut total = 0;
            let mut gt_cn_match = 0;
            for sv_match in sv2cnv_match_info.samples.iter().filter_map(|x| x.as_ref()) {
                total += 1;
                if sv_match.gt_cn_match {
                    gt_cn_match += 1;
                }
            }
            gt_cn_match as f64 / total as f64
        };

        if gt_cn_match_fraction < min_gt_cn_match_fraction {
            // Remove this SV from the CNV2SV data structure:
            for sv_match in sv2cnv_match_info.samples.iter().filter_map(|x| x.as_ref()) {
                let remove_key = if let Some(cnv_match) =
                    cnv_to_sv_matches.get_mut(&sv_match.refined_cnv_index)
                {
                    cnv_match.remove(sv_group_index);
                    cnv_match.is_empty()
                } else {
                    false
                };

                if remove_key {
                    cnv_to_sv_matches.remove(&sv_match.refined_cnv_index);
                }
            }

            false
        } else {
            true
        }
    });
}

/// Enumerate information on how a specific SV matches CNVs over all samples
///
/// This includes metrics to determine which SV to retain in case of multiple overlapping SV matches to the same CNV
///
struct SVCNVMatchScore {
    /// Total number of samples where the SV has a candiate CNV match with compatible CN shift direction for the SV type
    match_count: usize,

    /// Sum of the total breakend distances over all samples qualifying to be counted in match_count
    all_match_breakend_dist: i64,

    /// Index of the SV group within the full joint-call SVGroup list
    sv_group_index: usize,

    /// Information used to compose the unique VCF ID entry for the given SV
    id: SVUniqueIdData,
}

/// Enumerate information on how a specific SV matches CNVs over all samples
///
/// Return None if the SV has no CNV matches
///
fn get_sv_cnv_match_score(
    sv_to_cnv_matches: &SVToCNVMatchTracker,
    sv_group_index: usize,
    sv_groups: &[SVGroup],
) -> Option<SVCNVMatchScore> {
    if let Some(sv_match) = sv_to_cnv_matches.get(&sv_group_index) {
        let mut match_count = 0;
        let mut all_match_breakend_dist = 0;
        for sv_sample_match in sv_match
            .samples
            .iter()
            .filter_map(|x| x.as_ref())
            .filter(|&x| x.gt_cn_match)
        {
            match_count += 1;
            all_match_breakend_dist += sv_sample_match.total_breakend_dist;
        }

        let id = sv_groups[sv_group_index].refined_svs[0].id.clone();

        Some(SVCNVMatchScore {
            match_count,
            all_match_breakend_dist,
            sv_group_index,
            id,
        })
    } else {
        None
    }
}

/// Merge CNVs with similar large SVs
///
/// Assumed pre-conditions:
/// - Each sv_group has been scored and has expanded sample_score array
/// - Each refined_cnv has been scored and has an expanded sample_score array
///
fn merge_cnv_with_sv(
    settings: &JointCallSettings,
    sample_count: usize,
    sv_groups: &mut [SVGroup],
    refined_cnvs: &mut Vec<RefinedCNV>,
) {
    //
    // Hard code various parameters here for now
    //

    // Min SV size for CNV matching
    //
    // TODO: there's a lot of work to do better here for small SVs, this min size doesn't make sense with th max edge dist above for instance.
    //
    let min_sv_size = 10000;

    // Max dist between copy change boundary and breakend location
    //
    let max_sv_cnv_edge_dist = 8000;

    // Min reciprocal overlap between the sv and cnv spans. This filter is applied together
    // with the edge dist filter
    //
    let min_sv_cnv_recip_overlap = 0.9;

    let min_gt_cn_match_fraction = 0.5;

    // Start by sorting the sv_groups on their first sv's first breakend
    //
    sv_groups.sort_by_key(|x| x.refined_svs[0].bp.breakend1.segment.clone());

    // Next sort CNV segments by their region
    refined_cnvs.sort_by_key(|x| x.segment.clone());

    let refined_cnv_count = refined_cnvs.len();

    let mut head_refined_cnv_index = 0;

    // Track sv matches for each cnv segment. This structure is used to pick only one SV to be the winner which will
    // retain the CNV match in case of multiple matches
    //
    let mut cnv_to_sv_matches: CNVToSVMatchTracker = Default::default();

    // Track best match of each SV to CNVs, where there can be up to one best CNV per sample
    //
    // We potentially have multiple overlapping SVs (likely one of them is a FP duplicate record), overlapping with
    // multiple CNVs (one per sample, unlikely to represent FP duplicate records).
    //
    // Strategy to resolve this?
    // - If any two SVs share any CNV attachments, keep the SV record with (1) the most matches over all samples
    // (2) the lowest total_breakend_dist over all matching samples if match count is equal.
    //
    let mut sv_to_cnv_matches: SVToCNVMatchTracker = Default::default();

    for (sv_group_index, sv_group) in sv_groups
        .iter_mut()
        .enumerate()
        .filter(|(_, x)| !x.is_single_region())
    {
        // Filtering for non single-region SVs means this should always be true:
        assert_eq!(sv_group.refined_svs.len(), 1);

        let rsv = &mut sv_group.refined_svs[0];
        if !sv_eligible_for_cnv_matching(rsv, min_sv_size, settings.min_qual as f32) {
            continue;
        }

        let be1 = &rsv.bp.breakend1;
        let be2 = rsv.bp.breakend2.as_ref().unwrap();
        let sv_span = IntRange::from_pair(be1.segment.range.start, be2.segment.range.start);

        let mut is_first_cnv_be1_match = false;

        // Wind CNVs forward to see if we find matches to this SV
        //
        // Multiple matches are possible because we could have separate CNVs on each sample
        //
        for (refined_cnv_index, rcnv) in
            refined_cnvs.iter().enumerate().skip(head_refined_cnv_index)
        {
            match rcnv.segment.chrom_index.cmp(&be1.segment.chrom_index) {
                Ordering::Less => {
                    // The CNV chrom index is less than the SV chrom index, so wind CNV list forward
                    continue;
                }
                Ordering::Greater => {
                    // The CNV chrom index is greater than the SV chrom index, so no point in continuing search
                    break;
                }
                _ => (),
            };

            // SV and CNV chrom_index values must be equal at this point

            if rcnv.segment.range.start + max_sv_cnv_edge_dist < be1.segment.range.start {
                // The CNV start is less than the SV start + edge dist (ie. the CNV cannot match), so wind CNV list forward
                continue;
            } else if rcnv.segment.range.start > be1.segment.range.start + max_sv_cnv_edge_dist {
                // The CNV start is more than the SV start + edge dist (ie. the CNV cannot match), so no point in continuing search
                break;
            }

            // SV and CNV breakend1 positions must be in range at this point

            if !is_first_cnv_be1_match {
                // Mark the first CNV which comes close enough to this SV's breakend1. We'll restore to this CNV index when
                // the search for the next SV starts.
                head_refined_cnv_index = refined_cnv_index;
                is_first_cnv_be1_match = true;
            }

            if rcnv.segment.range.end + max_sv_cnv_edge_dist < be2.segment.range.start {
                // The CNV end is less than the SV start + edge dist (ie. the CNV cannot match), so wind the CNV list forward
                continue;
            } else if rcnv.segment.range.end > be2.segment.range.start + max_sv_cnv_edge_dist {
                // The CNV end is more than the SV start + edge dist (ie. the CNV cannot match). In this case we still move
                // the CNV list forward to check for matches from other samples.
                continue;
            }

            // SV and CNV have a full span match (within max_sv_cnv_edge_dist) at this point

            // Also check recip overlap, this helps filter incorrect smaller sv/cnv matches
            let recip_overlap = get_recip_overlap(&rcnv.segment.range, &sv_span);
            if recip_overlap < min_sv_cnv_recip_overlap {
                continue;
            }

            let total_breakend_dist = (rcnv.segment.range.start - be1.segment.range.start).abs()
                + (rcnv.segment.range.end - be2.segment.range.start).abs();

            // Match type to breakend direction to make sure DEL matches to a CN loss, DUP matches to a CN gain
            //
            // This match is tested separately for each sample defined in the refined CNV
            //
            for sample_index in 0..sample_count {
                let gt_cn_match = match is_gt_cn_match(sample_index, rcnv, rsv) {
                    Some(x) => x,
                    None => {
                        continue;
                    }
                };

                let s2c_entry =
                    sv_to_cnv_matches
                        .entry(sv_group_index)
                        .or_insert(SVToCNVMatchInfo {
                            samples: vec![None; sample_count],
                        });
                assert_eq!(s2c_entry.samples[sample_index], None);
                s2c_entry.samples[sample_index] = Some(SVToCNVSampleMatchInfo {
                    refined_cnv_index,
                    total_breakend_dist,
                    gt_cn_match,
                });

                let c2s_entry = cnv_to_sv_matches.entry(refined_cnv_index).or_default();
                c2s_entry.insert(sv_group_index);
            }
        }
    }

    // Now process the sv to cnv mapping structures to make sure that each CNV is matched to no more than one SV

    remove_poorly_matching_svs(
        &mut sv_to_cnv_matches,
        &mut cnv_to_sv_matches,
        min_gt_cn_match_fraction,
    );

    // Next, identify CNVs linked to more than one SV, for each such case, filter the SVs down to only one matching case
    for cnv_match in cnv_to_sv_matches.values().filter(|v| v.len() > 1) {
        let mut sv_index_scores = Vec::new();
        for &sv_group_index in cnv_match.iter() {
            match get_sv_cnv_match_score(&sv_to_cnv_matches, sv_group_index, sv_groups) {
                Some(x) => sv_index_scores.push(x),
                None => continue,
            }
        }

        if sv_index_scores.len() < 2 {
            // The CNV already matches no more than one SV, so the SV ranking and filtration steps below can be skipped
            continue;
        }

        // Use id as a final tie-break to keep the results deterministic, and to keep the CNV matches on
        // the lowest id in case of a tie, to keep it out of the dedup filter:
        sv_index_scores.sort_by(|a, b| {
            b.match_count
                .cmp(&a.match_count)
                .then(a.all_match_breakend_dist.cmp(&b.all_match_breakend_dist))
                .then(a.id.cmp(&b.id))
        });

        for sv_index_score in sv_index_scores.iter().skip(1) {
            sv_to_cnv_matches.remove(&sv_index_score.sv_group_index);
        }
    }

    // Now process the remaining items in sv_to_cnv_matches to update the corresponding sv and cnv data structures for merged vcf output
    let mut merged_cnv_indexes = BTreeSet::new();
    for (&sv_group_index, sv_to_cnv_match) in sv_to_cnv_matches.iter() {
        assert_eq!(sv_groups[sv_group_index].refined_svs.len(), 1);

        let rsv = &mut sv_groups[sv_group_index].refined_svs[0];

        rsv.ext.is_cnv_match = true;

        for (sample_index, sv_to_cnv_sample_match) in sv_to_cnv_match
            .samples
            .iter()
            .enumerate()
            .filter_map(|(i, x)| x.as_ref().map(|y| (i, y)))
            .filter(|(_, x)| x.gt_cn_match)
        {
            merged_cnv_indexes.insert(sv_to_cnv_sample_match.refined_cnv_index);
            let rcnv = &refined_cnvs[sv_to_cnv_sample_match.refined_cnv_index];
            let rcnv_sample_score = &rcnv.score.samples[sample_index];
            let rsv_sample_score = &mut rsv.score.samples[sample_index];
            rsv_sample_score.shared.copy_info = rcnv_sample_score.shared.copy_info.clone();
        }
        rsv.score.update_cnv_alt_score(settings.max_qscore);
    }

    info!(
        "Merged {} CNVs into corresponding breakpoint-based SVs",
        merged_cnv_indexes.len()
    );

    // Remove the CNVs that have been merged into similar SVs already:
    if !merged_cnv_indexes.is_empty() {
        remove_sorted_indices(refined_cnvs, merged_cnv_indexes.iter().cloned());
        assert_eq!(
            refined_cnv_count,
            refined_cnvs.len() + merged_cnv_indexes.len()
        );
    }
}

/// Mark any large SVs without CNV support for filtration into BND
///
fn mark_unsupported_large_deldup_for_bnd_conversion(
    min_sv_required_depth_support_size: usize,
    sv_groups: &mut [SVGroup],
) {
    for sv_group in sv_groups.iter_mut().filter(|x| !x.is_single_region()) {
        // Filtering for non single-region SVs means this should always be true:
        assert_eq!(sv_group.refined_svs.len(), 1);

        // Get SV breakends and filter out various types of non-CNV SVs:
        let rsv = &mut sv_group.refined_svs[0];
        if !rsv.bp.is_indel_or_dup() {
            // Filter out all types except indel/DUP:
            continue;
        }

        if rsv.bp.deldup_size().unwrap() < min_sv_required_depth_support_size {
            continue;
        }

        if rsv.ext.is_cnv_match {
            continue;
        }

        rsv.ext.force_breakpoint_representation = true;
    }
}

/// Reformat contents of all_sample_data to provide inputs needed for segmentation and CNV calling
///
fn get_segmentation_input<'a>(
    shared_data: &SharedJointCallData,
    all_sample_data: &'a [SampleJointCallData],
    transition_types: &'a [SampleTransitionTypes],
    debug: bool,
) -> Vec<SampleCopyNumberSegmentationInput<'a>> {
    let haploid_coverage_by_sample =
        estimate_haploid_coverage_by_sample(shared_data, all_sample_data, debug);

    let mut seg_input: Vec<SampleCopyNumberSegmentationInput<'_>> = Vec::new();
    for (sample_index, sample_data) in all_sample_data.iter().enumerate() {
        let haploid_coverage = match haploid_coverage_by_sample[sample_index].as_ref() {
            Some(x) => x,
            None => continue,
        };
        let sample_transition_types = Some(&transition_types[sample_index]);
        let sample_seg_input = SampleCopyNumberSegmentationInput {
            sample_index,
            discover_settings: &sample_data.discover_settings,
            gc_corrected_haploid_coverage: haploid_coverage.gc_corrected_depth,
            genome_depth_bins: &sample_data.genome_depth_bins,
            sample_gc_bias_data: sample_data.sample_gc_bias_data.as_ref().unwrap(),
            sample_transition_types,
            expected_copy_number_regions: sample_data.expected_copy_number_regions.as_ref(),
        };
        seg_input.push(sample_seg_input);
    }
    seg_input
}

fn estimate_high_state_depth_in_sample(
    shared_data: &SharedJointCallData,
    sample_seg_input: &SampleCopyNumberSegmentationInput,
    sample_cn_segments: &mut SampleCopyNumberSegments,
) {
    let sample_gc_bias_data = sample_seg_input.sample_gc_bias_data;
    let gc_corrected_haploid_depth_factor = 1. / sample_seg_input.gc_corrected_haploid_coverage;

    for (chrom_index, chrom_cn_segments) in sample_cn_segments.cn_segments.iter_mut().enumerate() {
        let chrom_depth_bins = &sample_seg_input.genome_depth_bins.depth_bins[chrom_index];
        let chrom_gc_levels = &shared_data.genome_gc_levels[chrom_index];

        for segment in chrom_cn_segments
            .iter_mut()
            .filter(|x| x.copy_number_info.state == CopyNumberState::High)
        {
            // In a single bin, solve for cn from:
            // haploid_coverage * gc_depth_reduction * copy_num_state = depth
            //
            // so...
            // copy_num_state = (depth / haploid_coverage * gc_depth_reduction)
            //
            // For a segment:
            // copy_num_state = 1/hap_coverage * sum(depth/gc_depth_reduction)/count
            //

            // Get GC-corrected bin average
            //
            // Excluded bins are skipped from the average. Set a min number & fraction of non-excluded observations where we bail out on depth
            // estimation for the segment, we shouldn't hit this if the segment got a High classification to begin with, so this is more of
            // a worst-case backstop assertion.
            //
            let min_nonexcluded_bin_count = 4;
            let min_nonexcluded_bin_fraction = 0.8;

            let mut nonexcluded_bin_count = 0;
            let mut total_gc_corrected_depth = 0f64;
            for bin_index in segment.begin_bin..segment.end_bin {
                let depth = match chrom_depth_bins[bin_index] {
                    DepthBin::Excluded => {
                        continue;
                    }
                    DepthBin::Depth(x) => x,
                };
                assert!(depth >= 0.0);
                nonexcluded_bin_count += 1;
                let gc_depth_reduction =
                    sample_gc_bias_data.gc_depth_reduction[chrom_gc_levels[bin_index]];

                total_gc_corrected_depth += depth / gc_depth_reduction;
            }

            if nonexcluded_bin_count < min_nonexcluded_bin_count {
                continue;
            }

            let segment_bin_count = segment.end_bin - segment.begin_bin;
            let nonexcluded_bin_fraction = nonexcluded_bin_count as f64 / segment_bin_count as f64;
            if nonexcluded_bin_fraction < min_nonexcluded_bin_fraction {
                continue;
            }

            let avg_gc_corrected_depth = total_gc_corrected_depth / nonexcluded_bin_count as f64;

            let estimated_copy_number = f64::max(
                gc_corrected_haploid_depth_factor * avg_gc_corrected_depth,
                CopyNumberState::High as usize as f64,
            );

            segment.copy_number_info.high_state_copy_number = Some(estimated_copy_number);
        }
    }
}

fn estimate_high_state_depth(
    shared_data: &SharedJointCallData,
    seg_input: &[SampleCopyNumberSegmentationInput],
    cn_segments: &mut [SampleCopyNumberSegments],
) {
    assert_eq!(seg_input.len(), cn_segments.len());
    for (sample_seg_input, sample_cn_segments) in seg_input.iter().zip(cn_segments.iter_mut()) {
        estimate_high_state_depth_in_sample(shared_data, sample_seg_input, sample_cn_segments);
    }
}

fn write_resegment_results_to_file(
    output_dir: &Utf8Path,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    seg_input: &[SampleCopyNumberSegmentationInput],
    cn_segments: &[SampleCopyNumberSegments],
) {
    use super::sample_output::get_sample_output_dir;

    assert_eq!(seg_input.len(), cn_segments.len());
    for (sample_seg_input, sample_cn_segments) in seg_input.iter().zip(cn_segments.iter()) {
        let sample_index = sample_seg_input.sample_index;
        let sample_data = &all_sample_data[sample_index];
        let sample_dir = get_sample_output_dir(output_dir, sample_index, &sample_data.sample_name);
        write_copy_number_segment_file(&sample_dir, &shared_data.chrom_list, sample_cn_segments);
    }
}

/// Merge joint-called SV output with CN segmentation information over all samples
///
/// The merge operation has many components, as describe in sawfish methods. Major components include:
/// 1. Resegment depth tracks in all samples based on SV breakpoint locations
/// 2. Complete 'high-depth' state depth estimation
/// 3. Write resegmented copy number to per-sample file outputs
/// 4. Complete multi-sample CNV boundary reconsiliation
/// 5. Merge CNVs and compute CN qualities
/// 6. Merge SV/CNV records
/// 7. Mark large unmerged SVs for BND conversion
///
/// It is assumed that all CNV parameters are identical in all input samples. This is not
/// actually checked anywhere yet.
///
pub fn merge_sv_with_depth_info(
    settings: &JointCallSettings,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    sv_groups: &mut [SVGroup],
) -> Vec<RefinedCNV> {
    info!("Starting joint CNV/SV refinement");

    let debug = false;

    let sample_count = all_sample_data.len();

    // Assert some input contracts
    {
        assert!(sample_count > 0);
        let first_sample_data = &all_sample_data[0];
        let bin_size = first_sample_data.genome_depth_bins.bin_size;

        // Check that bin_size is the same in all samples:
        assert!(
            all_sample_data
                .iter()
                .all(|x| x.genome_depth_bins.bin_size == bin_size),
            "All samples do not have the same bin_size"
        );
    }

    //
    // Hard coded parameters for this section
    //

    // Determines how much the copy number transition probability is relaxed at SV breakpoints during re-segmentation:
    //
    let breakpoint_transition_factor = 0.5;

    // Size above which deletions and duplications are reported as breakends if they lack correspondence
    // with depth segmentation.
    //
    let min_sv_required_depth_support_size = 50_000;

    // TO ADD LATER?: Loop over and identify all large overlapping SVs in each sample (should we do this as part of SV/CNV overlap detection instead?)

    // Find copy-number transition bonus locations from SV brekaends
    let transition_types =
        get_copynum_transition_types_from_breakends(shared_data, all_sample_data, sv_groups);

    // Rearrange sample data references to get all required CNV-calling elements for each sample
    let seg_input = get_segmentation_input(shared_data, all_sample_data, &transition_types, debug);

    // Resegment depth for every sample
    let mut cn_segments = get_single_pass_sample_copy_number_segments(
        &shared_data.genome_gc_levels,
        breakpoint_transition_factor,
        &seg_input,
    );

    // Estimate depth of "High" state segments
    estimate_high_state_depth(shared_data, &seg_input, &mut cn_segments);

    // Write new depth segments to file:
    write_resegment_results_to_file(
        &settings.output_dir,
        shared_data,
        all_sample_data,
        &seg_input,
        &cn_segments,
    );

    write_copy_number_info(
        &settings.output_dir,
        shared_data,
        all_sample_data,
        &seg_input,
        &cn_segments,
    );

    // Convert new segments to RefinedCNV type
    let mut refined_cnvs = get_refined_cnvs(
        settings.max_qscore,
        shared_data,
        all_sample_data,
        &seg_input,
        &cn_segments,
        debug,
    );

    // Run final SV <-> CNV matching operation, and update both SV and CNV records to reflect merging updates
    if !settings.disable_sv_cnv_merge {
        merge_cnv_with_sv(settings, sample_count, sv_groups, &mut refined_cnvs);
    }

    mark_unsupported_large_deldup_for_bnd_conversion(min_sv_required_depth_support_size, sv_groups);

    refined_cnvs
}
