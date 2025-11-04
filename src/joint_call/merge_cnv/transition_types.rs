use std::collections::BTreeMap;

use rust_vc_utils::int_range::{IntRange, get_recip_overlap};

use crate::breakpoint::BreakendDirection;
use crate::copy_number_segmentation::{CopyNumberState, SampleTransitionTypes};
use crate::joint_call::{SampleJointCallData, SharedJointCallData};
use crate::refine_sv::Genotype;
use crate::sv_group::SVGroup;

// This struct is used to maintain a pointer to the highest segment visited so far, enabling a per-sample
// merge of SVs to CN segments during the matching procedure:
#[derive(Clone, Default)]
struct CNSegmentID {
    chrom_index: usize,
    segment_index: usize,
}

/// Return true for large copy-changing SVs
fn is_large_deldup(sv_group: &SVGroup, min_copynum_checked_sv_size: usize) -> bool {
    if sv_group.is_single_region() {
        return false;
    }

    assert_eq!(sv_group.refined_svs.len(), 1);
    let rsv = &sv_group.refined_svs[0];
    rsv.bp
        .deldup_size()
        .is_some_and(|x| x >= min_copynum_checked_sv_size)
}

fn annotate_segmentation_power_boost_breakends(
    all_sample_data: &[SampleJointCallData],
    sv_groups: &[SVGroup],
    transition_types: &mut [SampleTransitionTypes],
    min_copynum_checked_sv_size: usize,
    debug: bool,
) {
    // For the copy-number track to be resegmented around a smaller SV with relaxed transition scores, the copy number
    // must be consistent already throughout the span of the SV AND through a left and right edge region of the
    // following size:
    let sv_edge_size = 5000;

    let sample_count = all_sample_data.len();
    assert!(sample_count > 0);

    let first_sample_data = &all_sample_data[0];
    let bin_size = first_sample_data.copy_number_segments.bin_size;

    // Iterate through large (multi-region SVs) to find those with the required CN segment overlap condition
    //
    let mut head_cn_segment = vec![CNSegmentID::default(); sample_count];
    for sv_group in sv_groups
        .iter()
        .filter(|&x| is_large_deldup(x, min_copynum_checked_sv_size))
    {
        assert_eq!(sv_group.refined_svs.len(), 1);
        let rsv = &sv_group.refined_svs[0];

        let be1 = &rsv.bp.breakend1;
        let be2 = rsv.bp.breakend2.as_ref().unwrap();
        let chrom_index = be1.segment.chrom_index;

        // Since CNV integration is still an essentially single-sample approach at this point, we'll set this up
        // as if each sample is independent:
        for (sample_index, sample_score) in rsv.score.samples.iter().enumerate() {
            if sample_score.is_max_scoring_depth {
                continue;
            }

            // Determine if the SV has a non-reference GT in this sample:
            if sample_score
                .shared
                .gt
                .as_ref()
                .is_none_or(|x| *x == Genotype::Ref)
            {
                continue;
            }

            let sample_transition_types = &mut transition_types[sample_index];

            // Wind CNV index forward from head segment to find the first segment intersecting this SV
            let sample_head_cn_segment = &mut head_cn_segment[sample_index];
            assert!(chrom_index >= sample_head_cn_segment.chrom_index);
            if chrom_index > sample_head_cn_segment.chrom_index {
                sample_head_cn_segment.chrom_index = chrom_index;
                sample_head_cn_segment.segment_index = 0;
            }

            let sample_cn_segments = &all_sample_data[sample_index].copy_number_segments;
            let chrom_cn_segments = &sample_cn_segments.cn_segments[chrom_index];
            let start_seg_index = sample_head_cn_segment.segment_index;
            for (cn_segment_index, cn_segment) in
                chrom_cn_segments.iter().enumerate().skip(start_seg_index)
            {
                let cn_segment_pos_range = cn_segment.to_range(bin_size);
                if cn_segment_pos_range.end > be1.segment.range.start {
                    // This assertion checks that we don't 'skip-over' the SV, which shouldn't be possible given that
                    // CN segments the entire genome
                    assert!(cn_segment_pos_range.start <= be2.segment.range.end);

                    sample_head_cn_segment.segment_index = cn_segment_index;

                    // Check if all criteria are met for the creation of a low-cost copy-number state transition on this
                    // SV's boundaries
                    let context_segment_large_enough = (cn_segment_pos_range.start + sv_edge_size)
                        <= be1.segment.range.start
                        && (be2.segment.range.end + sv_edge_size) <= cn_segment_pos_range.end;

                    let context_copy_number_consistent_with_sv = {
                        if be1.dir == BreakendDirection::LeftAnchor {
                            // DEL
                            cn_segment.copy_number_info.state != CopyNumberState::Zero
                        } else {
                            // DUP
                            cn_segment.copy_number_info.state != CopyNumberState::High
                        }
                    };

                    if context_copy_number_consistent_with_sv && context_segment_large_enough {
                        let start_pos = be1.segment.range.center();
                        let end_pos = be2.segment.range.center();
                        let is_del = be1.dir == BreakendDirection::LeftAnchor;
                        if is_del {
                            // DEL
                            sample_transition_types.set_loss_pos(chrom_index, start_pos);
                            sample_transition_types.set_gain_pos(chrom_index, end_pos);
                        } else {
                            // DUP
                            sample_transition_types.set_gain_pos(chrom_index, start_pos);
                            sample_transition_types.set_loss_pos(chrom_index, end_pos);
                        }
                        if debug {
                            let label = if is_del { "del" } else { "dup" };
                            eprintln!(
                                "Adding {label}-like-bnd power-boost segment annotation at: {chrom_index}:{start_pos}-{end_pos}"
                            );
                        }
                    }

                    // We're done marking this SV in this sample, so move to the next sample (or next eligible SV) without
                    // scanning any more copynum segments
                    break;
                }
            }
        }
    }
}

fn annotate_close_fit_breakends(
    all_sample_data: &[SampleJointCallData],
    sv_groups: &[SVGroup],
    transition_types: &mut [SampleTransitionTypes],
    min_copynum_checked_sv_size: usize,
    debug: bool,
) {
    // Hard-code parameters
    let min_sv_cnv_recip_overlap = 0.8;

    let sample_count = all_sample_data.len();
    assert!(sample_count > 0);

    let first_sample_data = &all_sample_data[0];
    let bin_size = first_sample_data.copy_number_segments.bin_size;

    let chrom_count = transition_types[0].sample_transition_types.len();

    #[derive(Default, Clone)]
    struct SVMatchInfo {
        /// If multiple SVs match the CN segment in one sample, then the annotation is not applied
        match_count: usize,

        start_pos: i64,
        end_pos: i64,

        /// State to be either DEL or DUP, so just use a bool for now
        is_del: bool,
    }

    // This data structure is used to keep us from associating 2 close fit SVs with one CN segment:
    //
    // The map is <usize,SVMatchInfo>. The key is the cn segment number, and the value is everything
    // we need to know about the matching sv to set transition probs correctly
    //
    let mut cn_segment_sv_match_tracker = vec![vec![BTreeMap::new(); chrom_count]; sample_count];

    // Iterate through large (multi-region SVs) to find those with the required CN segment overlap condition
    //
    let mut head_cn_segment = vec![CNSegmentID::default(); sample_count];
    for sv_group in sv_groups
        .iter()
        .filter(|&x| is_large_deldup(x, min_copynum_checked_sv_size))
    {
        assert_eq!(sv_group.refined_svs.len(), 1);
        let rsv = &sv_group.refined_svs[0];

        let be1 = &rsv.bp.breakend1;
        let be2 = rsv.bp.breakend2.as_ref().unwrap();
        let chrom_index = be1.segment.chrom_index;

        let sv_span = IntRange::from_pair(be1.segment.range.start, be2.segment.range.start);

        for (sample_index, rsv_sample_score) in rsv.score.samples.iter().enumerate() {
            if rsv_sample_score.is_max_scoring_depth {
                continue;
            }

            // Determine if the SV has a non-reference GT in this sample:
            if rsv_sample_score
                .shared
                .gt
                .as_ref()
                .is_none_or(|x| *x == Genotype::Ref)
            {
                continue;
            }

            let sv_locus_expected_copy_number =
                rsv_sample_score.expected_cn_info.expected_copy_number;

            let sample_cn_segment_sv_match_tracker = &mut cn_segment_sv_match_tracker[sample_index];

            // Wind CNV index forward from head segment to find the first segment intersecting this SV
            let sample_head_cn_segment = &mut head_cn_segment[sample_index];
            assert!(chrom_index >= sample_head_cn_segment.chrom_index);
            if chrom_index > sample_head_cn_segment.chrom_index {
                sample_head_cn_segment.chrom_index = chrom_index;
                sample_head_cn_segment.segment_index = 0;
            }

            let sample_cn_segments = &all_sample_data[sample_index].copy_number_segments;
            let chrom_cn_segments = &sample_cn_segments.cn_segments[chrom_index];
            let start_seg_index = sample_head_cn_segment.segment_index;

            let mut first_cnv_intersection_found = false;

            for (cn_segment_index, cn_segment) in
                chrom_cn_segments.iter().enumerate().skip(start_seg_index)
            {
                let cn_segment_pos_range = cn_segment.to_range(bin_size);
                if cn_segment_pos_range.end <= be1.segment.range.start {
                    continue;
                } else if !first_cnv_intersection_found {
                    sample_head_cn_segment.segment_index = cn_segment_index;
                    first_cnv_intersection_found = true;
                }

                if cn_segment_pos_range.start > be2.segment.range.end {
                    assert!(first_cnv_intersection_found);

                    // Past the end of the target SV, stop iterating through CNV segments at this point:
                    break;
                }

                // Check if all criteria are met for the creation of a low-cost copy-number state transition on this
                // SV's boundaries

                // 1. Check that CN is a variant state consistent with the SV type and GT in this sample
                //
                let copy_number = cn_segment.copy_number_info.state as u8;

                if be1.dir == BreakendDirection::LeftAnchor {
                    // Continue unless this is a CN loss matching the SV GT value
                    //
                    match rsv_sample_score.shared.gt {
                        Some(Genotype::Het) => {
                            if copy_number != sv_locus_expected_copy_number.saturating_sub(1) {
                                continue;
                            }
                        }
                        Some(Genotype::Hom) => {
                            if copy_number != 0 {
                                continue;
                            }
                        }
                        _ => continue,
                    }
                } else {
                    // Continue unless this is a CN gain
                    //
                    // Translating a GT to a CN gain match is a little too unreliable, so leave this for the DEL section only
                    if copy_number <= sv_locus_expected_copy_number {
                        continue;
                    }
                }

                // 2. check that the next and previous cn segments have expected cn
                if cn_segment_index > 0 {
                    let prior_cn_segment = &chrom_cn_segments[cn_segment_index - 1];
                    if prior_cn_segment.copy_number_info.state as u8
                        != sv_locus_expected_copy_number
                    {
                        continue;
                    }
                }

                if cn_segment_index + 1 < chrom_cn_segments.len() {
                    let next_cn_segment = &chrom_cn_segments[cn_segment_index + 1];
                    if next_cn_segment.copy_number_info.state as u8 != sv_locus_expected_copy_number
                    {
                        continue;
                    }
                }

                // 3. check that the SV and CN spans have the required recip overlap
                let recip_overlap = get_recip_overlap(&cn_segment_pos_range, &sv_span);
                if recip_overlap < min_sv_cnv_recip_overlap {
                    continue;
                }

                let entry = sample_cn_segment_sv_match_tracker[chrom_index]
                    .entry(cn_segment_index)
                    .or_insert_with(SVMatchInfo::default);
                if entry.match_count == 0 {
                    entry.start_pos = be1.segment.range.center();
                    entry.end_pos = be2.segment.range.center();
                    entry.is_del = be1.dir == BreakendDirection::LeftAnchor;
                }
                entry.match_count += 1;

                // We're done marking this SV in this sample, so move to the next sample (or next eligible SV) without
                // scanning any more copynum segments
                break;
            }
        }
    }

    // Transfer qualifying SV/CN segment matches into modified transition type scores:
    for sample_index in 0..sample_count {
        let sample_transition_types = &mut transition_types[sample_index];
        for chrom_index in 0..chrom_count {
            for sv_match_info in cn_segment_sv_match_tracker[sample_index][chrom_index]
                .values()
                .filter(|x| x.match_count == 1)
            {
                if sv_match_info.is_del {
                    // DEL
                    sample_transition_types.set_loss_pos(chrom_index, sv_match_info.start_pos);
                    sample_transition_types.set_gain_pos(chrom_index, sv_match_info.end_pos);
                } else {
                    // DUP
                    sample_transition_types.set_gain_pos(chrom_index, sv_match_info.start_pos);
                    sample_transition_types.set_loss_pos(chrom_index, sv_match_info.end_pos);
                }
                if debug {
                    let label = if sv_match_info.is_del { "del" } else { "dup" };
                    eprintln!(
                        "Adding {label}-like-bnd close-fit segment annotation at: {chrom_index}:{}-{}",
                        sv_match_info.start_pos, sv_match_info.end_pos,
                    );
                }
            }
        }
    }
}

/// Identify all changes to copy number transition scoring based on SV breakends
///
/// Right now this covers 2 different SV breakpoint identification schemes:
///
/// 1. "Segmentation power boost" - For SVs with no CN segmentation within or close to the SV boundaries, add the
///    SV breakends as enhanced CN transition points. The idea here is to boost the power of segmentation for
///    SVs that might be too small to segment from the depth track alone.
///
/// 2. "Close fit" - When there is a one-to-one match of an SV and CNV segment, but the spans of the two events
///    don't quite match up very well, add the SV breakends as enhanced CN transition points to see if the CN
///    segmentation to the SV will significantly improve with only a small nudge.
///
/// Return value is a vector with one SampleTransitionTypes for each sample, these contain information on which
/// bins should have boosted transition scores.
///
pub fn get_copynum_transition_types_from_breakends(
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    sv_groups: &mut [SVGroup],
) -> Vec<SampleTransitionTypes> {
    let debug = false;

    // Start by setting up transition_types data structure, which will contain all of the information produced by this
    // function to suggested alternate transition weights during resegmentation
    //
    let sample_count = all_sample_data.len();
    assert!(sample_count > 0);

    let first_sample_data = &all_sample_data[0];
    let bin_size = first_sample_data.copy_number_segments.bin_size;

    let mut transition_types =
        vec![SampleTransitionTypes::new(&shared_data.chrom_list, bin_size); sample_count];

    // Second step is to sort the sv_groups on their first sv's first breakend, this order is assumed in
    // the different annotation schemes below
    //
    sv_groups.sort_by_key(|x| x.refined_svs[0].bp.breakend1.segment.clone());

    // Parameters hard-coded here for now

    // Copy number changing SVs of this size or larger need to be evaluated against the depth track to determine
    // whether we retain their SV type or change them to BNDs:
    let min_copynum_checked_sv_size = 50000;

    annotate_segmentation_power_boost_breakends(
        all_sample_data,
        sv_groups,
        &mut transition_types,
        min_copynum_checked_sv_size,
        debug,
    );

    annotate_close_fit_breakends(
        all_sample_data,
        sv_groups,
        &mut transition_types,
        min_copynum_checked_sv_size,
        debug,
    );

    transition_types
}
