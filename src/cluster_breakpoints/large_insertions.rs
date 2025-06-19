use log::info;
use rust_vc_utils::ChromList;

use crate::breakpoint::{BreakendDirection, Breakpoint, BreakpointCluster, LargeInsertionInfo};
use crate::genome_segment::{GenomeSegment, IntRange, is_strict_order};
use crate::utils::drop_true;

/// Temp structure used to assemble the large insertions in place
struct LargeInsertionAssemblyData {
    segment: GenomeSegment,
    root_cluster_index: usize,
}

/// Return large insertion info for each anchor pairing
fn get_large_insertion_info(
    cluster_index: usize,
    bp: &Breakpoint,
    last_anchor_data: Option<&LargeInsertionAssemblyData>,
    max_pairing_dist: i64,
) -> Option<LargeInsertionInfo> {
    if let Some(large_insertion_data) = last_anchor_data {
        let this_seg = &bp.breakend1.segment;
        let root_seg = &large_insertion_data.segment;
        let expanded_root_seg = IntRange::from_pair(
            root_seg.range.start,
            root_seg.range.start + max_pairing_dist,
        );
        if root_seg.chrom_index == this_seg.chrom_index
            && expanded_root_seg.intersect_pos(this_seg.range.start)
        {
            let mut segment = large_insertion_data.segment.clone();
            segment.range.end = this_seg.range.end;
            Some(LargeInsertionInfo {
                segment,
                root_cluster_index: large_insertion_data.root_cluster_index,
                paired_cluster_index: cluster_index,
            })
        } else {
            None
        }
    } else {
        None
    }
}

/// Find soft-clipping clusters that should be candidates for the dedicated large-insertion calling
/// pathway
///
/// Use a simple pairing strategy to join candidate left and right breakends together when they form
/// a large-indel like-pattern.
///
fn discover_large_insertion_candidates(
    max_large_insertion_pairing_dist: i64,
    max_large_insertion_deletion_pairing_dist: i64,
    clusters: &[BreakpointCluster],
) -> Vec<LargeInsertionInfo> {
    let mut large_insertion_candidates = Vec::new();

    let mut last_left_anchor: Option<LargeInsertionAssemblyData> = None;
    let mut last_right_anchor: Option<LargeInsertionAssemblyData> = None;

    // Select for softclip clusters (breakend2 is undefined in softclip clusters)
    for (cluster_index, bpc) in clusters
        .iter()
        .enumerate()
        .filter(|(_, x)| x.breakpoint.breakend2.is_none())
    {
        let bp = &bpc.breakpoint;

        let (self_anchor, other_anchor, max_dist) =
            if bp.breakend1.dir == BreakendDirection::RightAnchor {
                (
                    &mut last_right_anchor,
                    &mut last_left_anchor,
                    max_large_insertion_deletion_pairing_dist,
                )
            } else {
                (
                    &mut last_left_anchor,
                    &mut last_right_anchor,
                    max_large_insertion_pairing_dist,
                )
            };

        let large_insertion_info =
            get_large_insertion_info(cluster_index, bp, other_anchor.as_ref(), max_dist);
        *other_anchor = None;
        if let Some(large_insertion_info) = large_insertion_info {
            large_insertion_candidates.push(large_insertion_info);
        } else {
            *self_anchor = Some(LargeInsertionAssemblyData {
                segment: bp.breakend1.segment.clone(),
                root_cluster_index: cluster_index,
            });
        }
    }

    large_insertion_candidates
}

/// Check through large insertion candidates and eliminate any that overlap with conventional
/// indel candidates
fn filter_large_insertion_candidates_intersecting_indels(
    chrom_list: &ChromList,
    assembly_read_flank_size: usize,
    mut large_insertion_candidates: Vec<LargeInsertionInfo>,
    clusters: &[BreakpointCluster],
) -> Vec<LargeInsertionInfo> {
    let mut large_insertion_head = 0;
    let mut filter_lii = vec![false; large_insertion_candidates.len()];

    // Select for non-softclip clusters (breakend2 is defined in non-softclip clusters)
    for bpc in clusters.iter().filter(|x| x.breakpoint.breakend2.is_some()) {
        let bp = &bpc.breakpoint;
        let is_indel_like = bpc.consolidated_assembly_segment.is_some()
            || (bp.is_indel() && bpc.assembly_segments.len() == 1);

        if !is_indel_like {
            continue;
        }

        let mut segment = if let Some(consolidated_region) = &bpc.consolidated_assembly_segment {
            consolidated_region.segment.clone()
        } else {
            bpc.assembly_segments[0].clone()
        };

        segment.expand_by(chrom_list, assembly_read_flank_size as i64);

        // Iterate through lii until we encounter any cases that could intersect segment:
        while large_insertion_head < large_insertion_candidates.len() {
            let lii_segment = &large_insertion_candidates[large_insertion_head].segment;
            if is_strict_order(&segment, lii_segment) {
                break;
            } else if !is_strict_order(lii_segment, &segment) {
                // The two segments must intersect if we've gotten here.
                //
                filter_lii[large_insertion_head] = true;
            }

            large_insertion_head += 1;
        }
    }

    drop_true(&mut large_insertion_candidates, &filter_lii);
    large_insertion_candidates
}

/// Modify breakpoint clusters with large insertion candidate information
///
/// Returns the total number of large insertion candidates
///
fn add_large_insertions_to_clusters(
    large_insertion_candidates: Vec<LargeInsertionInfo>,
    clusters: &mut [BreakpointCluster],
) -> usize {
    if large_insertion_candidates.is_empty() {
        return 0;
    }

    let mut large_insertion_head_index = 0;
    let mut root_cluster_index =
        large_insertion_candidates[large_insertion_head_index].root_cluster_index;
    for (cluster_index, bpc) in clusters.iter_mut().enumerate() {
        if cluster_index < root_cluster_index {
            continue;
        }
        assert_eq!(cluster_index, root_cluster_index);

        bpc.large_insertion_info =
            Some(large_insertion_candidates[large_insertion_head_index].clone());
        large_insertion_head_index += 1;
        if large_insertion_head_index >= large_insertion_candidates.len() {
            break;
        }
        root_cluster_index =
            large_insertion_candidates[large_insertion_head_index].root_cluster_index;
    }

    large_insertion_head_index
}

/// Find candidate regions for the dedicated large-insertion pathway
///
/// Modify clusters to incorporate large insertion candidate data
///
/// Returns the total number of large insertion candidates
///
pub fn find_large_insertion_candidates(
    chrom_list: &ChromList,
    assembly_read_flank_size: usize,
    max_large_insertion_pairing_dist: i64,
    max_large_insertion_deletion_pairing_dist: i64,
    clusters: &mut [BreakpointCluster],
) -> usize {
    info!("Finding large insertion candidates");

    let large_insertion_candidates = discover_large_insertion_candidates(
        max_large_insertion_pairing_dist,
        max_large_insertion_deletion_pairing_dist,
        clusters,
    );

    // Filter out large insertion candidates which intersect regular indel clusters:
    let large_insertion_candidates = filter_large_insertion_candidates_intersecting_indels(
        chrom_list,
        assembly_read_flank_size,
        large_insertion_candidates,
        clusters,
    );

    add_large_insertions_to_clusters(large_insertion_candidates, clusters)
}
