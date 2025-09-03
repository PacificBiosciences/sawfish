mod merge_sv_shared;
mod multi_region;
mod single_region;

use log::info;
use std::sync::mpsc::channel;

use self::merge_sv_shared::{ClusterMergeMap, MultiSampleClusterId};
use crate::breakpoint::BreakendNeighbor;
use crate::cli::{JointCallSettings, SharedSettings};
use crate::joint_call::SampleJointCallData;
use crate::run_stats::{MultiSampleCategoryMergeStats, MultiSampleMergeStats};
use crate::sv_group::SVGroup;

use super::read_sample_data::SharedJointCallData;

#[derive(Clone, Debug)]
struct CandidateSVGroupInfo {
    sample_index: usize,
    sv_group_index: usize,
}

/// Get duplicate stats after deduplicating an sv_group pool from the given input to output state
///
fn get_duplicate_stats(
    input_sv_groups: &[&SVGroup],
    output_sv_groups: &[&SVGroup],
) -> MultiSampleCategoryMergeStats {
    let input_sv_group_count = input_sv_groups.len();
    let input_hap_count = input_sv_groups
        .iter()
        .map(|x| x.group_haplotypes.len())
        .sum::<usize>();
    let input_sv_count = input_sv_groups
        .iter()
        .map(|x| x.refined_svs.len())
        .sum::<usize>();

    let output_sv_group_count = output_sv_groups.len();
    let output_hap_count = output_sv_groups
        .iter()
        .map(|x| x.group_haplotypes.len())
        .sum::<usize>();
    let output_sv_count = output_sv_groups
        .iter()
        .map(|x| x.refined_svs.len())
        .sum::<usize>();

    let mut duplicate_stats = MultiSampleCategoryMergeStats::default();
    duplicate_stats.sv_groups_merged_as_duplicates += input_sv_group_count - output_sv_group_count;
    duplicate_stats.sv_haplotypes_merged_as_duplicates += input_hap_count - output_hap_count;
    duplicate_stats.sv_candidates_merged_as_duplicates += input_sv_count - output_sv_count;
    duplicate_stats.sv_groups_after_merge += output_sv_group_count;
    duplicate_stats.sv_haplotypes_after_merge += output_hap_count;
    duplicate_stats.sv_candidates_after_merge += output_sv_count;
    duplicate_stats
}

/// Identify and merge similar haplotypes
///
/// This is primarily intended to merge haplotypes accross samples. We've found some cases where merging within samples can
/// help, but ideally this would be moved back to the discover step instead.
///
pub fn merge_haplotypes(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
) -> (Vec<SVGroup>, MultiSampleMergeStats) {
    info!("Merging SV haplotypes");
    let debug = false;

    let multi_region_duplicate_candidate_pools = multi_region::get_duplicate_candidate_pools(
        &shared_data.chrom_list,
        all_sample_data,
        debug,
    );
    let single_region_duplicate_candidate_pools =
        single_region::get_duplicate_candidate_pools(all_sample_data, debug);

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();

    let (tx, rx) = channel();

    worker_pool.scope(move |scope| {
        for duplicate_candidate_pool in multi_region_duplicate_candidate_pools {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let (sv_groups, cluster_merge_map, duplicate_stats) =
                    multi_region::process_duplicate_candidate_pool(
                        settings.treat_single_copy_as_haploid,
                        all_sample_data,
                        &duplicate_candidate_pool,
                        debug,
                    );
                tx.send((true, sv_groups, cluster_merge_map, duplicate_stats))
                    .unwrap();
            });
        }

        for duplicate_candidate_pool in single_region_duplicate_candidate_pools {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let (sv_groups, duplicate_stats) = single_region::process_duplicate_candidate_pool(
                    settings.treat_single_copy_as_haploid,
                    all_sample_data,
                    &duplicate_candidate_pool,
                    debug,
                );
                tx.send((
                    false,
                    sv_groups,
                    ClusterMergeMap::default(),
                    duplicate_stats,
                ))
                .unwrap();
            });
        }
    });

    let mut multi_region_merged_sv_groups = Vec::new();
    let mut multi_region_cluster_merge_map = ClusterMergeMap::default();
    let mut multi_region_duplicate_stats = MultiSampleCategoryMergeStats::default();
    let mut single_region_merged_sv_groups = Vec::new();
    let mut single_region_duplicate_stats = MultiSampleCategoryMergeStats::default();
    for (is_multi_region, sv_groups, cluster_merge_map, duplicate_stats) in rx {
        if is_multi_region {
            multi_region_cluster_merge_map.extend(cluster_merge_map);
            multi_region_duplicate_stats.merge(&duplicate_stats);
            multi_region_merged_sv_groups.extend(sv_groups);
        } else {
            single_region_duplicate_stats.merge(&duplicate_stats);
            single_region_merged_sv_groups.extend(sv_groups);
        }
    }

    // Final step for multi-region sv_groups is to update all breakend neighbor ids to their merged values:
    for x in multi_region_merged_sv_groups.iter_mut() {
        fn update_breakend_neighbor(
            breakend_neighbor: Option<&mut BreakendNeighbor>,
            cluster_merge_map: &ClusterMergeMap,
        ) {
            if let Some(x) = breakend_neighbor
                && let Some(y) = cluster_merge_map.data.get(&MultiSampleClusterId {
                    sample_index: x.sample_index,
                    cluster_index: x.cluster_index,
                })
            {
                x.sample_index = y.sample_index;
                x.cluster_index = y.cluster_index;
            }
        }

        let bp = &mut x.refined_svs[0].bp;
        update_breakend_neighbor(
            bp.breakend1_neighbor.as_mut(),
            &multi_region_cluster_merge_map,
        );
        update_breakend_neighbor(
            bp.breakend2_neighbor.as_mut(),
            &multi_region_cluster_merge_map,
        );
    }

    let mut merged_sv_groups = single_region_merged_sv_groups;
    merged_sv_groups.extend(multi_region_merged_sv_groups);

    if debug {
        eprintln!("Merged sv group count {}", merged_sv_groups.len());
    }

    let merge_stats = MultiSampleMergeStats {
        multi_region: multi_region_duplicate_stats,
        single_region: single_region_duplicate_stats,
    };

    info!("Finished merging SV haplotypes across samples");

    (merged_sv_groups, merge_stats)
}
