mod merge_sv_shared;
mod multi_region;
mod single_region;

use log::info;
use rust_vc_utils::ChromList;
use std::sync::mpsc::channel;

use crate::cli::SharedSettings;
use crate::joint_call::SampleJointCallData;
use crate::run_stats::{MultiSampleCategoryMergeStats, MultiSampleMergeStats};
use crate::sv_group::SVGroup;

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

/// Identify and merge similar haplotypes across samples
///
pub fn merge_haplotypes(
    shared_settings: &SharedSettings,
    chrom_list: &ChromList,
    all_sample_data: &[SampleJointCallData],
) -> (Vec<SVGroup>, MultiSampleMergeStats) {
    info!("Merging SV haplotypes across samples");
    let debug = false;

    let multi_region_duplicate_candidate_pools =
        multi_region::get_duplicate_candidate_pools(chrom_list, all_sample_data, debug);
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
                let (sv_groups, duplicate_stats) = multi_region::process_duplicate_candidate_pool(
                    all_sample_data,
                    &duplicate_candidate_pool,
                    debug,
                );
                tx.send((true, sv_groups, duplicate_stats)).unwrap();
            });
        }

        for duplicate_candidate_pool in single_region_duplicate_candidate_pools {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let (sv_groups, duplicate_stats) = single_region::process_duplicate_candidate_pool(
                    all_sample_data,
                    &duplicate_candidate_pool,
                    debug,
                );
                tx.send((false, sv_groups, duplicate_stats)).unwrap();
            });
        }
    });

    let mut merged_sv_groups = Vec::new();
    let mut multi_region_duplicate_stats = MultiSampleCategoryMergeStats::default();
    let mut single_region_duplicate_stats = MultiSampleCategoryMergeStats::default();
    for (is_multi_region, sv_groups, duplicate_stats) in rx {
        if is_multi_region {
            multi_region_duplicate_stats.merge(&duplicate_stats);
        } else {
            single_region_duplicate_stats.merge(&duplicate_stats);
        }
        merged_sv_groups.extend(sv_groups);
    }

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
