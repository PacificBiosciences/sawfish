use super::merge_sv_shared::{
    AnnoSVGroup, ClusterMergeMap, MultiSampleClusterId, clone_sv_group_as_multisample,
    set_sv_group_expected_cn_info,
};
use super::{CandidateSVGroupInfo, get_duplicate_stats};
use crate::bio_align_utils::{PairwiseAligner, print_pairwise_alignment};
use crate::joint_call::SampleJointCallData;
use crate::run_stats::MultiSampleCategoryMergeStats;
use crate::sv_group::SVGroup;
use crate::utils::print_fasta;

fn test_for_duplicate_haplotypes_with_alignment(
    sv_group1: &SVGroup,
    sv_group2: &SVGroup,
    debug_duptest: bool,
) -> bool {
    let min_merge_norm_score = 0.97;

    // All SV contig comparisons are normalized to be in the breakend1 frame of reference
    //

    // Determine which contig has a longer HQ range:
    let hq_range1 = {
        let sv_hap_index = sv_group1.sv_haplotype_map[0];
        let hap = &sv_group1.group_haplotypes[sv_hap_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        &breakend1_contig_alignment.high_quality_contig_range
    };
    let hq_range2 = {
        let sv_hap_index = sv_group2.sv_haplotype_map[0];
        let hap = &sv_group2.group_haplotypes[sv_hap_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        &breakend1_contig_alignment.high_quality_contig_range
    };

    let (base_group, query_group) = if hq_range2.size() > hq_range1.size() {
        (sv_group2, sv_group1)
    } else {
        (sv_group1, sv_group2)
    };

    let base_contig = {
        let sv_hap_index = base_group.sv_haplotype_map[0];
        let hap = &base_group.group_haplotypes[sv_hap_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        breakend1_contig_alignment.contig_seq.as_slice()
    };

    let query_contig = {
        let sv_hap_index = query_group.sv_haplotype_map[0];
        let hap = &query_group.group_haplotypes[sv_hap_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        let hq_range = &breakend1_contig_alignment.high_quality_contig_range;
        &breakend1_contig_alignment.contig_seq[hq_range.start as usize..hq_range.end as usize]
    };

    let mut aligner = PairwiseAligner::new_merge_aligner();
    let score = aligner.align(base_contig, query_contig);

    // normalize the score by the query_contig, which should be globally aligned.
    let norm_score = score as f64 / query_contig.len() as f64;

    if debug_duptest {
        eprintln!(
            "Multi-Merge: Contig base/query length: {}/{} norm_score: {norm_score}",
            base_contig.len(),
            query_contig.len()
        );

        print_fasta(150, &[base_contig, query_contig]);
        eprintln!("query to base alignment:");
        print_pairwise_alignment(&aligner, base_contig, query_contig);
    }
    norm_score >= min_merge_norm_score
}

/// Return true if the haplotypes from sv groups 1 and 2 should be treated as duplicates
///
fn test_for_duplicate_haplotypes(sv_group1: &SVGroup, sv_group2: &SVGroup) -> bool {
    let debug_duptest = false;

    // Use bp match as a quick shortcut test first:
    //
    // This relies on multi-region sv groups having exactly one refined_sv
    //
    // TODO: Note that extended information in the breakpoint, like neighbors can prevent this type
    // of equality, is that working the way we intend now?
    //
    if sv_group1.refined_svs[0].bp == sv_group2.refined_svs[0].bp {
        return true;
    }

    // Finally test for haplotype similarity via alignment
    test_for_duplicate_haplotypes_with_alignment(sv_group1, sv_group2, debug_duptest)
}

/// Merge candidates from the duplicate candidate pool
///
/// 1. Determine how many SVGroups to keep out of the pool members
/// 2. Reformat the single-sample SVGroup into multi-sample format
///
/// Returns a 2-tuple of:
/// 1. the deduplicated sv groups
/// 2. a map of deduplicated cluster ids to the cluster they were merged into
///
fn merge_duplicated_candidates(
    all_sample_data: &[SampleJointCallData],
    duplicate_candidate_pool: &Vec<CandidateSVGroupInfo>,
    debug: bool,
) -> (Vec<SVGroup>, ClusterMergeMap) {
    assert!(!duplicate_candidate_pool.is_empty());

    if debug {
        eprintln!("Multi-Merge: Processing duplicate candidate pool: {duplicate_candidate_pool:?}");
    }

    let sv_groups = duplicate_candidate_pool
        .iter()
        .map(|x| AnnoSVGroup {
            sv_group: &all_sample_data[x.sample_index].sv_groups[x.sv_group_index],
            sample_index: x.sample_index,
        })
        .collect::<Vec<_>>();

    if debug {
        eprintln!(
            "Multi-Merge: processing pool of {} candidates.",
            sv_groups.len()
        );
        for (sv_group_index, sv_group) in sv_groups.iter().enumerate() {
            let rsv = &sv_group.sv_group.refined_svs[0];
            eprintln!(
                "Candidate {sv_group_index}: cluster_id: {} bp: {:?}",
                rsv.id.cluster_index, rsv.bp
            );
        }
    }

    let mut cluster_merge_map = ClusterMergeMap::default();

    let sample_count = all_sample_data.len();

    // Step 1: Handle the single SV candidate case, which doesn't need any merging:
    if sv_groups.len() == 1 {
        // Reformat the only sv group for multi-sample output.
        let multisample_sv_group = clone_sv_group_as_multisample(sample_count, &sv_groups[0]);
        return (vec![multisample_sv_group], cluster_merge_map);
    }

    // Step 2: Run simple pairwise search to create similar sv group pools, then sort the best
    // sv group representative first in each pool
    //
    let sv_group_merge_pools = {
        let pre_merge_sv_group_count = sv_groups.len();

        // SV group pools are initialized with one sv group per pool
        let mut sv_group_merge_pools = (0..pre_merge_sv_group_count)
            .map(|x| vec![x])
            .collect::<Vec<_>>();

        for sv_group_index1 in 0..pre_merge_sv_group_count {
            if sv_group_merge_pools[sv_group_index1].is_empty() {
                continue;
            }
            for sv_group_index2 in (sv_group_index1 + 1)..pre_merge_sv_group_count {
                if sv_group_merge_pools[sv_group_index2].is_empty() {
                    continue;
                }

                if test_for_duplicate_haplotypes(
                    sv_groups[sv_group_index1].sv_group,
                    sv_groups[sv_group_index2].sv_group,
                ) {
                    if debug {
                        eprintln!("Multi-Merge: merging sv groups 1 and 2");
                        let rsv1 = &sv_groups[sv_group_index1].sv_group.refined_svs[0];
                        eprintln!(
                            "sv group 1: cluster_id: {} bp: {:?}",
                            rsv1.id.cluster_index, rsv1.bp
                        );
                        let rsv2 = &sv_groups[sv_group_index2].sv_group.refined_svs[0];
                        eprintln!(
                            "sv group 2: cluster_id: {} bp: {:?}",
                            rsv2.id.cluster_index, rsv2.bp
                        );
                    }
                    let mut tmp = Vec::new();
                    std::mem::swap(&mut sv_group_merge_pools[sv_group_index2], &mut tmp);
                    sv_group_merge_pools[sv_group_index1].append(&mut tmp);
                }
            }
        }

        // Now sort the pools so that the best representative is ordered first
        for pool in sv_group_merge_pools.iter_mut().filter(|x| x.len() > 1) {
            // Select the representative haplotype as the one with:
            // 1. Highest read support
            // 2. Longest core contig
            //
            // This is setup assuming that multi-region sv groups contain exactly one SV
            //
            pool.sort_by_key(|&sv_group_index| {
                let sv_group = sv_groups[sv_group_index].sv_group;
                let sv_hap_index = sv_group.sv_haplotype_map[0];
                let hap = &sv_groups[sv_group_index].sv_group.group_haplotypes[sv_hap_index];
                let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
                std::cmp::Reverse((
                    hap.contig_info.supporting_read_count,
                    breakend1_contig_alignment.contig_seq.len(),
                ))
            });
        }
        sv_group_merge_pools
    };

    // Step 3: Consolidate sv group pools
    let mut output_sv_groups = Vec::new();
    for pool in sv_group_merge_pools.into_iter().filter(|x| !x.is_empty()) {
        let mut base_multisample_sv_group =
            clone_sv_group_as_multisample(sample_count, &sv_groups[pool[0]]);

        // Merge other pool members into base:
        let base_sample_index = sv_groups[pool[0]].sample_index;
        let base_sample_hap_map =
            base_multisample_sv_group.sample_haplotype_list[base_sample_index].clone();
        let base_cluster_index = base_multisample_sv_group.refined_svs[0].id.cluster_index;

        for anno_sv_group in pool.into_iter().skip(1).map(|x| &sv_groups[x]) {
            let sample_index = anno_sv_group.sample_index;
            let rsv = &anno_sv_group.sv_group.refined_svs[0];
            let cluster_index = rsv.id.cluster_index;

            // Note that it turns out we can have repeated multi-region SVs from the same
            // sample. Although the starting target regions for the SVs are prevented from
            // duplicating, the assemblies from 2 nearby cases can end up forming the same
            // haplotype and SV. This is a rare case so just throw away the second copy below
            // and don't assert.
            //
            // assert_ne!(base_sample_index, sample_index);
            base_multisample_sv_group.sample_haplotype_list[sample_index] =
                base_sample_hap_map.clone();

            let from_cluster = MultiSampleClusterId {
                sample_index,
                cluster_index,
            };
            let to_cluster = MultiSampleClusterId {
                sample_index: base_sample_index,
                cluster_index: base_cluster_index,
            };

            let base_bp = &mut base_multisample_sv_group.refined_svs[0].bp;
            if base_bp.breakend1_neighbor.is_none() {
                base_bp.breakend1_neighbor = rsv.bp.breakend1_neighbor.clone();
            }
            if base_bp.breakend2_neighbor.is_none() {
                base_bp.breakend2_neighbor = rsv.bp.breakend2_neighbor.clone();
            }

            cluster_merge_map.data.insert(from_cluster, to_cluster);
        }
        output_sv_groups.push(base_multisample_sv_group);
    }

    if debug {
        eprintln!(
            "Multi-Merge: output_sv_group count: {}",
            output_sv_groups.len()
        );
    }

    (output_sv_groups, cluster_merge_map)
}

/// Return de-duplicated multi-sample SV groups
///
/// Return a 3-tuple of:
/// 1. A vector of ordered, de-duplicated candidate SV groups to evaluate over all samples
///    The groups in this case all have exactly 1 SV, because these are 2 region candidates we
///    don't deal with overlap processing.
/// 2. A map of how input clusters were consolidated to output clusters
/// 3. Multi-region SV merge stats
///
pub fn process_duplicate_candidate_pool(
    treat_single_copy_as_haploid: bool,
    all_sample_data: &[SampleJointCallData],
    duplicate_candidate_pool: &Vec<CandidateSVGroupInfo>,
    debug: bool,
) -> (Vec<SVGroup>, ClusterMergeMap, MultiSampleCategoryMergeStats) {
    assert!(!duplicate_candidate_pool.is_empty());

    let (mut sv_groups, cluster_merge_map) =
        merge_duplicated_candidates(all_sample_data, duplicate_candidate_pool, debug);

    set_sv_group_expected_cn_info(
        treat_single_copy_as_haploid,
        all_sample_data,
        &mut sv_groups,
    );

    // Reformat SV group info for stats routine:
    let input_sv_groups = duplicate_candidate_pool
        .iter()
        .map(|x| &all_sample_data[x.sample_index].sv_groups[x.sv_group_index])
        .collect::<Vec<_>>();
    let output_sv_groups = sv_groups.iter().collect::<Vec<_>>();

    let duplicate_stats = get_duplicate_stats(&input_sv_groups, &output_sv_groups);

    (sv_groups, cluster_merge_map, duplicate_stats)
}
