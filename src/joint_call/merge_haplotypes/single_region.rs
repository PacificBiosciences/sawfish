use std::collections::{BTreeSet, BinaryHeap};

use super::merge_sv_shared::{
    AnnoSVGroup, clone_sv_group_as_multisample, set_sv_group_expected_cn_info,
};
use super::{CandidateSVGroupInfo, get_duplicate_stats};
use crate::bio_align_utils::{PairwiseAligner, print_pairwise_alignment};
use crate::breakpoint::{Breakend, Breakpoint};
use crate::genome_segment::GenomeSegment;
use crate::joint_call::SampleJointCallData;
use crate::run_stats::MultiSampleCategoryMergeStats;
use crate::sv_group::SVGroup;
use crate::utils::{drop_true, print_fasta};

/// Merge input_group into base_group
///
/// This version doesn't do any duplicate haplotype consolidation, and just handles all the paperwork of merging
/// while keeping all the indexes valid.
///
fn merge_sv_group_into_base_group(
    input_sample_index: usize,
    input_group: &SVGroup,
    base_group: &mut SVGroup,
) {
    // Steps to merge one SV group onto another:
    // 0. Expand group_regions
    // 1. Concatenate the refined SVs (Relabel sample index in all refined SV id values)
    // 2. Concatenate the haplotypes
    // 3. Update sample_haplotype_list value for the input sv_group (remember to offset the
    // haplotype index values)
    // 4. extend sv_haplotype_map to cover new refined SVs (remember to offset the haplotype
    // index values)
    //

    // We'll need this to offset the imported haplotype index values:
    let start_haplotype_count = base_group.group_haplotypes.len();

    // 0. Expand group regions
    base_group.group_regions[0]
        .range
        .merge(&input_group.group_regions[0].range);

    // 1. Concat refined SVs
    base_group
        .refined_svs
        .extend(input_group.refined_svs.iter().map(|x| {
            let mut y = x.clone();
            y.id.sample_index = input_sample_index;
            y
        }));

    // 2. Copy input group haplotypes with updated sample_index:
    let input_group_haplotype_iter = input_group.group_haplotypes.iter().map(|x| {
        let mut y = x.clone();
        y.hap_id.sample_index = input_sample_index;
        y
    });

    // 3. Concat haplotypes
    base_group
        .group_haplotypes
        .extend(input_group_haplotype_iter);

    // 4. Update sample_haplotype_list

    // TODO bring this check back
    //if !base_group.sample_haplotype_list[input_sample_index].is_empty() {
    //    panic!("Sample conflict merging input_sample_index {input_sample_index} into sv_group: {:?}\ninput_sv_group: {:?}", base_group, input_group);
    //}
    base_group.sample_haplotype_list[input_sample_index] = input_group.sample_haplotype_list[0]
        .iter()
        .map(|x| x + start_haplotype_count)
        .collect();

    // 4. Update sv_haplotype_map
    base_group.sv_haplotype_map.extend(
        input_group
            .sv_haplotype_map
            .iter()
            .map(|x| x + start_haplotype_count),
    );
}

/// Return true if 2 breakend positions are within `max_breakend_distance`
///
fn are_breakends_close(be1: &Breakend, be2: &Breakend, max_breakend_distance: i64) -> bool {
    be1.segment.chrom_index == be2.segment.chrom_index
        && (be1.segment.range.start - be2.segment.range.start).abs() <= max_breakend_distance
}

/// Return true if 2 breakpoints have insertion sizes within `max_insert_size_diff`
///
/// * max_insert_size_diff - Max insertion size difference for 2 breakpoint insertions.
///   For insertion lengths a and b, the insertion size difference is 2*abs(a-b)/(a+b).
///
fn are_breakpoint_insertions_close(
    bp1: &Breakpoint,
    bp2: &Breakpoint,
    max_insert_size_diff: f64,
) -> bool {
    // This allows for small indel differences in non-indel breakends
    let non_indel_size = 100;

    let max_del_size = 500;

    let del_size1 = std::cmp::min(bp1.deletion_size().unwrap_or(non_indel_size), max_del_size);
    let del_size2 = std::cmp::min(bp2.deletion_size().unwrap_or(non_indel_size), max_del_size);
    let ins_size1 = bp1.insert_info.size();
    let ins_size2 = bp2.insert_info.size();
    let size1 = ins_size1 + del_size1;
    let size2 = ins_size2 + del_size2;
    if (size1 + size2) == 0 {
        return true;
    }
    let insert_size_diff =
        (2 * (size1 as i32 - size2 as i32).abs()) as f64 / (size1 + size2) as f64;
    insert_size_diff <= max_insert_size_diff
}

/// Return true if 2 breakpoints have insertion sizes within `max_insert_size_diff` and each breakend pair position is withi `max_breakend_distance`
///
fn are_breakpoints_close(
    bp1: &Breakpoint,
    bp2: &Breakpoint,
    max_breakend_distance: i64,
    max_insert_size_diff: f64,
) -> bool {
    if !are_breakends_close(&bp1.breakend1, &bp2.breakend1, max_breakend_distance) {
        return false;
    }

    if !are_breakends_close(
        bp1.breakend2.as_ref().unwrap(),
        bp2.breakend2.as_ref().unwrap(),
        max_breakend_distance,
    ) {
        return false;
    }

    are_breakpoint_insertions_close(bp1, bp2, max_insert_size_diff)
}

fn are_breakpoint_sets_very_close(
    haplotype1_breakpoints: &[&Breakpoint],
    haplotype2_breakpoints: &[&Breakpoint],
) -> bool {
    let max_breakend_distance = 10;
    let max_insert_size_diff = 0.01;

    haplotype1_breakpoints.len() == haplotype2_breakpoints.len()
        && haplotype1_breakpoints
            .iter()
            .zip(haplotype2_breakpoints.iter())
            .all(|(bp1, bp2)| {
                are_breakpoints_close(bp1, bp2, max_breakend_distance, max_insert_size_diff)
            })
}

/// Return true if the two sets of breakpoints are similar
///
/// Similarity is defined in the most minimal way to start -- we look for just one breakpoint
/// that is similar between sets 1 and 2
///
/// The breakpoint similarity criteria are that the two breakends are both in a similar area
/// and the insert length is similar, all as defined by parameters below.
///
/// This method assumes both breakpoint lists are sorted
///
#[allow(dead_code)]
fn are_similar_breakpoints_found(
    haplotype1_breakpoints: &Vec<&Breakpoint>,
    haplotype2_breakpoints: &Vec<&Breakpoint>,
) -> bool {
    let max_breakend_distance = 100;
    let max_insert_size_diff = 0.05;

    for bp1 in haplotype1_breakpoints {
        let bp1_be1_seg = &bp1.breakend1.segment;
        for bp2 in haplotype2_breakpoints {
            let bp2_be1_seg = &bp2.breakend1.segment;

            // Run the primary skip/stop filter based on breakend1 similarity:
            if bp2_be1_seg < bp1_be1_seg {
                if bp2_be1_seg.chrom_index != bp1_be1_seg.chrom_index
                    || ((bp1_be1_seg.range.start - bp2_be1_seg.range.start) > max_breakend_distance)
                {
                    continue;
                }
            } else {
                // Because the breakpoints are sorted we can break if bp2 is already
                // substantially greater than bp1
                if bp2_be1_seg.chrom_index != bp1_be1_seg.chrom_index
                    || ((bp2_be1_seg.range.start - bp1_be1_seg.range.start) > max_breakend_distance)
                {
                    break;
                }
            }

            // We've met breakend1 similarity so test breakend2 and insert size:
            if are_breakends_close(
                bp1.breakend2.as_ref().unwrap(),
                bp2.breakend2.as_ref().unwrap(),
                max_breakend_distance,
            ) && are_breakpoint_insertions_close(bp1, bp2, max_insert_size_diff)
            {
                return true;
            }
        }
    }
    false
}

fn test_for_duplicate_haplotypes_with_alignment(
    sv_group: &SVGroup,
    haplotype1_index: usize,
    haplotype2_index: usize,
    debug_duptest: bool,
) -> bool {
    let min_merge_norm_score = 0.97;

    // Determine which contig has a lorger HQ range:
    let hq_range1 = {
        let hap = &sv_group.group_haplotypes[haplotype1_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        &breakend1_contig_alignment.high_quality_contig_range
    };
    let hq_range2 = {
        let hap = &sv_group.group_haplotypes[haplotype2_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        &breakend1_contig_alignment.high_quality_contig_range
    };

    let (base_index, query_index) = if hq_range2.size() > hq_range1.size() {
        (haplotype2_index, haplotype1_index)
    } else {
        (haplotype1_index, haplotype2_index)
    };

    let base_contig = {
        let hap = &sv_group.group_haplotypes[base_index];
        let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
        breakend1_contig_alignment.contig_seq.as_slice()
    };

    let query_contig = {
        let hap = &sv_group.group_haplotypes[query_index];
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
            "Single-Merge: Contig base/query length: {}/{} score: {}",
            base_contig.len(),
            query_contig.len(),
            norm_score,
        );

        print_fasta(150, &[base_contig, query_contig]);
        eprintln!("query to base alignment:");
        print_pairwise_alignment(&aligner, base_contig, query_contig);
    }
    norm_score >= min_merge_norm_score
}

/// Return true if hap 1 and 2 should be treated as duplicates
///
fn test_for_duplicate_haplotypes(
    sv_group: &SVGroup,
    haplotype1_index: usize,
    haplotype1_breakpoints: &Vec<&Breakpoint>,
    haplotype2_index: usize,
    haplotype2_breakpoints: &Vec<&Breakpoint>,
) -> bool {
    let debug_duptest = false;
    if debug_duptest {
        eprintln!("Single-Merge: hap1_index: {haplotype1_index} bps: {haplotype1_breakpoints:?}");
        eprintln!("Single-Merge: hap2_index: {haplotype2_index} bps: {haplotype2_breakpoints:?}");
    }

    // Use bp match as a quick shortcut test
    if haplotype1_breakpoints == haplotype2_breakpoints {
        return true;
    }

    // Next check if all breakpoints are very close, if so merge the haplotypes as in the
    // exact match case above
    if are_breakpoint_sets_very_close(haplotype1_breakpoints, haplotype2_breakpoints) {
        return true;
    }

    /*
    Disabled this filtration step because it was preventing some true haplotype merges, and inflating false trio denovo rates, etc...
    It's possible that this was preventing some errors due to spurious bad haplotype merging, so we'll keep an eye out for
    this case and add a more focused filter if this case is found.

    // Next check if any of the breakpoints are roughly similar, if not don't bother with
    // any type of more expensive haplotype alignment below
    //
    if !are_similar_breakpoints_found(haplotype1_breakpoints, haplotype2_breakpoints) {
        if debug_duptest {
            eprintln!("Rejecting merge without alignment due to dissimilar breakpoints");
        }
        return false;
    }
    */

    // Finally test for haplotype similarity via alignment
    test_for_duplicate_haplotypes_with_alignment(
        sv_group,
        haplotype1_index,
        haplotype2_index,
        debug_duptest,
    )
}

/// For each haplotype index, provide a list of other haplotype indexes to merge with it
///
/// # Arguments
///
/// * `merged_sv_group` - all svs in a pool merged into a single sv group, but without yet having any
///   haplotypes consolidated
///
///
fn get_haplotype_merge_pools(
    merged_sv_group: &SVGroup,
    haplotype_to_sv_map: &[Vec<usize>],
    debug: bool,
) -> Vec<Vec<usize>> {
    let pre_merge_haplotype_count = merged_sv_group.group_haplotypes.len();

    // Haplotype pools are initialized with one haplotype per pool
    let mut haplotype_merge_pools = (0..pre_merge_haplotype_count)
        .map(|x| vec![x])
        .collect::<Vec<_>>();

    // Create a list of non-mergable haplotype indexes for each haplotype, all in pre-merge
    // haplotype index coordinates. The primary criteria for being non-mergable is that the
    // haplotypes occur in the same sample.
    //
    // Exclusion relationships are required to be symmetric
    //
    let mut excluded_haplotype_indexes = {
        let mut x = vec![BTreeSet::new(); pre_merge_haplotype_count];
        for sample_haplotype_list in &merged_sv_group.sample_haplotype_list {
            let sample_haplotype_count = sample_haplotype_list.len();
            if sample_haplotype_count > 1 {
                for &haplotype_index in sample_haplotype_list.iter() {
                    for &excluded_haplotype_index in sample_haplotype_list
                        .iter()
                        .filter(|&&x| x != haplotype_index)
                    {
                        x[haplotype_index].insert(excluded_haplotype_index);
                    }
                }
            }
        }
        x
    };

    for haplotype1_index in 0..pre_merge_haplotype_count {
        if haplotype_merge_pools[haplotype1_index].is_empty() {
            continue;
        }
        let haplotype1_breakpoints = haplotype_to_sv_map[haplotype1_index]
            .iter()
            .map(|&sv_index| &merged_sv_group.refined_svs[sv_index].bp)
            .collect::<Vec<_>>();
        if haplotype1_breakpoints.is_empty() {
            continue;
        }

        for haplotype2_index in (haplotype1_index + 1)..pre_merge_haplotype_count {
            if haplotype_merge_pools[haplotype2_index].is_empty() {
                continue;
            }
            let haplotype2_breakpoints = haplotype_to_sv_map[haplotype2_index]
                .iter()
                .map(|&sv_index| &merged_sv_group.refined_svs[sv_index].bp)
                .collect::<Vec<_>>();
            if haplotype2_breakpoints.is_empty() {
                continue;
            }

            // Test if haplotype2 is excluded from merging with the haplotype1 merge pool,
            // and skip consolidation test if so:
            //
            if excluded_haplotype_indexes[haplotype1_index].contains(&haplotype2_index) {
                continue;
            }

            if test_for_duplicate_haplotypes(
                merged_sv_group,
                haplotype1_index,
                &haplotype1_breakpoints,
                haplotype2_index,
                &haplotype2_breakpoints,
            ) {
                if debug {
                    eprintln!(
                        "Single-Merge: merging haplotype indexes: {haplotype1_index}/{haplotype2_index}"
                    );
                }

                // swaps here help to workaround the borrow checker:
                {
                    let mut tmp = Vec::new();
                    std::mem::swap(&mut haplotype_merge_pools[haplotype2_index], &mut tmp);
                    haplotype_merge_pools[haplotype1_index].append(&mut tmp);
                }

                if !excluded_haplotype_indexes[haplotype2_index].is_empty() {
                    let mut tmp = BTreeSet::new();
                    std::mem::swap(&mut excluded_haplotype_indexes[haplotype2_index], &mut tmp);
                    excluded_haplotype_indexes[haplotype1_index].append(&mut tmp);
                }
            }
        }
    }

    // Sort each pool so that the best representative haplotype is ordered first
    for pool in haplotype_merge_pools.iter_mut().filter(|x| x.len() > 1) {
        // Select the representative haplotype as the one with:
        // 1. Highest read support
        // 2. Longest core contig
        //
        pool.sort_by_key(|&haplotype_index| {
            let hap = &merged_sv_group.group_haplotypes[haplotype_index];
            let breakend1_contig_alignment = &hap.contig_info.contig_alignments[0];
            std::cmp::Reverse((
                hap.contig_info.supporting_read_count,
                breakend1_contig_alignment.contig_seq.len(),
            ))
        });
    }
    haplotype_merge_pools
}

/// Merge members of the duplicate candidate pool
///
/// Run more detailed tests on duplicate candidates to determine if they really will
/// be handled as duplicates, then merge the duplicates.
///
/// 1. Determine how many SVGroups to keep out of the pool members
/// 2. Reformat the single-sample SVGroup into multi-sample format
///
fn merge_duplicated_candidates(
    all_sample_data: &[SampleJointCallData],
    duplicate_candidate_pool: &Vec<CandidateSVGroupInfo>,
    debug: bool,
) -> Vec<SVGroup> {
    assert!(!duplicate_candidate_pool.is_empty());

    if debug {
        eprintln!(
            "Single-Merge: Processing duplicate candidate pool: {duplicate_candidate_pool:?}"
        );
    }

    let sv_groups = duplicate_candidate_pool
        .iter()
        .map(|x| AnnoSVGroup {
            sv_group: &all_sample_data[x.sample_index].sv_groups[x.sv_group_index],
            sample_index: x.sample_index,
        })
        .collect::<Vec<_>>();

    // Step 1: Arbitrarily pick the first SV group and format it for multi-sample output.
    //
    // Use this SV group as the base onto which any other SV groups will be merged.
    //
    let sample_count = all_sample_data.len();
    let mut merged_sv_group = clone_sv_group_as_multisample(sample_count, &sv_groups[0]);

    // Step 2: Handle the special N==1 case, where we don't have to worry about any merging:
    if sv_groups.len() == 1 {
        return vec![merged_sv_group];
    }

    // Step 3: Merge duplicate pool into one sv_group to start, then we'll merge haplotypes
    // within the unified sv_group in subsequent steps.
    //
    // This maybe isn't the best long-term strategy but for now this reduces complexity of
    // downstream steps, so start here and iterate.
    //
    for sv_group in sv_groups.iter().skip(1) {
        let input_sample_index = sv_group.sample_index;
        let input_group = sv_group.sv_group;
        //eprintln!("Merging input_sample_index: {input_sample_index} input_group: {:?}\n\tinto merged group: {:?}", input_group, merged_sv_group);
        merge_sv_group_into_base_group(input_sample_index, input_group, &mut merged_sv_group);
        //eprintln!("Finished merged group: {:?}", merged_sv_group);
    }

    if debug {
        eprintln!("Merged duplicate pool sv_group before consolidation: {merged_sv_group:?}");
    }

    // Step 4: Run simple pairwise search to create similar haplotype pools
    let haplotype_to_sv_map = merged_sv_group.haplotype_to_sv_map();
    let haplotype_merge_pools =
        get_haplotype_merge_pools(&merged_sv_group, &haplotype_to_sv_map, debug);

    // Step 5: Consolidate haplotype pools
    //
    // Consolidate with haplotype indexes locked, and then run a final step to eliminate certain
    // haplotype index numbers and reregister all data structures as we do. Same procedure to
    // update refined svs too.
    //
    let pre_merge_haplotype_count = merged_sv_group.group_haplotypes.len();

    // The remap list will ultimately be updated to reflect both:
    // (1) changes in haplotype index because a pool of similar haplotypes has been consolidated
    // down to a single representative copy.
    // (2) changes to haplotype index because a lower-index haplotype was removed due to the
    // above consolidation process.
    //
    let mut haplotype_remap_list = (0..pre_merge_haplotype_count).collect::<Vec<_>>();
    let mut haplotype_delete_list = vec![false; pre_merge_haplotype_count];
    let mut sv_delete_list = vec![false; merged_sv_group.refined_svs.len()];
    for pool in haplotype_merge_pools.into_iter().filter(|x| x.len() > 1) {
        let best_in_pool = pool[0];
        for haplotype_in_pool in pool.into_iter().skip(1) {
            // 1. Mark to replace all entries in sample_haplotype_list with 'best_in_pool'
            // 2. Mark haplotypes for deletion
            // 3. find all the refined SVs associated with this haplotype and mark then for deletion
            //
            haplotype_remap_list[haplotype_in_pool] = best_in_pool;
            haplotype_delete_list[haplotype_in_pool] = true;
            for delete_sv_index in haplotype_to_sv_map[haplotype_in_pool].iter() {
                sv_delete_list[*delete_sv_index] = true;
            }
        }
    }

    // Adjust haplotype_remap_list to account for deleted haplotypes
    {
        let mut shift = 0i32;
        let mut haplotype_delete_shift = vec![0; pre_merge_haplotype_count];
        for haplotype_index in 0..pre_merge_haplotype_count {
            if haplotype_delete_list[haplotype_index] {
                shift += 1;
                haplotype_delete_shift[haplotype_index] = -1;
            } else {
                haplotype_delete_shift[haplotype_index] = shift;
            }
        }

        for haplotype_index in 0..pre_merge_haplotype_count {
            let shift = haplotype_delete_shift[haplotype_remap_list[haplotype_index]];
            assert!(shift >= 0);
            assert!(haplotype_remap_list[haplotype_index] >= shift as usize);
            haplotype_remap_list[haplotype_index] -= shift as usize;
        }
    }

    // Remap haplotype ids
    {
        for sample_haps in merged_sv_group.sample_haplotype_list.iter_mut().flatten() {
            *sample_haps = haplotype_remap_list[*sample_haps];
        }
        for sv_haps in merged_sv_group.sv_haplotype_map.iter_mut() {
            *sv_haps = haplotype_remap_list[*sv_haps];
        }
    }

    // Drop removed haplotypes
    drop_true(
        &mut merged_sv_group.group_haplotypes,
        &haplotype_delete_list,
    );

    // Drop removed svs
    drop_true(&mut merged_sv_group.refined_svs, &sv_delete_list);

    // Drop sv_haplotype_map for removed svs
    drop_true(&mut merged_sv_group.sv_haplotype_map, &sv_delete_list);

    if debug {
        eprintln!("Merged duplicate pool sv_group after consolidation: {merged_sv_group:?}");
    }

    vec![merged_sv_group]
}

/// Merge SV group candidates across samples from the same locus
///
/// Once a set of SVGroups from multiple samples has been identified as a potentially overlapping group, this method
/// processes the candidates to identify duplicate haplotypes shared across samples and consolidates these down to
/// one or more multi-sample SVGroups ready to be processed for joint-genotyping downstream
///
/// # Arguments
///
/// * `treat_single_copy_as_haploid` - If true, treat expected single-copy regions as haploid
/// * `all_sample_data` - static collection of discovery data from all input samples
/// * `duplicate_candidate_pool` - The pool of candidate SV groups to process
/// * `debug` - If true, log some high-level process diagnostics
///
/// Return (1) the merged list of SV groups (2) duplicate merging stats
///
pub(super) fn process_duplicate_candidate_pool(
    treat_single_copy_as_haploid: bool,
    all_sample_data: &[SampleJointCallData],
    duplicate_candidate_pool: &Vec<CandidateSVGroupInfo>,
    debug: bool,
) -> (Vec<SVGroup>, MultiSampleCategoryMergeStats) {
    assert!(!duplicate_candidate_pool.is_empty());

    let mut sv_groups =
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

    (sv_groups, duplicate_stats)
}

/// The primary payload and sort key in the heap is segment.
///
/// Store insertion_index in the HeapKey to add deterministic tie-breaking
///
/// Store sample_index in the HeapKey but don't sort on it
#[derive(Eq)]
struct HeapKey {
    segment: GenomeSegment,

    /// Stores the order of heap insertions, so that otherwise equal segment candidates are
    /// pop'd first-in first-out.
    insertion_index: usize,
    sample_index: usize,
}

/// Sort the heap on the genome regions, since each single-region SV group will have only one
/// region it simplifies the meaning of the sort
///
/// Reverse the order of this sort so that we naturally pop segments off of the heap in
/// 'left-to-right' genomic order, and tie-break on low-to-high insertion order.
///
impl Ord for HeapKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .segment
            .cmp(&self.segment)
            .then(other.insertion_index.cmp(&self.insertion_index))
    }
}

impl PartialOrd for HeapKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapKey {
    fn eq(&self, other: &Self) -> bool {
        (self.segment == other.segment) && (self.insertion_index == other.insertion_index)
    }
}

fn add_next_sample_entry_to_heap(
    all_sample_data: &[SampleJointCallData],
    sample_svgroup_count: &[usize],
    sample_index: usize,
    insertion_index: &mut usize,
    sample_svgroup_head_indexes: &mut [usize],
    heap: &mut BinaryHeap<HeapKey>,
) {
    let sample_data = &all_sample_data[sample_index];
    while sample_svgroup_head_indexes[sample_index] < sample_svgroup_count[sample_index] {
        let sv_group = &sample_data.sv_groups[sample_svgroup_head_indexes[sample_index]];

        // Only consider single-region variants:
        if sv_group.is_single_region() {
            let segment = sv_group.group_regions[0].clone();
            heap.push(HeapKey {
                segment,
                insertion_index: *insertion_index,
                sample_index,
            });
            *insertion_index += 1;
            break;
        }
        sample_svgroup_head_indexes[sample_index] += 1;
    }
}

/// Cluster single-sample single-region SVs into overlapping groups to be processed into deduplicated multi-sample SV groups
///
/// Return a vector of overlapping sv_groups
///
pub(super) fn get_duplicate_candidate_pools(
    all_sample_data: &[SampleJointCallData],
    debug: bool,
) -> Vec<Vec<CandidateSVGroupInfo>> {
    let mut duplicate_candidate_pools = Vec::new();

    let sample_count = all_sample_data.len();
    let sample_svgroup_count = all_sample_data
        .iter()
        .map(|x| x.sv_groups.len())
        .collect::<Vec<_>>();

    // Setup heap merge data structures
    let mut sample_svgroup_head_indexes = vec![0; sample_count];
    let mut heap = BinaryHeap::new();
    let mut insertion_index = 0;

    // Populate heap with lowest position SVGroup from each sample to initialize it:
    for sample_index in 0..sample_count {
        add_next_sample_entry_to_heap(
            all_sample_data,
            &sample_svgroup_count,
            sample_index,
            &mut insertion_index,
            &mut sample_svgroup_head_indexes,
            &mut heap,
        );
    }

    // Collection of svgroups that should be evaluated to determine if any are duplicates
    //
    let mut duplicate_candidate_pool = Vec::new();

    // Scan through genome in breakpoint sort order, using the heap to merge breakpoints across all
    // samples in order, so that duplicates can be found:
    //
    let mut last_segment: Option<GenomeSegment> = None;
    while let Some(HeapKey {
        segment,
        sample_index,
        ..
    }) = heap.pop()
    {
        if debug {
            eprintln!("Single-Merge: Popped sample {sample_index} segment {segment:?}");
        }

        // Use simple duplicate pool criteria, any segment which intersects the first segment
        // in the pool is included.
        //
        let process_duplicate_pool = if let Some(last_segment) = &last_segment {
            !last_segment.intersect(&segment)
        } else {
            true
        };

        if process_duplicate_pool {
            if !duplicate_candidate_pool.is_empty() {
                duplicate_candidate_pools.push(duplicate_candidate_pool.clone());
            }
            duplicate_candidate_pool.clear();
            last_segment = Some(segment);
        }
        let sv_group_index = sample_svgroup_head_indexes[sample_index];
        duplicate_candidate_pool.push(CandidateSVGroupInfo {
            sample_index,
            sv_group_index,
        });

        sample_svgroup_head_indexes[sample_index] += 1;
        add_next_sample_entry_to_heap(
            all_sample_data,
            &sample_svgroup_count,
            sample_index,
            &mut insertion_index,
            &mut sample_svgroup_head_indexes,
            &mut heap,
        );
    }
    if !duplicate_candidate_pool.is_empty() {
        duplicate_candidate_pools.push(duplicate_candidate_pool.clone());
    }

    duplicate_candidate_pools
}
