use std::collections::BinaryHeap;

use rust_vc_utils::GenomeSegment;

use super::CandidateSVGroupInfo;
use super::merge_sv_shared::store_duplicate_pool;
use crate::joint_call::SampleJointCallData;
use crate::sv_group::SVGroup;

/// The primary payload and sort key in the heap is segment.
///
/// Store insertion_index in the HeapKey to add deterministic tie-breaking
///
/// Store sample_index and sv_group_index in the HeapKey but don't sort on them
#[derive(Debug, Eq)]
struct HeapKey {
    segment: GenomeSegment,

    /// Stores the order of heap insertions, so that otherwise equal segment candidates are
    /// pop'd first-in first-out.
    insertion_index: usize,
    sample_index: usize,
    sv_group_index: usize,
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

/// From an sv_group, return the region that will be used for sorting and merging single-region SVs across samples
///
/// The first group region is close to what we want for this purpose, but here we extend it with the actual first SV
/// candidate breakpoint. This helps sanitize merging behavior when the first breakpoint is discovered outside of
/// the target region.
///
pub fn get_sv_group_merge_region(sv_group: &SVGroup) -> GenomeSegment {
    assert!(!sv_group.group_regions.is_empty());
    assert!(!sv_group.refined_svs.is_empty());
    assert_eq!(
        sv_group.group_regions[0].chrom_index,
        sv_group.refined_svs[0].bp.breakend1.segment.chrom_index
    );
    let mut segment = sv_group.group_regions[0].clone();
    segment
        .range
        .merge(&sv_group.refined_svs[0].bp.breakend1.segment.range);
    segment
}

/// Add next entry from the indexed sample to the heap
///
/// # Parameters
/// * `all_sample_sv_group_indexes` - Ordered list of sv group indexes to select from
/// * `sample_index` - Index of the sample to add onto the heap
/// * `insertion_index`` - Counter of heap insertions used to make the heap pop order deterministic
/// * `all_sample_sv_group_indexes_head` - For each sample, track the current 'head' svgroup index
///
fn add_next_sample_entry_to_heap(
    all_sample_data: &[SampleJointCallData],
    all_sample_sv_group_indexes: &[Vec<usize>],
    sample_index: usize,
    insertion_index: &mut usize,
    all_sample_sv_group_indexes_head: &mut [usize],
    heap: &mut BinaryHeap<HeapKey>,
) {
    let sample_data = &all_sample_data[sample_index];
    let sample_sv_group_indexes = &all_sample_sv_group_indexes[sample_index];
    let sample_sv_group_indexes_head = &mut all_sample_sv_group_indexes_head[sample_index];

    if *sample_sv_group_indexes_head < sample_sv_group_indexes.len() {
        let sv_group_index = sample_sv_group_indexes[*sample_sv_group_indexes_head];
        let sv_group = &sample_data.sv_groups[sv_group_index];

        // SV group index selection should already be restricted to single-region variants:
        assert!(sv_group.is_single_region());

        let segment = get_sv_group_merge_region(sv_group);
        heap.push(HeapKey {
            segment,
            insertion_index: *insertion_index,
            sample_index,
            sv_group_index,
        });
        *insertion_index += 1;
        *sample_sv_group_indexes_head += 1;
    }
}

/// Find the subset of eligible single-region SV groups, and their sort order for the purpose of merging,
/// report this as a list of non-repeating sv_group indices per sample.
fn get_sv_group_index_order(all_sample_data: &[SampleJointCallData]) -> Vec<Vec<usize>> {
    let mut x = Vec::new();
    for sample_data in all_sample_data.iter() {
        // First restrict to just single-region sv_group indexes:
        let mut sv_group_indexes = sample_data
            .sv_groups
            .iter()
            .enumerate()
            .filter(|(_, x)| x.is_single_region())
            .map(|(i, _)| i)
            .collect::<Vec<_>>();

        // Second, sort indexes on the same crieria used for the heap in the pooling procedure below:
        sv_group_indexes.sort_by(|a, b| {
            // Add second comparison to make the sort stable:
            let region_a = get_sv_group_merge_region(&sample_data.sv_groups[*a]);
            let region_b = get_sv_group_merge_region(&sample_data.sv_groups[*b]);
            region_a.cmp(&region_b).then(a.cmp(b))
        });

        // eprintln!("AFTER");
        // for x in sv_group_indexes.iter() {
        //     let region = get_sv_group_merge_region(&sample_data.sv_groups[*x]);
        //     eprintln!("x: {} region: {region:?}", *x);
        // }

        x.push(sv_group_indexes);
    }
    x
}

/// Cluster single-sample single-region SVs into overlapping groups to be processed into deduplicated multi-sample SV groups
///
/// Return a vector of overlapping sv_groups
///
pub fn get_duplicate_candidate_pools(
    all_sample_data: &[SampleJointCallData],
    debug: bool,
) -> Vec<Vec<CandidateSVGroupInfo>> {
    let mut duplicate_candidate_pools = Vec::new();

    let sample_count = all_sample_data.len();
    let all_sample_sv_group_indexes = get_sv_group_index_order(all_sample_data);

    // Setup heap merge data structures
    let mut all_sample_sv_group_indexes_head = vec![0; sample_count];
    let mut heap = BinaryHeap::new();
    let mut insertion_index = 0;

    // Populate heap with lowest position SVGroup from each sample to initialize it:
    for sample_index in 0..sample_count {
        add_next_sample_entry_to_heap(
            all_sample_data,
            &all_sample_sv_group_indexes,
            sample_index,
            &mut insertion_index,
            &mut all_sample_sv_group_indexes_head,
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
        sv_group_index,
        ..
    }) = heap.pop()
    {
        if debug {
            eprintln!(
                "Single-Merge: Popped sample {sample_index} sv_group_index {sv_group_index} segment {segment:?}"
            );
        }

        // Use simple duplicate pool criteria, any segment which intersects the first segment
        // in the pool is included.
        //
        let store_current_duplicate_pool = if let Some(last_segment) = &last_segment {
            !last_segment.intersect(&segment)
        } else {
            true
        };

        if store_current_duplicate_pool {
            store_duplicate_pool(
                &mut duplicate_candidate_pool,
                &mut duplicate_candidate_pools,
                debug,
            );
            last_segment = Some(segment);
        }
        duplicate_candidate_pool.push(CandidateSVGroupInfo {
            sample_index,
            sv_group_index,
        });

        add_next_sample_entry_to_heap(
            all_sample_data,
            &all_sample_sv_group_indexes,
            sample_index,
            &mut insertion_index,
            &mut all_sample_sv_group_indexes_head,
            &mut heap,
        );
    }

    store_duplicate_pool(
        &mut duplicate_candidate_pool,
        &mut duplicate_candidate_pools,
        debug,
    );

    duplicate_candidate_pools
}
