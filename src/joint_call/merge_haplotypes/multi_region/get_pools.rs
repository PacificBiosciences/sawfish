use std::collections::BinaryHeap;

use rust_vc_utils::{ChromList, GenomeSegment};

use super::CandidateSVGroupInfo;
use super::merge_sv_shared::store_duplicate_pool;
use crate::breakpoint::Breakpoint;
use crate::joint_call::SampleJointCallData;

/// The primary payload and sort key in the heap is bp.
///
/// Store insertion_index in the HeapKey to add deterministic tie-breaking
///
/// Store sample_index and sv_group_index in the HeapKey but don't sort on them
#[derive(Debug, Eq)]
struct HeapKey {
    bp: Breakpoint,
    insertion_index: usize,
    sample_index: usize,
    sv_group_index: usize,
}

/// Sort the heap on the SV breakpoint, since each multi-region SV group will have only one SV
///
/// Reverse the order of this sort so that we naturally pop breakends off of the heap in
/// 'left-to-right' genomic order
impl Ord for HeapKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .bp
            .cmp(&self.bp)
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
        (self.bp == other.bp) && (self.insertion_index == other.insertion_index)
    }
}

#[allow(unused)]
fn debug_heap(heap: &BinaryHeap<HeapKey>) {
    eprintln!("Heap Status:");
    for x in heap.iter() {
        eprintln!("Heap key: {x:?}");
    }
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

        // SV group index selection should already be restricted to multi-region variants:
        assert!(!sv_group.is_single_region());

        let bp = sv_group.refined_svs[0].bp.clone();
        heap.push(HeapKey {
            bp,
            insertion_index: *insertion_index,
            sample_index,
            sv_group_index,
        });
        *insertion_index += 1;
        *sample_sv_group_indexes_head += 1;
    }
}

/// Find the subset of eligible multi-region SV groups, and their sort order for the purpose of merging,
/// report this as a list of non-repeating sv_group indices per sample.
fn get_sv_group_index_order(all_sample_data: &[SampleJointCallData]) -> Vec<Vec<usize>> {
    let mut x = Vec::new();
    for sample_data in all_sample_data.iter() {
        // First restrict to just multi-region sv_group indexes:
        let mut sv_group_indexes = sample_data
            .sv_groups
            .iter()
            .enumerate()
            .filter(|(_, x)| !x.is_single_region())
            .map(|(i, _)| i)
            .collect::<Vec<_>>();

        // Second, sort indexes on the same crieria used for the heap in the pooling procedure below:
        sv_group_indexes.sort_by(|a, b| {
            // Add second comparison to make the sort stable:
            let bp_a = &sample_data.sv_groups[*a].refined_svs[0].bp;
            let bp_b = &sample_data.sv_groups[*b].refined_svs[0].bp;
            bp_a.cmp(bp_b).then(a.cmp(b))
        });

        x.push(sv_group_indexes);
    }
    x
}

/// Cluster single-sample multi-region SVs into overlapping groups to be processed into deduplicated multi-sample SV groups
///
/// Return a vector of overlapping sv_groups
///
pub fn get_duplicate_candidate_pools(
    chrom_list: &ChromList,
    all_sample_data: &[SampleJointCallData],
    debug: bool,
) -> Vec<Vec<CandidateSVGroupInfo>> {
    let mut duplicate_candidate_pools = Vec::new();

    // Clustering strategy:
    // 1. Use a heap to pop sv_groups across all samples in chromosome order
    // 2. Rough clustering of potential duplicates into pools as they're popped off the heap
    //

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

    // Scan through genome, using the heap to merge sort breakpoints across all samples so that
    // duplicates can be found:
    //

    // Expand last breakpoint region by this amount to the left and right, before testing intersection with the current
    // breakpoint
    let duplicate_cluster_expansion = 100;

    // Collection of svgroups that should be evaluated to determine if any are duplicates:
    let mut duplicate_candidate_pool = Vec::new();
    let mut last_segments: Vec<GenomeSegment> = Vec::new();
    while let Some(HeapKey {
        bp,
        sample_index,
        sv_group_index,
        ..
    }) = heap.pop()
    {
        if debug {
            eprintln!(
                "Multi-Merge: Popped sample {sample_index} sv_group_index {sv_group_index} bp {bp:?}"
            );
        }

        // Use simple duplicate pool criteria, any segment set which intersects the
        // (expanded) first segment set in the pool is included.
        let segments = {
            let segment1 = bp.breakend1.segment.clone();
            let mut segments = vec![segment1];
            if let Some(be2) = &bp.breakend2 {
                let segment2 = be2.segment.clone();
                segments.push(segment2);
            }
            segments
        };
        let store_current_duplicate_pool = if last_segments.len() == segments.len() {
            last_segments
                .iter()
                .zip(segments.iter())
                .any(|(last_seg, new_seg)| !last_seg.intersect(new_seg))
        } else {
            true
        };

        if debug {
            eprintln!("Multi-Merge: last_segments: {last_segments:?} segments: {segments:?}");
            eprintln!("process_duplicate_pool: {store_current_duplicate_pool}");
        }

        if store_current_duplicate_pool {
            store_duplicate_pool(
                &mut duplicate_candidate_pool,
                &mut duplicate_candidate_pools,
                debug,
            );

            last_segments = segments
                .into_iter()
                .map(|mut x| {
                    x.expand_by(chrom_list, duplicate_cluster_expansion);
                    x
                })
                .collect::<Vec<_>>();
        }

        duplicate_candidate_pool.push(CandidateSVGroupInfo {
            sample_index,
            sv_group_index,
        });

        // Always add from the same sample we just popped, so that the heap contains (up to)
        // one sv_group from each sample
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
