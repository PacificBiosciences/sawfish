//! Support for reusing refinement results computed from a different breakpoint cluster list
//!
//! Breakpoint clusters can be refined speculatively, before all breakpoint observations are
//! available, from clusters built on the observations found so far. A speculative result can be
//! reused for a final cluster only if the speculative cluster is the same as the final cluster in
//! every way that refinement can see.
//!
//! Clusters refer to each other by their index in the cluster list (for consolidated indel
//! clusters, large insertion soft-clip pairs and neighboring breakends), and refinement results
//! carry their own cluster index (and any neighbor cluster indices) in SV and contig IDs. Indices
//! differ between cluster lists, so clusters are compared through a signature in which every
//! cluster index is replaced by the identity of the cluster it refers to, and reused results have
//! all cluster indices translated to the final cluster list.
//!

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rust_vc_utils::GenomeSegment;

use super::debug_output::ClusterDebugInfo;
use super::{RefineWorkerReturnType, RefinedSV};
use crate::breakpoint::{
    BreakendNeighbor, BreakendNeighborObservation, Breakpoint, BreakpointCluster,
    ConsolidatedAssemblySegmentInfo, InsertInfo, LargeInsertionInfo,
};
use crate::contig_output::ContigAlignmentInfo;
use crate::sv_group::{GroupHaplotypeId, SVGroup, SVGroupHaplotype};
use crate::sv_id::SVUniqueIdData;

/// Identity of a cluster for matching clusters between cluster lists
///
/// This is the cluster breakpoint without its neighbor references (which are cluster indices)
///
pub type ClusterIdentity = Breakpoint;

// The functions below destructure these structs without `..`, so that adding a field to any of
// them fails the build here. Whoever adds the field then has to decide whether it holds a cluster
// index (which has to be translated), or affects refinement (which has to be in the signature).

fn strip_neighbors(bp: &Breakpoint) -> Breakpoint {
    let Breakpoint {
        breakend1,
        breakend2,
        insert_info,
        is_precise,
        // Neighbors refer to clusters by index, so they are part of the normalized breakpoint
        breakend1_neighbor: _,
        breakend2_neighbor: _,
    } = bp.clone();
    Breakpoint {
        breakend1,
        breakend2,
        insert_info,
        is_precise,
        breakend1_neighbor: None,
        breakend2_neighbor: None,
    }
}

pub fn get_cluster_identity(bpc: &BreakpointCluster) -> ClusterIdentity {
    strip_neighbors(&bpc.breakpoint)
}

/// Breakpoint with neighbor references given as (neighbor cluster identity, breakend index)
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
struct NormalizedBreakpoint {
    breakpoint: Breakpoint,
    neighbors: [Option<(ClusterIdentity, usize, usize)>; 2],
}

/// Everything about a cluster which can affect its refinement, with cluster indices replaced by
/// cluster identities
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct ClusterSignature {
    breakpoint: NormalizedBreakpoint,
    evidence: Vec<usize>,
    evidence_qnames: BTreeSet<Vec<u8>>,
    assembly_segments: Vec<GenomeSegment>,
    consolidated_assembly_segment: Option<(GenomeSegment, InsertInfo, Vec<ClusterIdentity>)>,
    consolidated_cluster: Option<ClusterIdentity>,
    large_insertion: Option<(GenomeSegment, ClusterIdentity, ClusterIdentity)>,
    neighbor_observations: [Vec<BreakendNeighborObservation>; 2],
    neighbor_breakpoints: [Option<NormalizedBreakpoint>; 2],
}

fn normalize_breakpoint(clusters: &[BreakpointCluster], bp: &Breakpoint) -> NormalizedBreakpoint {
    let neighbor = |x: &Option<BreakendNeighbor>| {
        x.as_ref().map(|n| {
            let BreakendNeighbor {
                sample_index,
                cluster_index,
                breakend_index,
            } = n;
            (
                get_cluster_identity(&clusters[*cluster_index]),
                *breakend_index,
                *sample_index,
            )
        })
    };
    NormalizedBreakpoint {
        breakpoint: strip_neighbors(bp),
        neighbors: [
            neighbor(&bp.breakend1_neighbor),
            neighbor(&bp.breakend2_neighbor),
        ],
    }
}

pub fn get_cluster_signature(
    clusters: &[BreakpointCluster],
    cluster_index: usize,
) -> ClusterSignature {
    let BreakpointCluster {
        breakpoint,
        evidence,
        evidence_qnames,
        assembly_segments,
        consolidated_assembly_segment,
        consolidated_cluster_index,
        large_insertion_info,
        breakend1_neighbor_observations,
        breakend2_neighbor_observations,
        breakend1_neighbor_breakpoint,
        breakend2_neighbor_breakpoint,
    } = &clusters[cluster_index];
    let identity = |i: usize| get_cluster_identity(&clusters[i]);
    ClusterSignature {
        breakpoint: normalize_breakpoint(clusters, breakpoint),
        evidence: evidence.to_vec(),
        evidence_qnames: evidence_qnames.clone(),
        assembly_segments: assembly_segments.clone(),
        consolidated_assembly_segment: consolidated_assembly_segment.as_ref().map(|x| {
            let ConsolidatedAssemblySegmentInfo {
                segment,
                insert_info,
                cluster_ids,
            } = x;
            (
                segment.clone(),
                insert_info.clone(),
                cluster_ids.iter().map(|&i| identity(i)).collect(),
            )
        }),
        consolidated_cluster: consolidated_cluster_index.map(identity),
        large_insertion: large_insertion_info.as_ref().map(|x| {
            let LargeInsertionInfo {
                segment,
                root_cluster_index,
                paired_cluster_index,
            } = x;
            (
                segment.clone(),
                identity(*root_cluster_index),
                identity(*paired_cluster_index),
            )
        }),
        neighbor_observations: [
            breakend1_neighbor_observations.clone(),
            breakend2_neighbor_observations.clone(),
        ],
        neighbor_breakpoints: [
            breakend1_neighbor_breakpoint
                .as_ref()
                .map(|x| normalize_breakpoint(clusters, x)),
            breakend2_neighbor_breakpoint
                .as_ref()
                .map(|x| normalize_breakpoint(clusters, x)),
        ],
    }
}

/// Map from cluster identity to cluster index, for identities found in only one cluster
pub fn get_unique_identity_index_map(
    clusters: &[BreakpointCluster],
) -> BTreeMap<ClusterIdentity, usize> {
    let mut map = BTreeMap::new();
    let mut duplicates = BTreeSet::new();
    for (cluster_index, bpc) in clusters.iter().enumerate() {
        let identity = get_cluster_identity(bpc);
        if map.insert(identity.clone(), cluster_index).is_some() {
            duplicates.insert(identity);
        }
    }
    for identity in duplicates {
        map.remove(&identity);
    }
    map
}

/// All cluster indices found in a refinement result
pub fn get_result_cluster_indices(result: &RefineWorkerReturnType) -> BTreeSet<usize> {
    let mut indices = BTreeSet::new();
    let (debug_info, contig_alignments, sv_group) = result;
    let add_bp = |bp: &Breakpoint, indices: &mut BTreeSet<usize>| {
        for n in [&bp.breakend1_neighbor, &bp.breakend2_neighbor]
            .into_iter()
            .flatten()
        {
            let BreakendNeighbor {
                sample_index: _,
                cluster_index,
                breakend_index: _,
            } = n;
            indices.insert(*cluster_index);
        }
    };
    let ClusterDebugInfo { duration: _, bp } = debug_info;
    add_bp(bp, &mut indices);
    for x in contig_alignments.iter() {
        let ContigAlignmentInfo {
            sample_index: _,
            cluster_index,
            assembly_index: _,
            seq: _,
            segment_id: _,
            segments: _,
            supporting_read_count: _,
            high_quality_contig_range: _,
        } = x;
        indices.insert(*cluster_index);
    }
    if let Some(sv_group) = sv_group {
        let SVGroup {
            group_regions: _,
            group_haplotypes,
            sample_haplotype_list: _,
            sample_expected_cn_info: _,
            sv_haplotype_map: _,
            refined_svs,
        } = sv_group;
        for rsv in refined_svs.iter() {
            let RefinedSV {
                id,
                bp,
                breakend1_homology_seq: _,
                single_region_refinement: _,
                contig_pos_before_breakend1: _,
                assembly_regions: _,
                ext: _,
                score: _,
            } = rsv;
            let SVUniqueIdData {
                sample_index: _,
                cluster_index,
                assembly_index: _,
                alignment_index: _,
            } = id;
            indices.insert(*cluster_index);
            add_bp(bp, &mut indices);
        }
        for hap in group_haplotypes.iter() {
            let SVGroupHaplotype {
                hap_id,
                contig_info: _,
            } = hap;
            let GroupHaplotypeId {
                sample_index: _,
                cluster_index,
                assembly_index: _,
            } = hap_id;
            indices.insert(*cluster_index);
        }
    }
    indices
}

/// Translate all cluster indices in a refinement result
///
/// index_map must contain every index from get_result_cluster_indices
///
pub fn remap_result_cluster_indices(
    result: &mut RefineWorkerReturnType,
    index_map: &BTreeMap<usize, usize>,
) {
    let (debug_info, contig_alignments, sv_group) = result;
    let remap_bp = |bp: &mut Breakpoint| {
        for n in [&mut bp.breakend1_neighbor, &mut bp.breakend2_neighbor]
            .into_iter()
            .flatten()
        {
            let BreakendNeighbor {
                sample_index: _,
                cluster_index,
                breakend_index: _,
            } = n;
            *cluster_index = index_map[cluster_index];
        }
    };
    let ClusterDebugInfo { duration: _, bp } = debug_info;
    remap_bp(bp);
    for x in contig_alignments.iter_mut() {
        let ContigAlignmentInfo {
            sample_index: _,
            cluster_index,
            assembly_index: _,
            seq: _,
            segment_id: _,
            segments: _,
            supporting_read_count: _,
            high_quality_contig_range: _,
        } = x;
        *cluster_index = index_map[cluster_index];
    }
    if let Some(sv_group) = sv_group {
        let SVGroup {
            group_regions: _,
            group_haplotypes,
            sample_haplotype_list: _,
            sample_expected_cn_info: _,
            sv_haplotype_map: _,
            refined_svs,
        } = sv_group;
        for rsv in refined_svs.iter_mut() {
            let RefinedSV {
                id,
                bp,
                breakend1_homology_seq: _,
                single_region_refinement: _,
                contig_pos_before_breakend1: _,
                assembly_regions: _,
                ext: _,
                score: _,
            } = rsv;
            let SVUniqueIdData {
                sample_index: _,
                cluster_index,
                assembly_index: _,
                alignment_index: _,
            } = id;
            *cluster_index = index_map[cluster_index];
            remap_bp(bp);
        }
        for hap in group_haplotypes.iter_mut() {
            let SVGroupHaplotype {
                hap_id,
                contig_info: _,
            } = hap;
            let GroupHaplotypeId {
                sample_index: _,
                cluster_index,
                assembly_index: _,
            } = hap_id;
            *cluster_index = index_map[cluster_index];
        }
    }
}

/// All cluster indices referenced by a cluster, including its own index
///
/// Refinement results only refer to clusters from this set.
///
fn get_cluster_referenced_indices(
    bpc: &BreakpointCluster,
    cluster_index: usize,
) -> BTreeSet<usize> {
    let BreakpointCluster {
        breakpoint,
        evidence: _,
        evidence_qnames: _,
        assembly_segments: _,
        consolidated_assembly_segment,
        consolidated_cluster_index,
        large_insertion_info,
        breakend1_neighbor_observations: _,
        breakend2_neighbor_observations: _,
        breakend1_neighbor_breakpoint,
        breakend2_neighbor_breakpoint,
    } = bpc;
    let mut indices = BTreeSet::from([cluster_index]);
    let add_bp = |bp: &Breakpoint, indices: &mut BTreeSet<usize>| {
        for n in [&bp.breakend1_neighbor, &bp.breakend2_neighbor]
            .into_iter()
            .flatten()
        {
            indices.insert(n.cluster_index);
        }
    };
    add_bp(breakpoint, &mut indices);
    for bp in [breakend1_neighbor_breakpoint, breakend2_neighbor_breakpoint]
        .into_iter()
        .flatten()
    {
        add_bp(bp, &mut indices);
    }
    if let Some(x) = consolidated_assembly_segment {
        indices.extend(x.cluster_ids.iter().copied());
    }
    indices.extend(*consolidated_cluster_index);
    if let Some(x) = large_insertion_info {
        indices.insert(x.root_cluster_index);
        indices.insert(x.paired_cluster_index);
    }
    indices
}

/// Identity of each cluster index referenced by a cluster
fn get_cluster_referenced_identities(
    clusters: &[BreakpointCluster],
    cluster_index: usize,
) -> BTreeMap<usize, ClusterIdentity> {
    get_cluster_referenced_indices(&clusters[cluster_index], cluster_index)
        .into_iter()
        .map(|i| (i, get_cluster_identity(&clusters[i])))
        .collect()
}

/// A refinement result computed from a speculative cluster list
pub struct SpeculativeResult {
    /// Identity of each cluster index referenced in the result
    pub index_identities: BTreeMap<usize, ClusterIdentity>,
    pub result: RefineWorkerReturnType,
}

impl SpeculativeResult {
    /// Return None if the result refers to any cluster index without a known identity, so that the
    /// result can't be reused
    pub fn new(
        mut identities: BTreeMap<usize, ClusterIdentity>,
        result: RefineWorkerReturnType,
    ) -> Option<Self> {
        let mut index_identities = BTreeMap::new();
        for i in get_result_cluster_indices(&result) {
            index_identities.insert(i, identities.remove(&i)?);
        }
        Some(Self {
            index_identities,
            result,
        })
    }

    /// Convert the result to the final cluster list, or return None if any referenced cluster
    /// can't be found there
    pub fn into_final(
        mut self,
        final_identity_index_map: &BTreeMap<ClusterIdentity, usize>,
    ) -> Option<RefineWorkerReturnType> {
        let mut index_map = BTreeMap::new();
        for (index, identity) in self.index_identities.iter() {
            index_map.insert(*index, *final_identity_index_map.get(identity)?);
        }
        remap_result_cluster_indices(&mut self.result, &index_map);
        Some(self.result)
    }
}

/// Chromosomes of all regions used to refine a cluster
fn get_cluster_chrom_indices(bpc: &BreakpointCluster) -> BTreeSet<usize> {
    let mut chroms = BTreeSet::new();
    let add_bp = |bp: &Breakpoint, chroms: &mut BTreeSet<usize>| {
        chroms.insert(bp.breakend1.segment.chrom_index);
        if let Some(be2) = &bp.breakend2 {
            chroms.insert(be2.segment.chrom_index);
        }
    };
    add_bp(&bpc.breakpoint, &mut chroms);
    for bp in [
        &bpc.breakend1_neighbor_breakpoint,
        &bpc.breakend2_neighbor_breakpoint,
    ]
    .into_iter()
    .flatten()
    {
        add_bp(bp, &mut chroms);
    }
    for x in bpc.assembly_segments.iter() {
        chroms.insert(x.chrom_index);
    }
    if let Some(x) = &bpc.consolidated_assembly_segment {
        chroms.insert(x.segment.chrom_index);
    }
    if let Some(x) = &bpc.large_insertion_info {
        chroms.insert(x.segment.chrom_index);
    }
    chroms
}

/// Refine breakpoint clusters speculatively during the alignment scan
///
/// Each time chromosomes complete, the breakpoint observations found so far are clustered, and any
/// new clusters entirely on completed chromosomes are refined. The results can then be reused for
/// any final cluster with the same signature (see module documentation).
///
/// Runs until obs_rx is closed and all refinement jobs are complete. thread_count refinement jobs
/// run at once until obs_rx is closed, then all threads are used.
///
#[allow(clippy::too_many_arguments)]
pub fn run_speculative_refinement(
    obs_rx: std::sync::mpsc::Receiver<(usize, crate::cluster_breakpoints::BreakObservations)>,
    shared_settings: &crate::cli::SharedSettings,
    settings: &crate::cli::DiscoverSettings,
    static_settings: &crate::discover::StaticDiscoverSettings,
    worker_thread_dataset: &crate::worker_thread_data::BamReaderWorkerThreadDataSet,
    reference: &crate::genome_ref_utils::GenomeRef,
    chrom_list: &rust_vc_utils::ChromList,
    thread_count: usize,
) -> BTreeMap<Arc<ClusterSignature>, SpeculativeResult> {
    use std::sync::{Condvar, Mutex};

    // Limit on the number of refinement jobs running at once, which is raised to all threads once the
    // alignment scan is complete
    struct JobPermits {
        count: Mutex<usize>,
        cond: Condvar,
    }
    impl JobPermits {
        fn acquire(&self) {
            let mut count = self.count.lock().unwrap();
            while *count == 0 {
                count = self.cond.wait(count).unwrap();
            }
            *count -= 1;
        }
        fn release(&self, n: usize) {
            *self.count.lock().unwrap() += n;
            self.cond.notify_all();
        }
    }
    let permits = JobPermits {
        count: Mutex::new(thread_count),
        cond: Condvar::new(),
    };

    let refine_settings = &super::get_refine_settings(settings, static_settings);
    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();
    let results = Mutex::new(BTreeMap::new());

    // The scope runs on this thread, which only waits for observations and clusters them, while
    // refinement jobs run on the worker pool
    worker_pool.in_place_scope(|scope| {
        let results = &results;
        let permits = &permits;
        let mut submitted = BTreeSet::new();
        let mut observations = crate::cluster_breakpoints::BreakObservations::new();
        let mut completed_chroms = vec![false; chrom_list.data.len()];
        while let Ok(first) = obs_rx.recv() {
            // Take all completed chromosomes before clustering again
            let mut new_obs = vec![first];
            new_obs.extend(obs_rx.try_iter());
            for (chrom_index, mut chrom_obs) in new_obs {
                completed_chroms[chrom_index] = true;
                observations.merge(&mut chrom_obs);
            }

            let clusters = crate::cluster_breakpoints::get_breakpoint_clusters_quiet(
                settings,
                static_settings.assembly_read_flank_size,
                chrom_list,
                &mut observations,
            );
            for cluster_index in super::get_refine_job_indices(shared_settings, settings, &clusters)
            {
                if !get_cluster_chrom_indices(&clusters[cluster_index])
                    .into_iter()
                    .all(|x| completed_chroms[x])
                {
                    continue;
                }
                let signature = Arc::new(get_cluster_signature(&clusters, cluster_index));
                if !submitted.insert(signature.clone()) {
                    continue;
                }

                // Each job keeps only the data it needs from the cluster list, so that queued jobs
                // don't hold every cluster list built during the scan
                let cluster = clusters[cluster_index].clone();
                let identities = get_cluster_referenced_identities(&clusters, cluster_index);
                scope.spawn(move |_| {
                    permits.acquire();
                    let result = super::refine_cluster(
                        worker_thread_dataset,
                        refine_settings,
                        reference,
                        chrom_list,
                        &cluster,
                        cluster_index,
                    );
                    permits.release(1);
                    if let Some(result) = SpeculativeResult::new(identities, result) {
                        results.lock().unwrap().insert(signature, result);
                    }
                });
            }
        }

        // The alignment scan is complete, so all threads can be used
        permits.release(shared_settings.thread_count - thread_count);
    });

    results.into_inner().unwrap()
}

/// Convert speculative results into results for the final cluster list, for each refinement job
/// with a speculative result for the same cluster signature
pub fn get_reusable_results(
    shared_settings: &crate::cli::SharedSettings,
    settings: &crate::cli::DiscoverSettings,
    clusters: &[BreakpointCluster],
    mut speculative_results: BTreeMap<Arc<ClusterSignature>, SpeculativeResult>,
) -> BTreeMap<usize, RefineWorkerReturnType> {
    let identity_index_map = get_unique_identity_index_map(clusters);
    let mut reused_results = BTreeMap::new();
    if speculative_results.is_empty() {
        return reused_results;
    }
    let job_indices = super::get_refine_job_indices(shared_settings, settings, clusters);
    let job_count = job_indices.len();
    for cluster_index in job_indices {
        let signature = get_cluster_signature(clusters, cluster_index);
        if let Some(result) = speculative_results.remove(&signature)
            && let Some(result) = result.into_final(&identity_index_map)
        {
            reused_results.insert(cluster_index, result);
        }
    }
    log::info!(
        "Reusing speculative refinement results for {} of {} breakpoint clusters",
        reused_results.len(),
        job_count
    );
    reused_results
}

#[cfg(test)]
mod tests {
    use rust_vc_utils::IntRange;

    use super::*;
    use crate::refine_sv::RefinedSVExt;
    use crate::sv_group::ClusterAssembly;

    fn segment(chrom_index: usize, start: i64) -> GenomeSegment {
        GenomeSegment {
            chrom_index,
            range: IntRange::from_pair(start, start + 1),
        }
    }

    fn cluster(chrom_index: usize, pos: i64) -> BreakpointCluster {
        let mut bpc = BreakpointCluster::default();
        bpc.breakpoint.breakend1.segment = segment(chrom_index, pos);
        bpc.evidence[0] = 2;
        bpc
    }

    fn neighbor(cluster_index: usize) -> Option<BreakendNeighbor> {
        Some(BreakendNeighbor {
            sample_index: 0,
            cluster_index,
            breakend_index: 1,
        })
    }

    /// Four clusters which refer to each other in every way a cluster can: breakend neighbors
    /// (with the neighbor's breakpoint copied), consolidated indel clusters, and a large insertion
    fn linked_clusters() -> Vec<BreakpointCluster> {
        let mut clusters = (0..4)
            .map(|i| cluster(0, 100 * (i + 1)))
            .collect::<Vec<_>>();

        clusters[1].breakpoint.breakend1_neighbor = neighbor(2);
        clusters[1].breakend1_neighbor_breakpoint = Some(clusters[2].breakpoint.clone());
        clusters[1].breakend1_neighbor_observations = vec![BreakendNeighborObservation {
            breakend: clusters[2].breakpoint.breakend1.clone(),
            neighbor_dist: 50,
        }];
        clusters[2].breakpoint.breakend2_neighbor = neighbor(1);
        clusters[2].breakend2_neighbor_breakpoint = Some(clusters[1].breakpoint.clone());

        clusters[0].consolidated_cluster_index = Some(3);
        clusters[3].consolidated_assembly_segment = Some(ConsolidatedAssemblySegmentInfo {
            segment: segment(0, 90),
            insert_info: InsertInfo::Seq(b"ACGT".to_vec()),
            cluster_ids: vec![0, 3],
        });
        clusters[3].assembly_segments = vec![segment(0, 80), segment(0, 420)];
        clusters[3].evidence_qnames = BTreeSet::from([b"read1".to_vec()]);

        clusters[2].large_insertion_info = Some(LargeInsertionInfo {
            segment: segment(0, 150),
            root_cluster_index: 1,
            paired_cluster_index: 2,
        });
        clusters
    }

    /// Reorder clusters, translating all the indices they contain. new_index[i] is the new
    /// position of cluster i.
    fn permute_clusters(
        clusters: &[BreakpointCluster],
        new_index: &[usize],
    ) -> Vec<BreakpointCluster> {
        let remap_bp = |bp: &mut Breakpoint| {
            for n in [&mut bp.breakend1_neighbor, &mut bp.breakend2_neighbor]
                .into_iter()
                .flatten()
            {
                n.cluster_index = new_index[n.cluster_index];
            }
        };
        let mut result = vec![BreakpointCluster::default(); clusters.len()];
        for (i, bpc) in clusters.iter().enumerate() {
            let mut bpc = bpc.clone();
            remap_bp(&mut bpc.breakpoint);
            for bp in [
                &mut bpc.breakend1_neighbor_breakpoint,
                &mut bpc.breakend2_neighbor_breakpoint,
            ]
            .into_iter()
            .flatten()
            {
                remap_bp(bp);
            }
            if let Some(x) = bpc.consolidated_cluster_index.as_mut() {
                *x = new_index[*x];
            }
            if let Some(x) = bpc.consolidated_assembly_segment.as_mut() {
                x.cluster_ids.iter_mut().for_each(|y| *y = new_index[*y]);
            }
            if let Some(x) = bpc.large_insertion_info.as_mut() {
                x.root_cluster_index = new_index[x.root_cluster_index];
                x.paired_cluster_index = new_index[x.paired_cluster_index];
            }
            result[new_index[i]] = bpc;
        }
        result
    }

    #[test]
    fn test_signature_is_independent_of_cluster_order() {
        let clusters = linked_clusters();
        let new_index = [2, 0, 3, 1];
        let permuted = permute_clusters(&clusters, &new_index);
        for i in 0..clusters.len() {
            assert!(
                get_cluster_signature(&clusters, i)
                    == get_cluster_signature(&permuted, new_index[i]),
                "cluster {i}"
            );
        }
        // And the clusters really did move
        assert!(
            clusters[1]
                .breakpoint
                .breakend1_neighbor
                .as_ref()
                .unwrap()
                .cluster_index
                == 2
        );
        assert!(
            permuted[new_index[1]]
                .breakpoint
                .breakend1_neighbor
                .as_ref()
                .unwrap()
                .cluster_index
                == new_index[2]
        );
    }

    #[test]
    fn test_signature_changes_with_every_cluster_field() {
        let clusters = linked_clusters();
        let original = get_cluster_signature(&clusters, 1);

        type Edit = Box<dyn Fn(&mut Vec<BreakpointCluster>)>;
        let edits: Vec<(&str, Edit)> = vec![
            (
                "breakpoint",
                Box::new(|c| c[1].breakpoint.breakend1.segment = segment(0, 201)),
            ),
            (
                "is_precise",
                Box::new(|c| c[1].breakpoint.is_precise = true),
            ),
            (
                "insert_info",
                Box::new(|c| c[1].breakpoint.insert_info = InsertInfo::Seq(b"A".to_vec())),
            ),
            ("evidence", Box::new(|c| c[1].evidence[0] += 1)),
            (
                "evidence_qnames",
                Box::new(|c| {
                    c[1].evidence_qnames.insert(b"q".to_vec());
                }),
            ),
            (
                "assembly_segments",
                Box::new(|c| c[1].assembly_segments.push(segment(0, 5))),
            ),
            (
                "consolidated_assembly_segment",
                Box::new(|c| {
                    c[1].consolidated_assembly_segment = Some(ConsolidatedAssemblySegmentInfo {
                        segment: segment(0, 1),
                        insert_info: InsertInfo::Seq(Vec::new()),
                        cluster_ids: vec![1],
                    })
                }),
            ),
            (
                "consolidated_cluster_index",
                Box::new(|c| c[1].consolidated_cluster_index = Some(0)),
            ),
            (
                "large_insertion_info",
                Box::new(|c| {
                    c[1].large_insertion_info = Some(LargeInsertionInfo {
                        segment: segment(0, 1),
                        root_cluster_index: 0,
                        paired_cluster_index: 1,
                    })
                }),
            ),
            (
                "breakend1_neighbor_observations",
                Box::new(|c| c[1].breakend1_neighbor_observations[0].neighbor_dist += 1),
            ),
            (
                "breakend2_neighbor_observations",
                Box::new(|c| {
                    let x = c[1].breakend1_neighbor_observations[0].clone();
                    c[1].breakend2_neighbor_observations.push(x)
                }),
            ),
            (
                "breakend1_neighbor_breakpoint",
                Box::new(|c| {
                    c[1].breakend1_neighbor_breakpoint
                        .as_mut()
                        .unwrap()
                        .is_precise = true
                }),
            ),
            (
                "breakend2_neighbor_breakpoint",
                Box::new(|c| c[1].breakend2_neighbor_breakpoint = Some(Breakpoint::default())),
            ),
            // The neighbor is another cluster, with a different identity
            (
                "breakpoint neighbor",
                Box::new(|c| c[1].breakpoint.breakend1_neighbor = neighbor(3)),
            ),
            (
                "neighbor breakend index",
                Box::new(|c| {
                    c[1].breakpoint
                        .breakend1_neighbor
                        .as_mut()
                        .unwrap()
                        .breakend_index = 0
                }),
            ),
        ];
        for (label, edit) in edits {
            let mut edited = clusters.clone();
            edit(&mut edited);
            assert!(get_cluster_signature(&edited, 1) != original, "{label}");
        }

        // A change to a cluster that this one refers to changes the signature, through the
        // identity of that cluster
        let mut edited = clusters.clone();
        edited[2].breakpoint.breakend1.segment = segment(0, 301);
        assert!(get_cluster_signature(&edited, 1) != original);
    }

    fn refinement_result(
        neighbor_cluster_index: usize,
        cluster_index: usize,
    ) -> RefineWorkerReturnType {
        let bp = Breakpoint {
            breakend1_neighbor: neighbor(neighbor_cluster_index),
            ..Default::default()
        };
        let debug_info = ClusterDebugInfo {
            duration: Default::default(),
            bp: bp.clone(),
        };
        let contig_alignment = ContigAlignmentInfo {
            sample_index: 0,
            cluster_index,
            assembly_index: 0,
            seq: b"ACGT".to_vec(),
            segment_id: 0,
            segments: Vec::new(),
            supporting_read_count: 3,
            high_quality_contig_range: IntRange::from_int(0),
        };
        let refined_sv = RefinedSV {
            id: SVUniqueIdData {
                sample_index: 0,
                cluster_index,
                assembly_index: 0,
                alignment_index: 0,
            },
            bp,
            breakend1_homology_seq: Vec::new(),
            single_region_refinement: false,
            contig_pos_before_breakend1: 0,
            assembly_regions: Vec::new(),
            ext: RefinedSVExt::default(),
            score: Default::default(),
        };
        let sv_group = SVGroup {
            group_regions: vec![segment(0, 1)],
            group_haplotypes: vec![SVGroupHaplotype {
                hap_id: GroupHaplotypeId {
                    sample_index: 0,
                    cluster_index,
                    assembly_index: 0,
                },
                contig_info: ClusterAssembly::default(),
            }],
            sample_haplotype_list: Vec::new(),
            sample_expected_cn_info: Vec::new(),
            sv_haplotype_map: vec![0],
            refined_svs: vec![refined_sv],
        };
        (debug_info, vec![contig_alignment], Some(sv_group))
    }

    #[test]
    fn test_result_indices_are_translated() {
        let result = refinement_result(2, 1);
        assert_eq!(get_result_cluster_indices(&result), BTreeSet::from([1, 2]));

        let mut remapped = refinement_result(2, 1);
        remap_result_cluster_indices(&mut remapped, &BTreeMap::from([(1, 7), (2, 5)]));
        assert_eq!(
            get_result_cluster_indices(&remapped),
            BTreeSet::from([5, 7])
        );
        let (debug_info, contigs, group) = &remapped;
        let group = group.as_ref().unwrap();
        assert_eq!(
            debug_info
                .bp
                .breakend1_neighbor
                .as_ref()
                .unwrap()
                .cluster_index,
            5
        );
        assert_eq!(contigs[0].cluster_index, 7);
        assert_eq!(group.refined_svs[0].id.cluster_index, 7);
        assert_eq!(
            group.refined_svs[0]
                .bp
                .breakend1_neighbor
                .as_ref()
                .unwrap()
                .cluster_index,
            5
        );
        assert_eq!(group.group_haplotypes[0].hap_id.cluster_index, 7);
    }

    /// A result computed on one cluster list is reused on a permuted list, with its indices
    /// translated, and only if every cluster it refers to is found there
    #[test]
    fn test_speculative_result_reuse() {
        let clusters = linked_clusters();
        let new_index = [2, 0, 3, 1];
        let permuted = permute_clusters(&clusters, &new_index);

        let result = refinement_result(2, 1);
        let identities = get_cluster_referenced_identities(&clusters, 1);
        let speculative = SpeculativeResult::new(identities.clone(), result).unwrap();
        let final_map = get_unique_identity_index_map(&permuted);
        let (_, contigs, group) = speculative.into_final(&final_map).unwrap();
        assert_eq!(contigs[0].cluster_index, new_index[1]);
        let group = group.unwrap();
        assert_eq!(group.refined_svs[0].id.cluster_index, new_index[1]);
        assert_eq!(
            group.refined_svs[0]
                .bp
                .breakend1_neighbor
                .as_ref()
                .unwrap()
                .cluster_index,
            new_index[2]
        );

        // The referenced cluster is missing from the final list
        let mut missing = final_map.clone();
        missing.remove(&get_cluster_identity(&clusters[2]));
        let speculative =
            SpeculativeResult::new(identities.clone(), refinement_result(2, 1)).unwrap();
        assert!(speculative.into_final(&missing).is_none());

        // A result that refers to a cluster the speculative cluster didn't reference can't be kept
        assert!(SpeculativeResult::new(identities, refinement_result(3, 1)).is_none());

        // Duplicate identities are not matched
        let mut duplicated = permuted.clone();
        duplicated[0] = duplicated[1].clone();
        assert!(get_unique_identity_index_map(&duplicated).len() < permuted.len() - 1);
    }
}
