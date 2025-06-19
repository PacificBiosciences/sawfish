use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use camino::Utf8Path;
use log::info;
use rust_vc_utils::ChromList;
use unwrap::unwrap;

use super::large_insertions::find_large_insertion_candidates;
use super::{BreakObservations, BreakpointObservationMap};
use crate::breakpoint::*;
use crate::cli;
use crate::genome_segment::{GenomeSegment, get_segment_distance};
use crate::run_stats::ClusterStats;

fn debug_bpos(output_dir: &Utf8Path, bpos: &BreakpointObservationMap, num: usize) {
    let filename = output_dir.join("debug.bpos.txt");
    info!("Writing debug bpos to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create debug bpos file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    let mut reported = 0;
    for bpo_vec in bpos.data.values() {
        for bpo in bpo_vec {
            writeln!(f, "{bpo:?}",).unwrap();
            reported += 1;
            if reported >= num {
                return;
            }
        }
    }
}

fn debug_clusters(output_dir: &Utf8Path, breakpoint_clusters: &[BreakpointCluster], num: usize) {
    let filename = output_dir.join("debug.breakpoint_clusters.txt");
    info!("Writing breakpoint clusters debug file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create breakpoint clusters debug file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    for breakpoint_cluster in breakpoint_clusters.iter().take(num) {
        writeln!(f, "{breakpoint_cluster:?}",).unwrap();
    }
}

struct ClusterBedRecord {
    pub segment: GenomeSegment,
    pub cluster_id: usize,
    pub breakend_id: usize,
    pub dir: String,
    pub ins_size: usize,
    pub evidence_count: usize,
    pub assembly_segment_count: usize,
    pub consolidated_cluster_index: Option<usize>,
    pub breakend_neighbor_index: Option<usize>,
    pub unpaired_breakend: bool,
}

/// Produce bed file for debugging in IGV, this uses the IGV gffTags format for get a formatted
/// display of multiple key/value pairs in column #4
///
fn write_breakpoint_clusters_to_bed(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    clusters: &[BreakpointCluster],
) {
    // First pull all the records together for sorting:
    let mut brecs = Vec::new();

    for (cluster_index, bpc) in clusters.iter().enumerate() {
        let bp = &bpc.breakpoint;

        let assembly_segment_count = bpc.assembly_segments.len();

        // be1
        let unpaired_breakend = bp.breakend2.is_none();
        {
            let be1 = &bp.breakend1;
            let breakend_neighbor_index = bpc.breakend1_neighbor.as_ref().map(|x| x.cluster_index);
            brecs.push(ClusterBedRecord {
                segment: be1.segment.clone(),
                cluster_id: cluster_index,
                breakend_id: 1,
                dir: format!("{:?}", be1.dir),
                ins_size: bp.insert_info.size(),
                evidence_count: bpc.evidence_count(),
                assembly_segment_count,
                consolidated_cluster_index: bpc.consolidated_cluster_index,
                breakend_neighbor_index,
                unpaired_breakend,
            });
        }

        if let Some(be2) = bp.breakend2.as_ref() {
            // be2
            let breakend_neighbor_index = bpc.breakend2_neighbor.as_ref().map(|x| x.cluster_index);
            brecs.push(ClusterBedRecord {
                segment: be2.segment.clone(),
                cluster_id: cluster_index,
                breakend_id: 2,
                dir: format!("{:?}", be2.dir),
                ins_size: bp.insert_info.size(),
                evidence_count: bpc.evidence_count(),
                assembly_segment_count,
                consolidated_cluster_index: bpc.consolidated_cluster_index,
                breakend_neighbor_index,
                unpaired_breakend,
            });
        }
    }

    brecs.sort_by(|a, b| a.segment.partial_cmp(&b.segment).unwrap());

    let filename = output_dir.join("debug.breakpoint_clusters.bed");
    info!("Writing breakpoint_cluster debug bed file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create breakpoint cluster debug bed file: '{filename}'"
    );
    let mut f = BufWriter::new(f);
    writeln!(f, "#gffTags").unwrap();
    for brec in brecs {
        let chrom_label = &chrom_list.data[brec.segment.chrom_index].label;
        let range = &brec.segment.range;
        writeln!(
            f,
            "{chrom_label}\t{}\t{}\tName={};BreakendId={};BreakendDir={};InsertionSize={};EvidenceCount={}\
            ;AssemblySegmentCount={};ConsolidatedClusterIndex={};BreakendNeighborIndex={}{}",
            range.start,
            range.end,
            brec.cluster_id,
            brec.breakend_id,
            brec.dir,
            brec.ins_size,
            brec.evidence_count,
            brec.assembly_segment_count,
            if let Some(x) = brec.consolidated_cluster_index { x.to_string() } else { "None".to_string() },
            if let Some(x) = brec.breakend_neighbor_index { x.to_string() } else { "None".to_string() },
            if brec.unpaired_breakend { ";UnpairedBND" } else { "" },
        )
            .unwrap();
    }
}

mod cluster {
    use super::*;

    /// Clear a cache of breakpoint clusters of one dir-type
    ///
    /// Scan through clusters in the cache, and only save those with a minimum evidence count
    ///
    fn clear_cluster_cache(
        min_cluster_evidence_count: usize,
        dir_type_clusters: &mut Vec<BreakpointCluster>,
        saved_clusters: &mut Vec<BreakpointCluster>,
    ) {
        let mut tmp_clusters = Vec::new();
        std::mem::swap(&mut tmp_clusters, dir_type_clusters);
        for cluster in tmp_clusters {
            if cluster.evidence_count() >= min_cluster_evidence_count {
                saved_clusters.push(cluster);
            }
        }
    }

    type ClusterBuffer = [Vec<BreakpointCluster>; Breakpoint::DIR_TYPE_COUNT];

    fn add_breakpoint_observation_to_clusters(
        cluster_distance: usize,
        min_cluster_evidence_count: usize,
        bpo: &BreakpointObservation,
        cluster_buffer: &mut ClusterBuffer,
        breakpoint_clusters: &mut Vec<BreakpointCluster>,
    ) {
        let debug = false;

        if debug {
            eprintln!(
                "add_breakpoint_observation_to_clusters: New Breakpoint obs: {:?}",
                bpo
            );
        }
        let dir_type_index = bpo.breakpoint.get_dir_type_index();
        let dir_type_clusters = &mut cluster_buffer[dir_type_index];

        // First find the cluster that bpo is closest to:
        //
        let mut min_cluster_dist = None;
        let mut min_cluster_index = 0;
        for (cluster_index, cluster) in dir_type_clusters.iter().enumerate() {
            if debug {
                eprintln!("checking cluster_index: {cluster_index} cluster: {cluster:?}");
            }

            // check distance between this breakend and the cluster
            let dist = get_breakpoint_manhattan_distance(&cluster.breakpoint, &bpo.breakpoint);

            if debug {
                eprintln!("bp dist to cluster: {:?}", dist);
            }

            // Only keep the minimum cluster distance:
            if let Some(d) = dist {
                if let Some(mcd) = min_cluster_dist {
                    if d < mcd {
                        min_cluster_dist = dist;
                        min_cluster_index = cluster_index;
                    }
                } else {
                    min_cluster_dist = dist;
                    min_cluster_index = cluster_index;
                }
            }

            if debug {
                eprintln!(
                    "min_cluster_index: {min_cluster_index} min_cluster_dist: {min_cluster_dist:?}"
                );
            }
        }

        // Next test if the closest cluster is close enough to join bpo into it:
        let mut is_bpo_assigned = false;
        if let Some(mcd) = min_cluster_dist {
            if mcd <= cluster_distance {
                if debug {
                    eprintln!("assigning bpo to cluster index {min_cluster_index}");
                }
                dir_type_clusters[min_cluster_index].merge_breakpoint_observation(bpo);
                is_bpo_assigned = true;
            }
        }

        if !is_bpo_assigned {
            // If bpo was not assigned, then:
            // 1. determine whether we reset the dir-type cluster cache
            // 2. If so, reset the cache
            // 3. Use bpo to start new cluster

            // 1. Determine whether we reset the dir-type cluster cache
            //
            // Because the breakpoint obs are sorted on breakend1 only, we can only safely
            // clear the cache when the breakend1 distance exceeds cluster_dist
            let safe_to_clear_cache = {
                // get breakend1 distance
                let mut min_breakend1_dist = None;
                for cluster in dir_type_clusters.iter() {
                    // update min_breakend1_dist
                    let breakend1_dist = get_breakend_distance(
                        &cluster.breakpoint.breakend1,
                        &bpo.breakpoint.breakend1,
                    );
                    if let Some(min_breakend1_dist) = min_breakend1_dist.as_mut() {
                        if let Some(breakend1_dist) = breakend1_dist {
                            *min_breakend1_dist =
                                std::cmp::min(*min_breakend1_dist, breakend1_dist);
                        }
                    } else {
                        min_breakend1_dist = breakend1_dist;
                    }
                }
                if debug {
                    eprintln!("Min be1 dist: {min_breakend1_dist:?}");
                }

                if let Some(min_dist) = min_breakend1_dist {
                    min_dist > cluster_distance
                } else {
                    true
                }
            };

            if debug {
                eprintln!("safe to clear cache {safe_to_clear_cache}");
            }
            if safe_to_clear_cache {
                if debug {
                    eprintln!("Clearing tmp cluster cache");
                }
                clear_cluster_cache(
                    min_cluster_evidence_count,
                    dir_type_clusters,
                    breakpoint_clusters,
                );
            }

            // Use bpo to start new cluster
            if debug {
                eprintln!("Starting new tmp cluster with bpo");
            }
            dir_type_clusters.push(BreakpointCluster::from_breakpoint_observation(bpo));
        }
    }

    /// Cluster breakpoint observations into breakpoint clusters
    ///
    pub fn cluster_breakpoint_observations(bpos: &mut BreakObservations) -> Vec<BreakpointCluster> {
        info!("Completing breakpoint observation clustering");

        // cluster algo intention:
        // - cluster all breakpoint observations where both breakends are in the same direction and
        // within cluster_distance
        //
        // TODO A quick temp speedup is that we'll just test pairwise distance against the cluster.
        // This isn't strictly correct and means that you could get different results from different
        // input orderings. We can fix this later by having each cluster maintain a tmp list of
        // breakends during its build, if a breakend intersects the cluster range, then check for
        // intersection of an individual cluster point.
        //
        // TODO insert size plays no role in clustering right now

        // Temporarily we make the output deterministic by sorting the breakpoint evidence before
        // clustering. (Longer term the clustering should be insensitive to observation order)
        //
        bpos.singles.normalize();

        // A set of temporary clusters is built up for each type of breakend direction combination a
        // breakpoint could have
        //
        let mut cluster_buffer: ClusterBuffer = Default::default();

        let mut breakpoint_clusters = Vec::new();
        for bpo in bpos.singles.data.values().flatten() {
            add_breakpoint_observation_to_clusters(
                bpos.cluster_distance,
                bpos.min_cluster_evidence_count,
                bpo,
                &mut cluster_buffer,
                &mut breakpoint_clusters,
            );
        }

        // Clear the entire cluster buffer
        for dir_type_clusters in cluster_buffer.iter_mut() {
            clear_cluster_cache(
                bpos.min_cluster_evidence_count,
                dir_type_clusters,
                &mut breakpoint_clusters,
            );
        }

        // We've clustered everything at this point so get rid of the individual breakpoint
        // observations:
        bpos.singles.data.clear();

        // Sort the breakpoint clusters now that the we've finished finding them all
        breakpoint_clusters.sort_by_key(|x| x.breakpoint.breakend1.segment.clone());

        /*
        // Also get rid of the evidence names which we won't need anymore
        for bpc in breakpoint_clusters.iter_mut() {
            bpc.evidence_qnames.clear();
        }
        */

        info!("Finished breakpoint observation clustering");

        breakpoint_clusters
    }
}

/// Return the region or regions of the genome to scan to build assemblies for the given breakpoint
///
fn get_assembly_segments(assembly_read_flank_size: usize, bp: &Breakpoint) -> Vec<GenomeSegment> {
    assert!(bp.is_standardized());

    let seg1 = &bp.breakend1.segment;
    if bp.breakend2.is_none() {
        return vec![seg1.clone()];
    }

    let be2 = bp.breakend2.as_ref().unwrap();
    let seg2 = &be2.segment;

    if seg1.chrom_index != seg2.chrom_index {
        return vec![seg1.clone(), seg2.clone()];
    }

    if bp.is_indel()
        && (seg2.range.start - seg1.range.end <= ((assembly_read_flank_size as i64) * 2))
    {
        // Use a single region for indel-like candidates if the breakend ranges are close enough
        let mut merge_seg = seg1.clone();
        merge_seg.range.merge(&seg2.range);
        vec![merge_seg]
    } else {
        vec![seg1.clone(), seg2.clone()]
    }
}

/// Consolidation logic - single region clusters within N bases will be consolidated.
///
/// TODO for later - figure out how 2 region clusters should consolidate
///
/// For the single region cases, if we go through them in sorted order we can simply look back to
/// the right-most position observed so far to see if they join the last cluster.
///
/// Returns the total number of single region clusters remaining after consolidation
///
fn consolidate_indel_candidates(
    cluster_settings: &ClusterSettings,
    chrom_list: &ChromList,
    clusters: &mut [BreakpointCluster],
) -> usize {
    info!("Consolidating adjacent indel candidates");

    struct ConsolidatedAssemblyData<'a> {
        segment: GenomeSegment,
        insert_info: InsertInfo,
        expanded_segment: GenomeSegment,
        cluster_ids: Vec<usize>,
        root_cluster_index: usize,
        root_cluster: &'a mut BreakpointCluster,
    }

    let mut total_single_region_asmregions_after_consolidation = 0;
    let mut last_consolidated_assembly_data: Option<ConsolidatedAssemblyData> = None;
    for (cluster_index, bpc) in clusters.iter_mut().enumerate() {
        if bpc.assembly_segments.len() != 1 {
            continue;
        }

        // Don't consolidate single-sided breakends
        if bpc.breakpoint.breakend2.is_none() {
            continue;
        }

        let query_segment = bpc.assembly_segments[0].clone();
        let expanded_query_segment = {
            let mut x = query_segment.clone();
            x.expand_by(
                chrom_list,
                cluster_settings.indel_cluster_consolidation_range as i64,
            );
            x
        };
        let mut start_new_cluster = true;
        if let Some(consolidated_data) = last_consolidated_assembly_data.as_mut() {
            let no_extended_query_overlap = expanded_query_segment.chrom_index
                != consolidated_data.expanded_segment.chrom_index
                || expanded_query_segment.range.start
                    >= consolidated_data.expanded_segment.range.end;
            if no_extended_query_overlap
                || (query_segment.range.end - consolidated_data.segment.range.start)
                    > cluster_settings.max_consolidated_assembly_region_size
            {
                if consolidated_data.cluster_ids.len() > 1 {
                    // Store the consolidated cluster information back into the root cluster
                    consolidated_data.root_cluster.consolidated_assembly_segment =
                        Some(ConsolidatedAssemblySegmentInfo {
                            segment: consolidated_data.segment.clone(),
                            insert_info: consolidated_data.insert_info.clone(),
                            cluster_ids: consolidated_data.cluster_ids.clone(),
                        });
                }
            } else {
                // Join this cluster to the consolidated cluster:
                consolidated_data.segment.range.end =
                    std::cmp::max(consolidated_data.segment.range.end, query_segment.range.end);
                consolidated_data.expanded_segment.range.end = std::cmp::max(
                    consolidated_data.expanded_segment.range.end,
                    expanded_query_segment.range.end,
                );
                consolidated_data
                    .insert_info
                    .merge(&bpc.breakpoint.insert_info);
                consolidated_data.cluster_ids.push(cluster_index);
                bpc.consolidated_cluster_index = Some(consolidated_data.root_cluster_index);
                start_new_cluster = false;
            }
        }
        if start_new_cluster {
            last_consolidated_assembly_data = Some(ConsolidatedAssemblyData {
                segment: query_segment,
                insert_info: bpc.breakpoint.insert_info.clone(),
                expanded_segment: expanded_query_segment,
                cluster_ids: vec![cluster_index],
                root_cluster_index: cluster_index,
                root_cluster: bpc,
            });
            total_single_region_asmregions_after_consolidation += 1;
        }
    }

    total_single_region_asmregions_after_consolidation
}

/// Return a 2-tuple of intersect and union sizes
#[allow(unused)]
fn get_jaccard_index(s1: &BTreeSet<Vec<u8>>, s2: &BTreeSet<Vec<u8>>) -> (usize, usize) {
    (s1.intersection(s2).count(), s1.union(s2).count())
}

/// Search all 2-region clusters to find cases where breakends in a compatible direction are close, and test if these seem to be on the same haplotype.
///
/// This was primarly designed for small inversions, but is general enough to support many close breakend connections. Connections chaining deletions
/// together with short bridges are not included here, as these just tend to create noise.
///
fn annotate_close_breakends(
    max_close_breakend_distance: usize,
    clusters: &mut [BreakpointCluster],
) {
    let debug = false;

    // First pass through we make a new list of all breakends from the breakpoints
    #[derive(Debug)]
    struct BreakendInfo {
        be: Breakend,
        cluster_index: usize,

        /// Breakend index. This can only be 0 or 1, for the 1st and 2nd breakend in the corresponding breakpoint
        breakend_index: usize,
    }

    let mut breakends = Vec::new();
    for (cluster_index, bpc) in clusters.iter().enumerate() {
        // Filter single-sided (soft-clip) breakends
        if bpc.breakpoint.breakend2.is_none() {
            continue;
        }

        if bpc.assembly_segments.len() != 2 {
            continue;
        }

        let be1 = &bpc.breakpoint.breakend1;
        let be2 = bpc.breakpoint.breakend2.as_ref().unwrap();

        breakends.push(BreakendInfo {
            be: be1.clone(),
            cluster_index,
            breakend_index: 0,
        });

        breakends.push(BreakendInfo {
            be: be2.clone(),
            cluster_index,
            breakend_index: 1,
        });
    }

    breakends.sort_by_key(|x| (x.be.clone(), x.cluster_index));

    let min_close_annotations = 2;

    // Look at sorted breakends pairwise to find close neighbors:
    let be_count = breakends.len();
    for be_index in 1..be_count {
        let last_bei = &breakends[be_index - 1];

        let mut best_next_be_index = None;
        let mut best_next_be_dist = 0.0;

        // check if this breakend is facing the right way and has any close breakend annotations:
        {
            let last_bpc = &clusters[last_bei.cluster_index];
            let last_be = {
                if last_bei.breakend_index == 0 {
                    &last_bpc.breakpoint.breakend1
                } else {
                    last_bpc.breakpoint.breakend2.as_ref().unwrap()
                }
            };

            if last_be.dir != BreakendDirection::RightAnchor {
                continue;
            }

            let last_be_neighbors = if last_bei.breakend_index == 0 {
                &last_bpc.breakend1_neighbor_observations
            } else {
                &last_bpc.breakend2_neighbor_observations
            };
            if last_be_neighbors.len() < min_close_annotations {
                continue;
            }

            #[allow(clippy::needless_range_loop)]
            for next_be_index in be_index..be_count {
                let next_bei = &breakends[next_be_index];

                // Check if the breakends are close:
                let next_bpc = &clusters[next_bei.cluster_index];
                let next_be = {
                    if next_bei.breakend_index == 0 {
                        &next_bpc.breakpoint.breakend1
                    } else {
                        next_bpc.breakpoint.breakend2.as_ref().unwrap()
                    }
                };

                let is_close = match get_segment_distance(&last_be.segment, &next_be.segment) {
                    Some(x) => x <= max_close_breakend_distance,
                    None => false,
                };

                if !is_close {
                    break;
                }

                // check if this breakend is facing the right way and has any close breakend annotations:
                if next_be.dir != BreakendDirection::LeftAnchor {
                    continue;
                }

                // Add the exception against chaining deletions together
                if last_bpc.breakpoint.is_indel() && next_bpc.breakpoint.is_indel() {
                    continue;
                }

                let next_be_neighbors = if next_bei.breakend_index == 0 {
                    &next_bpc.breakend1_neighbor_observations
                } else {
                    &next_bpc.breakend2_neighbor_observations
                };
                if next_be_neighbors.len() < min_close_annotations {
                    continue;
                }

                // check the average distance of observation breakends to this breakend
                let mut dist = 0.0;
                for be in last_be_neighbors {
                    dist += get_segment_distance(&be.segment, &next_be.segment).unwrap() as f64;
                }

                if best_next_be_index.is_none() || dist < best_next_be_dist {
                    best_next_be_dist = dist;
                    best_next_be_index = Some(next_be_index);
                }
            }
        }

        if best_next_be_index.is_none() {
            continue;
        }

        // Passed the full screen so mark neighbor breakends now
        //
        let next_bei = &breakends[best_next_be_index.unwrap()];

        fn get_neighbor<'a>(
            bei: &BreakendInfo,
            clusters: &'a mut [BreakpointCluster],
        ) -> &'a mut Option<NeighborBreakend> {
            let bpc = &mut clusters[bei.cluster_index];

            if bei.breakend_index == 0 {
                &mut bpc.breakend1_neighbor
            } else {
                &mut bpc.breakend2_neighbor
            }
        }
        if debug {
            eprintln!(
                "Annotating close breakend connection between:\nlast: {:?}\nnext: {:?}",
                last_bei, next_bei
            );
        }

        let mut connect_breakends = |bei1: &BreakendInfo, bei2: &BreakendInfo| {
            *get_neighbor(bei1, clusters) = {
                let cluster_index = bei2.cluster_index;
                let breakend_index = bei2.breakend_index;
                let breakpoint = clusters[cluster_index].breakpoint.clone();
                Some(NeighborBreakend {
                    cluster_index,
                    breakend_index,
                    breakpoint,
                })
            };
        };

        connect_breakends(last_bei, next_bei);
        connect_breakends(next_bei, last_bei);
    }
}

/// Run additional operations on the initial cluster set to prepare for refinement
///
/// This includes:
/// 1. a secondary clustering of overlapping single-region clusters
/// 2. a large insertion candidate routine
///
fn process_clusters(
    settings: &cli::DiscoverSettings,
    cluster_settings: &ClusterSettings,
    assembly_read_flank_size: usize,
    chrom_list: &ChromList,
    clusters: &mut [BreakpointCluster],
    cluster_stats: &mut ClusterStats,
) {
    // Add assembly regions to clusters
    for bpc in clusters.iter_mut() {
        bpc.assembly_segments = get_assembly_segments(assembly_read_flank_size, &bpc.breakpoint);
    }

    // Get some simple cluster stats prior to consolidation:
    for bpc in clusters.iter() {
        if bpc.assembly_segments.len() == 1 {
            cluster_stats.single_region_candidate_count += 1;
        } else {
            cluster_stats.multi_region_candidate_count += 1;
        }
    }

    // Further cluster small indel candidates to keep refinement regions from overlapping
    cluster_stats.consolidated_single_region_candidate_count =
        consolidate_indel_candidates(cluster_settings, chrom_list, clusters);

    // Pair soft-clip breakpoints into large-insertion candidates
    cluster_stats.large_insertion_candidate_count = if !settings.disable_large_insertions {
        find_large_insertion_candidates(
            chrom_list,
            assembly_read_flank_size,
            cluster_settings.max_large_insertion_pairing_dist,
            cluster_settings.max_large_insertion_deletion_pairing_dist,
            clusters,
        )
    } else {
        0
    };

    annotate_close_breakends(settings.max_close_breakend_distance, clusters);

    write_breakpoint_clusters_to_bed(&settings.output_dir, chrom_list, clusters);
}

struct ClusterSettings {
    /// Expand all indel clusters by this size to determine if they overlap with an adjacent cluster
    ///
    /// Note both indel clusters will be expanded, so those within range*2 will consolidate
    ///
    indel_cluster_consolidation_range: usize,

    /// The maximum size that single assembly regions can be joined together into during the
    /// consolidation step
    max_consolidated_assembly_region_size: i64,

    /// Max distance between left and right-facing soft-clipping regions that will be converted into
    /// a large insertion candidate region
    ///
    /// Note this refers to the typical separation due to breakend homology, where the left-clipped
    /// bases are separated on the left-side and the right-clipped bases on the right:
    ///
    ///  Pairing Dist:         |-----|
    ///  Left-clipped: SSSSSSSSMMMMMMMMMMMMMMMMMMM
    /// Right-clipped: MMMMMMMMMMMMMMMSSSSSSSSSSSS
    ///
    max_large_insertion_pairing_dist: i64,

    /// Max distance between left and right-facing soft-clipping regions that will be converted into
    /// a large insertion candidate region
    ///
    /// Note this refers to the less common separation due to a breakend deletion, where the
    /// left-clipped bases are separated on the right-side and the right-clipped bases on the left:
    ///
    ///  Pairing Dist:         |-----|
    ///  Left-clipped: SSSSSSSSSSSSSSMMMMMMMMMMMMM
    /// Right-clipped: MMMMMMMMMSSSSSSSSSSSSSSSSSS
    ///
    max_large_insertion_deletion_pairing_dist: i64,
}

impl ClusterSettings {
    fn new() -> Self {
        Self {
            indel_cluster_consolidation_range: 300,
            max_consolidated_assembly_region_size: 8_000,
            max_large_insertion_pairing_dist: 500,
            max_large_insertion_deletion_pairing_dist: 200,
        }
    }
}

/// Cluster single breakpoint observations into canidate SV breakpoints for assembly
///
/// This includes an initial simple clustering step, and follow-on steps to indentify
/// small indel clusters and large insertions from soft-clip evidence.
///
/// Returns a 2-tuple of (1) breakpoint clusters and (2) clustering statistics
#[allow(clippy::field_reassign_with_default)]
pub fn process_breakpoint_clusters(
    settings: &cli::DiscoverSettings,
    assembly_read_flank_size: usize,
    chrom_list: &ChromList,
    genome_break_observations: &mut BreakObservations,
) -> (Vec<BreakpointCluster>, ClusterStats) {
    // Temp debug output for breakpoint scan:
    let debug = false;

    let cluster_settings = ClusterSettings::new();

    // Cluster breakpoint evidence:
    let mut cluster_stats = ClusterStats::default();
    cluster_stats.total_breakpoint_observation_count = genome_break_observations.singles.len();

    if debug {
        // Write first debug_output breakpoint observations to a debug output file:
        let debug_output = 100;
        let breaks = &genome_break_observations.singles;
        debug_bpos(&settings.output_dir, breaks, debug_output);
    }

    let mut breakpoint_clusters =
        cluster::cluster_breakpoint_observations(genome_break_observations);

    cluster_stats.total_breakpoint_cluster_count = breakpoint_clusters.len();

    if debug {
        // Write first debug_output clusters to a debug output file:
        let debug_output = 100;
        debug_clusters(&settings.output_dir, &breakpoint_clusters, debug_output);
    }

    process_clusters(
        settings,
        &cluster_settings,
        assembly_read_flank_size,
        chrom_list,
        &mut breakpoint_clusters,
        &mut cluster_stats,
    );

    (breakpoint_clusters, cluster_stats)
}
