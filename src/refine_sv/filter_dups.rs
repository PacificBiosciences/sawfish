//! Filter duplications already represented by an insertion assembly
//!
//! Only multi-region duplication SV groups under MAX_FILTERED_DUP_SIZE can be filtered, based on
//! their overlap with insertions in single-region SV groups. This allows the filter to run on
//! small keys extracted from each SV group as it is refined, so that all SV groups which can't be
//! filtered can be finalized and dropped immediately.
//!
//! Strategy:
//! 1. Record the segment of each duplication under MAX_FILTERED_DUP_SIZE
//! 2. Check each insertion for overlap with these duplications, marking intersecting
//!    duplications as filtered where the insertion is sufficiently large
//! 3. Remove filtered duplications
//!

use rust_vc_utils::util::drop_true;
use rust_vc_utils::{ChromList, GenomeSegment, IntRange};

use crate::breakpoint::{VcfSVType, get_breakpoint_span, get_breakpoint_vcf_sv_type};
use crate::genome_segment_utils::is_strict_order;
use crate::sv_group::SVGroup;

const MAX_FILTERED_DUP_SIZE: usize = 50_000;
const MAX_INSERTION_DUP_DISTANCE: i64 = 100;
const INSERTION_SIZE_EXPANSION_FACTOR: f64 = 1.1;

/// Key for a duplication SV group which could be filtered
pub struct DupFilterCandidateKey {
    cluster_index: usize,
    segment: GenomeSegment,
}

/// Key for an insertion which could cause a duplication to be filtered
pub struct InsertionKey {
    cluster_index: usize,

    /// Insertion breakend1 segment, expanded by MAX_INSERTION_DUP_DISTANCE
    segment: GenomeSegment,

    insert_size: usize,
}

/// If sv_group is a duplication which could be filtered, return its filter key
pub fn get_dup_filter_candidate_key(sv_group: &SVGroup) -> Option<DupFilterCandidateKey> {
    if sv_group.is_single_region() {
        return None;
    }

    // Multi region sv_groups should have a single sv candidate:
    assert_eq!(sv_group.refined_svs.len(), 1);
    let sv = &sv_group.refined_svs[0];
    let sv_type = get_breakpoint_vcf_sv_type(&sv.bp);
    if sv_type != VcfSVType::Duplication {
        return None;
    }

    let sv_span = get_breakpoint_span(&sv.bp).unwrap();
    if sv_span > MAX_FILTERED_DUP_SIZE {
        return None;
    }

    let start = sv.bp.breakend1.segment.range.start;
    let segment = GenomeSegment {
        chrom_index: sv.bp.breakend1.segment.chrom_index,
        range: IntRange {
            start,
            end: start + sv_span as i64,
        },
    };
    Some(DupFilterCandidateKey {
        cluster_index: sv_group.refined_svs[0].id.cluster_index,
        segment,
    })
}

/// Add keys for any insertions in sv_group which could cause a duplication to be filtered
pub fn add_insertion_keys(
    chrom_list: &ChromList,
    sv_group: &SVGroup,
    insertion_keys: &mut Vec<InsertionKey>,
) {
    if !sv_group.is_single_region() {
        return;
    }

    let cluster_index = sv_group.refined_svs.first().unwrap().id.cluster_index;
    for sv in sv_group.refined_svs.iter() {
        let sv_type = get_breakpoint_vcf_sv_type(&sv.bp);
        if sv_type != VcfSVType::Insertion {
            continue;
        }

        let mut segment = sv.bp.breakend1.segment.clone();
        segment.expand_by(chrom_list, MAX_INSERTION_DUP_DISTANCE);
        insertion_keys.push(InsertionKey {
            cluster_index,
            segment,
            insert_size: sv.bp.insert_info.size(),
        });
    }
}

/// Find all duplications made redundant by an insertion
///
/// Both key lists are sorted back into cluster index order, which approximates genomic order. The
/// sort is stable so that insertions from the same SV group stay in their original order.
///
/// Return a boolean vector set to true for each filtered duplication, in the sorted order of
/// dup_candidates
///
fn get_dup_filter_mask(
    dup_candidates: &mut [DupFilterCandidateKey],
    insertion_keys: &mut [InsertionKey],
) -> Vec<bool> {
    dup_candidates.sort_by_key(|x| x.cluster_index);
    insertion_keys.sort_by_key(|x| x.cluster_index);

    let mut filter_mask = vec![false; dup_candidates.len()];

    let mut dup_candidates_head = 0;
    for insertion in insertion_keys.iter() {
        // Iterate through dup_candidates until we encounter any cases that could intersect the
        // insertion:
        while dup_candidates_head < dup_candidates.len() {
            let candidate_segment = &dup_candidates[dup_candidates_head].segment;
            if is_strict_order(&insertion.segment, candidate_segment) {
                break;
            } else if !is_strict_order(candidate_segment, &insertion.segment) {
                // The two segments must intersect if we've gotten here. Now test if the insertion
                // is large enough. The insertion must be large enough to be at least the size of
                // the duplication, allowing for a small fudge factor.
                //
                if insertion.insert_size as f64 * INSERTION_SIZE_EXPANSION_FACTOR
                    >= candidate_segment.range.size() as f64
                {
                    filter_mask[dup_candidates_head] = true;
                }
            }

            dup_candidates_head += 1;
        }
    }
    filter_mask
}

/// Remove duplication SV groups made redundant by an insertion
///
/// # Arguments
/// * `dup_sv_groups` - All SV groups with a duplication filter candidate key, in the same order
///   as `dup_candidates`
///
/// Returns the remaining duplication SV groups and the number of duplications filtered
///
pub fn filter_redundant_duplications(
    mut dup_sv_groups: Vec<SVGroup>,
    mut dup_candidates: Vec<DupFilterCandidateKey>,
    mut insertion_keys: Vec<InsertionKey>,
) -> (Vec<SVGroup>, usize) {
    assert_eq!(dup_sv_groups.len(), dup_candidates.len());

    // Sort dup_sv_groups into the same order as dup_candidates will be sorted to:
    dup_sv_groups.sort_by_key(|x| x.refined_svs[0].id.cluster_index);

    let filter_mask = get_dup_filter_mask(&mut dup_candidates, &mut insertion_keys);

    let input_size = dup_sv_groups.len();
    drop_true(&mut dup_sv_groups, &filter_mask);
    let output_size = dup_sv_groups.len();

    (dup_sv_groups, input_size - output_size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::breakpoint::{Breakend, BreakendDirection, Breakpoint, InsertInfo};
    use crate::refine_sv::{RefinedSV, RefinedSVExt};
    use crate::sv_id::SVUniqueIdData;

    fn breakend(chrom_index: usize, pos: i64, dir: BreakendDirection) -> Breakend {
        Breakend {
            segment: GenomeSegment {
                chrom_index,
                range: IntRange::from_int(pos),
            },
            dir,
        }
    }

    fn refined_sv(cluster_index: usize, bp: Breakpoint) -> RefinedSV {
        RefinedSV {
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
        }
    }

    fn sv_group(region_count: usize, refined_svs: Vec<RefinedSV>) -> SVGroup {
        SVGroup {
            group_regions: vec![GenomeSegment::default(); region_count],
            group_haplotypes: Vec::new(),
            sample_haplotype_list: Vec::new(),
            sample_expected_cn_info: Vec::new(),
            sv_haplotype_map: Vec::new(),
            refined_svs,
        }
    }

    /// Multi-region duplication of [start, start + span)
    fn dup_group(cluster_index: usize, chrom_index: usize, start: i64, span: i64) -> SVGroup {
        let bp = Breakpoint {
            breakend1: breakend(chrom_index, start, BreakendDirection::RightAnchor),
            breakend2: Some(breakend(
                chrom_index,
                start + span,
                BreakendDirection::LeftAnchor,
            )),
            ..Default::default()
        };
        sv_group(2, vec![refined_sv(cluster_index, bp)])
    }

    /// Single-region group with one insertion, or a deletion if `size` is 0
    fn insertion_group(cluster_index: usize, chrom_index: usize, pos: i64, size: usize) -> SVGroup {
        let bp = Breakpoint {
            breakend1: breakend(chrom_index, pos, BreakendDirection::LeftAnchor),
            breakend2: Some(breakend(
                chrom_index,
                pos + if size == 0 { 50 } else { 0 },
                BreakendDirection::RightAnchor,
            )),
            insert_info: InsertInfo::Seq(vec![b'A'; size]),
            ..Default::default()
        };
        sv_group(1, vec![refined_sv(cluster_index, bp)])
    }

    fn chrom_list() -> ChromList {
        let mut chrom_list = ChromList::default();
        chrom_list.add_chrom("chr1", 1_000_000);
        chrom_list.add_chrom("chr2", 1_000_000);
        chrom_list
    }

    /// Cluster indices of the duplications that remain, and the number filtered
    fn run_filter(dups: Vec<SVGroup>, insertions: Vec<SVGroup>) -> (Vec<usize>, usize) {
        let chrom_list = chrom_list();
        let mut dup_candidates = Vec::new();
        let mut dup_sv_groups = Vec::new();
        for group in dups {
            if let Some(key) = get_dup_filter_candidate_key(&group) {
                dup_candidates.push(key);
                dup_sv_groups.push(group);
            }
        }
        let mut insertion_keys = Vec::new();
        for group in insertions.iter() {
            add_insertion_keys(&chrom_list, group, &mut insertion_keys);
        }
        let (remaining, filtered) =
            filter_redundant_duplications(dup_sv_groups, dup_candidates, insertion_keys);
        (
            remaining
                .iter()
                .map(|x| x.refined_svs[0].id.cluster_index)
                .collect(),
            filtered,
        )
    }

    #[test]
    fn test_dup_filter_candidate_and_insertion_keys() {
        let chrom_list = chrom_list();

        let key = get_dup_filter_candidate_key(&dup_group(7, 1, 1000, 500)).unwrap();
        assert_eq!(key.cluster_index, 7);
        assert_eq!(key.segment.chrom_index, 1);
        assert_eq!(key.segment.range.start, 1000);
        assert_eq!(key.segment.range.end, 1500);

        // Too large, single-region, and not a duplication
        assert!(get_dup_filter_candidate_key(&dup_group(1, 0, 1000, 50_001)).is_none());
        assert!(get_dup_filter_candidate_key(&insertion_group(1, 0, 1000, 80)).is_none());
        let mut deletion_like = dup_group(1, 0, 1000, 500);
        deletion_like.refined_svs[0].bp.breakend1.dir = BreakendDirection::LeftAnchor;
        deletion_like.refined_svs[0]
            .bp
            .breakend2
            .as_mut()
            .unwrap()
            .dir = BreakendDirection::RightAnchor;
        assert!(get_dup_filter_candidate_key(&deletion_like).is_none());

        // Insertions only come from single-region groups, and are expanded
        let mut keys = Vec::new();
        add_insertion_keys(&chrom_list, &insertion_group(3, 0, 1000, 80), &mut keys);
        add_insertion_keys(&chrom_list, &insertion_group(4, 0, 2000, 0), &mut keys);
        add_insertion_keys(&chrom_list, &dup_group(5, 0, 1000, 80), &mut keys);
        assert_eq!(keys.len(), 1);
        assert_eq!(keys[0].cluster_index, 3);
        assert_eq!(keys[0].insert_size, 80);
        assert_eq!(keys[0].segment.range.start, 900);
        assert_eq!(keys[0].segment.range.end, 1101);
    }

    #[test]
    fn test_filter_redundant_duplications() {
        let dups = || {
            vec![
                // Matched by a large enough insertion
                dup_group(10, 0, 10_000, 1000),
                // Insertion at the same place is too small
                dup_group(20, 0, 50_000, 1000),
                // Insertion at the same place is just large enough (1100 * 1.1 >= 1210)
                dup_group(30, 0, 70_000, 1210),
                // Same position on a different chromosome
                dup_group(40, 1, 10_000, 1000),
                // Insertion more than 100 bases away
                dup_group(50, 0, 90_000, 500),
                // Too large to be filtered at all, so it never reaches the filter
                dup_group(60, 0, 100_000, 60_000),
                // No insertion
                dup_group(70, 0, 300_000, 400),
            ]
        };
        let insertions = || {
            vec![
                insertion_group(11, 0, 10_050, 1100),
                insertion_group(21, 0, 50_500, 900),
                insertion_group(31, 0, 70_000, 1100),
                insertion_group(51, 0, 90_700, 5000),
                insertion_group(61, 0, 120_000, 60_000),
            ]
        };

        let (remaining, filtered) = run_filter(dups(), insertions());
        assert_eq!(filtered, 2);
        // Remaining duplications are in cluster index order
        assert_eq!(remaining, vec![20, 40, 50, 70]);

        // The result doesn't depend on the order the groups were found in
        let mut reversed_dups = dups();
        reversed_dups.reverse();
        let mut reversed_insertions = insertions();
        reversed_insertions.reverse();
        assert_eq!(
            run_filter(reversed_dups, reversed_insertions),
            (remaining, filtered)
        );

        // Without insertions nothing is filtered
        assert_eq!(
            run_filter(dups(), Vec::new()),
            (vec![10, 20, 30, 40, 50, 70], 0)
        );
    }
}
