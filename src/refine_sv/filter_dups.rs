use rust_vc_utils::ChromList;

use crate::breakpoint::{VcfSVType, get_breakpoint_span, get_breakpoint_vcf_sv_type};
use crate::genome_segment::{GenomeSegment, IntRange, is_strict_order};
use crate::sv_group::SVGroup;
use crate::utils::drop_true;

struct DupInfo {
    sv_group_index: usize,
    segment: GenomeSegment,
}

fn get_dup_filter_candandidates(sv_groups: &[SVGroup]) -> Vec<DupInfo> {
    let max_filtered_dup_size = 50_000;

    let mut filter_candidates = Vec::new();

    for (sv_group_index, sv_group) in sv_groups
        .iter()
        .enumerate()
        .filter(|(_, x)| !x.is_single_region())
    {
        // Multi region sv_groups should have a single sv candidate:
        assert_eq!(sv_group.refined_svs.len(), 1);
        let sv = &sv_group.refined_svs[0];
        let sv_type = get_breakpoint_vcf_sv_type(&sv.bp);
        if sv_type != VcfSVType::Duplication {
            continue;
        }

        let sv_span = get_breakpoint_span(&sv.bp).unwrap();
        if sv_span > max_filtered_dup_size {
            continue;
        }

        let start = sv.bp.breakend1.segment.range.start;
        let segment = GenomeSegment {
            chrom_index: sv.bp.breakend1.segment.chrom_index,
            range: IntRange {
                start,
                end: start + sv_span as i64,
            },
        };
        filter_candidates.push(DupInfo {
            sv_group_index,
            segment,
        });
    }
    filter_candidates
}

/// Find all sv groups corresponding to redundant duplications
///
/// Return a boolean vector set to true for each filtered entry
fn get_sv_group_filter_mask(
    chrom_list: &ChromList,
    sv_groups: &[SVGroup],
    filter_candidates: &[DupInfo],
) -> Vec<bool> {
    let mut sv_group_filter_mask = vec![false; sv_groups.len()];
    let max_insertion_dup_distance = 100;
    let insertion_size_expansion_factor = 1.1;

    let mut filter_candidates_head = 0;
    for sv_group in sv_groups.iter().filter(|x| x.is_single_region()) {
        for sv in sv_group.refined_svs.iter() {
            let sv_type = get_breakpoint_vcf_sv_type(&sv.bp);
            if sv_type != VcfSVType::Insertion {
                continue;
            }

            let mut segment = sv.bp.breakend1.segment.clone();
            segment.expand_by(chrom_list, max_insertion_dup_distance);

            // Iterate through filter_candidates until we encounter any cases that could intersect sv:
            while filter_candidates_head < filter_candidates.len() {
                let candidate_segment = &filter_candidates[filter_candidates_head].segment;
                if is_strict_order(&segment, candidate_segment) {
                    break;
                } else if !is_strict_order(candidate_segment, &segment) {
                    // The two segments must intersect if we've gotten here. Now test if the insertion is large enough. The insertion must be large enough
                    // to be at least the size of the deletion, allowing for a small fudge factor.
                    //
                    let ins_size = sv.bp.insert_info.size();
                    if ins_size as f64 * insertion_size_expansion_factor
                        >= candidate_segment.range.size() as f64
                    {
                        sv_group_filter_mask
                            [filter_candidates[filter_candidates_head].sv_group_index] = true;
                    }
                }

                filter_candidates_head += 1;
            }
        }
    }
    sv_group_filter_mask
}

/// Iterate through assembled sv candidates to find duplications already represented by an insertion assembly
///
/// Returns the number of duplications filtered
///
pub fn filter_redundant_duplications(
    chrom_list: &ChromList,
    sv_groups: &mut Vec<SVGroup>,
) -> usize {
    // Strategy
    // 1.  Iterate through candidates, and record all duplication segments less than 50kb
    // 2.  Iterate through candidates, check each insertion for duplication overlap.
    //       Mark intersections where insertion is sufficiently large as filtered.
    // 3.  Remove filtered duplications

    // Sort sv groups back into cluster index order, which approximates genomic order:
    //
    sv_groups.sort_by_key(|x| x.refined_svs.first().unwrap().id.cluster_index);

    let filter_candidates = get_dup_filter_candandidates(sv_groups);
    let sv_group_filter_mask = get_sv_group_filter_mask(chrom_list, sv_groups, &filter_candidates);

    let input_size = sv_groups.len();
    drop_true(sv_groups, &sv_group_filter_mask);
    let output_size = sv_groups.len();

    input_size - output_size
}
