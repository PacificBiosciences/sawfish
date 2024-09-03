use crate::breakpoint::{BreakendDirection, InsertInfo};
use crate::genome_segment::GenomeSegment;
use crate::int_range::IntRange;
use crate::refine_sv::Genotype;
use crate::sv_group::SVGroup;
use crate::sv_id::{get_sv_id_label, SVUniqueIdData};

fn get_recip_overlap(mut r1: IntRange, mut r2: IntRange) -> f64 {
    let min_span = 100;

    r1.end = std::cmp::max(r1.end, r1.start + min_span);
    r2.end = std::cmp::max(r2.end, r2.start + min_span);

    let olap = std::cmp::max(
        std::cmp::min(r1.end, r2.end) - std::cmp::max(r2.start, r1.start),
        0,
    );
    let span1 = r1.end - r1.start;
    let span2 = r2.end - r2.start;
    let span = std::cmp::min(span1, span2);

    olap as f64 / span as f64
}

struct CandBpInfo {
    /// Index of SV group in the current sorting of sv_group list
    sv_group_index: usize,
    id: SVUniqueIdData,
    dir: BreakendDirection,
    segment: GenomeSegment,
    gts: Vec<Option<Genotype>>,
}

/// Compare 2 inverted breakpoints to see if they form an inversion pair, return true if so
///
/// Require:
/// 1. On same chromosome
/// 2. Complementary anchor directions
/// 3. Reciprical overlap >= 0.8
///
fn is_inversion(c1: &CandBpInfo, c2: &CandBpInfo) -> bool {
    if c1.dir == c2.dir {
        return false;
    }
    if c1.segment.chrom_index != c2.segment.chrom_index {
        return false;
    }

    let min_inv_span_recip_overlap = 0.8;

    let recip_overlap = get_recip_overlap(c1.segment.range.clone(), c2.segment.range.clone());
    recip_overlap >= min_inv_span_recip_overlap
}

/// Duplicated breakpoints should be rare but can occur
fn is_duplicated_inv_breakpoint(c1: &CandBpInfo, c2: &CandBpInfo) -> bool {
    c1.segment == c2.segment && c1.dir == c2.dir
}

/// In a real inversion, the genotypes of the two breakpoints should be the same, check if this is the case for the
/// majority of samples
///
fn do_most_inversion_genotypes_match(c1: &CandBpInfo, c2: &CandBpInfo) -> bool {
    let mut same = 0;
    let mut total = 0;
    for (g1, g2) in c1.gts.iter().zip(c2.gts.iter()) {
        if g1.is_some() || g2.is_some() {
            total += 1;
            if g1 == g2 {
                same += 1;
            }
        }
    }
    let same_frac = if total > 0 {
        same as f64 / total as f64
    } else {
        1.0
    };

    same_frac >= 0.5
}

/// Given the SV group for a BND record marked as INV associated, update it with information needed to change the VCF output record per INV output policy.
///
fn mark_inv_bnd(sv_group: &mut SVGroup, inv_index: usize) {
    assert_eq!(sv_group.refined_svs.len(), 1);
    sv_group.refined_svs[0].ext.inversion_id = Some(inv_index);
}

/// Find likely simple inversions in the input breakend set
///
/// This routine does not include any logic to ensure that breakends are on the same haplotype, it just provides a simple pairing
///
/// Assumes sv_groups are not sorted
///
pub fn find_inversions(sv_groups: &mut Vec<SVGroup>) {
    // Store list of candidate breakpoints
    // Sort list
    // Walk sorted list to find simple pairs meeting all criteria.
    //
    let max_bp_span = 100_000;

    let mut inv_cand_bps = Vec::new();

    // Collect a list of inverted intrachromasomal breakpoints meeting eligibility criteria for inversion
    //
    for (sv_group_index, sv_group) in sv_groups.iter().enumerate() {
        // First filter down to just breakends found on a single chromosome with inverted orieintation:
        if sv_group.is_single_region() {
            continue;
        }
        assert_eq!(sv_group.refined_svs.len(), 1);
        let rsv = &sv_group.refined_svs[0];
        let bp = &rsv.bp;
        let be1 = &bp.breakend1;
        let be2 = bp.breakend2.as_ref().unwrap();

        if (be1.segment.chrom_index != be2.segment.chrom_index) || (be1.dir != be2.dir) {
            continue;
        }

        let span = be2.segment.range.start - be1.segment.range.start;

        if span > max_bp_span {
            continue;
        }

        // For each sample where this breakpoint is present, record the genotype
        let score_info = &rsv.score;
        let gts = score_info
            .samples
            .iter()
            .map(|x| {
                if x.on_sample_haplotype {
                    x.gt.clone()
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        // transform breakend segment into breakpoint span segment, removing the encoding of breakend homology range
        let mut segment = be1.segment.clone();
        segment.range.end = be2.segment.range.start + 1;

        inv_cand_bps.push(CandBpInfo {
            sv_group_index,
            id: rsv.id.clone(),
            dir: be1.dir,
            segment,
            gts,
        });
    }

    // Sort candidate inversion breakpoints
    //
    // Including direction in the sort ensures that we always start with the left-anchored end of the inversion given the (fairly common) case
    // that the starting left and right anchored inversion breakends are on the same position.
    //
    // Including id in the sort ensures determinism in the case of a duplicate breakpoint call. This should be very rare but does happen.
    //
    // In the VCF output deduping step the variant ID (as a string) is the final step to unambiguous sort duplicates as well, so this should make the
    // sort logic for inversions match to the VCF deduping routines, so that duplicates in the same order will be (1) skipped form inversion processing and
    // (2) removed from VCF output.
    //
    inv_cand_bps.sort_by_key(|x| {
        (
            x.segment.chrom_index,
            x.segment.range.start,
            x.dir,
            get_sv_id_label(&x.id),
        )
    });

    struct InvBpInfo {
        sv_group_index1: usize,
        sv_group_index2: usize,
        inv_index: usize,

        /// True if the two breakpoints have conflicting genotypes in most samples
        is_conflicting_gt: bool,
    }
    let mut inv_pairs = Vec::new();
    let mut inv_index = 0;

    // Find pairs in sorted candidates
    //
    // Keep overlap criteria super simple -- any cases of complex interactions
    // with close networks of interacting breakends shouldn't be in here anyway
    //
    let mut last_cand: Option<CandBpInfo> = None;
    for cand in inv_cand_bps {
        if let Some(base_cand) = &last_cand {
            if is_inversion(base_cand, &cand) {
                let is_conflicting_gt = !do_most_inversion_genotypes_match(base_cand, &cand);

                // Join cand and last cand into an inversion:
                inv_pairs.push(InvBpInfo {
                    sv_group_index1: base_cand.sv_group_index,
                    sv_group_index2: cand.sv_group_index,
                    inv_index,
                    is_conflicting_gt,
                });
                inv_index += 1;
                last_cand = None;
                continue;
            }

            // This should occur very rarely, but it is possible
            if is_duplicated_inv_breakpoint(base_cand, &cand) {
                continue;
            }
        }
        last_cand = Some(cand);
    }

    // Update SV groups to
    // 1. Mark both inversion BND record pairs as filtered (annotation at the RefinedSV level)
    // 2. Mark both inversion BND record pairs with the inversion event tag (annotation at the RefinedSV level)
    // 3. Add new inversion SV group, with the inversion event tag - How to fit this into the refinedSV scheme --  make up a fake breakpoint?
    //
    for ipair in inv_pairs.iter() {
        mark_inv_bnd(&mut sv_groups[ipair.sv_group_index1], ipair.inv_index);
        mark_inv_bnd(&mut sv_groups[ipair.sv_group_index2], ipair.inv_index);
        let mut inv_sv_group = sv_groups[ipair.sv_group_index1].clone();
        let rsv = &mut inv_sv_group.refined_svs[0];
        let be2 = rsv.bp.breakend2.as_mut().unwrap();
        rsv.ext.is_inversion = true;
        rsv.ext.is_conflicting_gt = ipair.is_conflicting_gt;

        // Set the end of the inversion to the furthest downstream breakend
        let sv_group2 = &sv_groups[ipair.sv_group_index2];
        let rsv2 = &sv_group2.refined_svs[0];
        let r2be2 = rsv2.bp.breakend2.as_ref().unwrap();
        if r2be2.segment.range.start > be2.segment.range.start {
            be2.segment.range = r2be2.segment.range.clone();
        }

        // Remove any breakend homology
        rsv.bp.breakend1.segment.range.end = rsv.bp.breakend1.segment.range.start + 1;
        be2.segment.range.end = be2.segment.range.start + 1;
        rsv.breakend1_homology_seq.clear();

        // Remove any breakpoint insertion
        rsv.bp.insert_info = InsertInfo::NoInsert;

        sv_groups.push(inv_sv_group);
    }
}
