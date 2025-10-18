use rust_vc_utils::GenomeSegment;
use rust_vc_utils::int_range::get_recip_overlap;

use crate::breakpoint::{BreakendDirection, BreakendNeighbor, InsertInfo};
use crate::refine_sv::Genotype;
use crate::sv_group::SVGroup;
use crate::sv_id::{SVUniqueIdData, get_inv_id_label, get_sv_id_label};

struct CandBpInfo {
    /// Index of SV group in the current sorting of sv_group list
    sv_group_index: usize,
    id: SVUniqueIdData,
    dir: BreakendDirection,
    segment: GenomeSegment,
    gts: Vec<Option<Genotype>>,
    breakend1_neighbor: Option<BreakendNeighbor>,
    breakend2_neighbor: Option<BreakendNeighbor>,
}

/// Check two inverted breakpoints to see if they are likely to comprise a balanced inversion, return true if so
///
/// Require:
/// 1. On same chromosome
/// 2. Complementary anchor directions
/// 3. Breakend partner dist is in acceptable range
/// 4. Reciprical overlap is in acceptable range
/// 5. If breakend phasing information exists, check that it doesn't conflict with an inversion pattern
///
fn is_inversion(inv_settings: &InversionSettings, c1: &CandBpInfo, c2: &CandBpInfo) -> bool {
    let debug = false;

    if c1.dir == c2.dir {
        return false;
    }
    if c1.segment.chrom_index != c2.segment.chrom_index {
        return false;
    }

    let start_be_dist = (c2.segment.range.start - c1.segment.range.start).abs();
    let end_be_dist = (c2.segment.range.end - c1.segment.range.end).abs();
    if std::cmp::max(start_be_dist, end_be_dist) > inv_settings.max_breakend_partner_dist {
        return false;
    }

    let recip_overlap = get_recip_overlap(&c1.segment.range, &c2.segment.range);
    if recip_overlap < inv_settings.min_breakpoint_span_recip_overlap {
        return false;
    }

    // Check if the downstream breakends are observed as adjacent neighbors on the same haplotype
    if let Some(x) = &c1.breakend1_neighbor
        && x.sample_index == c2.id.sample_index
        && x.cluster_index == c2.id.cluster_index
        && x.breakend_index == 0
    {
        if debug {
            eprintln!(
                "B1: c1 cluster {:?} and c2 cluster {:?} on the same haplotype, filtered from inversions",
                c1.id, c2.id
            );
        }
        return false;
    }

    // Check if the upstream breakends are observed as adjacent neighbors on the same haplotype
    if let Some(x) = &c1.breakend2_neighbor
        && x.sample_index == c2.id.sample_index
        && x.cluster_index == c2.id.cluster_index
        && x.breakend_index == 1
    {
        if debug {
            eprintln!(
                "B2: c1 cluster {:?} and c2 cluster {:?} on the same haplotype, filtered from inversions",
                c1.id, c2.id
            );
        }
        return false;
    }

    // We can also ask a less specific QC question, does the inversion breakend have a neighboring breakend on the same
    // haplotype that comes from a cluster other than the inversion? If so, this is a general indication that it's part
    // of a more complex SV, and not usually described as an inversion.
    //
    // This filtration idea is too non-specific in its current form to be applied to all inversions, but for larger
    // inversions we're even more conservative with inversion annotation, attempting to restrict to very obvious and
    // clean examples only.
    //
    let min_span_for_neighbor_breakend_filter = 100_000;
    let min_span = std::cmp::min(c1.segment.range.size(), c2.segment.range.size());
    if min_span >= min_span_for_neighbor_breakend_filter {
        fn has_non_inv_breakend_neighbor(ca: &CandBpInfo, cb: &CandBpInfo) -> bool {
            if let Some(x) = &ca.breakend1_neighbor
                && x.cluster_index != cb.id.cluster_index
            {
                true
            } else if let Some(x) = &ca.breakend2_neighbor
                && x.cluster_index != cb.id.cluster_index
            {
                true
            } else {
                false
            }
        }

        if has_non_inv_breakend_neighbor(c1, c2) || has_non_inv_breakend_neighbor(c2, c1) {
            return false;
        }
    }

    true
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
fn mark_inv_bnd(sv_group: &mut SVGroup, inversion_id: &str) {
    assert_eq!(sv_group.refined_svs.len(), 1);
    sv_group.refined_svs[0].ext.inversion_id = Some(inversion_id.to_string());
}

/// All parameters influencing inversion annotation
pub struct InversionSettings {
    /// Max span for either of the two breakpoints comprising the inversion
    ///
    /// If None, no max span is enforced
    max_breakpoint_span: Option<i64>,

    /// Max distance between the two partner breakends on either end of the inversion
    max_breakend_partner_dist: i64,

    /// The 'left-left' span and 'right-right' breakpoint spans need to have at least this recip overlap to be
    /// classified as an inversion event:
    min_breakpoint_span_recip_overlap: f64,
}

impl Default for InversionSettings {
    fn default() -> Self {
        Self {
            max_breakpoint_span: None,
            max_breakend_partner_dist: 10_000,
            min_breakpoint_span_recip_overlap: 0.6,
        }
    }
}

/// Annotate likely balanced inversions in the input breakpoint set
///
/// This routine does not include any logic to ensure that breakends are on the same haplotype, it just provides a
/// simple pairing inference.
///
/// Does not assume sv_groups are sorted. Will only edit BNDs determined to be part of inversions and add inversion
/// records to sv_groups.
///
pub fn annotate_inversions(sv_groups: &mut Vec<SVGroup>) {
    // Summary of ops:
    // 1. Store list of candidate breakpoints
    // 2. Sort list
    // 3. Walk sorted list to find simple pairs meeting all criteria.

    // Hard code parameters here for now:
    let inv_settings = InversionSettings::default();

    let mut inv_cand_bps = Vec::new();

    // Collect a list of inverted intrachromasomal breakpoints meeting eligibility criteria for inversion
    // Criteria:
    // 1. Not max scoring depth in any sample
    // 2. Inverted intrachromosomal breakpoint
    // 3. Optionally, breakpoint span is less than max_span
    // 4. At least one sample must have a variant genotype for the breakpoint
    //
    for (sv_group_index, sv_group) in sv_groups.iter().enumerate() {
        // First filter down to just breakends found on a single chromosome with inverted orieintation:
        if sv_group.is_single_region() {
            continue;
        }
        assert_eq!(sv_group.refined_svs.len(), 1);
        let rsv = &sv_group.refined_svs[0];
        let score_info = &rsv.score;

        let is_max_scoring_depth = score_info.samples.iter().any(|x| x.is_max_scoring_depth);
        if is_max_scoring_depth {
            continue;
        }

        let bp = &rsv.bp;
        let be1 = &bp.breakend1;
        let be2 = bp.breakend2.as_ref().unwrap();

        if (be1.segment.chrom_index != be2.segment.chrom_index) || (be1.dir != be2.dir) {
            continue;
        }

        // transform breakend segment into breakpoint span segment, removing the encoding of breakend homology range
        let mut segment = be1.segment.clone();
        segment.range.end = be2.segment.range.start + 1;

        if let Some(max_span) = inv_settings.max_breakpoint_span
            && segment.range.size() > max_span
        {
            continue;
        }

        // For each sample where this breakpoint is present, record the genotype
        let gts = score_info
            .samples
            .iter()
            .map(|x| {
                if x.on_sample_haplotype {
                    x.shared.gt.clone()
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        // Check that the breakpoint has a variant haplotype in at least one sample:
        let any_variants = gts
            .iter()
            .any(|x| x.as_ref().is_some_and(|y| *y != Genotype::Ref));
        if !any_variants {
            continue;
        }

        inv_cand_bps.push(CandBpInfo {
            sv_group_index,
            id: rsv.id.clone(),
            dir: be1.dir,
            segment,
            gts,
            breakend1_neighbor: bp.breakend1_neighbor.clone(),
            breakend2_neighbor: bp.breakend2_neighbor.clone(),
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

    /// All information stored for an inversion after a pair of inverted breakpoints has been found to
    /// meet all eligibility criteria
    ///
    struct InvBpInfo {
        sv_group_index1: usize,
        sv_group_index2: usize,
        inversion_id: String,

        /// True if the two breakpoints have conflicting genotypes in most samples
        is_conflicting_gt: bool,
    }
    let mut inv_pairs = Vec::new();

    // Find pairs in sorted candidates
    //
    // Keep overlap criteria super simple -- any cases of complex interactions
    // with close networks of interacting breakends shouldn't be in here anyway
    //
    let mut last_cand: Option<CandBpInfo> = None;
    for cand in inv_cand_bps {
        if let Some(base_cand) = &last_cand {
            if is_inversion(&inv_settings, base_cand, &cand) {
                let is_conflicting_gt = !do_most_inversion_genotypes_match(base_cand, &cand);
                let inversion_id = get_inv_id_label(&base_cand.id);

                // Join cand and last cand into an inversion:
                inv_pairs.push(InvBpInfo {
                    sv_group_index1: base_cand.sv_group_index,
                    sv_group_index2: cand.sv_group_index,
                    inversion_id,
                    is_conflicting_gt,
                });
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
        mark_inv_bnd(&mut sv_groups[ipair.sv_group_index1], &ipair.inversion_id);
        mark_inv_bnd(&mut sv_groups[ipair.sv_group_index2], &ipair.inversion_id);
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
