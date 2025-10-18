use rust_htslib::bam::record::Cigar;
use rust_vc_utils::bam_utils::cigar::{
    get_cigar_ref_offset, get_cigarseg_read_offset, get_cigarseg_ref_offset, update_read_pos,
    update_ref_pos,
};
use rust_vc_utils::indel_breakend_homology::get_indel_breakend_homology_info;
use rust_vc_utils::int_range::IntRange;

use super::TwoRegionAltHapInfo;
use crate::simple_alignment::SimpleAlignment;
use crate::utils::get_seq_pos_flanks;

/// Determine if the alignment to the alt haplotype supports the alt allele
///
/// Scan through cigar alignment and determine if (1) the alignment extends into both alt haplotype
/// regions and (2) if the alignment on both sides of the alt haplotype is of sufficient quality.
///
/// Sufficient quality currently means that there is a perfect match segment of at least size
/// min_assembly_edge_anchor on each side, with the assembly_high_quality_range region.
///
/// # Arguments
/// * `alt_hap_alignment` - Start position and cigar alignment of the contig to the alt hap sequence
///
/// Returns true if the alignment is anchored on both sides of the alt haplotype
///
fn does_alignment_support_alt_hap(
    min_assembly_edge_anchor: u32,
    alt_hap_alignment: &SimpleAlignment,
    alt_hap_segment_boundary1: i64,
    alt_hap_segment_boundary2: i64,
) -> bool {
    let mut left_segment_anchor = false;
    let mut right_segment_anchor = false;
    let mut alt_hap_pos = alt_hap_alignment.ref_offset;
    let mut contig_pos = 0;
    for c in alt_hap_alignment.cigar.iter() {
        use Cigar::*;
        // Test the anchor criteria
        if let Equal(len) = c
            && *len >= min_assembly_edge_anchor
        {
            // Determine if min size occurs within one region and if so which one(s)
            if (alt_hap_pos + min_assembly_edge_anchor as i64) < alt_hap_segment_boundary1 {
                left_segment_anchor = true;
            }

            if alt_hap_pos + *len as i64 - min_assembly_edge_anchor as i64
                >= alt_hap_segment_boundary2
            {
                right_segment_anchor = true;
            }
        }
        update_ref_pos(c, &mut alt_hap_pos);
        update_read_pos(c, &mut contig_pos, true);
    }

    left_segment_anchor && right_segment_anchor
}

/// Alignments to the left and right segments of the alt hap sequence
///
/// Each alignment is expressed in the coordinates of the left or right segment, so the reference
/// offset is wrt the left or right segment and not the alt-hap sequence as a whole. Similarly if
/// either segment has been revcomped compared to the reference this alignment is in the revcomped
/// orientation.
///
#[derive(Debug)]
pub struct AltHapLeftRightComponentAlignmentInfo {
    pub left_alignment: SimpleAlignment,
    pub right_alignment: SimpleAlignment,

    /// Range of homology offsets to each side of the left-to-right segment breakpoint
    ///
    /// Note these are not coordinates in ref or contig space, but rather offsets from the
    /// breakpoint locations implied in the left and right alignments
    pub left_right_breakend_homology_range: IntRange,

    pub left_right_breakend_homology_seq: Vec<u8>,
}

/// Determine if the alignment to the alt haplotype supports the alt allele, and if so process the
/// alignment information into two separate alignments to the left and right components of the alt
/// haplotype.
///
/// In addition to the left-right split, breakpoint homology and insertion sequence are found in this
/// routine.
///
/// Note the 'left' and 'right' components of the alt haplotype are still in a similar coordinate system to the
/// alt haplotype alignment coodinates, just with the spacer removed to yeild the corresponding left and right
/// sides of the alt haplotype, and with the alignments now expressed with respoect to those left and right
/// segments. Thus if reference region2 is rev-comped relative to reference in the alt-hap sequence it is still
/// rev-comped in the left-right result, and if reference region2 precedes region1 in the alt-hap sequence it
/// remains as the 'left' component in the left-right result.
///
/// Example:
///
/// # Input:
///
/// Reference:  TTACGTTTNNNTTACGTTT
/// Contig:        CG--------AC
///
/// AltHap SimpleAlignment:
/// Offset: 3 CIGAR: 2M8D2M
///
/// # Result:
///
/// Left SimpleAlignment
/// Offset: 3 CIGAR: 2M2S
///
/// Right SimpleAlignment
/// Offset: 2 CIGAR: 2S2M
///
/// The alt allele support checks here do not represent the complete set of filters applied to check
/// for alignment quality, these are rough checks applied early in the process to eliminate some of
/// the bigger problem cases.
///
/// # Arguments
///
/// * `assembly_contig_to_alt_hap_alignment` - Start position and cigar alignment of the contig to the
///   alt hap sequence
///
/// Return AltHapLeftRightComponentAlignmentInfo describing how the alignment is distributed over
/// the left and right segments of the alt haplotype if quality criteria are met, and None otherwise
///
pub fn transform_alt_hap_alignment_into_left_right_components(
    min_assembly_edge_anchor: u32,
    assembly_contig_to_alt_hap_alignment: &SimpleAlignment,
    alt_hap_info: &TwoRegionAltHapInfo,
    contig: &[u8],
) -> Option<AltHapLeftRightComponentAlignmentInfo> {
    let debug = false;
    let (alt_hap_segment_boundary1, alt_hap_segment_boundary2) = (
        alt_hap_info.spacer_range.start,
        alt_hap_info.spacer_range.end,
    );

    if debug {
        eprintln!("Starting transform_alt_hap_alignment_into_left_right_components");
        eprintln!(
            "alt_hap_segment_boundary 1/2 {alt_hap_segment_boundary1}/{alt_hap_segment_boundary2}"
        );
    }

    if !does_alignment_support_alt_hap(
        min_assembly_edge_anchor,
        assembly_contig_to_alt_hap_alignment,
        alt_hap_segment_boundary1,
        alt_hap_segment_boundary2,
    ) {
        if debug {
            eprintln!("failed does_alignment_support_alt_hap");
        }
        return None;
    }

    fn get_cigarseg_offsets(c: &Cigar) -> (i64, usize) {
        (
            get_cigarseg_ref_offset(c),
            get_cigarseg_read_offset(c, true),
        )
    }

    // The first pass through the contig alignment is to figure out the left segment cigar
    //
    let mut left_segment_cigar = Vec::new();
    let mut left_segment_soft_clip_len = 0;
    {
        let mut alt_hap_pos = assembly_contig_to_alt_hap_alignment.ref_offset;
        let mut in_segment = true;
        let mut last_cigar = None;
        for c in assembly_contig_to_alt_hap_alignment.cigar.iter() {
            let (ref_offset, read_offset) = get_cigarseg_offsets(c);

            if in_segment {
                // Test if we've transitioned out of the left segment, if so the remainder of the
                // read is soft-clipped
                //
                if (alt_hap_pos + ref_offset) > alt_hap_segment_boundary1 {
                    in_segment = false;

                    if let Some(last_cigar) = last_cigar {
                        let (last_ref_offset, last_read_offset) = get_cigarseg_offsets(last_cigar);
                        if last_ref_offset > 0 {
                            left_segment_cigar.push(*last_cigar);
                        } else {
                            left_segment_soft_clip_len += last_read_offset;
                        }
                    }
                } else if (alt_hap_pos + ref_offset) == alt_hap_segment_boundary1 {
                    // Only allow a cigar segment to end on the first position of the spacer if the
                    // corresponding ref segment is pinned:
                    if !alt_hap_info.left_boundary_pin {
                        if debug {
                            eprintln!("failed left_boundary_pin");
                        }
                        return None;
                    }
                    if let Some(last_cigar) = last_cigar {
                        left_segment_cigar.push(*last_cigar);
                    }
                } else if let Some(last_cigar) = last_cigar {
                    left_segment_cigar.push(*last_cigar);
                }
            } else {
                left_segment_soft_clip_len += read_offset;
            }

            alt_hap_pos += ref_offset;
            last_cigar = Some(c);
        }
        if left_segment_soft_clip_len > 0 {
            left_segment_cigar.push(Cigar::SoftClip(left_segment_soft_clip_len as u32));
        }
    }

    // Second pass through the contig alignment is to figure out the right side cigar:
    let mut right_segment_alt_hap_offset = None;
    let mut right_segment_cigar = Vec::new();
    let mut right_segment_soft_clip_len = 0;
    {
        let mut alt_hap_pos = assembly_contig_to_alt_hap_alignment.ref_offset;
        let mut in_segment = false;
        for c in assembly_contig_to_alt_hap_alignment.cigar.iter() {
            let (ref_offset, read_offset) = get_cigarseg_offsets(c);

            if !in_segment {
                // Test if we've transitioned into the right segment
                //
                if (alt_hap_pos + ref_offset) >= alt_hap_segment_boundary2 {
                    // This test is no longer valid with the extended_alt_hap_seq scheme
                    /*
                    // Only allow a cigar segment to end on the first position of the right segment
                    // if the corresponding ref segment is pinned:
                    if ((alt_hap_pos + ref_offset) == alt_hap_segment_boundary2)
                        && (!alt_hap_info.right_boundary_pin)
                    {
                        if debug {
                            eprintln!("failed right_boundary_pin");
                        }
                        return None;
                    }
                    */

                    // If spanning over from the left side, make sure this isn't a read-consuming
                    // segment. If it is, we assume a match block spanned the segment boundary,
                    // which shouldn't happen.
                    if (alt_hap_pos < alt_hap_segment_boundary1) && (read_offset > 0) {
                        if debug {
                            eprintln!("failed match-over-gap");
                        }
                        return None;
                    }

                    in_segment = true;
                } else {
                    // Don't allow cigar segments to end inside of the spacer
                    if (alt_hap_pos + ref_offset) > alt_hap_segment_boundary1 {
                        if debug {
                            eprintln!("failed segment-in-spacer");
                        }
                        return None;
                    }
                    right_segment_soft_clip_len += read_offset;
                }
            } else {
                // Set right_segment_alt_hap_offset at the first read-consuming cigar segment on the right
                // side of the alt hap:
                if right_segment_alt_hap_offset.is_none() && read_offset > 0 {
                    assert!(alt_hap_pos >= alt_hap_segment_boundary2);
                    right_segment_alt_hap_offset = Some(alt_hap_pos - alt_hap_segment_boundary2);
                }

                if right_segment_cigar.is_empty() {
                    if ref_offset == 0 {
                        right_segment_soft_clip_len += read_offset;
                    } else {
                        if right_segment_soft_clip_len > 0 {
                            right_segment_cigar
                                .push(Cigar::SoftClip(right_segment_soft_clip_len as u32));
                        }
                        right_segment_cigar.push(*c);
                    }
                } else {
                    right_segment_cigar.push(*c);
                }
            }

            alt_hap_pos += ref_offset;
        }
        assert!(!right_segment_cigar.is_empty());
    }

    let right_segment_alt_hap_offset = right_segment_alt_hap_offset.unwrap();

    // Range of the breakpoint breakend positions against the alt haplotype sequence. Due to the
    // arrangement of the alt haplotype seqeunce, all SV breakpoints will present as indels on this
    // sequence.
    let sv_alt_hap_range = {
        let start = assembly_contig_to_alt_hap_alignment.ref_offset
            + get_cigar_ref_offset(&left_segment_cigar);
        let end = alt_hap_segment_boundary2 + right_segment_alt_hap_offset;
        IntRange::from_pair(start, end)
    };

    // Range of the breakpoint breakend positions against the contig sequence
    let sv_contig_range = {
        let start = (contig.len() - left_segment_soft_clip_len) as i64;
        let end = right_segment_soft_clip_len as i64;
        IntRange::from_pair(start, end)
    };

    // Find breakpoint homology info between the left and right breakends:
    let (left_right_breakend_homology_range, left_right_breakend_homology_seq) = {
        // For the purpose of homology range calculation we use an extended version of the alt_hap_seq with a
        // full 'interior' sequence, and adjust sv_alt_hap_range to simulate an alignment to this
        // longer sequence:
        //
        assert!(alt_hap_info.extended_alt_hap_seq.len() >= alt_hap_info.alt_hap_seq.len());
        let extend_len = alt_hap_info.extended_alt_hap_seq.len() - alt_hap_info.alt_hap_seq.len();
        let extended_sv_alt_hap_range = IntRange::from_pair(
            sv_alt_hap_range.start,
            sv_alt_hap_range.end + extend_len as i64,
        );
        get_indel_breakend_homology_info(
            &alt_hap_info.extended_alt_hap_seq,
            &extended_sv_alt_hap_range,
            contig,
            &sv_contig_range,
        )
    };

    if debug {
        eprintln!("sv_alt_hap_range {sv_alt_hap_range:?} sv_contig_range {sv_contig_range:?}");
        let (before, after) =
            get_seq_pos_flanks(&alt_hap_info.alt_hap_seq, sv_alt_hap_range.start, 10);
        eprintln!("alt_hap_seq 10 bases below/above sv_alt_hap_range.start: {before}/{after}");
        let (before, after) =
            get_seq_pos_flanks(&alt_hap_info.alt_hap_seq, sv_alt_hap_range.end, 10);
        eprintln!("alt_hap_seq 10 bases below/above sv_alt_hap_range.end: {before}/{after}");
        let (before, after) = get_seq_pos_flanks(contig, sv_contig_range.start, 10);
        eprintln!("contig_seq 10 bases below/above contig_seq.start: {before}/{after}");
        eprintln!("left_right_hrange: {left_right_breakend_homology_range:?}");
    }

    Some(AltHapLeftRightComponentAlignmentInfo {
        left_alignment: SimpleAlignment {
            ref_offset: assembly_contig_to_alt_hap_alignment.ref_offset,
            cigar: left_segment_cigar,
        },
        right_alignment: SimpleAlignment {
            ref_offset: right_segment_alt_hap_offset,
            cigar: right_segment_cigar,
        },
        left_right_breakend_homology_range,
        left_right_breakend_homology_seq,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use rust_vc_utils::GenomeSegment;
    use rust_vc_utils::bam_utils::cigar::get_cigar_from_string;

    fn get_test_alt_hap_info() -> TwoRegionAltHapInfo {
        let alt_hap_seq = b"TTTTACTGACTGTTTNNNNTTTACTGACTGTTTT".to_vec();
        TwoRegionAltHapInfo {
            ref_segment1: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(100, 115),
            },
            ref_segment2: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(200, 215),
            },
            segment1_neighbor_extension_size: 0,
            segment2_neighbor_extension_size: 0,
            ref_segment2_revcomp: false,
            ref_segment2_after_ref_segment1: true,
            spacer_range: IntRange::from_pair(15, 15 + 4),
            left_boundary_pin: false,
            right_boundary_pin: false,
            extended_alt_hap_seq: alt_hap_seq.clone(),
            alt_hap_seq,
        }
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=10D8="),
        };
        let alt_hap_info = get_test_alt_hap_info();
        let contig = b"ACTGACTGACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            10,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_none());

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=8S"
        );
        assert_eq!(result.right_alignment.ref_offset, 3);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "8S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 0)
        );
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components_with_breakend_homology() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=10D8="),
        };
        let mut alt_hap_info = get_test_alt_hap_info();
        alt_hap_info.alt_hap_seq = b"TTTTACTGACTGACCNNNNTTTACTGACTGTTTT".to_vec();
        alt_hap_info.extended_alt_hap_seq = alt_hap_info.alt_hap_seq.clone();
        let contig = b"ACTGACTGACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=8S"
        );
        assert_eq!(result.right_alignment.ref_offset, 3);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "8S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 2)
        );
        assert_eq!(
            std::str::from_utf8(result.left_right_breakend_homology_seq.as_slice()).unwrap(),
            "AC"
        );
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components_with_right_breakend_insertion() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=10D4I8="),
        };
        let alt_hap_info = get_test_alt_hap_info();
        let contig = b"ACTGACTGCCCCACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=12S"
        );
        assert_eq!(result.right_alignment.ref_offset, 3);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "12S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 0)
        );
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components_with_left_breakend_insertion() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=4I10D8="),
        };
        let alt_hap_info = get_test_alt_hap_info();
        let contig = b"ACTGACTGCCCCACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=12S"
        );
        assert_eq!(result.right_alignment.ref_offset, 3);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "12S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 0)
        );
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components_with_left_pin() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=7D8="),
        };
        let alt_hap_seq = b"TTTTACTGACTGNNNNTTTACTGACTGTTTT".to_vec();
        let mut alt_hap_info = TwoRegionAltHapInfo {
            ref_segment1: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(100, 112),
            },
            ref_segment2: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(200, 215),
            },
            segment1_neighbor_extension_size: 0,
            segment2_neighbor_extension_size: 0,
            ref_segment2_revcomp: false,
            ref_segment2_after_ref_segment1: true,
            spacer_range: IntRange::from_pair(12, 12 + 4),
            left_boundary_pin: true,
            right_boundary_pin: false,
            extended_alt_hap_seq: alt_hap_seq.clone(),
            alt_hap_seq,
        };

        let contig = b"ACTGACTGACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=8S"
        );
        assert_eq!(result.right_alignment.ref_offset, 3);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "8S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 0)
        );

        // Test the invalid counterexample:
        alt_hap_info.left_boundary_pin = false;
        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_transform_alt_hap_alignment_into_left_right_components_with_right_pin() {
        let assembly_contig_to_alt_hap_alignment = SimpleAlignment {
            ref_offset: 4,
            cigar: get_cigar_from_string("8=7D8="),
        };
        let alt_hap_seq = b"TTTTACTGACTGTTTNNNNACTGACTGTTTT".to_vec();
        let alt_hap_info = TwoRegionAltHapInfo {
            ref_segment1: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(100, 115),
            },
            ref_segment2: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(200, 212),
            },
            segment1_neighbor_extension_size: 0,
            segment2_neighbor_extension_size: 0,
            ref_segment2_revcomp: false,
            ref_segment2_after_ref_segment1: true,
            spacer_range: IntRange::from_pair(15, 15 + 4),
            left_boundary_pin: false,
            right_boundary_pin: true,
            extended_alt_hap_seq: alt_hap_seq.clone(),
            alt_hap_seq,
        };

        let contig = b"ACTGACTGACTGACTG";

        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );

        assert!(result.is_some());
        let result = result.unwrap();
        assert_eq!(result.left_alignment.ref_offset, 4);
        assert_eq!(
            CigarString(result.left_alignment.cigar.clone()).to_string(),
            "8=8S"
        );
        assert_eq!(result.right_alignment.ref_offset, 0);
        assert_eq!(
            CigarString(result.right_alignment.cigar.clone()).to_string(),
            "8S8="
        );
        assert_eq!(
            result.left_right_breakend_homology_range,
            IntRange::from_pair(0, 0)
        );

        // Test the invalid counterexample (no longer valie with the extended_alt_hap_seq scheme)
        /*
        alt_hap_info.right_boundary_pin = false;
        let result = transform_alt_hap_alignment_into_left_right_components(
            3,
            &assembly_contig_to_alt_hap_alignment,
            &alt_hap_info,
            contig,
        );
        assert!(result.is_none());
        */
    }
}
