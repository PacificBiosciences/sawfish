mod segment_alignment;

use std::fmt;

use rust_htslib::bam::record::Cigar;
use rust_vc_utils::cigar::{
    get_cigar_ref_offset, get_cigarseg_read_offset, get_cigarseg_ref_offset,
    get_gap_compressed_identity_from_alignment, get_read_clip_positions, has_aligned_segments,
    update_read_pos, update_ref_and_read_pos, update_ref_pos,
};
use rust_vc_utils::indel_breakend_homology::get_indel_breakend_homology_info;
use rust_vc_utils::int_range::{IntRange, get_overlap_range};
use rust_vc_utils::{ChromList, GenomeRef, GenomeSegment, rev_comp_in_place};

use self::segment_alignment::{
    AltHapLeftRightComponentAlignmentInfo, transform_alt_hap_alignment_into_left_right_components,
};
use super::assemble::AssemblyResultContig;
use super::{AssemblyResult, RefineSVSettings, RefinedSV};
use crate::bam_utils::get_gap_compressed_identity_from_cigar_segment_range;
use crate::breakpoint::{
    Breakend, BreakendDirection, Breakpoint, BreakpointCluster, FullBreakendNeighborInfo,
    InsertInfo,
};
use crate::contig_output::{ContigAlignmentInfo, ContigAlignmentSegment};
use crate::expected_ploidy::SVLocusExpectedCNInfo;
use crate::genome_ref_utils::get_ref_segment_seq;
use crate::log_utils::debug_msg;
use crate::refine_sv::AnnotatedOverlappingHaplotype;
use crate::simple_alignment::{SimpleAlignment, clip_alignment_ref_edges};
use crate::sv_group::{
    ClusterAssembly, ClusterAssemblyAlignment, GroupHaplotypeId, SVGroup, SVGroupHaplotype,
};
use crate::sv_id::SVUniqueIdData;
use crate::utils::print_fasta;
use crate::wfa2_utils::PairwiseAligner;

/// Expand the target segment and extract the corresponding reference chromosome region
///
/// Returns a 4-tuple of
/// 1. Actual expanded target segment
/// 2. Sequence corresponding to the actual target segment in (1)
/// 3. Left shift - actual ref flank size expansion on the left side
/// 4. Right shift - actual ref flank size expansion on the right side
///
fn get_chrom_info_from_target_segment(
    left_ref_flank_size: i64,
    right_ref_flank_size: i64,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segment: &GenomeSegment,
) -> (GenomeSegment, Vec<u8>, i64, i64) {
    let mut ref_segment = target_segment.clone();
    let (actual_left_flank_size, actual_right_flank_size) =
        ref_segment.asymmetric_expand_by(chrom_list, left_ref_flank_size, right_ref_flank_size);

    let ref_seq = get_ref_segment_seq(chrom_list, reference, &ref_segment).to_vec();
    (
        ref_segment,
        ref_seq,
        actual_left_flank_size,
        actual_right_flank_size,
    )
}

/// Convert contig alignment indel into RefinedSV format so that it can be added to the
/// candidate SV set
///
/// # Arguments
/// * `contig_pos` - position in contig coordinates immediately after the start of the indel
///
#[allow(clippy::too_many_arguments)]
fn create_indel_refined_sv_candidate(
    id: SVUniqueIdData,
    target_chrom_index: usize,
    chrom_pos: i64,
    del_len: u32,
    ins_len: u32,
    contig_pos: i64,
    contig_seq: &[u8],
    target_segments: &[GenomeSegment],
    chrom_segment_start: i64,
    chrom_segment_seq: &[u8],
) -> RefinedSV {
    // Determine full breakpoint homology range for the indel candidate:
    //
    let indel_contig_range = IntRange::from_pair(contig_pos, contig_pos + ins_len as i64);
    assert!(chrom_segment_start <= chrom_pos);
    let chrom_segment_pos = chrom_pos - chrom_segment_start;
    let indel_chrom_segment_range =
        IntRange::from_pair(chrom_segment_pos, chrom_segment_pos + del_len as i64);
    let (homology_range, homology_seq) = get_indel_breakend_homology_info(
        chrom_segment_seq,
        &indel_chrom_segment_range,
        contig_seq,
        &indel_contig_range,
    );

    // Left shift the indel start positions if needed
    let chrom_pos = chrom_pos + homology_range.start;
    let contig_pos = (contig_pos + homology_range.start) as usize;
    let breakend_homology_length = homology_range.size();

    let breakend1 = {
        let start_pos = chrom_pos - 1;
        let end_pos = start_pos + breakend_homology_length + 1;
        Breakend {
            segment: GenomeSegment {
                chrom_index: target_chrom_index,
                range: IntRange::from_pair(start_pos, end_pos),
            },
            dir: BreakendDirection::LeftAnchor,
        }
    };
    let breakend2 = {
        let start_pos = chrom_pos - 1 + del_len as i64;
        let end_pos = start_pos + breakend_homology_length + 1;
        Breakend {
            segment: GenomeSegment {
                chrom_index: target_chrom_index,
                range: IntRange::from_pair(start_pos, end_pos),
            },
            dir: BreakendDirection::RightAnchor,
        }
    };
    let insert_info = if ins_len > 0 {
        let insert_seq = contig_seq[contig_pos..(contig_pos + ins_len as usize)].to_vec();
        InsertInfo::Seq(insert_seq)
    } else {
        InsertInfo::NoInsert
    };

    assert!(contig_pos >= 1);
    let contig_pos_before_breakend1 = contig_pos - 1;
    RefinedSV {
        id,
        bp: Breakpoint {
            breakend1,
            breakend2: Some(breakend2),
            insert_info,
            is_precise: true,
            ..Default::default()
        },
        breakend1_homology_seq: homology_seq,
        single_region_refinement: true,
        contig_pos_before_breakend1,
        assembly_regions: target_segments.to_vec(),
        ext: Default::default(),
        score: Default::default(),
    }
}

fn get_next_sv_id(
    sample_index: usize,
    cluster_index: usize,
    assembly_index: usize,
    alignment_index: &mut usize,
) -> SVUniqueIdData {
    let result = SVUniqueIdData {
        sample_index,
        cluster_index,
        assembly_index,
        alignment_index: *alignment_index,
    };
    *alignment_index += 1;
    result
}

#[allow(dead_code)]
fn get_contig_hq_overlap_size(
    contig_pos: usize,
    len: u32,
    assembly_high_quality_range: &IntRange,
) -> u32 {
    let contig_region = IntRange::from_pair(contig_pos as i64, (contig_pos + len as usize) as i64);
    let contig_hq_overlap = get_overlap_range(&contig_region, assembly_high_quality_range);
    if let Some(contig_hq_overlap) = contig_hq_overlap {
        contig_hq_overlap.size() as u32
    } else {
        0
    }
}

/// Check that enough of the anchor occurs in the high quality region of the read
#[allow(dead_code)]
fn is_good_indel_anchor(
    min_assembly_edge_anchor: u32,
    anchor_len: u32,
    read_pos: usize,
    read_hq_range: &IntRange,
) -> bool {
    let anchor_hq_len = get_contig_hq_overlap_size(read_pos, anchor_len, read_hq_range);
    anchor_hq_len >= min_assembly_edge_anchor
}

/// Filter contigs when either the left or right contig flank alignment quality is poor
///
/// Here alignment quality is measured as gap-compressed identity
///
/// Design note: It has been challanging to manage TP filtration with this filter. Two common
/// challange cases are (1) region is just diverged from the reference with a high SNP rate
/// (2) contig is not long enough for adjacent downstrem/upstream SV and partially extends into
/// it, creating high alignment noise.
///
/// See further discussion and examples on CR-332
///
#[allow(dead_code)]
fn flanking_gap_compressed_identity_filter(
    refined_svs: &mut Vec<(RefinedSV, usize)>,
    assembly_contig_to_chrom_segment_alignment: &SimpleAlignment,
    cluster_index: usize,
    assembly_index: usize,
    assembly_contig: &AssemblyResultContig,
    chrom_segment_seq: &[u8],
) {
    let debug = false;

    if refined_svs.is_empty() {
        return;
    }
    let firstsv_segment_index = refined_svs[0].1;
    let prefix_gci = get_gap_compressed_identity_from_cigar_segment_range(
        assembly_contig_to_chrom_segment_alignment.ref_offset,
        &assembly_contig.seq,
        &assembly_contig_to_chrom_segment_alignment.cigar,
        chrom_segment_seq,
        0,
        firstsv_segment_index,
    );

    let lastsv_segment_index = refined_svs[refined_svs.len() - 1].1;
    let suffix_gci = get_gap_compressed_identity_from_cigar_segment_range(
        assembly_contig_to_chrom_segment_alignment.ref_offset,
        &assembly_contig.seq,
        &assembly_contig_to_chrom_segment_alignment.cigar,
        chrom_segment_seq,
        lastsv_segment_index + 1,
        assembly_contig_to_chrom_segment_alignment.cigar.len(),
    );

    if debug {
        eprintln!(
            "clusterID:assemblyID {cluster_index}:{assembly_index}, prefix/suffix gci: {prefix_gci}/{suffix_gci}"
        );
    }

    let min_flank_gci = 0.92;
    if prefix_gci < min_flank_gci || suffix_gci < min_flank_gci {
        refined_svs.clear();
    }
}

/// Find 0 to many SV candidates from an assembly contig alignment to a single reference genome
/// region.
///
/// Multiple SV candidates may be found when they are in phase on the same haplotype, and thus
/// appear as part of the same contig alignment.
///
#[allow(clippy::too_many_arguments)]
fn get_single_region_refined_sv_candidates(
    refine_settings: &RefineSVSettings,
    mut chrom_pos: i64,
    assembly_contig_to_chrom_segment_alignment: &SimpleAlignment,
    cluster_index: usize,
    assembly_index: usize,
    assembly_contig: &AssemblyResultContig,
    target_segments: &[GenomeSegment],
    target_chrom_index: usize,
    chrom_segment_start: i64,
    chrom_segment_seq: &[u8],
    is_large_insertion: bool,
) -> Vec<RefinedSV> {
    // Find SV sized events in cigar string and locate their reference positions
    let mut contig_pos = 0;
    let mut found_anchor = false;
    let mut tmp_sv_buffer = Vec::new();
    let mut alignment_index = 0;
    let mut refined_svs = Vec::new();

    // This is a single sample method, so hard-code sample index to zero
    let sample_index = 0;

    for (cigar_segment_index, cigar_segment) in assembly_contig_to_chrom_segment_alignment
        .cigar
        .iter()
        .enumerate()
    {
        let mut get_indel_refined_sv_candidate = |del_len: u32, ins_len: u32| {
            let id = get_next_sv_id(
                sample_index,
                cluster_index,
                assembly_index,
                &mut alignment_index,
            );
            create_indel_refined_sv_candidate(
                id,
                target_chrom_index,
                chrom_pos,
                del_len,
                ins_len,
                contig_pos as i64,
                &assembly_contig.seq,
                target_segments,
                chrom_segment_start,
                chrom_segment_seq,
            )
        };

        let in_hiqual_contig_region = assembly_contig
            .high_quality_range
            .intersect_pos(contig_pos as i64);

        use Cigar::*;
        match cigar_segment {
            Equal(len) => {
                if *len >= refine_settings.min_assembly_edge_anchor {
                    if found_anchor {
                        refined_svs.append(&mut tmp_sv_buffer);
                    }
                    tmp_sv_buffer.clear();
                    found_anchor = true;
                }
            }
            Del(len) => {
                if *len >= refine_settings.min_indel_size && in_hiqual_contig_region {
                    tmp_sv_buffer
                        .push((get_indel_refined_sv_candidate(*len, 0), cigar_segment_index));
                }
            }
            Ins(len) => {
                if *len >= refine_settings.min_indel_size && in_hiqual_contig_region {
                    tmp_sv_buffer
                        .push((get_indel_refined_sv_candidate(0, *len), cigar_segment_index));
                }
            }
            _ => {}
        };
        update_ref_and_read_pos(cigar_segment, &mut chrom_pos, &mut contig_pos, true);
    }

    /*
    // See CR-332 for discussion of tradeoff for this filter
    flanking_gap_compressed_identity_filter(
        &mut refined_svs,
        assembly_contig_to_chrom_segment_alignment,
        cluster_index,
        assembly_index,
        assembly_contig,
        chrom_segment_seq,
    );
    */

    if is_large_insertion {
        const MIN_LARGE_INSERTION_LEN: usize = 3_000;

        // In this mode only output a large insertion
        refined_svs
            .into_iter()
            .filter(|(x, _)| x.bp.insert_info.size() > MIN_LARGE_INSERTION_LEN)
            .map(|(x, _)| x)
            .collect()
    } else {
        refined_svs.into_iter().map(|(x, _)| x).collect()
    }
}

/// Align assembly contig(s) to a single reference region
///
/// This process looks for smaller colinear SVs in a given region. Multiple overlapping SVs may be
/// discovered.
///
#[allow(clippy::too_many_arguments)]
fn align_single_ref_region_assemblies(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segments: &[GenomeSegment],
    assemblies: &[AssemblyResult],
    expected_cn_info: SVLocusExpectedCNInfo,
    cluster_index: usize,
    is_large_insertion: bool,
) -> (Vec<ContigAlignmentInfo>, Option<SVGroup>) {
    assert_eq!(target_segments.len(), 1);
    let target_segment = &target_segments[0];

    // Contigs of at least this size will be handled using a more efficient (if slower) memory model during alignment to the reference
    let min_long_contig_mode_len = 8_000;
    let is_long_contig = assemblies
        .iter()
        .flat_map(|x| &x.contigs)
        .any(|x| x.seq.len() > min_long_contig_mode_len);
    let mut aligner = PairwiseAligner::new_large_indel_aligner(is_long_contig);

    let (ref_segment, ref_segment_seq, _left_pin, _right_pin) = get_chrom_info_from_target_segment(
        refine_settings.hap_alignment_ref_flank_size,
        refine_settings.hap_alignment_ref_flank_size,
        reference,
        chrom_list,
        target_segment,
    );

    let target_chrom_index = target_segment.chrom_index;

    aligner.set_text(&ref_segment_seq);

    // output
    let mut contig_alignments = Vec::new();
    let mut refined_svs = Vec::new();
    let mut group_haplotypes = Vec::new();

    // Only accept one contig's worth of assembly output if in large insertion mode:
    let filter_all_overlapping_svs = is_large_insertion;
    let filter_close_overlapping_svs =
        refine_settings.reduce_overlapping_sv_alleles && (!filter_all_overlapping_svs);
    let is_any_filter = filter_all_overlapping_svs || filter_close_overlapping_svs;

    let debug = false;
    debug_msg!(
        debug,
        "SingleRef Cluster {cluster_index}: filter_all_overlapping_svs: {filter_all_overlapping_svs} filter_close_overlapping_svs: {filter_close_overlapping_svs}"
    );

    let mut overlapping_haplotypes = Vec::new();

    for (assembly_index, assembly) in assemblies.iter().enumerate() {
        if assembly.supporting_read_count() < refine_settings.min_assembly_read_support {
            debug_msg!(
                debug,
                "SingleRef cluster/contig: {cluster_index}/{assembly_index} - filtered for low read support: {}",
                assembly.supporting_read_count()
            );
            continue;
        }

        // Single ref region mode assumes only one assembly contig
        assert_eq!(assembly.contigs.len(), 1);
        let assembly_contig = &assembly.contigs[0];

        let (alignment_score, assembly_contig_to_chrom_segment_alignment) = match aligner
            .align(&assembly_contig.seq)
        {
            (score, Some(x)) => (score, x),
            (_, None) => {
                debug_msg!(
                    debug,
                    "SingleRef cluster/contig: {cluster_index}/{assembly_index} - filtered for poor contig to ref alignment"
                );
                // TODO add warning or statistics to track these filtration cases
                continue;
            }
        };
        debug_msg!(
            debug,
            "SingleRef cluster/contig: {cluster_index}/{assembly_index} - assembly_contig_to_chrom_segment_alignment: {:?}",
            assembly_contig_to_chrom_segment_alignment
        );
        debug_msg!(
            debug,
            "SingleRef cluster/contig: {cluster_index}/{assembly_index} - alignment_score {}",
            alignment_score
        );

        // Starting reference position of the assembly alignment
        let chrom_pos =
            ref_segment.range.start + assembly_contig_to_chrom_segment_alignment.ref_offset;

        // Output contigs to bam
        {
            let sample_index = 0;
            let segment = ContigAlignmentSegment {
                tid: target_chrom_index as i32,
                pos: chrom_pos,
                cigar: assembly_contig_to_chrom_segment_alignment.cigar.clone(),
                is_fwd_strand: true,
            };
            contig_alignments.push(ContigAlignmentInfo {
                sample_index,
                cluster_index,
                assembly_index,
                seq: assembly_contig.seq.clone(),
                segments: vec![segment],
                segment_id: 0,
                supporting_read_count: assembly.supporting_read_count(),
                high_quality_contig_range: assembly_contig.high_quality_range.clone(),
            });
        }

        if refined_svs.is_empty() || (!filter_all_overlapping_svs) {
            let new_refined_svs = get_single_region_refined_sv_candidates(
                refine_settings,
                chrom_pos,
                &assembly_contig_to_chrom_segment_alignment,
                cluster_index,
                assembly_index,
                assembly_contig,
                target_segments,
                target_chrom_index,
                ref_segment.range.start,
                &ref_segment_seq,
                is_large_insertion,
            );

            let start_sv_count = refined_svs.len();
            for new_sv in new_refined_svs {
                let mut accept_new_sv = true;
                if filter_close_overlapping_svs {
                    fn get_sv_dist(sv1: &RefinedSV, sv2: &RefinedSV) -> i64 {
                        (sv1.bp.breakend1.segment.range.start
                            - sv2.bp.breakend1.segment.range.start)
                            .abs()
                    }

                    for current_sv in refined_svs.iter().take(start_sv_count) {
                        if get_sv_dist(&new_sv, current_sv) < refine_settings.min_close_sv_range {
                            accept_new_sv = false;
                            break;
                        }
                    }
                }
                if accept_new_sv {
                    refined_svs.push(new_sv);
                }
            }

            let contig_alignment = SimpleAlignment {
                ref_offset: (assembly_contig_to_chrom_segment_alignment.ref_offset
                    + ref_segment.range.start),
                cigar: assembly_contig_to_chrom_segment_alignment.cigar.clone(),
            };

            // Extend group contigs list:
            {
                let hap_id = GroupHaplotypeId {
                    sample_index: 0,
                    cluster_index,
                    assembly_index,
                };
                let contig_alignment = ClusterAssemblyAlignment {
                    contig_seq: assembly_contig.seq.clone(),
                    chrom_index: target_segment.chrom_index,
                    contig_alignment: contig_alignment.clone(),
                    high_quality_contig_range: assembly_contig.high_quality_range.clone(),
                    is_fwd_strand: true,
                };
                let group_haplotype = SVGroupHaplotype {
                    hap_id,
                    contig_info: ClusterAssembly {
                        contig_alignments: vec![contig_alignment],
                        supporting_read_count: assembly.supporting_read_count(),
                    },
                };
                group_haplotypes.push(group_haplotype);
            }

            let sv_candidate_count = refined_svs.len() - start_sv_count;
            if (!is_any_filter) || (sv_candidate_count > 0) {
                // Add any overlapping haplotypes to this structure, unless SVs on the overlapping
                // haplotype could have been filtered out.
                overlapping_haplotypes.push(AnnotatedOverlappingHaplotype {
                    assembly_index,
                    supporting_read_count: assembly.supporting_read_count(),
                    sv_candidate_count,
                });
            }
        }
    }

    // Filter unlikely contigs from the overlapping haplotype set:
    let overlapping_haplotypes = {
        let min_contig_supporting_read_fraction = 0.1;
        let total_supporting_read_count = overlapping_haplotypes
            .iter()
            .map(|x| x.supporting_read_count)
            .sum::<usize>() as f64;
        overlapping_haplotypes
            .into_iter()
            .filter(|x| {
                x.sv_candidate_count != 0
                    || ((x.supporting_read_count as f64 / total_supporting_read_count)
                        >= min_contig_supporting_read_fraction)
            })
            .collect::<Vec<_>>()
    };

    let sv_group = if refined_svs.is_empty() {
        None
    } else {
        let sample_haplotype_list = vec![
            overlapping_haplotypes
                .iter()
                .map(|x| x.assembly_index)
                .collect(),
        ];
        let sample_expected_cn_info = vec![expected_cn_info];
        let sv_haplotype_map = refined_svs.iter().map(|x| x.id.assembly_index).collect();
        let sv_group = SVGroup {
            group_regions: refined_svs[0].assembly_regions.clone(),
            group_haplotypes,
            sample_haplotype_list,
            sample_expected_cn_info,
            sv_haplotype_map,
            refined_svs,
        };
        Some(sv_group)
    };

    (contig_alignments, sv_group)
}

/// Enumerate all details on the creation of the 2-region alt-haplotype 'reference' sequence, to
/// which the 2-region sample assemblies will be aligned.
///
pub struct TwoRegionAltHapInfo {
    /// Reference segment extracted from the reference region around breakend1
    pub ref_segment1: GenomeSegment,

    /// Reference segment extracted from the reference region around breakend2
    pub ref_segment2: GenomeSegment,

    /// Amount of additional sequence extended from the segment1 'neighbor' breakpoint, this is not part of ref_segment1
    pub segment1_neighbor_extension_size: usize,

    /// Amount of additional sequence extended from the segment2 'neighbor' breakpoint, this is not part of ref_segment2
    pub segment2_neighbor_extension_size: usize,

    /// If true, the ref segment2 region is reverse complemented in the alt hap derived chromosome
    /// model
    ///
    /// The alt hap sequence is constructed such that ref_segment1 never needs to be reversed
    ///
    pub ref_segment2_revcomp: bool,

    /// If true, the ref segment2 region is appended to the right side of ref segment1 in the alt
    /// hap derived chromosome model
    ///
    /// For a simple deletion this will always be true, given that ref segment1 is assigned to
    /// have a genomic sort order lower than ref segment2.
    ///
    /// An example where this would be false would be a tandem duplication breakpoint
    ///
    pub ref_segment2_after_ref_segment1: bool,

    /// Range of the 'N' sequence space inserted between the reference-derived segments of the alt
    /// hap sequence
    ///
    /// Coordinates are zero-indexed and expressed in the coordinate system of the alt haplotype
    /// sequence, rather than any source reference contig from which the alt haplotype was built.
    ///
    pub spacer_range: IntRange,

    /// If true the left side of the interior alt-hap segment boundary is from a segment of the
    /// reference that breaks at the end of the chromosome.
    pub left_boundary_pin: bool,

    /// If true the right side of the interior alt-hap segment boundary is from a segment of the
    /// reference that breaks at the start of the chromosome.
    pub right_boundary_pin: bool,

    /// The fully constructed alt haplotype alignment target sequence
    pub alt_hap_seq: Vec<u8>,

    /// A variation on alt_hap_seq with longer flank sequences. This is not intended for contig alingment,
    /// but can be useful to more accurately assess breakpoint homology.
    ///
    pub extended_alt_hap_seq: Vec<u8>,
}

impl TwoRegionAltHapInfo {
    pub fn get_ref_segment(&self, segment_index: usize) -> &GenomeSegment {
        assert!(segment_index < 2);
        if segment_index == 0 {
            &self.ref_segment1
        } else {
            &self.ref_segment2
        }
    }

    pub fn get_segment_neighbor_extension_size(&self, segment_index: usize) -> usize {
        assert!(segment_index < 2);
        if segment_index == 0 {
            self.segment1_neighbor_extension_size
        } else {
            self.segment2_neighbor_extension_size
        }
    }
}

impl fmt::Debug for TwoRegionAltHapInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ref_seg1/size/ext_size: {:?}/{}/{} ref_seg2/size/ext_size: {:?}/{}/{} ref_seg2_revcomp?: {} ref_seg2_after_ref_seg1?: {} \
            spacer_range: {:?} left_pin?: {} right_pin?: {} alt_hap_seq_size: {}",
            self.ref_segment1,
            self.ref_segment1.range.size(),
            self.segment1_neighbor_extension_size,
            self.ref_segment2,
            self.ref_segment2.range.size(),
            self.segment2_neighbor_extension_size,
            self.ref_segment2_revcomp,
            self.ref_segment2_after_ref_segment1,
            self.spacer_range,
            self.left_boundary_pin,
            self.right_boundary_pin,
            self.alt_hap_seq.len(),
        )
    }
}

/// Return a 2-tuple of (1) the alt-hap sequence and (2) the range of the spacer sequence within it
///
fn construct_alt_hap_seq(
    ref_segment1_seq: Vec<u8>,
    mut ref_segment2_seq: Vec<u8>,
    spacer_size: usize,
    ref_segment2_revcomp: bool,
    ref_segment2_after_ref_segment1: bool,
) -> (Vec<u8>, IntRange) {
    let spacer_seq = vec![b'N'; spacer_size];
    if ref_segment2_revcomp {
        rev_comp_in_place(&mut ref_segment2_seq);
    }

    let create_alt_hap_seq = |mut a: Vec<u8>, b: Vec<u8>| -> (Vec<u8>, IntRange) {
        let spacer_start = a.len() as i64;
        let spacer_range = IntRange::from_pair(spacer_start, spacer_start + spacer_size as i64);
        a.extend(spacer_seq);
        a.extend(b);
        (a, spacer_range)
    };

    if ref_segment2_after_ref_segment1 {
        create_alt_hap_seq(ref_segment1_seq, ref_segment2_seq)
    } else {
        create_alt_hap_seq(ref_segment2_seq, ref_segment1_seq)
    }
}

/// Find the max reference flanks size for the extension between two breakend regions on the same chromosome
///
/// For certain cases the reference flank size needs to be restricted to keep the two reference segments
/// from overlapping (but also keep the analysis as a two separate regions to limit false positives):
///
/// This limit is applied when we have noninverted breakpoints on a single chromosome.
///
/// The following example shows r1 and r2 regions corresponding to breakends 1 and 2:
///
///  |----r1----|    |----r2----|
///              zzzz
///
/// In this case the 'middle_flank_size' applied to each of r1 and r2 is restricted to half
/// of the 'z' region size
///
fn get_middle_flank_size(
    max_ref_flank_size: i64,
    target_segments: &[GenomeSegment],
    bp: &Breakpoint,
) -> i64 {
    if (!bp.same_orientation())
        && (target_segments[0].chrom_index == target_segments[1].chrom_index)
    {
        let ref_flank_size = std::cmp::max(
            0,
            target_segments[1].range.start - target_segments[0].range.end,
        ) / 2;
        std::cmp::min(max_ref_flank_size, ref_flank_size)
    } else {
        max_ref_flank_size
    }
}

/// Get the flank sizes used when extracting alt haplotype flank sequences
///
/// Return (left size, right size)
///
fn get_left_right_ref_flank_size(
    max_ref_flank_size: i64,
    is_segment_revcomped: bool,
    is_segment_first: bool,
) -> (i64, i64) {
    // Restrict the 'inside' of the contig, adjacent to the spacer is limited to a maximum size to
    //
    //   |---outer-flank--|-bp-|--inner-flank--|--spacer--|--inner-flank--|-bp-|---outer-flank---|
    //
    let max_interior_flank_size = 2000;
    let max_interior_flank_size = std::cmp::min(max_interior_flank_size, max_ref_flank_size);

    if is_segment_revcomped ^ is_segment_first {
        (max_ref_flank_size, max_interior_flank_size)
    } else {
        (max_interior_flank_size, max_ref_flank_size)
    }
}

struct RemoteBreakendInfo {
    /// left ref_flank_size of the target breakpoint after shortening to correctly place the remote_seq
    left_ref_flank_size: i64,

    /// right ref_flank_size of the target breakpoint after shortening to correctly place the remote_seq
    right_ref_flank_size: i64,

    remote_seq: Vec<u8>,
}

/// Get info on second breakpoint neighboring the current target breakend
///
/// If the 'target' breakend (the one comprising the breakpoint being assembled/scored/called) has a close neighboring breakpoint
/// on the same haplotype, then extract the additional information we need to represent this additional neighbor breakend in the
/// alt haplotype model of the event.
///
/// Diagram for terminology. Example shows 'target breakend' with 'neighbor breakend' downstream of it:
///
///                                           local          remote
///                      target            neighbor          neighbor
///                      breakend          breakend          breakend
/// ---------|          |--------------------------|        |---------------------------
///
/// # Arguments
/// * `full_ref_flank_size` - default ref flank size around each breakend prior to imposing any other length constraints
///
fn get_remote_breakend_info(
    full_ref_flank_size: i64,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segment: &GenomeSegment,
    target_breakend: &Breakend,
    breakend_neighbor_info: Option<FullBreakendNeighborInfo>,
    debug: bool,
) -> Option<RemoteBreakendInfo> {
    let breakend_neighbor_info = breakend_neighbor_info?;

    let mut left_ref_flank_size = full_ref_flank_size;
    let mut right_ref_flank_size = full_ref_flank_size;

    let breakend_neighbor_index = breakend_neighbor_info.breakend_neighbor.breakend_index;
    let local_breakend_neighbor = {
        let local_breakend_index = breakend_neighbor_index;
        breakend_neighbor_info
            .breakpoint
            .get_breakend(local_breakend_index)
    };

    let remote_breakend_neighbor = {
        let remote_breakend_index = (breakend_neighbor_index + 1) % 2;
        breakend_neighbor_info
            .breakpoint
            .get_breakend(remote_breakend_index)
    };

    // Verify the expected relationship between target and local neighbor breakends
    assert!(local_breakend_neighbor.dir != target_breakend.dir);
    assert!(local_breakend_neighbor.segment.chrom_index == target_breakend.segment.chrom_index);

    if debug {
        eprintln!(
            "target breakend {target_breakend:?} local neighbor breakend: {local_breakend_neighbor:?} remote neighbor breakend {remote_breakend_neighbor:?}"
        );
    }

    // Limit the flank size we use in the first extraction region based on the local breakpoint neighbor,
    // and find the length of sequence to extract from the remote breakend at the same time
    //
    // Initialize remote extraction to zero, because any breakend without a neighbor will naturally have zero
    // sequence extracted from a neighbor:
    let mut left_remote_flank_size = 0;
    let mut right_remote_flank_size = 0;
    match target_breakend.dir {
        BreakendDirection::LeftAnchor => {
            let max_left_flank =
                target_segment.range.start - local_breakend_neighbor.segment.range.start;

            // Give up on using neighbor info for odd cases until we can more carefully understand and handle these
            if max_left_flank <= 0 {
                return None;
            }

            left_ref_flank_size = std::cmp::min(left_ref_flank_size, max_left_flank);
            left_remote_flank_size = full_ref_flank_size - left_ref_flank_size;
        }
        BreakendDirection::RightAnchor => {
            let max_right_flank =
                local_breakend_neighbor.segment.range.end - target_segment.range.end;

            // Give up on using neighbor info for odd cases until we can more carefully understand and handle these
            if max_right_flank <= 0 {
                return None;
            }

            right_ref_flank_size = std::cmp::min(right_ref_flank_size, max_right_flank);
            right_remote_flank_size = full_ref_flank_size - right_ref_flank_size;
        }
    };

    if debug {
        eprintln!(
            "left_full_ref_flank_size {left_ref_flank_size} left_remote_flank_size {left_remote_flank_size}"
        );
        eprintln!(
            "right_full_ref_flank_size {right_ref_flank_size} right_remote_flank_size {right_remote_flank_size}"
        );
    }

    // Verify remote flank size constraints
    assert!(left_remote_flank_size >= 0 && right_remote_flank_size >= 0);
    assert!(left_remote_flank_size == 0 || right_remote_flank_size == 0);

    // If there's no need to add remote breakend sequence to the alt haplotype model
    if left_remote_flank_size <= 0 && right_remote_flank_size <= 0 {
        return None;
    }

    // Extract the corresponding amount we just trimmed off from the remote breakpoint region
    //

    // Determine if the remote is reversed relative to local
    let is_remote_reversed = local_breakend_neighbor.dir == remote_breakend_neighbor.dir;

    // Up to this point left/right have always been with respect to the target breakend, now these
    // values need to make sense with respect to the remote breakend
    if is_remote_reversed {
        std::mem::swap(&mut left_remote_flank_size, &mut right_remote_flank_size);
    }

    // Verify that we're extracting from the correct side of the remote breakend
    match remote_breakend_neighbor.dir {
        BreakendDirection::LeftAnchor => {
            assert!(left_remote_flank_size > 0);
        }
        BreakendDirection::RightAnchor => {
            assert!(right_remote_flank_size > 0);
        }
    }

    let mut remote_target_segment = remote_breakend_neighbor.segment.clone();

    // The homology range segment has already been included in the target flank sizes, so the remote segment range needs to be shrunk down to zero,
    // and here we select which side of the segment range to keep.
    if left_remote_flank_size > 0 {
        remote_target_segment.range.end = remote_target_segment.range.start;
    } else {
        remote_target_segment.range.start = remote_target_segment.range.end;
    }
    let (_, mut remote_seq, _, _) = get_chrom_info_from_target_segment(
        left_remote_flank_size,
        right_remote_flank_size,
        reference,
        chrom_list,
        &remote_target_segment,
    );

    if is_remote_reversed {
        rev_comp_in_place(&mut remote_seq);
    }

    if debug {
        eprintln!("is_remote_reversed {is_remote_reversed}");
        //eprintln!("remote_seq:");
        //print_fasta(250, &[remote_seq.as_slice()]);
    }

    Some(RemoteBreakendInfo {
        left_ref_flank_size,
        right_ref_flank_size,
        remote_seq,
    })
}

/// Extract extended and standard alternate haplotype segment info for one breakend
///
/// Workflow here is to get everything set up to the get the 'extended' version of the segment info
/// and then copy that and trim it down to the standard length.
///
/// The standard length is used for alignment. The extended version is used only for homology length analysis.
///
/// # Arguments
/// * `left_ref_flank_size` - the amount of reference flank to add on the left side of the breakend
/// * `right_ref_flank_size` - the amount of reference flank to add on the right side of the breakend
/// * `target_breakend` - this is only needed if breakend_neighbor_info is defined
/// * `breakend_neighbor_info` - details of a possible neighboring breakpoint candidate on the same haplotype.
///   This can be used to improve representation of the alt haplotype when breakends from different breakpoints
///   are close.
/// * `breakend_index` - index of target breakend, only used for debug messaging
///
#[allow(clippy::too_many_arguments, clippy::single_match)]
fn extract_reference_segment_info(
    full_ref_flank_size: i64,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segment: &GenomeSegment,
    target_breakend: &Breakend,
    left_ref_flank_size: i64,
    right_ref_flank_size: i64,
    breakend_neighbor_info: Option<FullBreakendNeighborInfo>,
    breakend_index: usize,
    debug: bool,
) -> (GenomeSegment, Vec<u8>, Vec<u8>, i64, i64) {
    if debug {
        eprintln!("Starting extract_reference_segment_info");
        eprintln!(
            "breakend{breakend_index} target_segment {target_segment:?} left_ref_flank_size/right_ref_flank_size {left_ref_flank_size}/{right_ref_flank_size}"
        );
        if let Some(breakend_neighbor_info) = &breakend_neighbor_info {
            eprintln!(
                "Breakend neighbor defined. Index: {} bp: {:?}",
                breakend_neighbor_info.breakend_neighbor.breakend_index,
                breakend_neighbor_info.breakpoint
            );
        }
        eprintln!("Breakend neighbor target breakend: {target_breakend:?}");
    }

    let remote_breakend_info = get_remote_breakend_info(
        full_ref_flank_size,
        reference,
        chrom_list,
        target_segment,
        target_breakend,
        breakend_neighbor_info,
        debug,
    );

    // Extracting the maximum flank size on each side of the primary target location, unless part of this needs to be taken out to leave room for the remote seqeunce:
    //
    let (
        extended_ref_segment,
        mut extended_ref_segment_seq,
        actual_left_flank_range_size,
        actual_right_flank_range_size,
    ) = {
        let (left_full_ref_flank_size, right_full_ref_flank_size) =
            if let Some(remote_breakend_info) = &remote_breakend_info {
                (
                    remote_breakend_info.left_ref_flank_size,
                    remote_breakend_info.right_ref_flank_size,
                )
            } else {
                (full_ref_flank_size, full_ref_flank_size)
            };

        get_chrom_info_from_target_segment(
            left_full_ref_flank_size,
            right_full_ref_flank_size,
            reference,
            chrom_list,
            target_segment,
        )
    };

    let mut actual_left_flank_seq_size = actual_left_flank_range_size;
    let mut actual_right_flank_seq_size = actual_right_flank_range_size;

    // Modify extended_ref_segment_seq to include the remote breakend neighbor sequence
    if let Some(remote_breakend_info) = remote_breakend_info {
        let remote_seq = remote_breakend_info.remote_seq;
        let remote_seq_len = remote_seq.len();
        let mut tmp = remote_seq;
        match target_breakend.dir {
            BreakendDirection::LeftAnchor => {
                tmp.append(&mut extended_ref_segment_seq);
                extended_ref_segment_seq = tmp;
                actual_left_flank_seq_size += remote_seq_len as i64;
            }
            BreakendDirection::RightAnchor => {
                extended_ref_segment_seq.append(&mut tmp);
                actual_right_flank_seq_size += remote_seq_len as i64;
            }
        }
    }

    // Now make a second version of the segment with the shortened left/right flanks
    let left_seq_trim = std::cmp::max(actual_left_flank_seq_size - left_ref_flank_size, 0) as usize;
    let right_seq_trim =
        std::cmp::max(actual_right_flank_seq_size - right_ref_flank_size, 0) as usize;

    let left_range_trim =
        std::cmp::max(actual_left_flank_range_size - left_ref_flank_size, 0) as usize;
    let right_range_trim =
        std::cmp::max(actual_right_flank_range_size - right_ref_flank_size, 0) as usize;

    if debug {
        eprintln!("Finished extract_reference_segment_info");
        eprintln!("extended_ref_segment: {:?}", &extended_ref_segment);
        eprintln!(
            "actual_left_flank_seq_size/actual_left_flank_range_size/left_ref_flank/left_seq_trim/left_range_trim: {actual_left_flank_seq_size}/{actual_left_flank_range_size}/{left_ref_flank_size}/{left_seq_trim}/{left_range_trim}"
        );
        eprintln!(
            "actual_right_flank_seq_size/actual_right_flank_range_size/right_ref_flank/right_seq_trim/right_range_trim: {actual_right_flank_seq_size}/{actual_right_flank_range_size}/{right_ref_flank_size}/{right_seq_trim}/{right_range_trim}"
        );
    }

    let seq_len = extended_ref_segment_seq.len();
    assert!((left_seq_trim + right_seq_trim) <= seq_len);

    let range_len = extended_ref_segment.range.size() as usize;
    assert!((left_range_trim + right_range_trim) <= range_len);

    let mut ref_segment = extended_ref_segment;
    ref_segment.range.start += left_range_trim as i64;
    ref_segment.range.end -= right_range_trim as i64;
    let ref_segment_seq =
        extended_ref_segment_seq[left_seq_trim..(seq_len - right_seq_trim)].to_vec();

    if debug {
        eprintln!("ref_segment: {:?}", &ref_segment);
    }

    (
        ref_segment,
        ref_segment_seq,
        extended_ref_segment_seq,
        actual_left_flank_range_size,
        actual_right_flank_range_size,
    )
}

/// - `cluster_index` - only used for generating debug messages
/// - `segment_index` - only used for generating debug messages
///
fn get_segment_neighbor_extension_size(
    segment_seq: &[u8],
    segment: &GenomeSegment,
    cluster_index: usize,
    segment_index: usize,
) -> usize {
    let slen = segment_seq.len();
    let rsize = segment.range.size() as usize;
    if slen < rsize {
        panic!(
            "get_segment_neighbor_extension_size failed for cluster: {cluster_index} segment{segment_index}. seq len: {slen} seg: {segment:?}"
        );
    } else {
        slen - rsize
    }
}

/// Extract reference segments corresponding to the two candidate breakpoint breakends, then orient and join them, with a
/// small 'N' spacer sequence, to create a low-resolution alt haplotype that putative SV haplotype contigs
/// can be aligned against.
///
fn get_two_region_alt_hap_info(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segments: &[GenomeSegment],
    cluster_index: usize,
    bpc: &BreakpointCluster,
) -> TwoRegionAltHapInfo {
    let debug = false;
    if debug {
        eprintln!(
            "starting get_two_region_alt_hap_info for cluster_index {cluster_index} target_segments {target_segments:?}"
        );
    }

    // The strategy to build the alternate haplotype model is to extract reference compositions corresponding to the
    // two breakend regions, then orient and fuse the two segments together with a spacer. The 'reference composition'
    // on each side of the breakpoint is almost always a single reference segment, but can be multiple segments when breakpoints
    // come close together.

    // The alt haplotype model and all reads are being assembled in the reference frame of breakend1, ie the first
    // target segment, so orienting everything to match up to the breakpoint is set by two binary choices:
    // 1) Whether to revcomp segment2
    // 2) Append segment2 before or after segment1
    //
    // These choices are critical to most of the operation so set these two values here first:
    let bp = &bpc.breakpoint;
    let breakend1 = &bp.breakend1;
    let breakend2 = bp.breakend2.as_ref().unwrap();
    let ref_segment2_revcomp = breakend1.dir == breakend2.dir;
    let ref_segment2_after_ref_segment1 = breakend1.dir == BreakendDirection::LeftAnchor;

    // Note that middle flank size is applied between DEL/DUP regions based on the *reference* space concept of middle
    //
    // In the discussion below about the 'interior' of the alt hap contig model, this is refering to the *alt-hap* space concept of middle.
    //
    // For deletions these two concepts refer to the same thing, but in many other cases they are different.
    //
    let middle_flank_size = get_middle_flank_size(
        refine_settings.hap_alignment_ref_flank_size,
        target_segments,
        bp,
    );

    if debug {
        eprintln!("bp {bp:?}");
        eprintln!(
            "ref_segment2_revcomp {ref_segment2_revcomp} ref_segment2_after_ref_segment1 {ref_segment2_after_ref_segment1} middle_flank_size {middle_flank_size}"
        );
    }

    // Extract the two standard and extended reference segments:
    let (
        ref_segment1,
        ref_segment1_seq,
        extended_ref_segment1_seq,
        actual_left_flank_range_size1,
        actual_right_flank_range_size1,
    ) = {
        // Find any modifications which will require shortening the left or right flank size:
        let (left_ref_flank_size, right_ref_flank_size) = get_left_right_ref_flank_size(
            refine_settings.hap_alignment_ref_flank_size,
            false,
            ref_segment2_after_ref_segment1,
        );
        let right_ref_flank_size = std::cmp::min(right_ref_flank_size, middle_flank_size);

        // Extract extended and standard reference segment info and return
        let breakend_index = 0;
        extract_reference_segment_info(
            refine_settings.hap_alignment_ref_flank_size,
            reference,
            chrom_list,
            &target_segments[0],
            breakend1,
            left_ref_flank_size,
            right_ref_flank_size,
            bpc.get_full_breakend_neighbor_info(breakend_index == 0),
            breakend_index,
            debug,
        )
    };

    let (
        ref_segment2,
        ref_segment2_seq,
        extended_ref_segment2_seq,
        actual_left_flank_range_size2,
        actual_right_flank_range_size2,
    ) = {
        // Find any modifications which will require shortening the left or right flank size:
        let (left_ref_flank_size, right_ref_flank_size) = get_left_right_ref_flank_size(
            refine_settings.hap_alignment_ref_flank_size,
            ref_segment2_revcomp,
            !ref_segment2_after_ref_segment1,
        );
        let left_ref_flank_size = std::cmp::min(left_ref_flank_size, middle_flank_size);

        // Extract extended and standard reference segment info and return
        let breakend_index = 1;
        extract_reference_segment_info(
            refine_settings.hap_alignment_ref_flank_size,
            reference,
            chrom_list,
            &target_segments[1],
            breakend2,
            left_ref_flank_size,
            right_ref_flank_size,
            bpc.get_full_breakend_neighbor_info(breakend_index == 0),
            breakend_index,
            debug,
        )
    };

    let segment1_neighbor_extension_size =
        get_segment_neighbor_extension_size(&ref_segment1_seq, &ref_segment1, cluster_index, 0);
    let segment2_neighbor_extension_size =
        get_segment_neighbor_extension_size(&ref_segment2_seq, &ref_segment2, cluster_index, 1);

    if debug {
        eprintln!(
            "ref_segment1_seq_len {} ref_segment1 {:?} ref_segment1_size {} segment1_extension_size {}",
            ref_segment1_seq.len(),
            &ref_segment1,
            ref_segment1.range.size(),
            segment1_neighbor_extension_size
        );
    }

    let spacer_size = refine_settings.two_region_alt_haplotype_spacer_size;
    let (alt_hap_seq, spacer_range) = construct_alt_hap_seq(
        ref_segment1_seq,
        ref_segment2_seq,
        spacer_size,
        ref_segment2_revcomp,
        ref_segment2_after_ref_segment1,
    );
    if debug {
        eprintln!("alt_hap_seq");
        print_fasta(250, &[alt_hap_seq.as_slice()]);
    }

    let (extended_alt_hap_seq, _) = construct_alt_hap_seq(
        extended_ref_segment1_seq,
        extended_ref_segment2_seq,
        spacer_size,
        ref_segment2_revcomp,
        ref_segment2_after_ref_segment1,
    );

    let (left_boundary_pin, right_boundary_pin) = {
        let left_pin1 = actual_left_flank_range_size1 == 0;
        let left_pin2 = actual_left_flank_range_size2 == 0;
        let right_pin1 = actual_right_flank_range_size1 == 0;
        let right_pin2 = actual_right_flank_range_size2 == 0;
        if ref_segment2_after_ref_segment1 {
            (
                right_pin1,
                if ref_segment2_revcomp {
                    right_pin2
                } else {
                    left_pin2
                },
            )
        } else {
            (
                if ref_segment2_revcomp {
                    left_pin2
                } else {
                    right_pin2
                },
                left_pin1,
            )
        }
    };

    TwoRegionAltHapInfo {
        ref_segment1,
        ref_segment2,
        segment1_neighbor_extension_size,
        segment2_neighbor_extension_size,
        ref_segment2_revcomp,
        ref_segment2_after_ref_segment1,
        spacer_range,
        left_boundary_pin,
        right_boundary_pin,
        alt_hap_seq,
        extended_alt_hap_seq,
    }
}

/// Alignment of the contig to the alt hap sequence after being translated all the way back to the
/// original reference segment regions
///
#[derive(Debug)]
struct AltHapRefSegmentAlignmentInfo {
    ref_segment1_alignment: SimpleAlignment,
    ref_segment2_alignment: SimpleAlignment,

    /// Range of homology from the perspective of the ref segment1 breakpoint
    ref_segment1_breakend_homology_range: IntRange,

    /// Breakend homology sequence from the perspective of the ref segment1 breakpoint
    ref_segment1_breakend_homology_seq: Vec<u8>,
}

/// Alt hap alignment info is reported for the left and right portions of the alt hap sequence. Here
/// these are transformed back to the order and orientation of the original two reference segments
/// used to construct the alt hap sequence
///
/// Returns the alignment information reformatted for reference segments 1 and 2
///
fn transform_left_right_alignments_back_to_ref_segments(
    alt_hap_info: &TwoRegionAltHapInfo,
    alt_hap_left_right_component_alignment_info: AltHapLeftRightComponentAlignmentInfo,
) -> AltHapRefSegmentAlignmentInfo {
    let optionally_reversed_seg2_alignment = |alignment: SimpleAlignment| -> SimpleAlignment {
        if alt_hap_info.ref_segment2_revcomp {
            let segment2_flank_size = if alt_hap_info.ref_segment2_after_ref_segment1 {
                alt_hap_info.alt_hap_seq.len() as i64 - alt_hap_info.spacer_range.end
            } else {
                alt_hap_info.spacer_range.start
            };
            alignment.get_reverse(segment2_flank_size)
        } else {
            alignment
        }
    };

    // Reformat alignments from the alt haplotype left-right format back to the original reference
    // segments used to contruct the alt haplotype sequence:
    //
    let (ref_segment1_alignment, ref_segment2_alignment) =
        if alt_hap_info.ref_segment2_after_ref_segment1 {
            (
                alt_hap_left_right_component_alignment_info.left_alignment,
                optionally_reversed_seg2_alignment(
                    alt_hap_left_right_component_alignment_info.right_alignment,
                ),
            )
        } else {
            (
                alt_hap_left_right_component_alignment_info.right_alignment,
                optionally_reversed_seg2_alignment(
                    alt_hap_left_right_component_alignment_info.left_alignment,
                ),
            )
        };

    // The alt haplotype has been arranged so that ref segment1 is never revcomped, so the breakend
    // homology details computed from the contig alignment in the left-right configuration are valid
    // from the perspective of ref segment1. For this reason, the left-right homology results are
    // just moved over to ref1 here.
    //
    let ref_segment1_breakend_homology_range =
        alt_hap_left_right_component_alignment_info.left_right_breakend_homology_range;
    let ref_segment1_breakend_homology_seq =
        alt_hap_left_right_component_alignment_info.left_right_breakend_homology_seq;

    AltHapRefSegmentAlignmentInfo {
        ref_segment1_alignment,
        ref_segment2_alignment,
        ref_segment1_breakend_homology_range,
        ref_segment1_breakend_homology_seq,
    }
}

/// Collapse small variants near the breakend of a two-region contig alignment
///
/// This transformation ensures that we get consistent representation of breakpoint insertions,
/// even as methods and scoring functions are adjusted for the pairwise contig to alt-hap model
/// alignment.
///
/// The method will collapse any insertion or mismatch within `min_assembly_collapse_anchor` bases
/// of the breakpoint into the breakpoint model.
///
fn collapse_indels_at_breakpoint(
    alt_hap_info: &TwoRegionAltHapInfo,
    assembly_contig_to_alt_hap_alignment: &mut SimpleAlignment,
) {
    let min_assembly_collapse_anchor = 4;
    let (alt_hap_segment_boundary1, alt_hap_segment_boundary2) = (
        alt_hap_info.spacer_range.start,
        alt_hap_info.spacer_range.end,
    );

    let mut segment1_is_anchored = false;
    let mut segment2_is_anchored = false;
    let mut first_edit_after_segment1_anchor_index = None;
    let mut last_edit_before_segment2_anchor_index = None;

    // Scan through alignment and identify the first and last block with a base match at least
    // min_assembly_collapse_anchor in length
    let mut alt_hap_pos = assembly_contig_to_alt_hap_alignment.ref_offset;
    let mut contig_pos = 0;
    for (index, c) in assembly_contig_to_alt_hap_alignment
        .cigar
        .iter()
        .enumerate()
    {
        // Test the anchor criteria
        if let Cigar::Equal(len) = c {
            if *len >= min_assembly_collapse_anchor {
                // Determine if min size occurs within the first region or second region of the alt
                // hap sequence, and record the appropriate cigar index
                if (alt_hap_pos + min_assembly_collapse_anchor as i64) < alt_hap_segment_boundary1 {
                    segment1_is_anchored = true;
                    first_edit_after_segment1_anchor_index = None;
                }

                if alt_hap_pos + *len as i64 - min_assembly_collapse_anchor as i64
                    >= alt_hap_segment_boundary2
                {
                    segment2_is_anchored = true;
                    break;
                }
            }
        } else if let Cigar::Del(_) | Cigar::Diff(_) | Cigar::Ins(_) = c {
            if segment1_is_anchored && first_edit_after_segment1_anchor_index.is_none() {
                first_edit_after_segment1_anchor_index = Some(index);
            }
            last_edit_before_segment2_anchor_index = Some(index);
        }

        update_ref_pos(c, &mut alt_hap_pos);
        update_read_pos(c, &mut contig_pos, true);
    }

    if !(segment1_is_anchored && segment2_is_anchored)
        || first_edit_after_segment1_anchor_index.is_none()
        || last_edit_before_segment2_anchor_index.is_none()
    {
        return;
    }

    let segment1_anchor_index = first_edit_after_segment1_anchor_index.unwrap();
    let segment2_anchor_index = last_edit_before_segment2_anchor_index.unwrap();

    assert!(segment2_anchor_index >= segment1_anchor_index);
    if segment1_anchor_index == segment2_anchor_index
        && let Cigar::Del(_) = assembly_contig_to_alt_hap_alignment.cigar[segment1_anchor_index]
    {
        // In this case no breakpoint adjustment is needed
        return;
    }

    let mut new_cigar = Vec::new();
    let mut update_bp = false;
    let mut total_bp_del = 0;
    let mut total_bp_ins = 0;
    for (index, c) in assembly_contig_to_alt_hap_alignment
        .cigar
        .iter()
        .enumerate()
    {
        if index < segment1_anchor_index || index > segment2_anchor_index {
            if update_bp {
                if total_bp_ins > 0 {
                    new_cigar.push(Cigar::Ins(total_bp_ins as u32));
                }
                if total_bp_del > 0 {
                    new_cigar.push(Cigar::Del(total_bp_del as u32));
                }
                update_bp = false;
                total_bp_del = 0;
                total_bp_ins = 0;
            }
            new_cigar.push(*c);
        } else {
            update_bp = true;
            total_bp_del += get_cigarseg_ref_offset(c);
            total_bp_ins += get_cigarseg_read_offset(c, true);
        }
    }

    std::mem::swap(
        &mut assembly_contig_to_alt_hap_alignment.cigar,
        &mut new_cigar,
    );
}

/// Check alignment quality of one of the two contig segments of a two region SV candidate
///
fn is_good_contig_segment_alignment_quality(
    anchor_min_gap_compressed_identity: f64,
    contig: &[u8],
    alt_hap_segment: &[u8],
    segment_alignment: &SimpleAlignment,
) -> bool {
    let gci = get_gap_compressed_identity_from_alignment(
        segment_alignment.ref_offset,
        &segment_alignment.cigar,
        contig,
        alt_hap_segment,
        true,
    );
    gci >= anchor_min_gap_compressed_identity
}

/// Check alignment quality separately for the left and right sides of the contig alignment
///
fn is_good_left_right_segment_alignment_quality(
    anchor_min_gap_compressed_identity: f64,
    contig: &[u8],
    alt_hap_info: &TwoRegionAltHapInfo,
    left_alignment: &SimpleAlignment,
    right_alignment: &SimpleAlignment,
) -> bool {
    let left_good = is_good_contig_segment_alignment_quality(
        anchor_min_gap_compressed_identity,
        contig,
        &alt_hap_info.alt_hap_seq[..(alt_hap_info.spacer_range.start as usize)],
        left_alignment,
    );
    let right_good = is_good_contig_segment_alignment_quality(
        anchor_min_gap_compressed_identity,
        contig,
        &alt_hap_info.alt_hap_seq[(alt_hap_info.spacer_range.end as usize)..],
        right_alignment,
    );
    left_good && right_good
}

/// If the segment has extention sequence, then clip the extended regeion off of the
/// segment alignment so that it corresponds to the true local segment alignment
///
fn get_clipped_contig_alignment_segment(
    alt_hap_info: &TwoRegionAltHapInfo,
    bp: &Breakpoint,
    breakend_index: usize,
    ref_segment_alignment: &SimpleAlignment,
    debug: bool,
) -> ContigAlignmentSegment {
    let ref_segment = alt_hap_info.get_ref_segment(breakend_index);
    let extension_size = alt_hap_info.get_segment_neighbor_extension_size(breakend_index) as i64;
    let clipped_alignment = if extension_size == 0 {
        ref_segment_alignment
    } else {
        let (left_extension_size, right_extension_size) = {
            if bp.get_breakend(breakend_index).dir == BreakendDirection::LeftAnchor {
                (extension_size, 0)
            } else {
                (0, extension_size)
            }
        };

        let ref_size = (ref_segment.range.size() + extension_size) as usize;
        let mut clipped_alignment = clip_alignment_ref_edges(
            ref_segment_alignment,
            ref_size,
            left_extension_size,
            right_extension_size,
        );
        clipped_alignment.ref_offset -= left_extension_size;

        if debug {
            eprintln!(
                "segment{breakend_index} left_ext/right_ext {left_extension_size}/{right_extension_size}"
            );
            eprintln!("segment{breakend_index} input alignment {ref_segment_alignment:?}");
            eprintln!("segment{breakend_index} output alignment {clipped_alignment:?}");
        }

        &Box::new(clipped_alignment)
    };

    let is_fwd_strand = if breakend_index == 0 {
        true
    } else {
        !alt_hap_info.ref_segment2_revcomp
    };

    ContigAlignmentSegment {
        tid: ref_segment.chrom_index as i32,
        pos: ref_segment.range.start + clipped_alignment.ref_offset,
        cigar: clipped_alignment.cigar.clone(),
        is_fwd_strand,
    }
}

/// Update contig alignments used for bam output to debug the SV refinement routine
///
#[allow(clippy::too_many_arguments)]
fn update_two_ref_segment_contig_alignments(
    cluster_index: usize,
    assembly_index: usize,
    bp: &Breakpoint,
    alt_hap_info: &TwoRegionAltHapInfo,
    supporting_read_count: usize,
    assembly_contig: &AssemblyResultContig,
    ref_segment1_alignment: &SimpleAlignment,
    ref_segment2_alignment: &SimpleAlignment,
    contig_alignments: &mut Vec<ContigAlignmentInfo>,
) {
    let debug = false;
    let sample_index = 0;

    if debug {
        eprintln!("Creating bam record for cluster/assembly {cluster_index}/{assembly_index}");
    }

    let segment1 =
        get_clipped_contig_alignment_segment(alt_hap_info, bp, 0, ref_segment1_alignment, debug);
    let segment2 =
        get_clipped_contig_alignment_segment(alt_hap_info, bp, 1, ref_segment2_alignment, debug);

    /*
    eprintln!(
        "MultiRef cluster/contig ids {cluster_index}/{assembly_index} - segment1: {:?} segment2: {:?}",
        segment1, segment2
    );
     */

    for segment_id in 0..2 {
        let mut seq = assembly_contig.seq.clone();
        let mut high_quality_contig_range = assembly_contig.high_quality_range.clone();

        // Determine if segment2 needs to be revcomp'ed
        //
        if segment_id == 1 && alt_hap_info.ref_segment2_revcomp {
            rev_comp_in_place(&mut seq);

            // Also reverse the HQ contig range to match the sequence:
            let read_size = seq.len() as i64;
            let rev_start = read_size - high_quality_contig_range.end;
            let rev_end = read_size - high_quality_contig_range.start;
            high_quality_contig_range = IntRange::from_pair(rev_start, rev_end);
        }

        contig_alignments.push(ContigAlignmentInfo {
            sample_index,
            cluster_index,
            assembly_index,
            seq,
            segment_id,
            segments: vec![segment1.clone(), segment2.clone()],
            supporting_read_count,
            high_quality_contig_range,
        });
    }
}

fn check_segment_for_invalid_extension_clipping(
    bp: &Breakpoint,
    alt_hap_info: &TwoRegionAltHapInfo,
    breakend_index: usize,
    ref_segment_alignment: &SimpleAlignment,
) -> bool {
    if alt_hap_info.get_segment_neighbor_extension_size(breakend_index) > 0 {
        let segment = get_clipped_contig_alignment_segment(
            alt_hap_info,
            bp,
            breakend_index,
            ref_segment_alignment,
            false,
        );
        !has_aligned_segments(&segment.cigar)
    } else {
        false
    }
}

/// If neighbor extension clipping occurs for this breakpoint, then check if the clipping creates invalid alignment output
///
/// Ideally this should never occur but the check allows us to guard against some rarely occuring complications.
///
fn check_for_invalid_extension_clipping(
    bp: &Breakpoint,
    alt_hap_info: &TwoRegionAltHapInfo,
    ref_segment1_alignment: &SimpleAlignment,
    ref_segment2_alignment: &SimpleAlignment,
) -> bool {
    check_segment_for_invalid_extension_clipping(bp, alt_hap_info, 0, ref_segment1_alignment)
        || check_segment_for_invalid_extension_clipping(bp, alt_hap_info, 1, ref_segment2_alignment)
}

/// Get the zero-indexed ref offset of the position immediately before the breakend
///
/// This assumes that alignment is confined to the region of the breakend and is soft-clipped
/// starting at the breakpoint
///
/// The returned ref offset is in the same coordinate system as alignment
///
fn get_position_before_breakend(
    dir: BreakendDirection,
    alignment: &SimpleAlignment,
    neighbor_extension_size: usize,
) -> i64 {
    alignment.ref_offset - 1
        + match dir {
            BreakendDirection::LeftAnchor => {
                get_cigar_ref_offset(&alignment.cigar) - neighbor_extension_size as i64
            }
            BreakendDirection::RightAnchor => 0,
        }
}

/// Return a 3-tuple of:
/// 1. The position immediately before the ref segment1 breakend in ref segment1 coordinates
/// 2. The position immediately before the ref segment2 breakend in ref segment2 coordinates
/// 3. The range of the insert sequence in read (contig) coordinates.
///
fn get_ref_segment_breakend_offsets_and_insert_range(
    alt_hap_info: &TwoRegionAltHapInfo,
    ref_segment1_alignment: &SimpleAlignment,
    ref_segment2_alignment: &SimpleAlignment,
    bp: &Breakpoint,
) -> (i64, i64, std::ops::Range<usize>) {
    let ref_segment1_breakend_offset = get_position_before_breakend(
        bp.breakend1.dir,
        ref_segment1_alignment,
        alt_hap_info.segment1_neighbor_extension_size,
    );

    let breakend2 = bp.breakend2.as_ref().unwrap();
    let ref_segment2_breakend_offset = get_position_before_breakend(
        breakend2.dir,
        ref_segment2_alignment,
        alt_hap_info.segment2_neighbor_extension_size,
    );

    let insert_range = {
        let ignore_hard_clip = false;
        let (seg1_left_clip_contig_pos, seg1_right_clip_contig_pos, contig_size) =
            get_read_clip_positions(&ref_segment1_alignment.cigar, ignore_hard_clip);
        let (seg2_left_clip_contig_pos, seg2_right_clip_contig_pos, _) =
            get_read_clip_positions(&ref_segment2_alignment.cigar, ignore_hard_clip);
        if alt_hap_info.ref_segment2_after_ref_segment1 {
            let insert_end_pos = if alt_hap_info.ref_segment2_revcomp {
                contig_size - seg2_right_clip_contig_pos
            } else {
                seg2_left_clip_contig_pos
            };
            seg1_right_clip_contig_pos..insert_end_pos
        } else {
            let insert_start_pos = if alt_hap_info.ref_segment2_revcomp {
                contig_size - seg2_left_clip_contig_pos
            } else {
                seg2_right_clip_contig_pos
            };
            insert_start_pos..seg1_left_clip_contig_pos
        }
    };

    (
        ref_segment1_breakend_offset,
        ref_segment2_breakend_offset,
        insert_range,
    )
}

/// Given the contig alignment already translated back to the original 2 reference regions,
/// transform this into a refined SV candidate
///
/// # Arguments
/// `alt_hap_ref_segment_alignment_info` - contig alignment transformed back onto reference chromosomes
/// `bp` - Low resolution breakpoint candidate
///
#[allow(clippy::too_many_arguments)]
fn get_two_region_refined_sv_candidate(
    chrom_list: &ChromList,
    cluster_index: usize,
    assembly_index: usize,
    alt_hap_info: &TwoRegionAltHapInfo,
    alt_hap_ref_segment_alignment_info: &AltHapRefSegmentAlignmentInfo,
    contig_seq: &[u8],
    target_segments: &[GenomeSegment],
    bp: &Breakpoint,
) -> Option<RefinedSV> {
    let debug = false;

    // This is a single sample method, so hard-code sample to zero
    let sample_index = 0;

    let id = get_next_sv_id(sample_index, cluster_index, assembly_index, &mut 0);

    let (ref_segment1_offset, ref_segment2_offset, insert_range) =
        get_ref_segment_breakend_offsets_and_insert_range(
            alt_hap_info,
            &alt_hap_ref_segment_alignment_info.ref_segment1_alignment,
            &alt_hap_ref_segment_alignment_info.ref_segment2_alignment,
            bp,
        );

    let ref1_homology_range =
        &alt_hap_ref_segment_alignment_info.ref_segment1_breakend_homology_range;
    let breakend_homology_length = ref1_homology_range.size();

    let breakend1 = {
        let mut be = bp.breakend1.clone();
        let start =
            alt_hap_info.ref_segment1.range.start + ref_segment1_offset + ref1_homology_range.start;
        let end = start + breakend_homology_length + 1;
        be.segment.range = IntRange::from_pair(start, end);
        be
    };
    let breakend2 = {
        let mut be = bp.breakend2.clone().unwrap();
        let hom_range_adjustment = if bp.same_orientation() {
            -ref1_homology_range.end
        } else {
            ref1_homology_range.start
        };
        let start =
            alt_hap_info.ref_segment2.range.start + ref_segment2_offset + hom_range_adjustment;
        let end = start + breakend_homology_length + 1;
        be.segment.range = IntRange::from_pair(start, end);
        be
    };

    // Sanity check that each breakend has a valid range on its chromosome
    //
    // Anything skipped here should represent an upstream bug in the contig alignment
    // procedure, this final filtration process helps to protect overall stability of
    // the SV caller.
    if !(breakend1.is_valid_range(chrom_list) && breakend2.is_valid_range(chrom_list)) {
        if debug {
            eprintln!("Cluster/Assembly ID: {cluster_index}/{assembly_index} - bp: {bp:?}");
            eprintln!("Filtered for invalid breakend extending off contig");
        }
        return None;
    }

    let contig_pos_before_breakend1 = {
        let contig_pos = (insert_range.start as i64) + ref1_homology_range.start - 1;
        if contig_pos < 0 {
            if debug {
                eprintln!("Cluster/Assembly ID: {cluster_index}/{assembly_index} - bp: {bp:?}");
                eprintln!(
                    "Invalid contig_pos value {contig_pos}. ir {insert_range:?} r1hr {ref1_homology_range:?}"
                );
            }
            return None;
        }
        contig_pos as usize
    };

    let insert_info = if insert_range.is_empty() {
        InsertInfo::NoInsert
    } else {
        InsertInfo::Seq(contig_seq[insert_range].to_vec())
    };

    let breakend2 = Some(breakend2);
    let mut bp = Breakpoint {
        breakend1,
        breakend2,
        insert_info,
        is_precise: true,
        breakend1_neighbor: bp.breakend1_neighbor.clone(),
        breakend2_neighbor: bp.breakend2_neighbor.clone(),
    };
    bp.standardize();

    let rsv = RefinedSV {
        id,
        bp,
        breakend1_homology_seq: alt_hap_ref_segment_alignment_info
            .ref_segment1_breakend_homology_seq
            .clone(),
        single_region_refinement: false,
        contig_pos_before_breakend1,
        assembly_regions: target_segments.to_vec(),
        ext: Default::default(),
        score: Default::default(),
    };

    Some(rsv)
}

fn get_group_haplotype(
    cluster_index: usize,
    assembly_index: usize,
    alt_hap_info: &TwoRegionAltHapInfo,
    supporting_read_count: usize,
    assembly_contig: &AssemblyResultContig,
    alt_hap_ref_segment_alignment_info: &AltHapRefSegmentAlignmentInfo,
) -> SVGroupHaplotype {
    let hap_id = GroupHaplotypeId {
        sample_index: 0,
        cluster_index,
        assembly_index,
    };

    let segment1_alignement = {
        let mut contig_alignment = alt_hap_ref_segment_alignment_info
            .ref_segment1_alignment
            .clone();
        contig_alignment.ref_offset += alt_hap_info.ref_segment1.range.start;
        ClusterAssemblyAlignment {
            contig_seq: assembly_contig.seq.clone(),
            chrom_index: alt_hap_info.ref_segment1.chrom_index,
            contig_alignment,
            high_quality_contig_range: assembly_contig.high_quality_range.clone(),
            is_fwd_strand: true,
        }
    };
    let segment2_alignement = {
        let mut contig_seq = assembly_contig.seq.clone();
        let mut high_quality_contig_range = assembly_contig.high_quality_range.clone();
        if alt_hap_info.ref_segment2_revcomp {
            rev_comp_in_place(&mut contig_seq);
            high_quality_contig_range.reverse(contig_seq.len() as i64);
        }
        let mut contig_alignment = alt_hap_ref_segment_alignment_info
            .ref_segment2_alignment
            .clone();
        contig_alignment.ref_offset += alt_hap_info.ref_segment2.range.start;

        ClusterAssemblyAlignment {
            contig_seq,
            chrom_index: alt_hap_info.ref_segment2.chrom_index,
            contig_alignment,
            high_quality_contig_range,
            is_fwd_strand: !alt_hap_info.ref_segment2_revcomp,
        }
    };

    SVGroupHaplotype {
        hap_id,
        contig_info: ClusterAssembly {
            contig_alignments: vec![segment1_alignement, segment2_alignement],
            supporting_read_count,
        },
    }
}

/// # Arguments
/// `bp` - Low resolution breakpoint candidate
///
/// Output is a 2-tuple of:
/// 1. Contig alignment information, which will be written out to a bam file
/// 2. Refined SVs produced from contig alignments
///
#[allow(clippy::too_many_arguments)]
fn align_multi_ref_region_assemblies(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segments: &[GenomeSegment],
    assemblies: &[AssemblyResult],
    expected_cn_info: SVLocusExpectedCNInfo,
    cluster_index: usize,
    bpc: &BreakpointCluster,
) -> (Vec<ContigAlignmentInfo>, Option<SVGroup>) {
    let bp = &bpc.breakpoint;

    assert_eq!(target_segments.len(), 2);
    assert!(bp.breakend2.is_some());

    // Multi-region contigs are never exceptionally long
    let is_long_contig = false;
    let mut aligner = PairwiseAligner::new_large_indel_aligner(is_long_contig);

    let alt_hap_info = get_two_region_alt_hap_info(
        refine_settings,
        reference,
        chrom_list,
        target_segments,
        cluster_index,
        bpc,
    );

    let debug = false;
    if debug {
        eprintln!("MultiRef cluster id {cluster_index} - alt_hap_info: {alt_hap_info:?}");
    }

    aligner.set_text(&alt_hap_info.alt_hap_seq);

    // output structures
    let mut contig_alignments = Vec::new();
    let mut refined_svs = Vec::new();
    let mut group_haplotypes = Vec::new();

    for (assembly_index, assembly) in assemblies.iter().enumerate() {
        if assembly.supporting_read_count() < refine_settings.min_assembly_read_support {
            if debug {
                eprintln!(
                    "MultiRef cluster_id/assembly_id {cluster_index}/{assembly_index} - filtered for low read support: {}",
                    assembly.supporting_read_count()
                );
            }
            continue;
        }

        // Loop over multiple contigs for each assembly -- typically this means there's a contig with longer reference flanks to try out if the smaller one fails
        for (contig_index, assembly_contig) in assembly.contigs.iter().enumerate() {
            if debug {
                eprintln!(
                    "MultiRef cluster_id/assembly_id/contig_id {cluster_index}/{assembly_index}/{contig_index} - START"
                );
            }

            let (_aligment_score, mut assembly_contig_to_alt_hap_alignment) =
                match aligner.align(&assembly_contig.seq) {
                    (score, Some(x)) => {
                        if debug {
                            eprintln!("assembly_contig_to_alt_hap_alignment: {x:?}");
                        }
                        (score, x)
                    }
                    (_, None) => {
                        if debug {
                            eprintln!("Filtered for poor contig to alt hap alignment");
                        }
                        // TODO add warning or statistics to track these filtration cases
                        continue;
                    }
                };

            if refine_settings.collapse_indels_at_breakpoint {
                collapse_indels_at_breakpoint(
                    &alt_hap_info,
                    &mut assembly_contig_to_alt_hap_alignment,
                );
            }

            let alt_hap_left_right_component_alignment_info =
                transform_alt_hap_alignment_into_left_right_components(
                    refine_settings.min_assembly_edge_anchor,
                    &assembly_contig_to_alt_hap_alignment,
                    &alt_hap_info,
                    &assembly_contig.seq,
                );

            let alt_hap_left_right_component_alignment_info =
                match alt_hap_left_right_component_alignment_info {
                    Some(x) => {
                        if debug {
                            eprintln!("alt_hap_left_right_component_alignment_info: {x:?}");
                        }
                        x
                    }
                    None => {
                        if debug {
                            eprintln!("Filtered for poor alt_hap left-right segmentation");
                        }
                        // TODO add warning or statistics to track these filtration cases
                        continue;
                    }
                };

            // Check the alignment quality separately for the left and right sides of the contig aligment, and require both
            // to reach a minimum quality level to proceed
            if !is_good_left_right_segment_alignment_quality(
                refine_settings.anchor_min_gap_compressed_identity,
                &assembly_contig.seq,
                &alt_hap_info,
                &alt_hap_left_right_component_alignment_info.left_alignment,
                &alt_hap_left_right_component_alignment_info.right_alignment,
            ) {
                if debug {
                    eprintln!("Filtered for poor alt_hap left or right side alignment quality");
                }
                continue;
            }

            let alt_hap_ref_segment_alignment_info =
                transform_left_right_alignments_back_to_ref_segments(
                    &alt_hap_info,
                    alt_hap_left_right_component_alignment_info,
                );

            if debug {
                eprintln!(
                    "alt_hap_ref_segment_alignment_info: {alt_hap_ref_segment_alignment_info:?}"
                );
            }

            // Add additional check here to make sure we still meet the minimum alignment size after any
            // segment1 or segment2 extension sequence is accounted for
            //
            if check_for_invalid_extension_clipping(
                bp,
                &alt_hap_info,
                &alt_hap_ref_segment_alignment_info.ref_segment1_alignment,
                &alt_hap_ref_segment_alignment_info.ref_segment2_alignment,
            ) {
                if debug {
                    eprintln!("Filtered for invalid alignment clipping of neighbor extension");
                }
                continue;
            }

            // Convert alignment into refined SV candidate -- note that for multi-region analyses,
            // we're assuming/enforcing that only one assembly will contain an SV, so will
            // filter out once a suitable candidate is found, but we can still produce additional
            // contig alignments for debugging so don't break out of the loop after the first SV
            // is found.
            //
            if refined_svs.is_empty() {
                let refined_sv = get_two_region_refined_sv_candidate(
                    chrom_list,
                    cluster_index,
                    assembly_index,
                    &alt_hap_info,
                    &alt_hap_ref_segment_alignment_info,
                    &assembly_contig.seq,
                    target_segments,
                    bp,
                );
                let refined_sv = match refined_sv {
                    Some(x) => x,
                    None => continue,
                };

                refined_svs.push(refined_sv);

                // Note we don't need to capture a correct sv_group_haplotype here since this will just be spilled to disk and read back in
                // for genotyping in the next step, but completing this for now until we can entirely remove it.
                //
                let group_haplotype = get_group_haplotype(
                    cluster_index,
                    assembly_index,
                    &alt_hap_info,
                    assembly.supporting_read_count(),
                    assembly_contig,
                    &alt_hap_ref_segment_alignment_info,
                );

                group_haplotypes.push(group_haplotype);
            }

            update_two_ref_segment_contig_alignments(
                cluster_index,
                assembly_index,
                bp,
                &alt_hap_info,
                assembly.supporting_read_count(),
                assembly_contig,
                &alt_hap_ref_segment_alignment_info.ref_segment1_alignment,
                &alt_hap_ref_segment_alignment_info.ref_segment2_alignment,
                &mut contig_alignments,
            );

            // This contig length was not filtered so we break out of the contig length loop, but not the assembly loop
            break;
        }
    }

    let sv_group = if refined_svs.is_empty() {
        None
    } else {
        let sample_haplotype_list = vec![vec![0]];
        let sample_expected_cn_info = vec![expected_cn_info];
        let sv_haplotype_map = refined_svs.iter().map(|x| x.id.assembly_index).collect();
        let sv_group = SVGroup {
            group_regions: refined_svs[0].assembly_regions.clone(),
            group_haplotypes,
            sample_haplotype_list,
            sample_expected_cn_info,
            sv_haplotype_map,
            refined_svs,
        };
        Some(sv_group)
    };

    (contig_alignments, sv_group)
}

/// Extract the appropriate sections of the reference sequence and align assemblies to these
///
/// # Arguments
/// `bp` - Low resolution breakpoint candidate
///
/// Output is a 2-tuple of:
/// 1. Contig alignment information, which will be written out to a bam file
/// 2. Refined SVs produced from contig alignments
///
#[allow(clippy::too_many_arguments)]
pub(super) fn align_assemblies(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    target_segments: &[GenomeSegment],
    assemblies: &[AssemblyResult],
    expected_cn_info: SVLocusExpectedCNInfo,
    bpc: &BreakpointCluster,
    cluster_index: usize,
    is_large_insertion: bool,
) -> (Vec<ContigAlignmentInfo>, Option<SVGroup>) {
    // 1. Determine type of low-res breakpoint
    //   - Close & conlinear breakpoints can be represented with a single segment of the reference sequence,
    //   - More distant breakpoints have to be represented by the reference sequence around each
    //   breakend, plus a low-res model of the breakpoint allele obtained by concatenating the two
    //   regions together, possibly with one of the regions reverse complemented.
    //
    // 2. Build the corresponding reference target(s)
    //
    // 3. Align assemblies to reference target
    //
    // 4. Analyze alignments for evidence of SV-scale variants and convert these into high-res
    // candidates

    // 1. Use the same target region criteria that were used during assembly, so assembly of a single
    // non bp haplotype also means we're aligning back to a single reference region.
    //
    let use_single_ref_region = target_segments.len() == 1;

    if use_single_ref_region {
        align_single_ref_region_assemblies(
            refine_settings,
            reference,
            chrom_list,
            target_segments,
            assemblies,
            expected_cn_info,
            cluster_index,
            is_large_insertion,
        )
    } else {
        align_multi_ref_region_assemblies(
            refine_settings,
            reference,
            chrom_list,
            target_segments,
            assemblies,
            expected_cn_info,
            cluster_index,
            bpc,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use rust_vc_utils::bam_utils::cigar::get_cigar_from_string;

    #[test]
    fn test_collapse_indels_at_breakpoint() {
        let alt_hap_info = TwoRegionAltHapInfo {
            ref_segment1: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(0, 300),
            },
            ref_segment2: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(1000, 1300),
            },
            segment1_neighbor_extension_size: 0,
            segment2_neighbor_extension_size: 0,
            ref_segment2_revcomp: false,
            ref_segment2_after_ref_segment1: true,
            spacer_range: IntRange::from_pair(300, 300 + 100),
            left_boundary_pin: false,
            right_boundary_pin: false,
            alt_hap_seq: Vec::new(),
            extended_alt_hap_seq: Vec::new(),
        };

        // Test noop case:
        let in_cigar = "100=300D100=";
        let mut align = SimpleAlignment {
            ref_offset: 100,
            cigar: get_cigar_from_string(in_cigar),
        };

        collapse_indels_at_breakpoint(&alt_hap_info, &mut align);
        assert_eq!(CigarString(align.cigar.clone()).to_string(), in_cigar);

        // Test preceding insertion:
        let in_cigar = "87=10I3=300D100=";
        let mut align = SimpleAlignment {
            ref_offset: 100,
            cigar: get_cigar_from_string(in_cigar),
        };

        collapse_indels_at_breakpoint(&alt_hap_info, &mut align);

        assert_eq!(
            CigarString(align.cigar.clone()).to_string(),
            "87=13I303D100="
        );

        // Test following insertion:
        let in_cigar = "100=300D3=10I87=";
        let mut align = SimpleAlignment {
            ref_offset: 100,
            cigar: get_cigar_from_string(in_cigar),
        };

        collapse_indels_at_breakpoint(&alt_hap_info, &mut align);

        assert_eq!(
            CigarString(align.cigar.clone()).to_string(),
            "100=13I303D87="
        );

        // Test leading and following mismatch:
        let in_cigar = "96=1X3=300D3=1X96=";
        let mut align = SimpleAlignment {
            ref_offset: 100,
            cigar: get_cigar_from_string(in_cigar),
        };

        collapse_indels_at_breakpoint(&alt_hap_info, &mut align);

        assert_eq!(CigarString(align.cigar.clone()).to_string(), "96=8I308D96=");
    }
}
