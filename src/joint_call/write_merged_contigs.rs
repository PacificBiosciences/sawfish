use rust_vc_utils::cigar::clip_alignment_read_edges;

use crate::contig_output::{ContigAlignmentInfo, ContigAlignmentSegment, write_contig_alignments};
use crate::joint_call::read_sample_data::SharedJointCallData;
use crate::joint_call::{JointCallSettings, SharedSettings};
use crate::sv_group::SVGroup;

/// Convert the multi-sample merged SVGroups haplotypes into the contig
/// alignment format used to write the alingments out to BAM
///
/// If 'cleanup_mode' is true, this operation is arranged to simplify some of the internal
/// contig complexity for end-users to reason about their SV calls. Otherwise, the full
/// haplotype information is converted as directly as possible to facilitate debugging.
///
fn convert_merged_sv_groups_to_contig_alignments(
    merged_sv_groups: &[SVGroup],
    cleanup_mode: bool,
) -> Vec<ContigAlignmentInfo> {
    let mut contig_alignments = Vec::new();

    for sv_group in merged_sv_groups.iter() {
        let haplotype_to_sv_map = sv_group.haplotype_to_sv_map();
        for (haplotype_index, sv_group_haplotype) in sv_group.group_haplotypes.iter().enumerate() {
            // Check if this group haplotype is one that doesn't have an associated SV, and skip the
            // haplotype if this is the case.
            //
            // Note that non-SV haplotypes are maintained in some cases to provide an accurate alignment
            // target for the second haplotype in a region, which can improve genotype accuracy. For the
            // user-facing contig alignment output we want to simplify this complexity out and only show
            // the contigs corresponding to the SV call ouput.
            //
            if cleanup_mode && haplotype_to_sv_map[haplotype_index].is_empty() {
                continue;
            }

            let hap_id = &sv_group_haplotype.hap_id;
            let sample_index = hap_id.sample_index;
            let cluster_index = hap_id.cluster_index;
            let assembly_index = hap_id.assembly_index;
            let supporting_read_count = sv_group_haplotype.contig_info.supporting_read_count;

            let segment_count = sv_group_haplotype.contig_info.contig_alignments.len();

            // Build all segments:
            let segments = {
                let mut x = Vec::new();
                for segment_id in 0..segment_count {
                    let segment_alignment =
                        &sv_group_haplotype.contig_info.contig_alignments[segment_id];

                    // Trim the first and last N bases of the alignment
                    //
                    // This is done to remove some of the alignemnt noise created by they way we're trying to use WFA for
                    // local contig to reference alignment, which tends to create small alignment segments at the ends of
                    // each read. Ideally we should clean up the artifact instead of trying to remove it here
                    //
                    let (trimmed_cigar, trimmed_ref_offset) = if cleanup_mode {
                        let contig_alignment_edge_trim_size = 20;
                        clip_alignment_read_edges(
                            &segment_alignment.contig_alignment.cigar,
                            contig_alignment_edge_trim_size,
                            contig_alignment_edge_trim_size,
                        )
                    } else {
                        (segment_alignment.contig_alignment.cigar.clone(), 0)
                    };

                    let contig_alignment_segment = ContigAlignmentSegment {
                        tid: segment_alignment.chrom_index as i32,
                        pos: segment_alignment.contig_alignment.ref_offset + trimmed_ref_offset,
                        cigar: trimmed_cigar,
                        is_fwd_strand: segment_alignment.is_fwd_strand,
                    };
                    x.push(contig_alignment_segment);
                }
                x
            };

            // Push contig_alignment for each segment:
            for segment_id in 0..segment_count {
                let segment_alignment =
                    &sv_group_haplotype.contig_info.contig_alignments[segment_id];

                let segment_contig_alignment = ContigAlignmentInfo {
                    sample_index,
                    cluster_index,
                    assembly_index,
                    seq: segment_alignment.contig_seq.clone(),
                    segment_id,
                    segments: segments.clone(),
                    supporting_read_count,
                    high_quality_contig_range: segment_alignment.high_quality_contig_range.clone(),
                };
                contig_alignments.push(segment_contig_alignment);
            }
        }
    }

    contig_alignments
}

pub fn write_merged_contig_alignments(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    shared_data: &SharedJointCallData,
    merged_sv_groups: &[SVGroup],
) {
    let cleanup_mode = !settings.no_contig_cleanup;

    let contig_alignments =
        convert_merged_sv_groups_to_contig_alignments(merged_sv_groups, cleanup_mode);
    write_contig_alignments(
        &settings.output_dir,
        shared_settings.thread_count,
        &shared_data.chrom_list,
        "merged",
        contig_alignments,
    );
}
