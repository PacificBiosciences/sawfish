use rust_vc_utils::cigar::{get_cigar_ref_offset, get_read_clip_positions};
use rust_vc_utils::int_range::IntRange;

use super::{AssemblyResult, AssemblyResultContig, RefineSVSettings};
use crate::refine_sv::trimmed_reads::TrimmedReadInfo;
use crate::simple_alignment::SimpleAlignment;

/// Tracks consensus sequences at various stages of processing, together with their source reads
pub(super) struct ConsensusSequence {
    pub source_reads: Vec<usize>,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug)]
pub(super) struct TrimmedReadConsensusAlignmentInfo {
    /// Distance from the left edge of the ref sequence where alignment starts
    pub ref_left_clip: usize,

    /// Distance from the right edge of the ref sequence where alignment ends
    pub ref_right_clip: usize,

    /// Distance from the left edge of the read sequence where alignment starts
    pub read_left_clip: usize,

    /// Distance from the right edge of the read sequence where alignment ends
    pub read_right_clip: usize,
}

pub(super) struct ConsensusSequenceInfo {
    pub consensus_seq: ConsensusSequence,
    pub source_read_alignment_info: Vec<TrimmedReadConsensusAlignmentInfo>,
}

pub fn get_consensus_alignment_info(
    consensus_len: usize,
    alignment: &SimpleAlignment,
) -> TrimmedReadConsensusAlignmentInfo {
    let ignore_hard_clip = false;
    let (read_start, read_end, read_size) =
        get_read_clip_positions(&alignment.cigar, ignore_hard_clip);
    assert!(read_end <= read_size);

    let ref_offset = get_cigar_ref_offset(&alignment.cigar);
    let total_offset = (alignment.ref_offset + ref_offset) as usize;
    assert!(total_offset <= consensus_len);

    TrimmedReadConsensusAlignmentInfo {
        ref_left_clip: alignment.ref_offset as usize,
        ref_right_clip: consensus_len - total_offset,
        read_left_clip: read_start,
        read_right_clip: read_size - read_end,
    }
}

/// Select the best read to use for the fake flank prefix sequence
///
/// The best read is defined as the one with the longest prefix length
///
/// Returns index of the best prefix read from the source_reads list, or None
///
fn get_best_prefix_index(
    trimmed_reads: &[TrimmedReadInfo],
    source_reads: &[usize],
    source_read_alignment_info: &[TrimmedReadConsensusAlignmentInfo],
    max_ref_clip: usize,
) -> Option<usize> {
    let debug_prefix_selection = false;

    let mut best_source_read_index = None;
    for (source_read_index, trimmed_read_index) in source_reads.iter().copied().enumerate() {
        let read_info = &trimmed_reads[trimmed_read_index];
        let align_info = &source_read_alignment_info[source_read_index];

        if debug_prefix_selection {
            eprintln!(
                "Prefix qname: {} source_read_index {source_read_index} trimmmed_read_index {trimmed_read_index}",
                read_info.qname
            );
            eprintln!(
                "prefix_len {} align_info {:?}",
                read_info.prefix_len(),
                align_info
            );
        }

        // Filter to ensure that the left end of the trimmed read completely aligns to the start of the
        // assembly contig:
        if align_info.read_left_clip > 0 {
            continue;
        }

        // Filter reads that don't align close enough to the left end of the contig
        if align_info.ref_left_clip > max_ref_clip {
            continue;
        }

        let mut update = false;
        if let Some(best_source_read_index) = best_source_read_index {
            let best_trimmed_read_index = source_reads[best_source_read_index];
            let best_read_info: &TrimmedReadInfo = &trimmed_reads[best_trimmed_read_index];
            if read_info.prefix_len() > best_read_info.prefix_len() {
                update = true;
            }
        } else {
            update = true;
        }

        if update {
            if debug_prefix_selection {
                eprintln!("selected as prefix read");
            }
            best_source_read_index = Some(source_read_index);
        }
    }
    best_source_read_index
}

/// Select the best read to use for the fake flank suffix sequence
///
/// The best read is defined as the one with the longest suffix length
///
/// Returns index of the best suffix read from the source_reads list, or None
///
fn get_best_suffix_index(
    trimmed_reads: &[TrimmedReadInfo],
    source_reads: &[usize],
    source_read_alignment_info: &[TrimmedReadConsensusAlignmentInfo],
    max_ref_clip: usize,
) -> Option<usize> {
    let debug_suffix_selection = false;

    let mut best_source_read_index = None;
    for (source_read_index, trimmed_read_index) in source_reads.iter().copied().enumerate() {
        let read_info = &trimmed_reads[trimmed_read_index];
        let align_info = &source_read_alignment_info[source_read_index];

        if debug_suffix_selection {
            eprintln!(
                "Suffix qname: {} source_read_index {source_read_index} trimmmed_read_index {trimmed_read_index}",
                read_info.qname
            );
            eprintln!(
                "suffix_len {} align_info {:?}",
                read_info.suffix_len(),
                align_info
            );
        }

        // Filter to ensure that the right end of the read completely aligns to the end of the
        // assembly contig:
        if align_info.read_right_clip > 0 {
            continue;
        }

        // Filter reads that don't align close enough to the right end of the contig
        if align_info.ref_right_clip > max_ref_clip {
            continue;
        }

        let mut update = false;
        if let Some(best_source_read_index) = best_source_read_index {
            let best_trimmed_read_index = source_reads[best_source_read_index];
            let best_read_info: &TrimmedReadInfo = &trimmed_reads[best_trimmed_read_index];
            if read_info.suffix_len() > best_read_info.suffix_len() {
                update = true;
            }
        } else {
            update = true;
        }

        if update {
            if debug_suffix_selection {
                eprintln!("selected as suffix read");
            }
            best_source_read_index = Some(source_read_index);
        }
    }
    best_source_read_index
}

/// Produce a final assembly result
///
/// One of the core operations here is to extend the flanks of the actual consensus assembly by
/// using the appropriate section of single prefix and suffix reads. The flanks can be less accurate
/// and just help correctly align the assembly in the context of high breakend homology.
///
/// This function should use component computations from the assembly such that it is largely
/// independent of the exact assembly method being used.
///
/// # Arguments
///
/// * `max_contig_fake_flank_size` - maximum amount of fake flank that can be added on to each end of the
///   asseembled contig segment
///
/// * `trimmed_reads` - full list of read segments cut out around the target assembly region
///
/// * `source_reads` - Indexes of the trimmed reads used to assemble the contig in question
///
/// * `core_consensus` - The assembled contig sequence. It is the 'core' here because approximated
///   left and right flanks will be added onto it.
///
/// * `source_read_alignment_info` - Summary information from the source read to consensus alignment
///
#[allow(clippy::too_many_arguments)]
fn process_consensus_info_to_assembly_result_contig(
    max_contig_fake_flank_size: i64,
    max_ref_clip_for_fake_flank_selection: usize,
    cluster_index: usize,
    assembly_index: usize,
    trimmed_reads: &[TrimmedReadInfo],
    source_reads: &[usize],
    core_consensus: Vec<u8>,
    source_read_alignment_info: Vec<TrimmedReadConsensusAlignmentInfo>,
) -> AssemblyResultContig {
    // TODO Instead of a global max_ref_clip_for_fake_flank_selection, we could also take the mode
    // of the ref distance from the end, and allow up to that amount -- this would guard against the
    // case of 1 strange extension alignment messing things up and being selected as the
    // prefix/suffix even though it's a terrible fake flank choice.
    //
    let best_prefix_source_read_index = get_best_prefix_index(
        trimmed_reads,
        source_reads,
        &source_read_alignment_info,
        max_ref_clip_for_fake_flank_selection,
    );

    let best_suffix_source_read_index = get_best_suffix_index(
        trimmed_reads,
        source_reads,
        &source_read_alignment_info,
        max_ref_clip_for_fake_flank_selection,
    );

    let prefix_suffix_debug = false;
    if prefix_suffix_debug {
        eprintln!(
            "ClusterIndex/AssemblyIndex {cluster_index}/{assembly_index} contig prefix/suffix info:"
        );
        if let Some(index) = best_prefix_source_read_index {
            let trimmed_read_index = source_reads[index];
            let read_info = &trimmed_reads[trimmed_read_index];
            eprintln!(
                "Best prefix read: {} prefix_size: {} suffix_size: {}",
                read_info.qname,
                read_info.prefix_len(),
                read_info.suffix_len(),
            );
            let align_info = &source_read_alignment_info[index];
            eprintln!(
                "Best prefix read: read left/right: {}/{}",
                align_info.read_left_clip, align_info.read_right_clip,
            );
            eprintln!(
                "Best prefix read: ref left/right: {}/{}",
                align_info.ref_left_clip, align_info.ref_right_clip,
            );
        } else {
            eprintln!("No prefix read identified");
        }
        if let Some(index) = best_suffix_source_read_index {
            let trimmed_read_index = source_reads[index];
            let read_info = &trimmed_reads[trimmed_read_index];
            eprintln!(
                "Best suffix read: {} prefix_size: {} suffix_size: {}",
                read_info.qname,
                read_info.prefix_len(),
                read_info.suffix_len(),
            );
            let align_info = &source_read_alignment_info[index];
            eprintln!(
                "Best suffix read: read left/right: {}/{}",
                align_info.read_left_clip, align_info.read_right_clip,
            );
            eprintln!(
                "Best suffix read: ref left/right: {}/{}",
                align_info.ref_left_clip, align_info.ref_right_clip,
            );
        } else {
            eprintln!("No suffix read identified");
        }
    }

    let mut prefix_read_flank = {
        if let Some(index) = best_prefix_source_read_index {
            let trimmed_read_index = source_reads[index];
            let read_info = &trimmed_reads[trimmed_read_index];
            let prefix_len = read_info.prefix_len() as i64;
            let prefix_start = std::cmp::max(prefix_len - max_contig_fake_flank_size, 0) as usize;
            read_info.read_prefix.as_slice()[prefix_start..].to_vec()
        } else {
            Vec::new()
        }
    };

    let mut suffix_read_flank = {
        if let Some(index) = best_suffix_source_read_index {
            let trimmed_read_index = source_reads[index];
            let read_info = &trimmed_reads[trimmed_read_index];
            let suffix_len = read_info.suffix_len() as i64;
            let suffix_end = std::cmp::min(max_contig_fake_flank_size, suffix_len) as usize;
            read_info.read_suffix.as_slice()[..suffix_end].to_vec()
        } else {
            Vec::new()
        }
    };

    // Trim core consensus to reflect end points of prefix and suffix reads if they didn't extend
    // all the way to the end
    //
    let mut trim_core_start = {
        if let Some(index) = best_prefix_source_read_index {
            source_read_alignment_info[index].ref_left_clip
        } else {
            0
        }
    };

    let mut trim_core_end_offset = {
        if let Some(index) = best_suffix_source_read_index {
            source_read_alignment_info[index].ref_right_clip
        } else {
            0
        }
    };

    if prefix_suffix_debug {
        eprintln!(
            "trim_core_start: {trim_core_start} trim_core_end_offset: {trim_core_end_offset} core_consensus_len: {}",
            core_consensus.len()
        );
    }

    // Sanity check to bail out if this gets messy.
    if (trim_core_start + trim_core_end_offset) > (core_consensus.len() / 2) {
        trim_core_start = 0;
        trim_core_end_offset = 0;
        prefix_read_flank = Vec::new();
        suffix_read_flank = Vec::new();
        if prefix_suffix_debug {
            eprintln!("bailing out on prefix/suffix flanks");
        }
    }

    let trim_core_end = core_consensus.len() - trim_core_end_offset;

    let trimmed_core_consensus = &core_consensus[trim_core_start..trim_core_end];

    let core_start = prefix_read_flank.len() as i64;
    let core_end = core_start + trimmed_core_consensus.len() as i64;
    let high_quality_range = IntRange::from_pair(core_start, core_end);
    let consensus = {
        prefix_read_flank.extend(trimmed_core_consensus.iter());
        prefix_read_flank.extend(suffix_read_flank.iter());
        prefix_read_flank
    };

    AssemblyResultContig {
        high_quality_range,
        seq: consensus,
        fake_flank_size: max_contig_fake_flank_size,
    }
}

fn process_consensus_sequence_to_assembly_result(
    fake_flank_sizes: &[i64],
    max_ref_clip_for_fake_flank_selection: usize,
    cluster_index: usize,
    assembly_index: usize,
    trimmed_reads: &[TrimmedReadInfo],
    consensus_sequence_info: ConsensusSequenceInfo,
) -> AssemblyResult {
    let mut assembly_result_contigs = Vec::new();
    for fake_flank_size in fake_flank_sizes.iter() {
        let assembly_contig = process_consensus_info_to_assembly_result_contig(
            *fake_flank_size,
            max_ref_clip_for_fake_flank_selection,
            cluster_index,
            assembly_index,
            trimmed_reads,
            &consensus_sequence_info.consensus_seq.source_reads,
            consensus_sequence_info.consensus_seq.seq.clone(),
            consensus_sequence_info.source_read_alignment_info.clone(),
        );
        assembly_result_contigs.push(assembly_contig);
    }

    AssemblyResult {
        source_reads: consensus_sequence_info.consensus_seq.source_reads.clone(),
        contigs: assembly_result_contigs,
    }
}

pub(super) fn add_flanks(
    refine_settings: &RefineSVSettings,
    cluster_index: usize,
    max_assembly_count: usize,
    trimmed_reads: &[TrimmedReadInfo],
    is_multi_region: bool,
    consensus_seq_infos: Vec<ConsensusSequenceInfo>,
) -> Vec<AssemblyResult> {
    let mut fake_flanks_sizes = vec![refine_settings.max_contig_fake_flank_size];
    if is_multi_region {
        fake_flanks_sizes.push(refine_settings.breakpoint_retry_max_contig_fake_flank_size);
    }

    consensus_seq_infos
        .into_iter()
        .take(max_assembly_count)
        .enumerate()
        .map(|(assembly_index, x)| {
            process_consensus_sequence_to_assembly_result(
                &fake_flanks_sizes,
                refine_settings.max_ref_clip_for_fake_flank_selection,
                cluster_index,
                assembly_index,
                trimmed_reads,
                x,
            )
        })
        .collect()
}
