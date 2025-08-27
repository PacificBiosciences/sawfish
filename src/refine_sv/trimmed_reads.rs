use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use rust_htslib::bam::{self, Read, record::CigarString};
use rust_vc_utils::bam_utils::aux::{get_optional_float_aux_tag, is_aux_tag_found};
use rust_vc_utils::bam_utils::cigar::{
    get_cigar_ref_offset, is_hard_clipped, update_ref_and_hard_clipped_read_pos,
};
use rust_vc_utils::bam_utils::filter_out_alignment_record;
use rust_vc_utils::{ChromList, GenomeRef};
use unwrap::unwrap;

use crate::bam_sa_parser::{SeqOrderSplitReadSegment, get_seq_order_read_split_segments};
use crate::bam_utils::{
    LargeInsertionSoftClipState, bam_fetch_segment, get_gap_compressed_identity,
    get_simplified_dna_seq, test_read_for_large_insertion_soft_clip,
    translate_ref_range_to_hardclipped_read_range,
};
use crate::breakpoint::{Breakend, BreakendDirection, Breakpoint, BreakpointCluster, InsertInfo};
use crate::genome_segment::{GenomeSegment, IntRange, get_int_range_distance};
use crate::refine_sv::RefineSVSettings;
use crate::utils::remove_sorted_indices;

/// Expand read range by the expansion factor and extract the corresponding read sequence
///
///
fn get_trimmed_read_from_read_range(
    expand_by: usize,
    record: &bam::Record,
    read_range: &IntRange,
    large_insertion_clip_state: LargeInsertionSoftClipState,
) -> TrimmedReadInfo {
    let read_range = IntRange::from_pair(
        std::cmp::max(0, read_range.start - expand_by as i64),
        std::cmp::min(record.seq_len() as i64, read_range.end + expand_by as i64),
    );

    let trim_start = read_range.start as usize;
    let trim_end = read_range.end as usize;
    let read = get_simplified_dna_seq(record);
    let read_prefix = read[..trim_start].to_vec();
    let trimmed_read = read[trim_start..trim_end].to_vec();
    let read_suffix = read[trim_end..].to_vec();

    // Read quality has experimental uses by one of the alternate assemblers but isn't critical for performance, so allow this tag
    // to be optional. The most common reason an rq tag will be missing is read mapping from fastqs.
    let default_read_quality = 0.99;
    let read_quality = get_optional_float_aux_tag(record, b"rq").unwrap_or(default_read_quality);

    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    TrimmedReadInfo {
        qname,
        is_supplementary: record.is_supplementary(),
        read_range,
        trimmed_read,
        read_prefix,
        read_suffix,
        large_insertion_clip_state,
        is_trimmed_read_nul_terminated: false,
        read_quality,
    }
}

/// Sort large insertion trimmed reads so that an extension-style POA is more likely to succeed
///
/// 1. Primary sort order:
///    1. Right side soft-clip reads
///    2. Non soft-clip reads
///    3. Left-side soft-clip reads
/// 2. Right side soft-clip reads should be sorted from longest to shortest
/// 3. Left side soft-clip reads should be sorted from longest to shortest
///
fn sort_large_insertion_reads(a: &TrimmedReadInfo, b: &TrimmedReadInfo) -> Ordering {
    let order1 = a
        .large_insertion_clip_state
        .cmp(&b.large_insertion_clip_state);
    if let Ordering::Equal = order1 {
        use LargeInsertionSoftClipState::*;
        match a.large_insertion_clip_state {
            Right => b.trimmed_read.len().cmp(&a.trimmed_read.len()),
            Left => b.trimmed_read.len().cmp(&a.trimmed_read.len()),
            _ => order1,
        }
    } else {
        order1
    }
}

fn get_large_del_candidate_read_range(record: &bam::Record, bp: &Breakpoint) -> Option<IntRange> {
    let debug = false;

    // As a pre-condition of this method, we already know that bp has an indel-type pattern, so see
    // if (a) the read supports the bp in the cigar string and (b) if the corresponding ref
    // coordinates can be parsed out
    //

    // The criteria to match the expected del will be very loose, the idea being that a FP inclusion
    // will just be pushed back out at the assembly stage.
    //
    // Require that the size and start position is within 50% of the target del size, up to a max
    // of MAX_DEL_DIST bases
    //
    const DEL_MATCH_FACTOR: f32 = 0.5;
    const MAX_DEL_DIST: i64 = 100;
    let del_start_range = &bp.breakend1.segment.range;
    let del_size = bp.breakend2.as_ref().unwrap().segment.range.center() - del_start_range.center();
    let del_match_dist = std::cmp::min(MAX_DEL_DIST, (del_size as f32 * DEL_MATCH_FACTOR) as i64);

    if debug {
        eprintln!(
            "del_start_range: {del_start_range:?} del_size: {del_size} del_match_dist: {del_match_dist}"
        );
    }

    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    for c in record.cigar().iter() {
        use bam::record::Cigar::*;
        if let Del(len) = c {
            let del_start_dist =
                get_int_range_distance(del_start_range, &IntRange::from_int(ref_pos)) as i64;
            // Is this del close to the expected breakpoint position and size?
            let is_del_close = (del_start_dist <= del_match_dist)
                && ((*len as i64 - del_size).abs() <= del_match_dist);

            if debug {
                eprintln!(
                    "del_start_dist; {del_start_dist} len: {len} is_del_close: {is_del_close}"
                );
            }

            if is_del_close {
                return Some(IntRange::from_pair(read_pos as i64 - 1, read_pos as i64));
            }
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }

    None
}

fn get_breakend2_split_segment(
    segments: &[SeqOrderSplitReadSegment],
    is_higher: bool,
) -> Option<&SeqOrderSplitReadSegment> {
    for (i, s) in segments.iter().enumerate() {
        if s.from_primary_bam_record {
            if is_higher {
                if i + 1 < segments.len() {
                    return Some(&segments[i + 1]);
                }
            } else if i > 0 {
                return Some(&segments[i - 1]);
            }
            break;
        }
    }
    None
}

fn is_breakend_close(
    max_dist: i64,
    be: &Breakend,
    seg_chrom_index: usize,
    seg_ref_pos: i64,
    seg_cigar: &CigarString,
) -> bool {
    if be.segment.chrom_index != seg_chrom_index {
        false
    } else {
        let query_be_match_ref_pos = seg_ref_pos
            + if be.dir == BreakendDirection::LeftAnchor {
                get_cigar_ref_offset(seg_cigar) - 1
            } else {
                0
            };
        (query_be_match_ref_pos - be.segment.range.start).abs() <= max_dist
    }
}

/// Process split alignments consistent with the breakpoint, so that we can extract the
/// corresponding read locations
///
/// Returns the corresponding read range for the split alignment if a split matching bp is found,
/// otherwise returns None
///
/// read_range is returned in a read coordinate system consistent with the alignment direction of
/// record
///
fn get_split_candidate_read_range(
    chrom_list: &ChromList,
    record: &bam::Record,
    bp: &Breakpoint,
) -> Option<IntRange> {
    // Max distance between split read position and the corresponding input bp breakends
    const MAX_BREAKEND_DIST: i64 = 200;

    const SA_AUX_TAG: &[u8] = b"SA";
    if !is_aux_tag_found(record, SA_AUX_TAG) {
        return None;
    }

    // Try to use either primary or supplementary reads here, but note that this won't work if the
    // supplementary reads have hard clipping, so we'll be sure to assert that case.
    //
    // This can be solved later to support minimap hard-clipped supplementary alignments, but we'll
    // have to do some read buffering or similar to guarantee the same result as a no-hard-clipping
    // bam.
    //
    if is_hard_clipped(&record.cigar()) {
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        eprintln!("ERROR: Unsupported hard-clipped split alignment in read: {qname}");
        std::process::exit(exitcode::DATAERR);
    }

    // Test whether breakend1 maps to the current primary bam record:
    if !is_breakend_close(
        MAX_BREAKEND_DIST,
        &bp.breakend1,
        record.tid() as usize,
        record.pos(),
        &record.cigar(),
    ) {
        return None;
    }

    // For the breakend2 match, get the map of all split reads segments for this read
    //
    // We will use this to extract the 1 segment adjacent to the primary alignment which we just
    // matched to breakend1
    //
    let seq_order_read_split_segments = get_seq_order_read_split_segments(chrom_list, record);

    // Determine if breakend2 has a higher sequencing order in the read compared to breakend1
    let is_breakend2_higher_seq_order_read_pos =
        (bp.breakend1.dir == BreakendDirection::LeftAnchor) ^ record.is_reverse();
    let breakend2_split_segment = get_breakend2_split_segment(
        &seq_order_read_split_segments,
        is_breakend2_higher_seq_order_read_pos,
    );

    // Test whether breakend2 maps to the expected SA alignment record:
    let is_breakend2_close = if let Some(breakend2_split_segment) = breakend2_split_segment {
        is_breakend_close(
            MAX_BREAKEND_DIST,
            bp.breakend2.as_ref().unwrap(),
            breakend2_split_segment.chrom_index,
            breakend2_split_segment.pos,
            &breakend2_split_segment.cigar,
        )
    } else {
        false
    };

    if !is_breakend2_close {
        return None;
    }

    // The breakend matches the read's split pattern so now report the corresponding read_range
    let breakend2_split_segment = breakend2_split_segment.unwrap();
    let seq_order_read_start = if is_breakend2_higher_seq_order_read_pos {
        breakend2_split_segment.seq_order_read_start - 1
    } else {
        breakend2_split_segment.seq_order_read_end
    };

    // Orient read coordinates to match the primary alignment's frame of reference:
    let read_start = if record.is_reverse() {
        let seq_order_read_end = seq_order_read_start + 1;
        record.seq_len() - (seq_order_read_end + 1)
    } else {
        seq_order_read_start
    } as i64;

    Some(IntRange::from_int(read_start))
}

#[derive(Clone)]
pub struct TrimmedReadInfo {
    /// QNAME of the bam record from which the read is trimmed
    ///
    /// Used to uniquely identify the read and prevent the same read from being counted
    /// multiple times, it is also used for many debugging outputs
    ///
    pub qname: String,

    /// True if the bam records used to create the trimmed read is supplementary
    ///
    /// This flag helps to distinguish cases where we would expect to see repeated qnames from
    /// potential read trimming or bam formatting bugs
    ///
    pub is_supplementary: bool,

    /// Range of trimmed read section in read coordinates
    pub read_range: IntRange,

    /// The trimmed read itself to be used for assembly
    pub trimmed_read: Vec<u8>,

    /// Segment of read before trimmed_read
    pub read_prefix: Vec<u8>,

    /// Segment of read after trimmed_read
    pub read_suffix: Vec<u8>,

    /// Record whether there is a large soft-clip on one end of the read, consistent with a large
    /// insertion
    pub large_insertion_clip_state: LargeInsertionSoftClipState,

    /// True if a nul-char has been added to trimmed_read
    ///
    /// This addition is used to interact with certain c libraries.
    ///
    /// Together with methods on this object which check the bool to operate on the trimmed
    /// read correctly, this will help make a safe transition if we switch libraries and
    /// stop adding the nul-term in future.
    ///
    pub is_trimmed_read_nul_terminated: bool,

    /// From the rq tag
    pub read_quality: f32,
}

impl TrimmedReadInfo {
    pub fn add_trimmed_read_nul_term(&mut self) {
        self.trimmed_read.push(0u8);
        self.is_trimmed_read_nul_terminated = true;
    }

    pub fn trimmed_read_len(&self) -> usize {
        self.trimmed_read.len()
            - if self.is_trimmed_read_nul_terminated {
                1
            } else {
                0
            }
    }

    pub fn prefix_len(&self) -> usize {
        self.read_prefix.len()
    }

    pub fn suffix_len(&self) -> usize {
        self.read_suffix.len()
    }
}

/// Return true if indel SV evidence is found in the targeted region of the reference
///
pub fn does_alignment_region_have_sv_indel_evidence(
    record: &bam::Record,
    ref_range: &IntRange,
    min_evidence_indel_size: u32,
) -> bool {
    use bam::record::Cigar::*;

    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    for c in record.cigar().iter() {
        match c {
            Del(len) => {
                if *len > min_evidence_indel_size
                    && ref_range
                        .intersect_range(&IntRange::from_pair(ref_pos, ref_pos + *len as i64))
                {
                    return true;
                }
            }
            Ins(len) => {
                if *len > min_evidence_indel_size
                    && ref_range.intersect_range(&IntRange::from_int(ref_pos))
                {
                    return true;
                }
            }
            _ => {}
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    false
}

/// Return true if softclip SV evidence is found in the targeted region of the reference
///
pub fn does_alignment_region_have_sv_softclip_evidence(
    record: &bam::Record,
    ref_range: &IntRange,
    min_evidence_indel_size: u32,
) -> bool {
    use bam::record::Cigar::*;

    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    for c in record.cigar().iter() {
        if let SoftClip(len) = c
            && *len > min_evidence_indel_size
            && ref_range.intersect_range(&IntRange::from_int(ref_pos))
        {
            return true;
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    false
}

#[allow(clippy::too_many_arguments)]
fn get_bam_record_trimmed_read(
    assembly_read_flank_size: usize,
    min_large_insertion_soft_clip_len: usize,
    min_evidence_indel_size: u32,
    chrom_list: &ChromList,
    target_segment: &GenomeSegment,
    is_large_insertion_inferred_from_clipping: bool,
    is_large_insertion_inferred_from_alignment: bool,
    insert_size_estimate: i64,
    record: &bam::record::Record,
) -> Option<TrimmedReadInfo> {
    let is_large_insertion_type =
        is_large_insertion_inferred_from_clipping || is_large_insertion_inferred_from_alignment;

    // Filter out reads lacking any SV support, by checking the region around target segment
    // range (in a reference coordinate window) for SV support evidence.
    //
    let mut check_segment = target_segment.clone();
    check_segment.expand_by(chrom_list, assembly_read_flank_size as i64);
    let indel_evidence_in_region = does_alignment_region_have_sv_indel_evidence(
        record,
        &check_segment.range,
        min_evidence_indel_size,
    );
    let softclip_evidence_in_region = does_alignment_region_have_sv_softclip_evidence(
        record,
        &check_segment.range,
        min_evidence_indel_size,
    );

    let filter_read_for_standard_sv_candidate =
        || -> bool { !(indel_evidence_in_region || softclip_evidence_in_region) };

    // Filter read out if no evidence is found in the target region, except for large insertion
    // types where large soft clipping is searched for even outside of the target region:
    //
    if !is_large_insertion_type && filter_read_for_standard_sv_candidate() {
        return None;
    }

    // Map from the breakend locations to locations in the hardclipped version of the read
    let mut read_range =
        translate_ref_range_to_hardclipped_read_range(record, &target_segment.range)?;

    use LargeInsertionSoftClipState::*;
    let large_insertion_clip_state = if is_large_insertion_type {
        test_read_for_large_insertion_soft_clip(record, min_large_insertion_soft_clip_len)
    } else {
        Null
    };

    // Modify read_range for large insertion soft-clipped reads
    if is_large_insertion_inferred_from_clipping {
        // For the clipping-inference case, we don't have an estimate of the insertion size, so for
        // any soft clipped reads we extract the entire clipped side of the read to the end.
        //
        // Any reads besides those with large insertion clipped-read support are filtered out
        //
        // TODO check that reads clip in the expected region
        match large_insertion_clip_state {
            Null => {
                return None;
            }
            Left => {
                read_range.start = 0;
            }
            Right => {
                read_range.end = record.seq_len() as i64;
            }
        }
    } else if is_large_insertion_inferred_from_alignment {
        // For the alignment-inference case, we have an estimate of the insert size, so limit the
        // size of the trimmed reads that we extract.
        //
        // Any reads besides those with large insertion clipped-read support are handled according
        // to the standard case logic.
        //
        match large_insertion_clip_state {
            Null => {
                if filter_read_for_standard_sv_candidate() {
                    return None;
                }
            }
            Left => {
                read_range.start = std::cmp::max(0, read_range.start - insert_size_estimate);
            }
            Right => {
                read_range.end = std::cmp::min(
                    record.seq_len() as i64,
                    read_range.end + insert_size_estimate,
                );
            }
        }
    }

    Some(get_trimmed_read_from_read_range(
        assembly_read_flank_size,
        record,
        &read_range,
        large_insertion_clip_state,
    ))
}

/// Dedup trimmed reads
///
/// A locus can have both the primary and supplementary alignments of the same read intersect the
/// trimmed read locus (mainly in large duplicated insertions). This method helps to remove the
/// duplicate cases so they aren't counted as double supporting evidence.
///
fn dedup_trimmed_reads(trimmed_reads: &mut Vec<TrimmedReadInfo>) {
    // Get duplicate qnames
    let mut duplicate_qnames = HashSet::new();
    {
        let mut qnames = HashSet::new();
        for trimmed_read in trimmed_reads.iter() {
            if !qnames.insert(trimmed_read.qname.clone()) {
                duplicate_qnames.insert(trimmed_read.qname.clone());
            }
        }
    }

    if duplicate_qnames.is_empty() {
        return;
    }

    // Extract the duplicated reads:
    let mut duplicate_reads = HashMap::new();
    for (tr_index, trimmed_read) in trimmed_reads.iter().enumerate() {
        if duplicate_qnames.contains(&trimmed_read.qname) {
            duplicate_reads
                .entry(&trimmed_read.qname)
                .or_insert_with(Vec::new)
                .push((tr_index, trimmed_read));
        }
    }

    // Choose the best read to keep from each duplicated read set and mark the others for removal
    let mut tr_remove_indices = Vec::new();
    for duplicate_read_set in duplicate_reads.values_mut() {
        // Choose the best representative trimmed read among the duplicated set:
        // (1) longest read
        // (2) for same read size choose the primary of supp reads
        duplicate_read_set.sort_by(|a, b| {
            b.1.trimmed_read
                .len()
                .cmp(&a.1.trimmed_read.len())
                .then(a.1.is_supplementary.cmp(&b.1.is_supplementary))
        });
        tr_remove_indices.extend(duplicate_read_set.iter().skip(1).map(|x| x.0));
    }
    tr_remove_indices.sort();

    remove_sorted_indices(trimmed_reads, tr_remove_indices);
}

/// Trim reads to include no more than 'assembly_read_flank_size' outside of the target_segment
/// range.
///
/// The method only gathers trimmed reads which contain evidence of an SV-like signature within
/// the flank size.
///
/// This method assumes simple indel like SVs where both breakend regions are nearby and colinear.
/// It supports a special mode targeting very large insertions.
///
/// Returns a vector of trimmed reads for assembly
///
#[allow(clippy::too_many_arguments)]
pub fn get_sv_candidate_region_trimmed_reads(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    target_segment: &GenomeSegment,
    bpc: &BreakpointCluster,
) -> Vec<TrimmedReadInfo> {
    // True for large insertions inferred from pairs of soft-clipped breakpoints, where no cigar
    // based large insertion candidate was found
    let is_large_insertion_inferred_from_clipping = bpc.large_insertion_info.is_some();

    let insert_size_estimate = if let InsertInfo::SizeRange(size) = &bpc.breakpoint.insert_info {
        size.end
    } else {
        0
    };

    // True for large insertions discovered via conventional cigar string alignments, where we seek
    // to combine cigar alignment based evidence with soft-clipping evidence
    let is_large_insertion_inferred_from_alignment = if !is_large_insertion_inferred_from_clipping {
        let min_insertion_soft_clip_recovery_size = 1000;
        insert_size_estimate >= min_insertion_soft_clip_recovery_size
    } else {
        false
    };

    // Get trimmed reads
    let mut trimmed_reads = Vec::new();
    {
        // For insertions with left-side soft-clipped read support, we need to go at least one base past the
        // breakpoint to detect them. The expansize size of 10 is designed to give a little more wiggle room
        // for read noise and other factors that could push the right-anchored breakpoint back.
        //
        let trimmed_region_search_expansion_size = 10;
        let target_segment = {
            let mut x = target_segment.clone();
            x.expand_by(chrom_list, trimmed_region_search_expansion_size);
            x
        };
        let assembly_read_flank_size = refine_settings
            .assembly_read_flank_size
            .saturating_sub(trimmed_region_search_expansion_size as usize);

        let chrom_sequence = {
            let chrom_label = &chrom_list.data[target_segment.chrom_index].label;
            reference.chroms.get(chrom_label).unwrap()
        };
        bam_fetch_segment(bam_reader, &target_segment);
        let mut record = bam::Record::new();
        while let Some(r) = bam_reader.read(&mut record) {
            unwrap!(r, "Failed to parse alignment record");

            if filter_out_alignment_record(&record) {
                continue;
            }

            if record.mapq() < refine_settings.min_sv_mapq as u8 {
                continue;
            }

            if get_gap_compressed_identity(&record, chrom_sequence)
                < refine_settings.min_gap_compressed_identity
            {
                continue;
            }

            if let Some(trimmed_read_info) = get_bam_record_trimmed_read(
                assembly_read_flank_size,
                refine_settings.min_large_insertion_soft_clip_len,
                refine_settings.min_evidence_indel_size,
                chrom_list,
                &target_segment,
                is_large_insertion_inferred_from_clipping,
                is_large_insertion_inferred_from_alignment,
                insert_size_estimate,
                &record,
            ) {
                trimmed_reads.push(trimmed_read_info);
            }
        }
    }

    dedup_trimmed_reads(&mut trimmed_reads);

    if is_large_insertion_inferred_from_clipping || is_large_insertion_inferred_from_alignment {
        trimmed_reads.sort_by(sort_large_insertion_reads);
    }

    trimmed_reads
}

/// Trim bp-supporting reads to include no more than 'refine_flank_size' outside of the
/// target_segment ranges.
///
/// This is a generalization of the `get_sv_candidate_region_trimmed_reads` to handle arbitrary
/// SVs expressed as split reads
///
pub fn get_sv_breakpoint_candidate_trimmed_reads(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    target_segments: &[GenomeSegment],
    bp: &Breakpoint,
) -> Vec<TrimmedReadInfo> {
    let debug = false;

    assert_eq!(target_segments.len(), 2);
    assert!(bp.is_standardized());
    assert!(bp.breakend2.is_some());

    // The breakpoint region is expanded by this amount when searching for supporting split reads
    //
    // Note that for right-anchored breakends at least one base of window is required to pick up
    // the split reads starting at the next position after the breakpoint.
    //
    // For both left and right anchored breakpoints, some extra wiggle room should help to pick up
    // additional supporting reads that may have some wobble in the breakpoint, although in theory
    // this should be picked up already in the diffention of the breakpoint region. This wiggle room
    // idea hasn't been thoroughly tested or optimized, but probably doesn't have a big up or down
    // side.
    //
    let breakpoint_candidate_trimmed_region_search_expansion_size = 10;
    let mut trimmed_reads_search_segment = target_segments[0].clone();
    trimmed_reads_search_segment.expand_by(
        chrom_list,
        breakpoint_candidate_trimmed_region_search_expansion_size,
    );

    if debug {
        eprintln!("Starting trimmed read search for bp {bp:?}");
        eprintln!("Expanded trimmed read search segment: {trimmed_reads_search_segment:?}",);
    }

    // for a first version just read everything from the breakend 1 region, since we aren't including
    // non-SA soft clipped read evidence at this point.

    // If the breakpoint does not have an indel-like pattern, then don't bother looking at the alignment
    // CIGAR strings
    let is_indel_bp = bp.is_indel();

    let mut trimmed_reads = Vec::new();

    let chrom_label = &chrom_list.data[trimmed_reads_search_segment.chrom_index].label;
    let chrom_sequence = reference.chroms.get(chrom_label).unwrap();

    bam_fetch_segment(bam_reader, &trimmed_reads_search_segment);
    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        if filter_out_alignment_record(&record) {
            continue;
        }

        if record.mapq() < refine_settings.min_sv_mapq as u8 {
            continue;
        }

        let gci = get_gap_compressed_identity(&record, chrom_sequence);
        if gci < refine_settings.min_gap_compressed_identity {
            if debug {
                eprintln!("Filtering read gci: {gci}");
            }
            continue;
        }

        // process cigar alignments consistent with the bp (this should only apply to large
        // deletions
        let read_range = if is_indel_bp {
            get_large_del_candidate_read_range(&record, bp)
        } else {
            None
        };

        let read_range = if read_range.is_none() {
            get_split_candidate_read_range(chrom_list, &record, bp)
        } else {
            read_range
        };

        if debug {
            eprintln!("read range {read_range:?}");
        }
        if let Some(read_range) = read_range {
            trimmed_reads.push(get_trimmed_read_from_read_range(
                refine_settings.assembly_read_flank_size,
                &record,
                &read_range,
                LargeInsertionSoftClipState::Null,
            ));
        }
    }

    dedup_trimmed_reads(&mut trimmed_reads);

    trimmed_reads
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trimmed_read_sort() {
        use LargeInsertionSoftClipState::*;

        let tr_right_len6 = TrimmedReadInfo {
            qname: "foo".to_string(),
            is_supplementary: false,
            read_range: IntRange::from_int(0),
            trimmed_read: b"ACTCTA".to_vec(),
            read_prefix: b"".to_vec(),
            read_suffix: b"GGT".to_vec(),
            large_insertion_clip_state: Right,
            is_trimmed_read_nul_terminated: false,
            read_quality: 0.99,
        };

        let tr_right_len3 = {
            let mut x = tr_right_len6.clone();
            x.trimmed_read = b"ACT".to_vec();
            x
        };

        let tr_left_len6 = {
            let mut x = tr_right_len6.clone();
            x.large_insertion_clip_state = Left;
            x
        };

        let tr_left_len3 = {
            let mut x = tr_left_len6.clone();
            x.trimmed_read = b"ACT".to_vec();
            x
        };

        let mut rval = vec![tr_left_len3, tr_left_len6, tr_right_len6, tr_right_len3];
        rval.sort_by(sort_large_insertion_reads);

        let expected_clip_state = vec![Right, Right, Left, Left];
        assert_eq!(
            rval.iter()
                .map(|x| x.large_insertion_clip_state)
                .collect::<Vec<_>>(),
            expected_clip_state
        );

        let expected_tread_len = vec![6, 3, 6, 3];
        assert_eq!(
            rval.iter()
                .map(|x| x.trimmed_read.len())
                .collect::<Vec<_>>(),
            expected_tread_len
        );
    }
}
