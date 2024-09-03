//! This module covers many utilities related to bam-format processing, including cigar-based alignment utilities that may not otherwise
//! be related to a bam record
//!
//! Many of these utilities are candidates for future integration into the rust-vc-utils bam_utils module as they stabilize or as use cases
//! are found outside of sawfish
//!

use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString},
};
use rust_vc_utils::bam_utils::cigar::{
    get_hard_clipped_read_clip_positions, update_ref_and_hard_clipped_read_pos,
};
use rust_vc_utils::{bam_utils::aux::is_aux_tag_found, cigar::compress_cigar};

use crate::genome_segment::GenomeSegment;
use crate::int_range::IntRange;

/// Get gap-compressed sequence identity from sequence alignment information
///
/// Metric is discussed here:
/// <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
///
/// This method is modified to allow for match (=), mismatch (X), and alignment match (M).
/// M cigars are compared to the reference to determine match/mismatch.
///
pub fn get_gap_compressed_identity_from_alignment(
    mut ref_pos: i64,
    read_seq: &[u8],
    cigar: &[Cigar],
    chrom_seq: &[u8],
) -> f64 {
    use Cigar::*;

    let mut mismatch_events = 0u32;
    let mut match_bases = 0u32;

    // set to usize since we need to index our sequences
    let mut read_pos = 0;

    for c in cigar.iter() {
        match c {
            Ins(_) | Del(_) => {
                mismatch_events += 1;
            }
            Diff(len) => {
                mismatch_events += len;
            }
            Equal(len) => {
                match_bases += len;
            }
            Match(len) => {
                let ref_pos = ref_pos as usize;
                for (offset, &ref_base) in chrom_seq[ref_pos..ref_pos + (*len as usize)]
                    .iter()
                    .enumerate()
                {
                    if ref_base == read_seq[read_pos + offset] {
                        match_bases += 1;
                    } else {
                        mismatch_events += 1;
                    }
                }
            }
            RefSkip(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }

    if (match_bases + mismatch_events) == 0 {
        1.0f64
    } else {
        match_bases as f64 / (match_bases + mismatch_events) as f64
    }
}

// Specify a subset of an alignment for gap_compressed identity metric
//
pub fn get_gap_compressed_identity_from_cigar_segment_range(
    mut ref_pos: i64,
    read_seq: &[u8],
    cigar: &[Cigar],
    chrom_seq: &[u8],
    cigar_start_index: usize,
    cigar_end_index: usize,
) -> f64 {
    let debug = false;

    assert!(cigar_start_index < cigar_end_index);
    assert!(cigar_end_index <= cigar.len());

    let sub_cigar = &cigar[cigar_start_index..cigar_end_index];

    if debug {
        eprintln!("full cigar: {}", CigarString(cigar.to_vec()));
        eprintln!("sub cigar: {}", CigarString(sub_cigar.to_vec()));
    }

    let mut read_offset = 0;
    for c in cigar[..cigar_start_index].iter() {
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_offset);
    }

    get_gap_compressed_identity_from_alignment(
        ref_pos,
        &read_seq[read_offset..],
        sub_cigar,
        chrom_seq,
    )
}

/// Get gap-compressed sequence identity from the bam record
///
pub fn get_gap_compressed_identity(record: &bam::Record, chrom_seq: &[u8]) -> f64 {
    get_gap_compressed_identity_from_alignment(
        record.pos(),
        &record.seq().as_bytes(),
        record.cigar().as_slice(),
        chrom_seq,
    )
}

/// Controls the closest match used by get_closest_read_to_ref_mapping
#[allow(dead_code)]
#[derive(Clone, Copy)]
pub enum TargetMatchType {
    Any,
    TargetOrLower,
    TargetOrHigher,
    Target,
}

#[derive(Debug, PartialEq)]
pub struct ReadToRefMapping {
    pub read_pos: usize,
    pub ref_pos: i64,
}

/// Given an alignment of a read sequence to a reference, and a target reference position, find the
/// read-to-ref mapped position that is closest to the target reference position
///
/// The solution can be additionally constrained by `target_match_type` to be greater or less than
/// the target, or to an exact match.
///
/// This method assumes the client has already checked chromosome match.
///
pub fn get_alignment_closest_to_target_ref_pos(
    alignment_start_ref_pos: i64,
    alignment_cigar: &[Cigar],
    target_ref_pos: i64,
    target_match_type: TargetMatchType,
) -> Option<ReadToRefMapping> {
    let mut closest = None;

    let mut ref_pos = alignment_start_ref_pos;
    let mut read_pos = 0usize;

    for cigar_element in alignment_cigar.iter() {
        use Cigar::*;
        match cigar_element {
            Diff(len) | Equal(len) | Match(len) => {
                let start_diff = target_ref_pos - ref_pos;
                let end_diff = start_diff - (*len as i64);
                if (start_diff >= 0) && (end_diff < 0) {
                    // Exact match:
                    closest = Some(ReadToRefMapping {
                        read_pos: read_pos + (start_diff as usize),
                        ref_pos: ref_pos + start_diff,
                    });
                    break;
                } else {
                    if let Some(closest) = &closest {
                        let low_diff = std::cmp::min(start_diff.abs(), end_diff.abs());
                        let cdiff = target_ref_pos - closest.ref_pos;
                        if cdiff <= low_diff {
                            break;
                        }
                    }
                    use TargetMatchType::*;
                    closest = if start_diff.abs() < end_diff.abs() {
                        // Past target position. Use start of alignment match region.
                        match target_match_type {
                            Any | TargetOrHigher => Some(ReadToRefMapping { read_pos, ref_pos }),
                            _ => None,
                        }
                    } else {
                        // Before target position. Use end of alignment match region.
                        match target_match_type {
                            Any | TargetOrLower => Some(ReadToRefMapping {
                                read_pos: read_pos + (*len as usize) - 1,
                                ref_pos: ref_pos + (*len as i64) - 1,
                            }),
                            _ => None,
                        }
                    };
                }
            }
            _ => {}
        }
        update_ref_and_hard_clipped_read_pos(cigar_element, &mut ref_pos, &mut read_pos);
    }

    closest
}

pub fn get_bam_alignment_closest_to_target_ref_pos(
    record: &bam::Record,
    target_ref_pos: i64,
    target_match_type: TargetMatchType,
) -> Option<ReadToRefMapping> {
    get_alignment_closest_to_target_ref_pos(
        record.pos(),
        &record.cigar(),
        target_ref_pos,
        target_match_type,
    )
}

/// Translate a range in ref coordinates into a range in read coordinates
///
/// All read coordinates reflect the position in the read after hard-clipping
///
/// Policy for read deletions spanning the ref start or end coordinates:
/// 1. The start ref pos maps to the right-most read pos not more than start ref pos.
/// 2. The end ref pos maps to the left-most read pos not less than end ref pos.
///
/// Policy for soft-clipping and alignment start/stop (see unit tests)
///
/// 1. If the read alignment is soft-clipped or starts between the ref range, the starting read
///    pos is the first local alignment position.
/// 2. If the alignment is soft-clipped or ends between the ref range, the ending read pos
///    is the last local alignment position
///
/// Returns None when the alignment is outside of ref_range
///
pub fn translate_ref_range_to_hardclipped_read_range(
    record: &bam::Record,
    ref_range: &IntRange,
) -> Option<IntRange> {
    use Cigar::*;

    let mut read_start = None;
    let mut read_end = None;

    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    // Read position after last alignment match to reference
    let mut read_pos_match_end = 0usize;
    for c in record.cigar().iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                if read_start.is_none() {
                    if ref_range.start >= ref_pos && ref_range.start < ref_pos + *len as i64 {
                        // If this case is true then we've exactly matched a read base to the start
                        // ref base at some point in the alignment:
                        //
                        let offset = ref_range.start - ref_pos;
                        read_start = Some(read_pos as i64 + offset);
                    } else if ref_range.start < ref_pos && ref_range.end > ref_pos {
                        // In this case, if read_start is still None, then the alignment has started
                        // after ref_range.start, so we pick it up were it starts.
                        read_start = Some(read_pos as i64);
                    }
                }
                if ref_range.end >= ref_pos && ref_range.end < ref_pos + *len as i64 {
                    let offset = ref_range.end - ref_pos;
                    read_end = Some(read_pos as i64 + offset);
                    break;
                }
                read_pos_match_end = read_pos + *len as usize;
            }
            Del(len) | RefSkip(len) => {
                if ref_range.start >= ref_pos && ref_range.start < ref_pos + *len as i64 {
                    read_start = Some(read_pos as i64 - 1);
                }
                if ref_range.end >= ref_pos && ref_range.end < ref_pos + *len as i64 {
                    read_end = Some(read_pos as i64);
                    break;
                }
            }
            _ => {}
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }

    if read_start.is_some() && read_end.is_none() && ref_range.end >= ref_pos {
        // This is the case where the alignment ended before reaching ref_range.end:
        read_end = Some(read_pos_match_end as i64);
    }

    if let (Some(read_start), Some(read_end)) = (read_start, read_end) {
        Some(IntRange::from_pair(read_start, read_end))
    } else {
        None
    }
}

pub fn bam_fetch_segment(bam_reader: &mut bam::IndexedReader, target_segment: &GenomeSegment) {
    bam_reader
        .fetch(bam::FetchDefinition::Region(
            target_segment.chrom_index as i32,
            target_segment.range.start,
            target_segment.range.end,
        ))
        .unwrap();
}

/// Report the read offset of the cigar alignment after hard-clipping is removed
///
#[allow(dead_code)]
pub fn get_cigar_hard_clipped_read_offset(cigar: &[Cigar]) -> usize {
    use rust_vc_utils::bam_utils::cigar::update_hard_clipped_read_pos;
    let mut read_pos = 0;
    for c in cigar.iter() {
        update_hard_clipped_read_pos(c, &mut read_pos);
    }
    read_pos
}

/// Get all edge soft and hard clipped sizes, ignoring any hard clipping present
///
/// Returns a 2-tuple of the left and right side combined soft-clip+insertion clipping
///
#[allow(dead_code)]
pub fn get_cigar_total_edge_clipping(cigar: &[Cigar]) -> (usize, usize) {
    use Cigar::*;

    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            SoftClip(len) | Ins(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            HardClip(_) => {}
            _ => {
                right_clip_size = 0;
                left_clip = false;
            }
        };
    }
    (left_clip_size, right_clip_size)
}

/// Note that enum order is non-arbitrary and used as part of read sort for assembly
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum LargeInsertionSoftClipState {
    Right,

    /// Null state represents everything besides the right or left side clipping states
    Null,

    Left,
}

/// Check bam record for soft-clipping of at least the specified length on one side of the read only
///
pub fn test_read_for_large_insertion_soft_clip(
    record: &bam::Record,
    min_soft_clip_len: usize,
) -> LargeInsertionSoftClipState {
    let (left_sclip_len, right_sclip_start, read_size) =
        get_hard_clipped_read_clip_positions(&record.cigar());
    let right_sclip_len = read_size.checked_sub(right_sclip_start).unwrap();

    let is_left_sclip = left_sclip_len >= min_soft_clip_len;
    let is_right_sclip = right_sclip_len >= min_soft_clip_len;

    if !(is_left_sclip || is_right_sclip) {
        return LargeInsertionSoftClipState::Null;
    }

    // Only one side should be clipped:
    if (is_left_sclip && right_sclip_len > 0) || (is_right_sclip && left_sclip_len > 0) {
        return LargeInsertionSoftClipState::Null;
    }

    if is_left_sclip {
        LargeInsertionSoftClipState::Left
    } else {
        LargeInsertionSoftClipState::Right
    }
}

pub fn is_split_read(record: &bam::Record) -> bool {
    const SA_AUX_TAG: &[u8] = b"SA";
    is_aux_tag_found(record, SA_AUX_TAG)
}

fn unexpected_aux_val_err(
    record: &bam::Record,
    aux_tag: &[u8],
    aux_val: bam::record::Aux<'_>,
) -> ! {
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    panic!(
        "Unexpected {} tag format in read {qname}: {:?}",
        std::str::from_utf8(aux_tag).unwrap(),
        aux_val,
    );
}

fn missing_aux_tag_err(record: &bam::Record, aux_tag: &[u8]) -> ! {
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    panic!(
        "Missing {} tag in read {qname}",
        std::str::from_utf8(aux_tag).unwrap(),
    );
}

/// Retrieve an aux tag float value from bam record, if the tag exists.
///
/// In this version the the tag itself is optional, but the function will still panic if the tag is present but has a non-float value
///
pub fn get_optional_float_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<f32> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::Float(val) => val,
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Retrieve an aux tag float value from bam record
///
/// Function will panic if the tag is missing or has a non-float value
///
#[allow(unused)]
pub fn get_float_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> f32 {
    get_optional_float_aux_tag(record, aux_tag)
        .unwrap_or_else(|| missing_aux_tag_err(record, aux_tag))
}

/// Modify a cigar string so that the read is soft-clipped on the left side to achieve at least the specified reference start position shift
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the actual downstream shift of the reference start position
///
pub fn clip_alignment_ref_start(cigar_in: &[Cigar], min_left_ref_clip: i64) -> (Vec<Cigar>, i64) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;

    let mut cigar_out = Vec::new();
    let mut left_ref_clip_shift = 0;
    for c in cigar_in.iter() {
        match c {
            Del(len) | RefSkip(len) => {
                if ref_pos <= min_left_ref_clip {
                    left_ref_clip_shift += *len as i64;
                } else {
                    cigar_out.push(*c);
                }
            }
            Ins(len) => {
                if ref_pos < min_left_ref_clip {
                    cigar_out.push(SoftClip(*len));
                } else {
                    cigar_out.push(*c);
                }
            }
            Match(len) | Diff(len) | Equal(len) => {
                if ref_pos < min_left_ref_clip {
                    let remaining_clip = min_left_ref_clip - left_ref_clip_shift;
                    let match_size = std::cmp::max(*len as i64 - remaining_clip, 0);
                    let clip_size = *len as i64 - match_size;
                    cigar_out.push(SoftClip(clip_size as u32));

                    if match_size > 0 {
                        let mut m = *c;
                        if let Match(ref mut len) | Diff(ref mut len) | Equal(ref mut len) = m {
                            *len = match_size as u32;
                        }
                        cigar_out.push(m);
                    }

                    left_ref_clip_shift += clip_size;
                } else {
                    cigar_out.push(*c);
                }
            }
            _ => {
                cigar_out.push(*c);
            }
        };
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (cigar_out, left_ref_clip_shift)
}

/// Modify a cigar string so that the read is soft-clipped on the left and right sides to create at least the specified reference start and end position shifts
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the actual downstream shift of the reference start position
///
pub fn clip_alignment_ref_edges(
    cigar_in: &[Cigar],
    min_left_ref_clip: i64,
    min_right_ref_clip: i64,
) -> (Vec<Cigar>, i64) {
    let rev_cigar = cigar_in.iter().rev().copied().collect::<Vec<_>>();
    let (mut right_clip_cigar, _) = clip_alignment_ref_start(&rev_cigar, min_right_ref_clip);

    right_clip_cigar.reverse();

    let (clip_cigar, ref_shift) = clip_alignment_ref_start(&right_clip_cigar, min_left_ref_clip);

    let cigar_out = compress_cigar(&clip_cigar);

    (cigar_out, ref_shift)
}

/// Modify a cigar string so that the read is soft-clipped on the left side to at least the specified length
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the downstream shift of the reference start position implied by the any new left-side clipping.
///
fn clip_alignment_read_start(cigar_in: &[Cigar], min_left_clip: usize) -> (Vec<Cigar>, i64) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;

    let mut cigar_out = Vec::new();
    let mut left_ref_clip_shift = 0;
    for c in cigar_in.iter() {
        match c {
            Del(len) | RefSkip(len) => {
                if read_pos <= min_left_clip {
                    left_ref_clip_shift += *len as i64;
                } else {
                    cigar_out.push(*c);
                }
            }
            Ins(len) => {
                if read_pos < min_left_clip {
                    cigar_out.push(SoftClip(*len));
                } else {
                    cigar_out.push(*c);
                }
            }
            Match(len) | Diff(len) | Equal(len) => {
                if read_pos < min_left_clip {
                    let remaining_clip = min_left_clip - read_pos;
                    let match_size = std::cmp::max(*len as i64 - remaining_clip as i64, 0);
                    let clip_size = *len as i64 - match_size;
                    cigar_out.push(SoftClip(clip_size as u32));

                    if match_size > 0 {
                        let mut m = *c;
                        if let Match(ref mut len) | Diff(ref mut len) | Equal(ref mut len) = m {
                            *len = match_size as u32;
                        }
                        cigar_out.push(m);
                    }

                    left_ref_clip_shift += clip_size;
                } else {
                    cigar_out.push(*c);
                }
            }
            _ => {
                cigar_out.push(*c);
            }
        };
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (cigar_out, left_ref_clip_shift)
}

/// Modify a cigar string so that the read is soft-clipped on the left and right sides to at least the specified length
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the downstream shift of the reference start position implied by the any new left-side clipping.
///
pub fn clip_alignment_read_edges(
    cigar_in: &[Cigar],
    min_left_clip: usize,
    min_right_clip: usize,
) -> (Vec<Cigar>, i64) {
    let rev_cigar = cigar_in.iter().rev().copied().collect::<Vec<_>>();
    let (mut right_clip_cigar, _) = clip_alignment_read_start(&rev_cigar, min_right_clip);

    right_clip_cigar.reverse();

    let (clip_cigar, ref_shift) = clip_alignment_read_start(&right_clip_cigar, min_left_clip);

    let cigar_out = compress_cigar(&clip_cigar);

    (cigar_out, ref_shift)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{header, Header, HeaderView};

    #[test]
    fn test_gap_compressed_identity() {
        use rust_htslib::bam::record::{Cigar, CigarString};
        use rust_vc_utils::genome_ref;

        // "contig_gap_test" => "AACGTGTTAACCCCT"
        let reference_filename = "test_data/test_reference.fa";
        let reference_genome = genome_ref::get_genome_ref_from_fasta(reference_filename);
        let reference_contig = reference_genome.chroms.get("contig_gap_test").unwrap();

        // (read_sequence, Cigar String, read POS, expected gap identity)
        let test_cases = [
            // exact reference
            (
                "AACGTGTTAACCCCT",
                CigarString(vec![Cigar::Equal(15)]),
                0,
                1.0,
            ),
            // exact reference, but using Match
            (
                "AACGTGTTAACCCCT",
                CigarString(vec![Cigar::Match(15)]),
                0,
                1.0,
            ),
            // exact reference, but using Match, and with some clipping at both ends
            (
                "AAAAACGTGTTAACCCCT",
                CigarString(vec![
                    Cigar::SoftClip(3),
                    Cigar::Match(15),
                    Cigar::HardClip(100),
                ]),
                0,
                1.0,
            ),
            // exact reference starting from the 5th position
            ("GTTAACCCCT", CigarString(vec![Cigar::Match(10)]), 5, 1.0),
            // two base changes reference starting from the 5th position
            ("GTCAACCGCT", CigarString(vec![Cigar::Match(10)]), 5, 0.8),
            // all the things at once, we have 2 modified bases, one per Match; an Insertion; both types of clipping; a position offset; and a deletion
            // AACGTGTTA--ACCCCT
            // ---GTATTANNACCA-T
            // ---SSX===II===XD= => S2M4I2M4D1M1+H100
            // total matching bases = 7
            // total mismatches = 2
            // total indels = 2
            // 7 / (7+4)
            (
                "GTATTANNACCAT",
                CigarString(vec![
                    Cigar::SoftClip(2),
                    Cigar::Match(4),
                    Cigar::Ins(2),
                    Cigar::Match(4),
                    Cigar::Del(1),
                    Cigar::Match(1),
                    Cigar::HardClip(100),
                ]),
                5,
                7.0 / 11.0,
            ),
        ];

        for (test_index, (test_sequence, test_cigar, test_pos, expected_score)) in
            test_cases.iter().enumerate()
        {
            let mut test_record = bam::Record::new();
            test_record.set(
                "test_qname".as_bytes(),
                Some(test_cigar),
                test_sequence.as_bytes(),
                &vec![20; test_sequence.len()],
            );
            test_record.set_pos(*test_pos);

            assert_eq!(
                *expected_score,
                get_gap_compressed_identity(&test_record, reference_contig),
                "Failed test case at index {test_index}, {test_sequence}"
            )
        }

        // assert_eq!(variant.convert_index(2), u8::MAX);
    }

    #[test]
    fn test_get_gap_compressed_identity_from_alignment() {
        let ref_pos = 2;
        let chrom_seq = b"ACGTACGTACGT";
        let read_seq = b"GTAATCTTAC";
        let cigar = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let gci = get_gap_compressed_identity_from_alignment(ref_pos, read_seq, &cigar, chrom_seq);
        approx::assert_ulps_eq!(gci, 6.0 / 9.0, max_ulps = 4);
    }

    #[test]
    fn test_get_gap_compressed_identity_from_cigar_segment_range() {
        let ref_pos = 2;
        let chrom_seq = b"ACGTACGTACGT";
        let read_seq = b"GTAATCTTAC";
        let cigar = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let gci = get_gap_compressed_identity_from_cigar_segment_range(
            ref_pos, read_seq, &cigar, chrom_seq, 2, 3,
        );
        approx::assert_ulps_eq!(gci, 3.0 / 4.0, max_ulps = 4);
    }

    fn get_test_header() -> HeaderView {
        let mut _header = Header::new();
        _header.push_record(
            header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 10000000),
        );
        HeaderView::from_header(&_header)
    }

    #[test]
    fn test_get_bam_alignment_closest_to_target_ref_pos() {
        use Cigar::*;
        use TargetMatchType::*;

        let target_ref_pos = 6;

        // Test exact match case:
        //
        let cigar = vec![Match(10)];
        let closest = get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, Any);
        assert_eq!(
            closest,
            Some(ReadToRefMapping {
                read_pos: 5,
                ref_pos: 6,
            })
        );

        // Test left-side clip case:
        //
        let expect = Some(ReadToRefMapping {
            read_pos: 1,
            ref_pos: 2,
        });
        let cigar = vec![Match(2), SoftClip(8)];
        let closest = get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, Any);
        assert_eq!(closest, expect);
        let closest =
            get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, TargetOrLower);
        assert_eq!(closest, expect);
        let closest =
            get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, TargetOrHigher);
        assert_eq!(closest, None);
        let closest = get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, Target);
        assert_eq!(closest, None);

        // Test spanning deletion case:
        //
        let expect = Some(ReadToRefMapping {
            read_pos: 2,
            ref_pos: 8,
        });
        let cigar = vec![Match(2), Del(5), Match(8)];
        let closest = get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, Any);
        assert_eq!(closest, expect);
        let closest =
            get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, TargetOrHigher);
        assert_eq!(closest, expect);
        let closest =
            get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, TargetOrLower);
        assert_eq!(closest, None);
        let closest = get_alignment_closest_to_target_ref_pos(1, &cigar, target_ref_pos, Target);
        assert_eq!(closest, None);
    }

    /// Test some basic ref to read range translation ops
    #[test]
    fn test_translate_ref_range_to_hardclipped_read_range() {
        let header = get_test_header();

        let ref_range = IntRange::from_pair(6, 8); // 1-indexed pos: [7-8]

        // Ungapped read test
        //
        // In this case read range is just an offset version of the ref range.
        //
        let sam_line = b"qname\t0\tchr1\t2\t60\t5H10M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(5, 7)));

        // Read insertion test
        //
        // In this case the read insertion is inside of target ref range, so the read range is
        // longer than the ref range.
        //
        let sam_line = b"qname\t0\tchr1\t2\t60\t6M1I3M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(5, 8)));

        // Read deletion test
        //
        // In this case the read deletion is inside of the target ref range, so the read range is
        // shorter than the ref range.
        //
        let sam_line = b"qname\t0\tchr1\t2\t60\t3M4D7M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(2, 3)));
    }

    /// Test more difficult ref to read range translation cases encountered from real data
    #[test]
    fn test_translate_ref_range_to_hardclipped_read_range_edge_cases() {
        let header = get_test_header();

        let ref_range = IntRange::from_pair(4, 8); // 1-indexed pos: [5-8]

        // Test case where read doesn't overlap with ref range at all:
        let sam_line = b"qname\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, None);

        // Test case where read alignment starts partway into the ref range.
        //
        // Here the read maps to 1-indexed ref positions [7-8], which translates to 1-indexed read
        // range [6-7] because the first 5 read bases are soft-clipped.
        let sam_line = b"qname\t0\tchr1\t7\t60\t5S5M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(5, 7)));

        // Test case where read alignment ends partway into the ref range.
        //
        // Here the read maps to 1-indexed ref positions [5-6] before starting a right-side
        // soft-clip
        let sam_line = b"qname\t0\tchr1\t1\t60\t6M4S\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(4, 6)));

        // Test case where read alignment has a deletion spanning the ref range.
        //
        // The resulting read map is the two read bases flanking the deletion, and therefor spanning
        // the target ref range.
        //
        let sam_line = b"qname\t0\tchr1\t1\t60\t4M10D6M\t*\t0\t0\tACGCCGTATC\tDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let read_range = translate_ref_range_to_hardclipped_read_range(&rec, &ref_range);
        assert_eq!(read_range, Some(IntRange::from_pair(3, 4)));
    }

    #[test]
    fn test_clip_alignment_ref_edges() {
        {
            let cigar = vec![Cigar::SoftClip(3), Cigar::Match(15)];
            let (cigar_out, ref_shift) = clip_alignment_ref_edges(&cigar, 5, 2);
            assert_eq!(
                cigar_out,
                vec![Cigar::SoftClip(8), Cigar::Match(8), Cigar::SoftClip(2)]
            );
            assert_eq!(ref_shift, 5);
        }

        {
            let cigar = vec![
                Cigar::SoftClip(3),
                Cigar::Match(2),
                Cigar::Del(3),
                Cigar::Match(13),
            ];
            let (cigar_out, ref_shift) = clip_alignment_ref_edges(&cigar, 5, 2);
            assert_eq!(
                cigar_out,
                vec![Cigar::SoftClip(5), Cigar::Match(11), Cigar::SoftClip(2)]
            );
            assert_eq!(ref_shift, 5);
        }
    }

    #[test]
    fn test_clip_alignment_read_edges() {
        {
            let cigar = vec![Cigar::SoftClip(3), Cigar::Match(15)];
            let (cigar_out, ref_shift) = clip_alignment_read_edges(&cigar, 5, 2);
            assert_eq!(
                cigar_out,
                vec![Cigar::SoftClip(5), Cigar::Match(11), Cigar::SoftClip(2)]
            );
            assert_eq!(ref_shift, 2);
        }

        {
            let cigar = vec![
                Cigar::SoftClip(3),
                Cigar::Match(2),
                Cigar::Del(3),
                Cigar::Match(13),
            ];
            let (cigar_out, ref_shift) = clip_alignment_read_edges(&cigar, 5, 2);
            assert_eq!(
                cigar_out,
                vec![Cigar::SoftClip(5), Cigar::Match(11), Cigar::SoftClip(2)]
            );
            assert_eq!(ref_shift, 5);
        }

        {
            let cigar = vec![Cigar::SoftClip(3), Cigar::Ins(3), Cigar::Match(12)];
            let (cigar_out, ref_shift) = clip_alignment_read_edges(&cigar, 5, 2);
            assert_eq!(
                cigar_out,
                vec![Cigar::SoftClip(6), Cigar::Match(10), Cigar::SoftClip(2)]
            );
            assert_eq!(ref_shift, 0);
        }
    }
}
