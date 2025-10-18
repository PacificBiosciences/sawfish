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
use rust_vc_utils::bam_utils::aux::is_aux_tag_found;
use rust_vc_utils::bam_utils::cigar::{
    get_gap_compressed_identity_from_alignment, get_read_clip_positions, update_ref_and_read_pos,
    update_ref_pos,
};
use rust_vc_utils::{GenomeSegment, IntRange};

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

    let ignore_hard_clip = true;

    assert!(cigar_start_index < cigar_end_index);
    assert!(cigar_end_index <= cigar.len());

    let sub_cigar = &cigar[cigar_start_index..cigar_end_index];

    if debug {
        eprintln!("full cigar: {}", CigarString(cigar.to_vec()));
        eprintln!("sub cigar: {}", CigarString(sub_cigar.to_vec()));
    }

    let mut read_offset = 0;
    for c in cigar[..cigar_start_index].iter() {
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_offset, ignore_hard_clip);
    }

    get_gap_compressed_identity_from_alignment(
        ref_pos,
        sub_cigar,
        &read_seq[read_offset..],
        chrom_seq,
        ignore_hard_clip,
    )
}

/// Get gap-compressed sequence identity from the bam record
///
pub fn get_gap_compressed_identity(record: &bam::Record, chrom_seq: &[u8]) -> f64 {
    let ignore_hard_clip = true;

    get_gap_compressed_identity_from_alignment(
        record.pos(),
        record.cigar().as_slice(),
        &get_simplified_dna_seq(record),
        chrom_seq,
        ignore_hard_clip,
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
        update_ref_and_read_pos(cigar_element, &mut ref_pos, &mut read_pos, true);
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
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, true);
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
    use rust_vc_utils::bam_utils::cigar::update_read_pos;
    let mut read_pos = 0;
    for c in cigar.iter() {
        update_read_pos(c, &mut read_pos, true);
    }
    read_pos
}

/// Note that enum order is non-arbitrary and used as part of read sort for assembly
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum LargeInsertionSoftClipState {
    /// Large soft-clipped region on the right side of the read alignment
    Right,

    /// Null state represents everything besides the right or left side clipping states
    Null,

    /// Large soft-clipped region on the left side of the read alignment
    Left,
}

/// Check bam record for soft-clipping of at least the specified length on one side of the read only
///
pub fn test_read_for_large_insertion_soft_clip(
    record: &bam::Record,
    min_soft_clip_len: usize,
) -> LargeInsertionSoftClipState {
    let (left_sclip_len, right_sclip_start, read_size) =
        get_read_clip_positions(&record.cigar(), true);
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

pub fn get_allowed_u8_lut(allowed: &[u8]) -> [bool; 256] {
    let mut x = [false; 256];
    for &c in allowed.iter() {
        x[c as usize] = true;
    }
    x
}

/// Get sequence from bam record, but convert any non-ACGT bases to N
///
pub fn get_simplified_dna_seq(record: &bam::Record) -> Vec<u8> {
    let allowed_lut = get_allowed_u8_lut(b"ACGTN");
    let seq = record.seq();
    (0..seq.len())
        .map(|i| {
            let b = seq[i];
            if allowed_lut[b as usize] { b } else { b'N' }
        })
        .collect()
}

/// Parse all ref-ranges from cigar alignment with gaps no larger than min_gap_size
///
pub fn get_alignment_ref_segments(
    mut ref_pos: i64,
    cigar: &[Cigar],
    min_gap_size: u32,
) -> Vec<IntRange> {
    use Cigar::*;

    let mut ref_segments = Vec::new();
    let mut begin_pos = ref_pos;

    for c in cigar.iter() {
        match c {
            Del(d) | RefSkip(d) => {
                if *d >= min_gap_size {
                    let ref_range = IntRange::from_pair(begin_pos, ref_pos);
                    if ref_range.size() > 0 {
                        ref_segments.push(ref_range)
                    }
                    begin_pos = ref_pos + *d as i64;
                }
            }
            _ => {}
        }
        update_ref_pos(c, &mut ref_pos);
    }

    let ref_range = IntRange::from_pair(begin_pos, ref_pos);
    if ref_range.size() > 0 {
        ref_segments.push(ref_range)
    }

    ref_segments
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Header, HeaderView, header};

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
                b"test_qname",
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
    fn test_get_alignment_ref_segments() {
        let cigar = vec![
            Cigar::Match(10),
            Cigar::Del(2),
            Cigar::Match(10),
            Cigar::Del(10),
            Cigar::Match(10),
        ];
        let ref_segs = get_alignment_ref_segments(100, &cigar, 10);
        assert_eq!(ref_segs.len(), 2);
        assert_eq!(ref_segs[0], IntRange::from_pair(100, 122));
        assert_eq!(ref_segs[1], IntRange::from_pair(132, 142));
    }
}
