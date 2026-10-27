//! This module covers many utilities related to bam-format processing, including cigar-based alignment utilities that may not otherwise
//! be related to a bam record
//!
//! Many of these utilities are candidates for future integration into the rust-vc-utils bam_utils module as they stabilize or as use cases
//! are found outside of sawfish
//!

use std::borrow::Cow;

use rust_htslib::bam::{
    self,
    record::{Cigar, CigarString, CigarStringView},
};
use rust_vc_utils::bam_utils::aux::is_aux_tag_found;
use rust_vc_utils::bam_utils::cigar::{
    get_gap_compressed_identity_from_alignment, get_read_clip_positions, update_ref_and_read_pos,
    update_ref_pos,
};
use rust_vc_utils::{GenomeSegment, IntRange};

/// Get gap-compressed sequence identity directly from the raw (packed BAM) CIGAR words, without
/// the read or reference sequence
///
/// Returns None if the CIGAR contains any alignment match (M) operation, because the identity
/// then depends on the sequences. Otherwise the result is exactly the value of
/// `get_gap_compressed_identity_from_alignment`: reference skips, clips and pads are ignored.
///
pub fn get_gap_compressed_identity_from_raw_cigar_no_align_match(raw_cigar: &[u32]) -> Option<f64> {
    let mut mismatch_events = 0u32;
    let mut match_bases = 0u32;

    for &c in raw_cigar.iter() {
        let len = c >> 4;
        match c & 0xf {
            // M
            0 => return None,
            // I, D
            1 | 2 => mismatch_events += 1,
            // =
            7 => match_bases += len,
            // X
            8 => mismatch_events += len,
            // N, S, H, P and anything else
            _ => {}
        }
    }

    // Same final step as get_gap_compressed_identity_from_alignment
    Some(if (match_bases + mismatch_events) == 0 {
        1.0f64
    } else {
        match_bases as f64 / (match_bases + mismatch_events) as f64
    })
}

/// Pack one CIGAR operation into BAM's raw 32-bit form (length in the high 28 bits, operation in
/// the low 4)
fn get_raw_cigar_word(cigar: &Cigar) -> u32 {
    use Cigar::*;
    let (len, op) = match cigar {
        Match(l) => (*l, 0),
        Ins(l) => (*l, 1),
        Del(l) => (*l, 2),
        RefSkip(l) => (*l, 3),
        SoftClip(l) => (*l, 4),
        HardClip(l) => (*l, 5),
        Pad(l) => (*l, 6),
        Equal(l) => (*l, 7),
        Diff(l) => (*l, 8),
    };
    (len << 4) | op
}

/// Inverse of `get_raw_cigar_word`
fn get_cigar_from_raw_word(raw: u32) -> Cigar {
    use Cigar::*;
    let len = raw >> 4;
    match raw & 0xf {
        0 => Match(len),
        1 => Ins(len),
        2 => Del(len),
        3 => RefSkip(len),
        4 => SoftClip(len),
        5 => HardClip(len),
        6 => Pad(len),
        7 => Equal(len),
        8 => Diff(len),
        op => panic!("Invalid CIGAR operation code {op}"),
    }
}

/// Get the record's CIGAR, borrowing the one cached by `Record::cache_cigar` when there is one
///
/// `Record::cigar` copies the cached CIGAR on every call, so a caller that reads the CIGAR more
/// than once per record should call `cache_cigar` once and then use this function.
///
/// The cache is not cleared when a new alignment is read into the record, so a cached CIGAR is
/// only used if it is identical to the record's current CIGAR (compared with the raw CIGAR, which
/// is cheap next to decoding it) and the record's position hasn't moved. Otherwise the CIGAR is
/// decoded from the record. This makes the function correct for any cache state, but
/// `Record::cigar` is not: code that can run on a record with a stale cache has to call this
/// function instead.
///
pub fn get_record_cigar(record: &bam::Record) -> Cow<'_, CigarStringView> {
    match record.cigar_cached() {
        Some(cigar)
            if cigar.pos() == record.pos()
                && cigar.len() == record.cigar_len()
                && cigar
                    .iter()
                    .zip(record.raw_cigar())
                    .all(|(c, &raw)| get_raw_cigar_word(c) == raw) =>
        {
            Cow::Borrowed(cigar)
        }
        // Record::cigar would return the stale cached CIGAR, so decode from the raw CIGAR
        _ => Cow::Owned(
            CigarString(
                record
                    .raw_cigar()
                    .iter()
                    .map(|&raw| get_cigar_from_raw_word(raw))
                    .collect(),
            )
            .into_view(record.pos()),
        ),
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
/// See genome_ref_utils::get_gap_compressed_identity for the version used with the reference
///
#[cfg(test)]
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
    for c in get_record_cigar(record).iter() {
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
        get_read_clip_positions(&get_record_cigar(record), true);
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

/// Get sequence from bam record, but convert any non-ACGT bases to N
///
pub fn get_simplified_dna_seq(record: &bam::Record) -> Vec<u8> {
    decode_simplified_dna_seq(record.seq().encoded, record.seq_len())
}

/// Decode a sequence from BAM's 4-bit encoding, converting any non-ACGT bases to N
///
/// `encoded` is the packed sequence of a bam record and `len` is its number of bases.
///
pub fn decode_simplified_dna_seq(encoded: &[u8], len: usize) -> Vec<u8> {
    // Decode both bases of each byte of BAM's 4-bit sequence encoding ("=ACMGRSVTWYHKDBN") at once
    const DECODE: &[u8; 16] = b"NACNGNNNTNNNNNNN";
    const PAIR_DECODE: [[u8; 2]; 256] = {
        let mut x = [[0u8; 2]; 256];
        let mut i = 0;
        while i < 256 {
            x[i] = [DECODE[i >> 4], DECODE[i & 0xf]];
            i += 1;
        }
        x
    };
    let mut seq = Vec::with_capacity(len + 1);
    for &b in encoded {
        seq.extend_from_slice(&PAIR_DECODE[b as usize]);
    }
    seq.truncate(len);
    seq
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
            Del(d) | RefSkip(d) if *d >= min_gap_size => {
                let ref_range = IntRange::from_pair(begin_pos, ref_pos);
                if ref_range.size() > 0 {
                    ref_segments.push(ref_range)
                }
                begin_pos = ref_pos + *d as i64;
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

    /// A cigar cached for one record must never be used for the next record read into the buffer
    #[test]
    fn test_get_record_cigar_ignores_stale_cache() {
        use bam::{Read, Writer};
        use rust_htslib::bam::record::CigarString;

        let cigars = [
            vec![Cigar::Match(10), Cigar::Del(2), Cigar::Match(8)],
            // Same number of operations as the first, to defeat a length-only check
            vec![Cigar::Match(7), Cigar::Ins(3), Cigar::Match(10)],
            vec![Cigar::Equal(20)],
            // Same CIGAR ops as the previous record, but at a different position
            vec![Cigar::Equal(20)],
        ];
        let tmp_path = std::env::temp_dir().join(format!(
            "sawfish_test_stale_cigar_{}.bam",
            std::process::id()
        ));
        {
            let mut header = bam::Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", 1000),
            );
            let mut writer = Writer::from_path(&tmp_path, &header, bam::Format::Bam).unwrap();
            for (i, cigar) in cigars.iter().enumerate() {
                let mut record = bam::Record::new();
                let len = cigar.iter().map(|c| c.len() as usize).sum::<usize>();
                let seq = vec![b'A'; if i == 1 { 20 } else { len }];
                record.set(
                    format!("r{i}").as_bytes(),
                    Some(&CigarString(cigar.clone())),
                    &seq,
                    &vec![30u8; seq.len()],
                );
                record.set_tid(0);
                record.set_pos(100 * (i as i64 + 1));
                writer.write(&record).unwrap();
            }
        }

        let mut reader = bam::Reader::from_path(&tmp_path).unwrap();
        let mut record = bam::Record::new();
        for (i, cigar) in cigars.iter().enumerate() {
            reader.read(&mut record).unwrap().unwrap();
            let observed = get_record_cigar(&record);
            assert_eq!(observed.as_slice(), cigar.as_slice(), "record {i}");
            assert_eq!(observed.pos(), record.pos());
            // Stale cache left over from this record
            record.cache_cigar();
            assert_eq!(get_record_cigar(&record).as_slice(), cigar.as_slice());
            assert!(matches!(get_record_cigar(&record), Cow::Borrowed(_)));
        }
        std::fs::remove_file(&tmp_path).unwrap();
    }

    /// The CIGAR-only identity must be bit-identical to the sequence-based one for =/X cigars
    #[test]
    fn test_gci_from_raw_cigar_matches_alignment_version() {
        use Cigar::*;
        let mut state = 12345u64;
        let mut next = |n: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % n
        };
        for _ in 0..2000 {
            let mut cigar = Vec::new();
            for _ in 0..(1 + next(12)) {
                let len = 1 + next(50) as u32;
                cigar.push(match next(8) {
                    0 => Ins(len),
                    1 => Del(len),
                    2 => RefSkip(len),
                    3 => SoftClip(len),
                    4 => HardClip(len),
                    5 => Pad(len),
                    6 => Diff(len),
                    _ => Equal(len),
                });
            }
            let raw = cigar.iter().map(get_raw_cigar_word).collect::<Vec<u32>>();
            // The sequence-based version never reads the sequences for =/X cigars
            let expected = get_gap_compressed_identity_from_alignment(0, &cigar, &[], &[], true);
            let observed = get_gap_compressed_identity_from_raw_cigar_no_align_match(&raw).unwrap();
            assert_eq!(expected.to_bits(), observed.to_bits(), "{cigar:?}");
        }
        assert!(
            get_gap_compressed_identity_from_raw_cigar_no_align_match(&[
                (10 << 4) | 7,
                (10 << 4) | 8
            ]) == Some(0.5)
        );
        assert!(
            get_gap_compressed_identity_from_raw_cigar_no_align_match(&[10 << 4]).is_none(),
            "an M operation needs the sequences"
        );
    }

    /// Compare with decoding each base, then converting bases other than ACGTN to N
    #[test]
    fn test_get_simplified_dna_seq() {
        let bases = b"=ACMGRSVTWYHKDBN";
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for len in 0..50 {
            let seq = (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    bases[(state % 16) as usize]
                })
                .collect::<Vec<_>>();
            let mut record = bam::Record::new();
            record.set(b"r", None, &seq, &vec![30; len]);
            let mut expected = record.seq().as_bytes();
            for b in expected.iter_mut() {
                if !matches!(*b, b'A' | b'C' | b'G' | b'T' | b'N') {
                    *b = b'N';
                }
            }
            assert_eq!(get_simplified_dna_seq(&record), expected);
        }
    }
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
