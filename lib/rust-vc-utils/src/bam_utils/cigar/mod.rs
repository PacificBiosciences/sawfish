//! BAM record cigar-processing utilities
//!

mod clip_alignment;
mod score_alignment;
mod shift_indels;

pub use clip_alignment::*;
pub use score_alignment::*;
pub use shift_indels::*;

use rust_htslib::bam::record::{self, Cigar};

/// Is the cigar element any clip type?
///
pub fn is_clip(c: &Cigar) -> bool {
    matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_))
}

/// Is the cigar element any of the alignment match types?
///
pub fn is_alignment_match(c: &Cigar) -> bool {
    matches!(c, Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_))
}

pub fn get_cigarseg_read_offset(c: &Cigar, ignore_hard_clip: bool) -> usize {
    use Cigar::*;
    match c {
        Ins(len) | SoftClip(len) | Diff(len) | Equal(len) | Match(len) => *len as usize,
        HardClip(len) => {
            if ignore_hard_clip {
                0
            } else {
                *len as usize
            }
        }
        _ => 0,
    }
}

pub fn get_cigarseg_ref_offset(c: &Cigar) -> i64 {
    use Cigar::*;
    match c {
        Del(len) | RefSkip(len) | Diff(len) | Equal(len) | Match(len) => *len as i64,
        _ => 0,
    }
}

/// A utility method to track read positions while iterating through a cigar string
///
pub fn update_read_pos(c: &Cigar, read_pos: &mut usize, ignore_hard_clip: bool) {
    *read_pos += get_cigarseg_read_offset(c, ignore_hard_clip);
}

/// A utility method to track ref positions while iterating through a cigar string
pub fn update_ref_pos(c: &Cigar, ref_pos: &mut i64) {
    *ref_pos += get_cigarseg_ref_offset(c);
}

/// A utility method to track ref and read positions while iterating through a cigar string
///
/// # Example
/// ```ignore
/// let mut ref_pos = 100;
/// let mut read_pos = 100;
/// for (index, c) in record.cigar().iter().enumerate() {
///     update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, false);
/// }
/// ```
pub fn update_ref_and_read_pos(
    c: &Cigar,
    ref_pos: &mut i64,
    read_pos: &mut usize,
    ignore_hard_clip: bool,
) {
    update_read_pos(c, read_pos, ignore_hard_clip);
    update_ref_pos(c, ref_pos);
}

/// Report the following positions in read coordinates:
/// 1. The first position after all left-side clipping (in read coordinates)
/// 2. The first position of all right-side clipping (in read coordaintes)
/// 3. The read length
///
pub fn get_read_clip_positions(cigar: &[Cigar], ignore_hard_clip: bool) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            SoftClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            HardClip(len) => {
                if !ignore_hard_clip {
                    if left_clip {
                        left_clip_size += *len as usize;
                    } else {
                        right_clip_size += *len as usize;
                    }
                }
            }
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// Report the following positions in read coordinates:
/// 1. The first position after all left-side hard clipping
/// 2. The first position of all right-side hard clipping
/// 3. The read length
///
pub fn get_read_hard_clip_positions(cigar: &[Cigar]) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            HardClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, false);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// Report the reference and read offset of the cigar alignment
///
pub fn get_cigar_ref_and_read_offset(cigar: &[Cigar], ignore_hard_clip: bool) -> (i64, usize) {
    let mut read_pos = 0;
    let mut ref_pos = 0;
    for c in cigar.iter() {
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    (ref_pos, read_pos)
}

/// Report the read offset of the cigar alignment
///
pub fn get_cigar_read_offset(cigar: &[Cigar], ignore_hard_clip: bool) -> usize {
    let mut read_pos = 0;
    for c in cigar.iter() {
        update_read_pos(c, &mut read_pos, ignore_hard_clip);
    }
    read_pos
}

/// Report the reference offset of the cigar alignment
///
pub fn get_cigar_ref_offset(cigar: &[Cigar]) -> i64 {
    let mut ref_pos = 0;
    for c in cigar.iter() {
        update_ref_pos(c, &mut ref_pos);
    }
    ref_pos
}

/// Return true if any part of the alignment is hard-clipped
///
pub fn is_hard_clipped(cigar: &[Cigar]) -> bool {
    cigar.iter().any(|x| matches!(x, Cigar::HardClip { .. }))
}

/// Convert CIGAR in string format into the format used in the this library
///
/// This is a convenience function so that client code doesn't need to add their own (potentially
/// conflicting) rust-htslib dependency
///
pub fn get_cigar_from_string(cigar_str: &str) -> Vec<Cigar> {
    record::CigarString::try_from(cigar_str.as_bytes())
        .unwrap()
        .into()
}

/// Compress CIGAR string down to canonical format:
///
/// 1. Convert any matching adjacent cigar elements into a single element
/// 2. Remove any zero-length elements
///
pub fn compress_cigar(cigar_in: &[Cigar]) -> Vec<Cigar> {
    let mut cigar_out = Vec::new();
    let mut last_elem = Cigar::Match(0);
    for new_elem in cigar_in.iter().filter(|x| !x.is_empty()) {
        if std::mem::discriminant(new_elem) == std::mem::discriminant(&last_elem) {
            use Cigar::*;
            if let Match(ref mut n) | Equal(ref mut n) | Diff(ref mut n) | Del(ref mut n)
            | Ins(ref mut n) | HardClip(ref mut n) | SoftClip(ref mut n)
            | RefSkip(ref mut n) = last_elem
            {
                *n += new_elem.len();
            }
        } else {
            if !last_elem.is_empty() {
                cigar_out.push(last_elem);
            }
            last_elem = *new_elem;
        }
    }
    if !last_elem.is_empty() {
        cigar_out.push(last_elem)
    }

    cigar_out
}

/// Any insertion segments found at either edge of the cigar string will be converted to soft-clip
///
/// Run compress_cigar after this to ensure there aren't unmerged adjacent soft-clip sequences
///
pub fn cigar_edge_insertion_to_softclip(cigar: &mut [Cigar]) {
    fn update_element(c: &mut Cigar) {
        if let Cigar::Ins(len) = c {
            *c = Cigar::SoftClip(*len);
        }
    }

    for c in cigar.iter_mut().take_while(|x| !is_alignment_match(x)) {
        update_element(c);
    }

    for c in cigar
        .iter_mut()
        .rev()
        .take_while(|x| !is_alignment_match(x))
    {
        update_element(c);
    }
}

/// Clean up all insertions and deletions at the edges of the alignment
///
/// This routine converts all insertions on the edge to soft-clip, and all deletions on the edge to zero-length
/// soft-clip. The alignment edge is defined here as the segment on each side of the alignment before the first match
/// occurs.
///
/// This routine may result in repeated cigar segments of the same type, and zero-length segments, so it is designed to
/// be followed by a call to compress_cigar.
///
/// Returns the shift in alignment start position, if any leading edge deletion length was removed
///
pub fn clean_up_cigar_edge_indels(cigar: &mut [Cigar]) -> usize {
    fn update_element(c: &mut Cigar) -> usize {
        let mut ret = 0;
        if let Cigar::Del(len) = c {
            ret = *len as usize;
            *c = Cigar::SoftClip(0);
        } else if let Cigar::Ins(len) = c {
            *c = Cigar::SoftClip(*len);
        }
        ret
    }

    let mut del_shift = 0;
    for c in cigar.iter_mut().take_while(|x| !is_alignment_match(x)) {
        del_shift += update_element(c);
    }

    for c in cigar
        .iter_mut()
        .rev()
        .take_while(|x| !is_alignment_match(x))
    {
        update_element(c);
    }

    del_shift
}

/// Return true if the CIGAR string contains any aligned (M/X/=) segments
///
pub fn has_aligned_segments(cigar: &[Cigar]) -> bool {
    cigar.iter().any(is_alignment_match)
}

/// Remove all leading soft/hard clipping from the cigar
pub fn strip_leading_clip(cigar: &mut Vec<Cigar>) {
    let mut non_clip_found = false;
    cigar.retain(|x| {
        if non_clip_found {
            true
        } else if is_clip(x) {
            false
        } else {
            non_clip_found = true;
            true
        }
    });
}

/// Remove all trailing soft/hard clipping from the cigar
pub fn strip_trailing_clip(cigar: &mut Vec<Cigar>) {
    let mut non_clip_found = false;
    cigar.retain(|x| {
        if non_clip_found {
            !is_clip(x)
        } else {
            if !is_clip(x) {
                non_clip_found = true;
            }
            true
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, Header, HeaderView, header};

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
    fn test_update_ref_and_read_pos() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5H5S5M5D5I5=5N5X5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGACTCTAGAGCT\tDDDDDEEEEEDDDDDEEEEEDDDDDEEEEE";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let ref_pos_expected = [100, 100, 105, 110, 110, 115, 120, 125, 125];
        let read_pos_expected = [5, 10, 15, 15, 20, 25, 25, 30, 35];

        let mut ref_pos = 100;
        let mut read_pos = 0;
        for (index, c) in record.cigar().iter().enumerate() {
            update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, false);
            assert_eq!(ref_pos, ref_pos_expected[index]);
            assert_eq!(read_pos, read_pos_expected[index]);
        }
        assert_eq!(record.cigar().iter().count(), ref_pos_expected.len());
    }

    #[test]
    fn test_update_read_pos() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5H5S5M5D5I5=5N5X5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGACTCTAGAGCT\tDDDDDEEEEEDDDDDEEEEEDDDDDEEEEE";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let read_pos_expected = [0, 5, 10, 10, 15, 20, 20, 25, 30];

        let mut read_pos = 0;
        for (index, c) in record.cigar().iter().enumerate() {
            update_read_pos(c, &mut read_pos, true);
            assert_eq!(read_pos, read_pos_expected[index]);
        }
    }

    #[test]
    fn test_get_read_clip_positions() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();

        let result = get_read_clip_positions(&cigar, true);
        assert_eq!(result, (10, 20, 30));

        let result = get_read_clip_positions(&cigar, false);
        assert_eq!(result, (20, 30, 50));
    }

    #[test]
    fn test_get_read_hard_clip_positions() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        let result = get_read_hard_clip_positions(&cigar);
        assert_eq!(result, (10, 40, 50));
    }

    #[test]
    fn test_is_hard_clipped() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        assert!(is_hard_clipped(&cigar));

        let cigar = record::CigarString::try_from("10S10M10S".as_bytes()).unwrap();
        assert!(!is_hard_clipped(&cigar));
    }

    #[test]
    fn test_compress_cigar() {
        use Cigar::*;
        let cigar = vec![
            HardClip(1),
            HardClip(1),
            SoftClip(1),
            SoftClip(1),
            Match(1),
            Match(1),
            Diff(1),
            Diff(0),
            Diff(1),
            Equal(1),
            Equal(1),
            Ins(1),
            Ins(1),
            Del(1),
            Del(1),
            Match(1),
            Match(1),
        ];
        assert_eq!(
            compress_cigar(&cigar),
            vec![
                HardClip(2),
                SoftClip(2),
                Match(2),
                Diff(2),
                Equal(2),
                Ins(2),
                Del(2),
                Match(2)
            ]
        );
    }

    #[test]
    fn test_cigar_edge_insertion_to_softclip() {
        use Cigar::*;
        let mut cigar = vec![
            HardClip(1),
            SoftClip(1),
            Ins(1),
            Match(1),
            Ins(1),
            Match(1),
            Ins(1),
            SoftClip(1),
        ];

        cigar_edge_insertion_to_softclip(&mut cigar);

        assert_eq!(
            cigar,
            vec![
                HardClip(1),
                SoftClip(1),
                SoftClip(1),
                Match(1),
                Ins(1),
                Match(1),
                SoftClip(1),
                SoftClip(1),
            ]
        );
    }

    #[test]
    fn test_clean_up_cigar_edge_indels() {
        use Cigar::*;
        let mut cigar = vec![
            HardClip(1),
            SoftClip(1),
            Ins(1),
            Del(2),
            Match(1),
            Ins(1),
            Del(1),
            Match(1),
            Ins(1),
            Del(1),
            SoftClip(1),
        ];

        let ret = clean_up_cigar_edge_indels(&mut cigar);

        assert_eq!(ret, 2);
        assert_eq!(
            cigar,
            vec![
                HardClip(1),
                SoftClip(1),
                SoftClip(1),
                SoftClip(0),
                Match(1),
                Ins(1),
                Del(1),
                Match(1),
                SoftClip(1),
                SoftClip(0),
                SoftClip(1),
            ]
        );
    }

    #[test]
    fn test_has_aligned_segments() {
        use Cigar::*;
        let cigar = vec![HardClip(2), SoftClip(2)];
        assert!(!has_aligned_segments(&cigar));

        let cigar = vec![Match(2)];
        assert!(has_aligned_segments(&cigar));
    }

    #[test]
    fn test_strip_leading_clip() {
        use Cigar::*;
        let mut cigar = vec![
            HardClip(2),
            SoftClip(2),
            Match(2),
            Ins(2),
            Match(2),
            SoftClip(2),
            HardClip(2),
        ];

        strip_leading_clip(&mut cigar);
        assert_eq!(
            cigar,
            vec![Match(2), Ins(2), Match(2), SoftClip(2), HardClip(2),]
        );
    }

    #[test]
    fn test_strip_trailing_clip() {
        use Cigar::*;
        let mut cigar = vec![
            HardClip(2),
            SoftClip(2),
            Match(2),
            Ins(2),
            Match(2),
            SoftClip(2),
            HardClip(2),
        ];

        strip_trailing_clip(&mut cigar);
        assert_eq!(
            cigar,
            vec![HardClip(2), SoftClip(2), Match(2), Ins(2), Match(2),]
        );
    }
}
