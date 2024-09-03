use std::fmt;

use rust_htslib::bam::record::{Cigar, CigarString};
use rust_vc_utils::bam_utils::cigar::{get_cigar_hard_clipped_read_offset, get_cigar_ref_offset};

use crate::bam_utils;

#[derive(Clone, Default)]
pub struct SimpleAlignment {
    pub ref_offset: i64,
    pub cigar: Vec<Cigar>,
}

impl SimpleAlignment {
    pub fn new() -> Self {
        Self {
            ref_offset: 0,
            cigar: Vec::new(),
        }
    }

    pub fn get_reverse(&self, ref_size: i64) -> Self {
        let ref_end_offset = self.ref_offset + get_cigar_ref_offset(&self.cigar);
        assert!(ref_end_offset <= ref_size);
        Self {
            ref_offset: (ref_size - ref_end_offset),
            cigar: self.cigar.iter().rev().copied().collect(),
        }
    }
}

impl fmt::Debug for SimpleAlignment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let cigar_ref_len = get_cigar_ref_offset(&self.cigar);
        let cigar_read_len = get_cigar_hard_clipped_read_offset(&self.cigar);
        write!(
            f,
            "offset: {} cigar_ref_len: {cigar_ref_len} cigar_read_len: {cigar_read_len} cigar: {}",
            self.ref_offset,
            CigarString(self.cigar.clone())
        )
    }
}

/// Modify an alignment so that it is soft-clipped on the left and right sides to create at least the specified reference start and end position shifts
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a new alignment with cigar and ref offset modified to add the required clipping length
///
pub fn clip_alignment_ref_edges(
    alignment: &SimpleAlignment,
    ref_size: usize,
    min_left_ref_clip: i64,
    min_right_ref_clip: i64,
) -> SimpleAlignment {
    let ref_end_offset = alignment.ref_offset + get_cigar_ref_offset(&alignment.cigar);
    let right_ref_offset = ref_size as i64 - ref_end_offset;

    let min_left_cigar_ref_clip = min_left_ref_clip - alignment.ref_offset;
    let min_right_cigar_ref_clip = min_right_ref_clip - right_ref_offset;

    let (clipped_cigar, clipped_ref_shift) = bam_utils::clip_alignment_ref_edges(
        &alignment.cigar,
        min_left_cigar_ref_clip,
        min_right_cigar_ref_clip,
    );

    SimpleAlignment {
        ref_offset: alignment.ref_offset + clipped_ref_shift,
        cigar: clipped_cigar,
    }
}

/// Modify an alignment so that it is soft-clipped on the left and right sides to at least the specified length
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out.
/// If the required clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a new alignment with cigar and ref offset modified to add the required clipping length
///
#[allow(unused)]
pub fn clip_alignment_read_edges(
    alignment: &SimpleAlignment,
    min_left_clip: usize,
    min_right_clip: usize,
) -> SimpleAlignment {
    let (clipped_cigar, clipped_ref_shift) =
        bam_utils::clip_alignment_read_edges(&alignment.cigar, min_left_clip, min_right_clip);

    SimpleAlignment {
        ref_offset: alignment.ref_offset + clipped_ref_shift,
        cigar: clipped_cigar,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_reverse() {
        let al = SimpleAlignment {
            ref_offset: 1,
            cigar: vec![Cigar::Match(1)],
        };
        let alr = al.get_reverse(4);

        assert_eq!(alr.ref_offset, 2);
        assert_eq!(alr.cigar, vec![Cigar::Match(1)]);

        let al = SimpleAlignment {
            ref_offset: 1,
            cigar: vec![Cigar::Match(2)],
        };
        let alr = al.get_reverse(4);

        assert_eq!(alr.ref_offset, 1);
        assert_eq!(alr.cigar, vec![Cigar::Match(2)]);

        let al = SimpleAlignment {
            ref_offset: 2,
            cigar: vec![Cigar::Match(2), Cigar::Del(2), Cigar::Diff(1)],
        };
        let alr = al.get_reverse(20);

        assert_eq!(alr.ref_offset, 13);
        assert_eq!(
            alr.cigar,
            vec![Cigar::Diff(1), Cigar::Del(2), Cigar::Match(2)]
        );
    }
}
