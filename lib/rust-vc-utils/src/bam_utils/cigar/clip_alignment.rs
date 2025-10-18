//! Alignment clipping facilities for read or ref based clip sizes
//!

use super::{compress_cigar, update_ref_and_read_pos, update_ref_pos};
use rust_htslib::bam::record::Cigar;

/// Modify a cigar string so that the read is soft-clipped on the left side to achieve at least the specified reference
/// start position shift
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out. If the required
/// clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the actual downstream shift of the reference start position
///
pub fn clip_alignment_ref_start(cigar_in: &[Cigar], min_left_ref_clip: i64) -> (Vec<Cigar>, i64) {
    use Cigar::*;

    let mut ref_pos = 0;

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
        update_ref_pos(c, &mut ref_pos);
    }
    (cigar_out, left_ref_clip_shift)
}

/// Modify a cigar string so that the read is soft-clipped on the left and right sides to create at least the specified
/// reference start and end position shifts
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out. If the required
/// clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
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
/// If the required clip length would end within an insertion, the entire insertion will be clipped out. If the required
/// clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the downstream shift of the reference start position implied
/// by the any new left-side clipping.
///
fn clip_alignment_read_start(cigar_in: &[Cigar], min_left_clip: usize) -> (Vec<Cigar>, i64) {
    use Cigar::*;

    let ignore_hard_clip = false;

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
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    (cigar_out, left_ref_clip_shift)
}

/// Modify a cigar string so that the read is soft-clipped on the left and right sides to at least the specified length
///
/// If the required clip length would end within an insertion, the entire insertion will be clipped out. If the required
/// clip length would create an unanchored deletion on the end of the read, the deletion will be removed.
///
/// Return a 2-tuple of (1) the modified cigar string (2) the downstream shift of the reference start position implied
/// by the any new left-side clipping.
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
