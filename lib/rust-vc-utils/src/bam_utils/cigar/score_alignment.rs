use rust_htslib::bam::record::Cigar;
use simple_error::{SimpleResult, bail};

use super::update_ref_and_read_pos;

/// Get edit distance from sequence alignment
///
/// This method will work with cigar strings using either the =/X or M match state.
///
pub fn get_edit_distance(
    mut ref_pos: i64,
    cigar: &[Cigar],
    read_seq: &[u8],
    ref_seq: &[u8],
) -> u32 {
    use Cigar::*;

    let ignore_hard_clip = false;

    let mut dist = 0u32;
    let mut read_pos = 0;
    for c in cigar.iter() {
        match c {
            Ins(len) | Del(len) | RefSkip(len) | Diff(len) => {
                dist += len;
            }
            Match(len) => {
                let ref_pos = ref_pos as usize;
                dist += ref_seq[ref_pos..ref_pos + (*len as usize)]
                    .iter()
                    .zip(read_seq[read_pos..].iter())
                    .filter(|(a, b)| a != b)
                    .count() as u32;
            }
            Equal(_) | SoftClip(_) | HardClip(_) | Pad(_) => (),
        }
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }

    dist
}

/// Get edit distance from sequence alignment
///
/// This method will only work with cigar strings using the =/X match state
///
pub fn get_edit_distance_no_align_match(cigar: &[Cigar]) -> SimpleResult<u32> {
    use Cigar::*;

    let mut dist = 0u32;
    for c in cigar.iter() {
        match c {
            Ins(len) | Del(len) | RefSkip(len) | Diff(len) => {
                dist += len;
            }
            Match(_) => {
                bail!(
                    "Method assumes alignment CIGAR strings use seq match/mismatch (=/X) instead of alignment match (M)"
                );
            }
            Equal(_) | SoftClip(_) | HardClip(_) | Pad(_) => (),
        }
    }

    Ok(dist)
}

fn get_final_gci(match_bases: u32, mismatch_events: u32) -> f64 {
    if (match_bases + mismatch_events) == 0 {
        1.0f64
    } else {
        match_bases as f64 / (match_bases + mismatch_events) as f64
    }
}

/// Get gap-compressed sequence identity from sequence alignment information
///
/// Metric is discussed here:
/// <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
///
/// This method will work with cigar strings using either the =/X or M match state.
///
pub fn get_gap_compressed_identity_from_alignment(
    mut ref_pos: i64,
    cigar: &[Cigar],
    read_seq: &[u8],
    ref_seq: &[u8],
    ignore_hard_clip: bool,
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
                for (offset, &ref_base) in ref_seq[ref_pos..ref_pos + (*len as usize)]
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
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }

    get_final_gci(match_bases, mismatch_events)
}

/// Get gap-compressed sequence identity from sequence alignment information
///
/// This is a deprecated/legacy function name which keeps the previous ignore_hard_clip=false default
///
pub fn get_gap_compressed_identity(ref_pos: i64, cigar: &[Cigar], read_seq: &[u8], ref_seq: &[u8]) {
    get_gap_compressed_identity_from_alignment(ref_pos, cigar, read_seq, ref_seq, false);
}

/// Get gap-compressed sequence identity from sequence alignment information
///
/// Metric is discussed here:
/// <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
///
/// This method will only work with cigar strings using the =/X match state.
///
pub fn get_gap_compressed_identity_no_align_match(cigar: &[Cigar]) -> SimpleResult<f64> {
    use Cigar::*;

    let mut mismatch_events = 0u32;
    let mut match_bases = 0u32;

    for c in cigar.iter() {
        match c {
            Ins(_) | Del(_) | RefSkip(_) => {
                mismatch_events += 1;
            }
            Diff(len) => {
                mismatch_events += len;
            }
            Equal(len) => {
                match_bases += len;
            }
            Match(_) => {
                bail!(
                    "Method assumes alignment CIGAR strings use seq match/mismatch (=/X) instead of alignment match (M)"
                );
            }
            SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }

    Ok(get_final_gci(match_bases, mismatch_events))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_edit_distance() {
        let ref_pos = 2;
        let ref_seq = b"ACGTACGTACGT";
        let read_seq = b"GTAATCTTAC";
        let cigar = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let dist = get_edit_distance(ref_pos, &cigar, read_seq, ref_seq);
        assert_eq!(dist, 4);
    }

    #[test]
    fn test_get_gap_compressed_identity_from_alignment() {
        let ref_pos = 2;
        let ref_seq = b"ACGTACGTACGT";
        let read_seq = b"GTAATCTTAC";
        let cigar = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let gci =
            get_gap_compressed_identity_from_alignment(ref_pos, &cigar, read_seq, ref_seq, false);
        approx::assert_ulps_eq!(gci, 6.0 / 9.0, max_ulps = 4);
    }
}
