use crate::int_range::IntRange;

/// Find the range over which an indel can vary without changing the alignment's edit distance
///
/// # Arguments
///
/// * `ref_range` - ref span of the indel in its current represented location
/// * `read_range` - read span of the indel in its current represented location
///
/// Range coordinates are zero-indexed and start at the first position affected by the indel. For instance:
/// - The deletion 2M1D2M would have refRange(2,3), readRange(2,2)
/// - The insertion 2M1I2M would have refRange(2,2), readRange(2,3)
///
/// Note that this method will push indels all the way to the edge of a read if homology supports it, which may
/// translate into an alignment like 1D4M. Client code can check and apply alternate logic for such a case if required.
///
/// Additionally note that this method does not account for adjacent indels as part of a more complex alignment, the
/// target indel is assessed assuming perfect mapping along the left and right flanks.
///
/// Return a 2-tuple of
/// 1. The range of indel offsets relative to the current pos that maintain the edit distance,
/// 2. The corresponding homology sequence
///
pub fn get_indel_breakend_homology_info(
    ref_seq: &[u8],
    ref_range: &IntRange,
    read_seq: &[u8],
    read_range: &IntRange,
) -> (IntRange, Vec<u8>) {
    let mut hom_seq = Vec::new();

    // Test how far the indel can be translated to the left of its current position:
    let max_left_offset = std::cmp::min(ref_range.start, read_range.start);
    let mut left_offset = 0;
    loop {
        if left_offset >= max_left_offset {
            break;
        }
        let ref_base = ref_seq[(ref_range.end - left_offset - 1) as usize];
        let read_base = read_seq[(read_range.end - left_offset - 1) as usize];

        if ref_base != read_base {
            break;
        }
        hom_seq.push(ref_base);
        left_offset += 1;
    }

    hom_seq.reverse();

    // Test how far the indel can be translated to the right of its current position:
    let max_right_offset = std::cmp::min(
        ref_seq.len() as i64 - ref_range.end,
        read_seq.len() as i64 - read_range.end,
    );
    let mut right_offset = 0;
    loop {
        if right_offset >= max_right_offset {
            break;
        }
        let ref_base = ref_seq[(ref_range.start + right_offset) as usize];
        let read_base = read_seq[(read_range.start + right_offset) as usize];
        if ref_base != read_base {
            break;
        }
        hom_seq.push(ref_base);
        right_offset += 1;
    }

    let hom_range = IntRange::from_pair(-left_offset, right_offset);

    (hom_range, hom_seq)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_indel_homology_range() {
        let seq1 = b"ABCDDABC";
        let seq2 = b"ABCDDDABC";

        {
            // left shifted case:
            let seq1_range = IntRange::from_pair(3, 3);
            let seq2_range = IntRange::from_pair(3, 4);

            // order reflects a deletion
            let (range, seq) =
                get_indel_breakend_homology_info(seq2, &seq2_range, seq1, &seq1_range);
            assert_eq!(range, IntRange::from_pair(0, 2));
            assert_eq!(seq, b"DD");

            // order reflects an insertion
            let (range, seq) =
                get_indel_breakend_homology_info(seq1, &seq1_range, seq2, &seq2_range);
            assert_eq!(range, IntRange::from_pair(0, 2));
            assert_eq!(seq, b"DD");
        }

        {
            // right shifted case:
            let seq1_range = IntRange::from_pair(5, 5);
            let seq2_range = IntRange::from_pair(5, 6);

            // order reflects a deletion
            let (range, seq) =
                get_indel_breakend_homology_info(seq2, &seq2_range, seq1, &seq1_range);
            assert_eq!(range, IntRange::from_pair(-2, 0));
            assert_eq!(seq, b"DD");

            // order reflects an insertion
            let (range, seq) =
                get_indel_breakend_homology_info(seq1, &seq1_range, seq2, &seq2_range);
            assert_eq!(range, IntRange::from_pair(-2, 0));
            assert_eq!(seq, b"DD");
        }
    }

    #[test]
    fn test_get_indel_homology_range_edge_checks() {
        {
            // Bump into left edge:
            let seq1 = b"DDDDABC";
            let seq2 = b"DDDDDDABC";

            let seq1_range = IntRange::from_pair(2, 2);
            let seq2_range = IntRange::from_pair(3, 4);
            let (range, seq) =
                get_indel_breakend_homology_info(seq2, &seq2_range, seq1, &seq1_range);
            assert_eq!(range, IntRange::from_pair(-2, 2));
            assert_eq!(seq, b"DDDD");
        }
        {
            // Bump into right edge:
            let seq1 = b"ABCDDDD";
            let seq2 = b"ABCDDDDDD";

            let seq1_range = IntRange::from_pair(3, 3);
            let seq2_range = IntRange::from_pair(3, 4);
            let (range, seq) =
                get_indel_breakend_homology_info(seq2, &seq2_range, seq1, &seq1_range);
            assert_eq!(range, IntRange::from_pair(-0, 4));
            assert_eq!(seq, b"DDDD");
        }
    }
}
