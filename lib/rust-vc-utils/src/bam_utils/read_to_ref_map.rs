use std::collections::BTreeMap;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use super::cigar::update_ref_and_read_pos;

/// Create an array of read length, mapping from the read to the reference position for the given read alignment
///
/// Allowing the ref_pos and cigar alignment input enables the processing of split alignment segments with this
/// function.
///
/// Read positions which are not mapped to a reference position are None. Read positions include hard-clipped segments.
///
/// Ref positions are 0-indexed
///
pub fn get_read_segment_to_ref_pos_map(
    seq_len: usize,
    mut ref_pos: i64,
    cigar: &[Cigar],
    ignore_hard_clip: bool,
) -> Vec<Option<i64>> {
    use rust_htslib::bam::record::Cigar::*;

    let mut read_to_ref = vec![None; seq_len];
    let mut read_pos = 0usize;

    for c in cigar.iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                let len = *len as usize;
                for i in 0..len {
                    read_to_ref[read_pos + i] = Some(ref_pos + i as i64);
                }
            }
            _ => {}
        }
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    read_to_ref
}

/// Create an array of read length, mapping from the read to the reference position for the read's primary alignment
/// segment
///
/// Read positions which are not mapped to a reference position are None. Read positions include hard-clipped segments.
///
/// Ref positions are 0-indexed
///
pub fn get_read_to_ref_pos_map(record: &bam::Record, ignore_hard_clip: bool) -> Vec<Option<i64>> {
    get_read_segment_to_ref_pos_map(
        record.seq_len(),
        record.pos(),
        &record.cigar(),
        ignore_hard_clip,
    )
}

#[derive(Clone, Default)]
pub struct ReadToRefTreeMap {
    // Key is the read_pos starting a block
    // Val is None for read not mapping to ref, or the ref_pos at the start of a match block
    map: BTreeMap<usize, Option<i64>>,
}

impl ReadToRefTreeMap {
    pub fn get_ref_pos(&self, read_pos: usize) -> Option<i64> {
        match self.map.range(..=read_pos).next_back() {
            Some((k, v)) => v.map(|ref_pos| ref_pos + (read_pos - k) as i64),
            None => None,
        }
    }

    pub fn get_ref_range(
        &self,
        read_start_pos: usize,
        read_end_pos: usize,
    ) -> std::collections::btree_map::Range<'_, usize, Option<i64>> {
        let read_start_block_pos = match self.map.range(..=read_start_pos).next_back() {
            Some((k, _)) => *k,
            None => read_start_pos,
        };

        self.map.range(read_start_block_pos..read_end_pos)
    }

    pub fn get_map(&self) -> &BTreeMap<usize, Option<i64>> {
        &self.map
    }
}

/// Create a data structure enabling rapid lookup from the read to the reference position for the given read alignment
///
/// Allowing the ref_pos and cigar alignment input enables the processing of split alignment segments with this
/// function.
///
/// Read positions which are not mapped to a reference position are None. Read positions include hard-clipped segments.
///
/// Ref positions are 0-indexed
///
pub fn get_read_segment_to_ref_pos_tree_map(
    mut ref_pos: i64,
    cigar: &[Cigar],
    ignore_hard_clip: bool,
) -> ReadToRefTreeMap {
    use rust_htslib::bam::record::Cigar::*;

    let mut read_to_ref = ReadToRefTreeMap::default();
    let mut read_pos = 0usize;

    let mut update_map = |ref_pos: i64, read_pos: usize, match_len: &mut usize| {
        if *match_len > 0 {
            read_to_ref
                .map
                .insert(read_pos - *match_len, Some(ref_pos - *match_len as i64));
            read_to_ref.map.insert(read_pos, None);
            *match_len = 0;
        }
    };

    let mut match_len = 0;

    for c in cigar.iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                match_len += *len as usize;
            }
            _ => {
                update_map(ref_pos, read_pos, &mut match_len);
            }
        }
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    update_map(ref_pos, read_pos, &mut match_len);

    read_to_ref
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Header, HeaderView, header};

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
    fn test_get_read_to_ref_pos_map() {
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M\t*\t0\t0\tACTG\tDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let rval = get_read_to_ref_pos_map(&record, true);
        assert_eq!(rval, vec![Some(9), Some(10), None, Some(11)]);
    }

    #[test]
    fn test_get_read_segment_to_ref_pos_tree_map() {
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M\t*\t0\t0\tACTG\tDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let rval_tree = get_read_segment_to_ref_pos_tree_map(record.pos(), &record.cigar(), true);
        let rval = (0..4).map(|x| rval_tree.get_ref_pos(x)).collect::<Vec<_>>();
        assert_eq!(rval, vec![Some(9), Some(10), None, Some(11)]);

        let rrange = rval_tree.get_ref_range(0, 2).collect::<Vec<_>>();
        assert_eq!(rrange, vec![(&0, &Some(9))]);
    }
}
