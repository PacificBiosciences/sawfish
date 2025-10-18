use rust_htslib::bam::record::{CigarString, Record};

use super::aux::get_optional_string_aux_tag;
use super::aux::sa_tag_parser::parse_sa_aux_val;
use super::cigar::{get_cigar_ref_offset, get_read_clip_positions, has_aligned_segments};
use crate::ChromList;

/// Object summarizing information for a single segment of a split read alignment
///
/// This is an extension of the simpler SplitReadSegment structure directly parsed from the bam SA tag. It includes read
/// position information consistently translated in sequencing order, and is designed to represent both primary and
/// supplementary alignments.
///
#[derive(Debug, PartialEq)]
pub struct SeqOrderSplitReadSegment {
    /// Sequencer-order read position of left-most position of the alignment segment
    pub seq_order_read_start: usize,

    /// Sequencer-order read position one base after the right-most position of the alignment segment
    pub seq_order_read_end: usize,
    pub chrom_index: usize,
    pub pos: i64,
    pub is_fwd_strand: bool,
    pub cigar: CigarString,
    pub mapq: u8,

    /// This is set to true for the split segment derived directly from the bam record presented to the
    /// fwd_read_split_segments parsing routine.
    ///
    /// All split segments found from SA tags in this bam record should be false
    pub from_primary_bam_record: bool,
}

impl SeqOrderSplitReadSegment {
    /// The derived Debug output can be dominated by the Cigar output, so provide a more compact display option
    pub fn short_display(&self) -> String {
        let end = self.pos + get_cigar_ref_offset(&self.cigar);
        format!(
            "seq_order_read_start/end: {}/{} ref_segment: {}:{}-{} fwd: {} mapq: {}",
            self.seq_order_read_start,
            self.seq_order_read_end,
            self.chrom_index,
            self.pos,
            end,
            self.is_fwd_strand,
            self.mapq
        )
    }
}

/// Parse all primary- and split-read segments bam record, and order segments in the read sequencing order
///
/// All segments are ordered by their read start position, where the read position is always consistently expressed in
/// the sequencing order of the read.
///
pub fn get_seq_order_read_split_segments(
    chrom_list: &ChromList,
    record: &Record,
) -> Vec<SeqOrderSplitReadSegment> {
    let ignore_hard_clip = false;

    /// Get the start and end positions of a split read segment in read coordinates oriented in original sequencing
    /// order
    ///
    /// Given read_start and read_end computed for the read segment in the current alignment direction, we need to
    /// reverse the computation of start and end to get the read coordinates back to sequencing order.
    ///
    /// Note that start and end here correspond to the bed-style range coordinates (zero-indexed, half-closed). This
    /// means that the reverse stand flip in the logic below does not preserve exact read positions, but also flips the
    /// start and end to the respective right side of the breakend in each order.
    ///
    /// Example: For a 100 base read segment mapped in reverse orientation to the reference with CIGAR string 80S15M5S,
    /// implying start and end read positions of 80 to 95 in segment mapping orieintation, this method reports the
    /// position tuple (5,20)
    ///
    fn get_seq_order_read_pos(
        read_start: usize,
        read_end: usize,
        read_size: usize,
        is_fwd_strand: bool,
    ) -> (usize, usize) {
        if is_fwd_strand {
            (read_start, read_end)
        } else {
            (read_size - read_end, read_size - read_start)
        }
    }

    let (primary_read_size, mut seq_order_read_split_segments) = {
        let (read_start, read_end, read_size) =
            get_read_clip_positions(&record.cigar(), ignore_hard_clip);
        let (seq_order_read_start, seq_order_read_end) =
            get_seq_order_read_pos(read_start, read_end, read_size, !record.is_reverse());
        let primary_segment = SeqOrderSplitReadSegment {
            seq_order_read_start,
            seq_order_read_end,
            chrom_index: record.tid() as usize,
            pos: record.pos(),
            is_fwd_strand: !record.is_reverse(),
            cigar: record.cigar().take(),
            mapq: record.mapq(),
            from_primary_bam_record: true,
        };
        (read_size, vec![primary_segment])
    };

    if let Some(sa_aux_val) = get_optional_string_aux_tag(record, b"SA") {
        let sa_segments = parse_sa_aux_val(&sa_aux_val);

        for (sa_segment_index, sa_segment) in sa_segments.iter().enumerate() {
            if !has_aligned_segments(&sa_segment.cigar) {
                let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                panic!("Bam record split segment id unaligned in read {qname}");
            }
            let (read_start, read_end, read_size) =
                get_read_clip_positions(&sa_segment.cigar, ignore_hard_clip);
            assert_eq!(primary_read_size, read_size);
            let (seq_order_read_start, seq_order_read_end) =
                get_seq_order_read_pos(read_start, read_end, read_size, sa_segment.is_fwd_strand);
            let chrom_index = match chrom_list.label_to_index.get(&sa_segment.rname) {
                Some(&x) => x,
                None => {
                    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                    panic!(
                        "In read '{qname}', the SA aux tag describes a split read mapped to {}:{} (in segment {}), which is not found in the input reference fasta",
                        sa_segment.rname, sa_segment.pos, sa_segment_index
                    );
                }
            };
            seq_order_read_split_segments.push({
                SeqOrderSplitReadSegment {
                    seq_order_read_start,
                    seq_order_read_end,
                    chrom_index,
                    pos: sa_segment.pos,
                    is_fwd_strand: sa_segment.is_fwd_strand,
                    cigar: sa_segment.cigar.clone(),
                    mapq: sa_segment.mapq,
                    from_primary_bam_record: false,
                }
            });
        }

        seq_order_read_split_segments.sort_by_key(|x| x.seq_order_read_start);
    }

    // Sanity check the final split read set to ensure that start and end coordinates define a non-empty range in each
    // segment:
    //
    for s in seq_order_read_split_segments.iter() {
        if s.seq_order_read_start >= s.seq_order_read_end {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            panic!(
                "Can't parse consistent split read information from SA tag format in read: {qname}"
            );
        };
    }

    seq_order_read_split_segments
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;
    use rust_htslib::bam::{self, Header, HeaderView, header};

    fn get_test_header() -> HeaderView {
        let mut _header = Header::new();
        for i in 0..3 {
            let label = format!("chr{i}");
            _header.push_record(
                header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", label.as_str())
                    .push_tag(b"LN", 1000),
            );
        }
        HeaderView::from_header(&_header)
    }

    fn get_expected_seg(
        seq_order_read_start: usize,
        seq_order_read_end: usize,
        chrom_index: usize,
        pos: i64,
        is_fwd_strand: bool,
        cigar_str: &str,
        from_primary_bam_record: bool,
    ) -> SeqOrderSplitReadSegment {
        SeqOrderSplitReadSegment {
            seq_order_read_start,
            seq_order_read_end,
            chrom_index,
            pos,
            is_fwd_strand,
            cigar: CigarString::try_from(cigar_str).unwrap(),
            mapq: 60,
            from_primary_bam_record,
        }
    }

    #[test]
    fn test_get_seq_order_read_split_segments() {
        let header = get_test_header();
        let chrom_list = &ChromList::from_bam_header(&header);

        // Test degenerate case of no split reads
        let sam_line =
            "qname\t0\tchr2\t10\t60\t10S5M5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line.as_bytes()).unwrap();

        let result = get_seq_order_read_split_segments(chrom_list, &rec);

        assert_eq!(
            result,
            vec![get_expected_seg(10, 15, 2, 9, true, "10S5M5S", true),]
        );

        // Test case of multiple split-reads and primary reads in a random order wrt sequencing order
        // Test sequence is 20 bases, and split into 4 5-base alignments:
        // |chr1:200-| |chr0:100+| |chr2:10+| |chr0:20-|
        let sam_line = "qname\t0\tchr2\t10\t60\t10S5M5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\t\
            SA:Z:chr0,20,-,5M15S,60,0;chr0,100,+,5S5M10S,60,0;chr1,200,-,15S5M,60,0;";
        let rec = bam::Record::from_sam(&header, sam_line.as_bytes()).unwrap();

        let result = get_seq_order_read_split_segments(chrom_list, &rec);

        assert_eq!(
            result,
            vec![
                get_expected_seg(0, 5, 1, 199, false, "15S5M", false),
                get_expected_seg(5, 10, 0, 99, true, "5S5M10S", false),
                get_expected_seg(10, 15, 2, 9, true, "10S5M5S", true),
                get_expected_seg(15, 20, 0, 19, false, "5M15S", false),
            ]
        );
    }
}
