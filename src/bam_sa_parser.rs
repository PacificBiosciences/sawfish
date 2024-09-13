use rust_htslib::bam::record::{CigarString, Record};
use rust_vc_utils::bam_utils::aux::get_string_aux_tag;
use rust_vc_utils::bam_utils::cigar::get_complete_read_clip_positions;
use rust_vc_utils::ChromList;

use crate::bam_utils::has_aligned_segments;

/// Object to directly represent one segment from a BAM split alignment
pub struct SplitReadSegment {
    /// reference sequence name
    pub rname: String,

    /// reference zero-indexed alignment start position
    pub pos: i64,

    /// Alignment using rust_htslib::bam::record::Cigar object
    pub cigar: CigarString,

    pub is_fwd_strand: bool,

    /// mapping quality
    pub mapq: u8,

    /// alignment edit distance
    pub _nm: i32,
}

/// Parse one segment from the bam SA aux tag string into a split alignment object
///
pub fn parse_sa_segment(seg: &str) -> SplitReadSegment {
    let sa_fields = seg.split_terminator(',').collect::<Vec<_>>();
    assert_eq!(
        sa_fields.len(),
        6,
        "Unexpected segment in bam SA tag: {seg}"
    );
    let rname = sa_fields[0].to_string();
    let pos = sa_fields[1].parse::<i64>().unwrap() - 1;
    let is_fwd_strand = sa_fields[2] == "+";
    let cigar = CigarString::try_from(sa_fields[3].as_bytes()).unwrap();
    let mapq = sa_fields[4].parse::<u8>().unwrap();
    let _nm = sa_fields[5].parse::<i32>().unwrap();
    SplitReadSegment {
        rname,
        pos,
        is_fwd_strand,
        cigar,
        mapq,
        _nm,
    }
}

/// Split the bam SA aux tag into each supplementary alignment, and parse each into a split alignment object
///
pub fn parse_sa_aux_val(sa_aux_val: &str) -> Vec<SplitReadSegment> {
    sa_aux_val
        .split_terminator(';')
        .map(parse_sa_segment)
        .collect::<Vec<SplitReadSegment>>()
}

/// Object summarizing information for a single segment of a split read alignment
///
/// This is an extension of the simpler SplitReadSegment structure directly parsed from the bam SA tag.
/// It includes read position information consistenty translated in sequencing order, and is designed
/// to represent both primary and supplemenary alignments.
///
#[derive(Debug)]
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

    /// This is set to true for the split segment derived directly from the bam record presented
    /// to the fwd_read_split_segments parsing routine.
    ///
    /// All split segments found from SA tags in this bam record should be false
    pub from_primary_bam_record: bool,
}

/// Parse all read segments from split-read bam record, and order segments in the read sequencing order
///
/// All segments are ordered by their read start position, where the read position is always consistently
/// expressed in the sequencing order of the read.
///
pub fn get_seq_order_read_split_segments(
    chrom_list: &ChromList,
    record: &Record,
) -> Vec<SeqOrderSplitReadSegment> {
    const SA_AUX_TAG: &[u8] = b"SA";
    let sa_aux_val = get_string_aux_tag(record, SA_AUX_TAG);
    let sa_segments = parse_sa_aux_val(&sa_aux_val);

    let mut seq_order_read_split_segments = Vec::new();

    /// Get the start and end positions of a split read segment in read coordinates oriented in original sequencing order
    ///
    /// Given read_start and read_end computed for the read segment in the current alignment direction, we need to reverse the
    /// computation of start and end to get the read coordinates back to sequencing order.
    ///
    /// Note that start and end here correspond to the bed-style range coordinates (zero-indexed, half-closed). This means that
    /// the reverse stand flip in the logic below does not preserve exact read positions, but also flips the start and end to the
    /// respective right side of the breakend in each order.
    ///
    /// Example:
    /// For a 100 base read segment mapped in reverse orientation to the reference with CIGAR string 80S15M5S, implying start and end
    /// read positions of 80 to 95 in segment mapping orieintation, this method reports the position tuple (5,20)
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

    // Add the primary alignment first:
    let primary_read_size = {
        let (read_start, read_end, read_size) = get_complete_read_clip_positions(&record.cigar());
        let (seq_order_read_start, seq_order_read_end) =
            get_seq_order_read_pos(read_start, read_end, read_size, !record.is_reverse());
        seq_order_read_split_segments.push(SeqOrderSplitReadSegment {
            seq_order_read_start,
            seq_order_read_end,
            chrom_index: record.tid() as usize,
            pos: record.pos(),
            is_fwd_strand: !record.is_reverse(),
            cigar: record.cigar().take(),
            mapq: record.mapq(),
            from_primary_bam_record: true,
        });
        read_size
    };

    for sa_segment in sa_segments.iter() {
        if !has_aligned_segments(&sa_segment.cigar) {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            panic!("Bam record split segment id unaligned in read {qname}");
        }
        let (read_start, read_end, read_size) = get_complete_read_clip_positions(&sa_segment.cigar);
        assert_eq!(primary_read_size, read_size);
        let (seq_order_read_start, seq_order_read_end) =
            get_seq_order_read_pos(read_start, read_end, read_size, sa_segment.is_fwd_strand);
        let chrom_index = *chrom_list.label_to_index.get(&sa_segment.rname).unwrap();
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

    // Sanity check the final split read set to ensure that start and end coordinates define a non-empty range in each segment:
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

    #[test]
    fn test_parse_sa_aux_val() {
        let test_val = "chr3,10001,+,5535S10=1D39=2X11438S,60,192;\
        chr3,10001,+,3073S15=2D20=2X11=1X5=1I23=1X5=14798S,22,44;\
        chr4,106872270,-,23=1I226=1I195=1X147=1D1021=7362S,60,19;";

        let result = parse_sa_aux_val(test_val);

        assert_eq!(result.len(), 3);
        assert_eq!(result[2].rname, "chr4");
        assert_eq!(result[1].pos, 10_000);
        assert!(!result[2].is_fwd_strand);
    }
}
