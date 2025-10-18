//! Miscelanious BAM record processing utilities
//!
use rust_htslib::{bam, htslib};

use super::cigar::update_ref_pos;

/// Check if the alignment record should be filtered from consideration in any part of the variant
/// calling pipeline
///
pub fn filter_out_alignment_record(record: &bam::Record) -> bool {
    static FLAG_FILTER: u32 =
        htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY | htslib::BAM_FQCFAIL | htslib::BAM_FDUP;

    ((record.flags() as u32) & FLAG_FILTER) != 0
}

/// Report the end reference position of a bam record
///
/// The end position is the zero-indexed right-most mapped position + 1
///
pub fn get_alignment_end(record: &bam::Record) -> i64 {
    let mut ref_pos = record.pos();
    for c in record.cigar().iter() {
        update_ref_pos(c, &mut ref_pos);
    }
    ref_pos
}

/// Translate the (zero-indexed) read position into the corresponding read position in the reverse orientation
///
pub fn get_reverse_read_position(record: &bam::Record, read_pos: usize) -> usize {
    let read_len = record.seq_len();
    if read_pos >= read_len {
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        panic!(
            "Invalid read position {read_pos}, exceeds the read_length {read_len}, in read {qname}"
        );
    }
    read_len - (read_pos + 1)
}

/// Translate the (zero-indexed) read position in fwd-aligned read order to the read position for the
/// read in sequencer order
///
pub fn get_seq_order_read_position(record: &bam::Record, read_pos: usize) -> usize {
    if record.is_reverse() {
        get_reverse_read_position(record, read_pos)
    } else {
        read_pos
    }
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
    fn test_filter_out_alignment_record() {
        let header = get_test_header();

        // Unmapped read:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert!(filter_out_alignment_record(&rec));

        // Mapped read:
        let sam_line =
            b"qname\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert!(!filter_out_alignment_record(&rec));
    }

    #[test]
    fn test_get_alignment_end() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5S5M10D5I5M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert_eq!(get_alignment_end(&rec), 29);
    }

    #[test]
    fn test_get_seq_order_read_position() {
        let header = get_test_header();
        let sam_line =
            b"qname\t16\tchr1\t10\t60\t20M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert_eq!(get_seq_order_read_position(&rec, 1), 18);
    }
}
