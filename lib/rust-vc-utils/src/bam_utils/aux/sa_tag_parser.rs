use rust_htslib::bam::record::CigarString;

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
/// This routine does check for the SA tag conforming to the expected syntax, but doesn't add any checking that contig
/// names are restricted to an expected set, or that the contig length is within an expected range. All such checks are
/// left to the client.
///
pub fn parse_sa_aux_val(sa_aux_val: &str) -> Vec<SplitReadSegment> {
    sa_aux_val
        .split_terminator(';')
        .map(parse_sa_segment)
        .collect::<Vec<SplitReadSegment>>()
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
