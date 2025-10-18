//! Miscelanious BAM-related utilities that don't fit any other modules
//!
use rust_htslib::bam;

/// Recreation of htslib hts_reg2bin in rust, I've copied this from noodles (MIT License)
/// with minor type adaptions.
///
/// begin and end should follow bed zero-based half-closed format
///
fn hts_reg2bin(begin: usize, end: usize, min_shift: u8, depth: u8) -> usize {
    let end = end - 1;
    let mut l = depth;
    let mut s = min_shift;
    let mut t = ((1 << (depth * 3)) - 1) / 7;

    while l > 0 {
        if begin >> s == end >> s {
            return t + (begin >> s);
        }

        l -= 1;
        s += 3;
        t -= 1 << (l * 3);
    }

    0
}

/// Recreation of htslib bam_reg2bin in rust
///
/// begin and end should follow bed zero-based half-closed format
///
pub fn bam_reg2bin(begin: usize, end: usize) -> u16 {
    hts_reg2bin(begin, end, 14, 5) as u16
}

/// Extract sample name from bam header
///
/// This uses the sample name from the first read group found in the header, and does not
/// check for additional read groups. If no read group is found, `default_name` will
/// be used.
///
pub fn get_sample_name(header: &bam::HeaderView, default_name: &str) -> String {
    for line in std::str::from_utf8(header.as_bytes()).unwrap().split('\n') {
        for (i, word) in line.split('\t').enumerate() {
            if i == 0 {
                if word != "@RG" {
                    break;
                }
            } else if let Some(sample_name) = word.strip_prefix("SM:") {
                return sample_name.to_string();
            }
        }
    }
    default_name.to_string()
}
