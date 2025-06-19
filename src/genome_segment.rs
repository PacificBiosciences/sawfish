use std::fmt;

use rust_vc_utils::ChromList;

use crate::int_range::get_int_range_dir_distance;
pub use crate::int_range::{IntRange, get_int_range_distance};

/// The structure represents a contiguous region of the genome on a single chromosome
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct GenomeSegment {
    /// chrom_index is defined by the indexing scheme used in the input bam file
    pub chrom_index: usize,
    pub range: IntRange,
}

impl GenomeSegment {
    pub fn new() -> Self {
        Self {
            chrom_index: 0,
            range: IntRange::new(),
        }
    }

    /// Convert from a string in 'samtools' region format (e.g. chr20:100-200)
    ///
    #[allow(dead_code)]
    pub fn from_region_str(chrom_list: &ChromList, str: &str) -> Self {
        let (chrom_index, start, end) = parse_samtools_region_string(chrom_list, str);
        Self {
            chrom_index,
            range: IntRange::from_pair(start, end),
        }
    }

    /// Convert to a string in 'samtools' region format (e.g. chr20:100-200)
    ///
    pub fn to_region_str(&self, chrom_list: &ChromList) -> String {
        let chrom = &chrom_list.data[self.chrom_index].label;
        format!("{chrom}:{}-{}", self.range.start + 1, self.range.end)
    }

    #[allow(dead_code)]
    pub fn intersect(&self, other: &Self) -> bool {
        self.chrom_index == other.chrom_index && self.range.intersect_range(&other.range)
    }

    /// Expand the genomic region, restricted by the bounds of the chromosome length
    pub fn expand_by(&mut self, chrom_list: &ChromList, size: i64) -> (i64, i64) {
        self.asymmetric_expand_by(chrom_list, size, size)
    }

    /// Expand the genomic region separately on left and right sides, restricted by the bounds of
    /// the chromosome length
    ///
    /// Note the expansion could technically be negative, but it's up to the client to ensure that start
    /// and end haven't flipped into an invalid range statement in this case
    ///
    /// Returns the actual left and right expansion after chromosome end clipping
    ///
    pub fn asymmetric_expand_by(
        &mut self,
        chrom_list: &ChromList,
        left_size: i64,
        right_size: i64,
    ) -> (i64, i64) {
        let chrom_size = chrom_list.data[self.chrom_index].length;
        let new_start = std::cmp::max(self.range.start - left_size, 0);
        let new_end = std::cmp::min(self.range.end + right_size, chrom_size as i64);
        let left_shift = self.range.start - new_start;
        let right_shift = new_end - self.range.end;
        self.range.start = new_start;
        self.range.end = new_end;
        (left_shift, right_shift)
    }
}

impl fmt::Debug for GenomeSegment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Segment: {}:{:?}", self.chrom_index, self.range)
    }
}

/// Srict order here means that `a` comes before `b` without intersection
///
pub fn is_strict_order(a: &GenomeSegment, b: &GenomeSegment) -> bool {
    a.chrom_index < b.chrom_index
        || (a.chrom_index == b.chrom_index && a.range.end <= b.range.start)
}

/// Returns None if on different chromosomes, Return 0 if the ranges intersect or are adjacent
///
pub fn get_segment_distance(gs1: &GenomeSegment, gs2: &GenomeSegment) -> Option<usize> {
    if gs1.chrom_index != gs2.chrom_index {
        None
    } else {
        Some(get_int_range_distance(&gs1.range, &gs2.range))
    }
}

/// Returns None if on different chromosomes, Return (direction, distance) between the segments
/// otherwise. The distance is 0 if the ranges intersect or are adjacent.
///
pub fn get_segment_dir_distance(gs1: &GenomeSegment, gs2: &GenomeSegment) -> Option<(bool, usize)> {
    if gs1.chrom_index != gs2.chrom_index {
        None
    } else {
        Some(get_int_range_dir_distance(&gs1.range, &gs2.range))
    }
}

/// Parse the chromosome string out of a samtools-style region string
///
/// Return the index of the chromosome from the expected chromosome list, and
/// an optional position string following the chromosome name
///
fn parse_chrom_index_from_samtools_region_string<'a>(
    chrom_list: &ChromList,
    str: &'a str,
) -> (usize, Option<&'a str>) {
    // Note that rsplitn orders words in reverse order compared to how they appear in the string:
    let s1 = str.rsplitn(2, ':').collect::<Vec<_>>();
    let s1l = s1.len();
    assert!(
        s1l > 0 && s1l < 3,
        "Unexpected format in genome region string '{str}'"
    );
    let chrom = *s1.last().unwrap();
    if let Some(&chrom_index) = chrom_list.label_to_index.get(chrom) {
        let pos_string = if s1l == 2 { Some(s1[0]) } else { None };
        (chrom_index, pos_string)
    } else if let Some(&chrom_index) = chrom_list.label_to_index.get(str) {
        (chrom_index, None)
    } else {
        let msg = if str != chrom {
            format!(
                "Unexpected format in genome region string '{str}': can't find chromosome name '{chrom}' or '{str}' in bam file header"
            )
        } else {
            format!(
                "Unexpected format in genome region string '{str}': can't find chromosome '{chrom}' in bam file header"
            )
        };
        panic!("{}", msg);
    }
}

/// Parse position range from samtools-style genomic interval string, return
/// start-end coordinate in bedtools zero-index half-open format.
///
/// In the samtools-style string, "100-300" would return (99,300). Just "100"
/// should retunr (99, chrom_length)
///
/// # Arguments
/// * `region_str` - Only used to improve error messages
///
fn parse_samtools_pos_range(
    region_str: &str,
    pos_range_str: Option<&str>,
    chrom_size: i64,
) -> (i64, i64) {
    if let Some(pos_range_str) = pos_range_str {
        let s2 = pos_range_str.split('-').collect::<Vec<_>>();
        let s2l = s2.len();
        assert!(
            s2l <= 2,
            "Unexpected format in position range '{pos_range_str}' from genome region string {region_str}"
        );

        // Strip any commas out of the number field (same as tabix cmdline behavior)
        let s2 = s2
            .into_iter()
            .map(|s| {
                let mut s = String::from(s);
                s.retain(|c| c != ',');
                s
            })
            .collect::<Vec<_>>();
        let start = s2[0].parse::<i64>().unwrap() - 1;
        let end = if s2l == 1 {
            chrom_size
        } else {
            s2[1].parse::<i64>().unwrap()
        };
        (start, end)
    } else {
        (0, chrom_size)
    }
}

/// Convert from a string in 'samtools' region format (e.g. chr20:100-200) to a tuple of
/// (chrom_index, chrom_label, start, end)
/// ...where start and end are converted to the zero-indexed half-open convention used for bed
///
/// Commas will be stripped out of coordinates if present
///
/// This parser makes a 'best-effort' to parse contig names with colons in them, such as HLA alleles
/// like "HLA-DRB1*10:01:01". Given that samtools region format already has an optinoal colon, it may
/// be impossible to resolve some cases.
///
pub fn parse_samtools_region_string(chrom_list: &ChromList, region_str: &str) -> (usize, i64, i64) {
    let (chrom_index, pos_str) =
        parse_chrom_index_from_samtools_region_string(chrom_list, region_str);
    let chrom_size = chrom_list.data[chrom_index].length as i64;
    let (start, end) = parse_samtools_pos_range(region_str, pos_str, chrom_size);
    (chrom_index, start, end)
}

#[allow(dead_code)]
pub fn get_target_segments(chrom_list: &ChromList, regions: &[String]) -> Vec<GenomeSegment> {
    let mut rval = Vec::new();
    for region in regions.iter() {
        rval.push(GenomeSegment::from_region_str(chrom_list, region));
    }
    rval
}

#[cfg(test)]
mod tests {
    use super::*;

    /// This test makes sure the auto-generated ordering for GenomeSegment is doing what we assume
    ///
    #[test]
    fn test_segment_order() {
        // Ensure chrom_index has priority over pos
        let segment1 = GenomeSegment {
            chrom_index: 0,
            range: IntRange::from_int(10),
        };
        let segment2 = GenomeSegment {
            chrom_index: 1,
            range: IntRange::from_int(1),
        };
        assert!(segment1 < segment2);

        // Ensure begin pos has priority over end pos
        let segment1 = GenomeSegment {
            chrom_index: 0,
            range: IntRange::from_pair(1, 20),
        };
        let segment2 = GenomeSegment {
            chrom_index: 0,
            range: IntRange::from_pair(10, 11),
        };
        assert!(segment1 < segment2);

        // Ensure that equal segments are not gt
        let segment1 = GenomeSegment {
            chrom_index: 1,
            range: IntRange::from_int(10),
        };
        let segment2 = GenomeSegment {
            chrom_index: 1,
            range: IntRange::from_int(10),
        };
        assert!(segment1 >= segment2);
    }

    #[test]
    fn test_to_region_string() {
        let mut chrom_list = ChromList::default();
        chrom_list.add_chrom("chr1", 100);
        chrom_list.add_chrom("chr2", 100);

        let segment1 = GenomeSegment {
            chrom_index: 1,
            range: IntRange::from_int(10),
        };

        assert_eq!(
            segment1.to_region_str(&chrom_list),
            "chr2:11-11".to_string()
        );
    }

    #[test]
    fn test_samtools_region_string_splitter() {
        let chrom_list = {
            let mut x = ChromList::default();
            x.add_chrom("chr1", 10000);
            x.add_chrom("chr2", 10000);
            x.add_chrom("chr3", 10000);
            x
        };

        // A simple case
        let s = "chr2:1000-2000";
        let (chrom_index, start, end) = parse_samtools_region_string(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(start, 999);
        assert_eq!(end, 2000);

        // Simple case with commas
        let s = "chr2:1,000-2,000";
        let (chrom_index, start, end) = parse_samtools_region_string(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(start, 999);
        assert_eq!(end, 2000);

        // No end
        let s = "chr2:1,000";
        let (chrom_index, start, end) = parse_samtools_region_string(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(start, 999);
        assert_eq!(end, 10000);
    }

    #[test]
    fn test_samtools_region_string_splitter_hla() {
        let chrom_list = {
            let mut x = ChromList::default();
            x.add_chrom("HLA-DRB1*10:01:01", 10000);
            x
        };

        let s = "HLA-DRB1*10:01:01:1000-2000";
        let (chrom_index, start, end) = parse_samtools_region_string(&chrom_list, s);
        assert_eq!(chrom_index, 0);
        assert_eq!(start, 999);
        assert_eq!(end, 2000);

        let s = "HLA-DRB1*10:01:01";
        let (chrom_index, start, end) = parse_samtools_region_string(&chrom_list, s);
        assert_eq!(chrom_index, 0);
        assert_eq!(start, 0);
        assert_eq!(end, 10000);
    }
}
