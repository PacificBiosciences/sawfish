use std::fmt;

use rust_vc_utils::ChromList;

use crate::int_range::get_int_range_dir_distance;
pub use crate::int_range::{get_int_range_distance, IntRange};

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
        let (chrom_index, _chrom_label, start, end) =
            samtools_region_string_splitter(chrom_list, str);
        Self {
            chrom_index,
            range: IntRange::from_pair(start, end),
        }
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

/// Convert from a string in 'samtools' region format (e.g. chr20:100-200) to a tuple of
/// (chrom_index, chrom_label, start, end)
/// ...where start and end are converted to the zero-indexed half-open convention used for bed
///
/// Commas will be stripped out of coordinates if present
///
pub fn samtools_region_string_splitter(
    chrom_list: &ChromList,
    str: &str,
) -> (usize, String, i64, i64) {
    let s1 = str.split(':').collect::<Vec<_>>();
    let s1l = s1.len();
    assert!(
        s1l > 0 && s1l < 3,
        "Unexpected format in genome region string {}",
        str
    );
    let chrom = s1[0].to_string();
    let chrom_index = match chrom_list.label_to_index.get(s1[0]) {
        Some(x) => *x,
        None => {
            panic!("Can't find chromosome '{}' in bam file header", s1[0]);
        }
    };
    let chrom_size = chrom_list.data[chrom_index].length as i64;
    let (start, end) = if s1l == 1 {
        (0, chrom_size)
    } else {
        let s2 = s1[1].split('-').collect::<Vec<_>>();
        let s2l = s2.len();
        assert!(
            s2l > 0 && s2l < 3,
            "Unexpected format in genome region string {}",
            str
        );
        // Strip any commas out of the number field (don't know if samtools does this but just a
        // nice ease of use bonus:
        let s2 = s2
            .into_iter()
            .map(|s| {
                let mut s = String::from(s);
                s.retain(|c| c != ',');
                s
            })
            .collect::<Vec<_>>();
        let start = s2[0].parse::<i64>().unwrap() - 1;
        if s2l == 1 {
            (start, chrom_size)
        } else {
            let end = s2[1].parse::<i64>().unwrap();
            (start, end)
        }
    };
    (chrom_index, chrom, start, end)
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
    fn test_samtools_region_string_splitter() {
        let mut chrom_list = ChromList::default();
        chrom_list.add_chrom("chr1", 10000);
        chrom_list.add_chrom("chr2", 10000);
        chrom_list.add_chrom("chr3", 10000);
        let chrom_list = chrom_list;

        // A simple case
        let s = "chr2:1000-2000";
        let (chrom_index, chrom_label, start, end) =
            samtools_region_string_splitter(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(chrom_label, "chr2");
        assert_eq!(start, 999);
        assert_eq!(end, 2000);

        // Simple case with commas
        let s = "chr2:1,000-2,000";
        let (chrom_index, chrom_label, start, end) =
            samtools_region_string_splitter(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(chrom_label, "chr2");
        assert_eq!(start, 999);
        assert_eq!(end, 2000);

        // No end
        let s = "chr2:1,000";
        let (chrom_index, chrom_label, start, end) =
            samtools_region_string_splitter(&chrom_list, s);
        assert_eq!(chrom_index, 1);
        assert_eq!(chrom_label, "chr2");
        assert_eq!(start, 999);
        assert_eq!(end, 10000);
    }
}
