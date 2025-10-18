use rust_vc_utils::{ChromList, GenomeSegment, get_int_range_dir_distance, get_int_range_distance};

/// Strict order here means that `a` comes before `b` without intersection
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

#[allow(dead_code)]
pub fn get_target_segments(chrom_list: &ChromList, regions: &[String]) -> Vec<GenomeSegment> {
    let mut rval = Vec::new();
    for region in regions.iter() {
        rval.push(GenomeSegment::from_region_str(chrom_list, region));
    }
    rval
}
