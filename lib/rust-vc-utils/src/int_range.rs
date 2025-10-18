use std::fmt;

/// A simple type for integer ranges
///
/// All ranges follow the bed file range convention: 0-indexed, half-closed, [start,end)
///
/// This struct is used instead of the native rust Range type just to focus on the specific goals of
/// primarily genomic region intervals.
///
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct IntRange {
    pub start: i64,
    pub end: i64,
}

impl Default for IntRange {
    fn default() -> Self {
        Self::new()
    }
}

impl IntRange {
    /// initialize to a detectably invalid state:
    pub fn new() -> Self {
        Self::from_int(-1)
    }

    pub fn from_int(start: i64) -> Self {
        Self {
            start,
            end: start + 1,
        }
    }

    pub fn from_pair(start: i64, end: i64) -> Self {
        Self { start, end }
    }

    pub fn size(&self) -> i64 {
        self.end - self.start
    }

    pub fn center(&self) -> i64 {
        (self.start + self.end) / 2
    }

    /// Return true if pos intersects range (adjacency does not count)
    ///
    pub fn intersect_pos(&self, pos: i64) -> bool {
        pos >= self.start && pos < self.end
    }

    /// Return true if the ranges intersect (adjacency does not count)
    ///
    #[allow(dead_code)]
    pub fn intersect_range(&self, other: &IntRange) -> bool {
        other.end >= self.start && other.start < self.end
    }

    pub fn merge(&mut self, other: &IntRange) {
        if other.start < self.start {
            self.start = other.start;
        }
        if other.end > self.end {
            self.end = other.end;
        }
    }

    /// Translate a range to its reversed version, given a certain region size
    ///
    /// # Example
    ///
    /// reverse_range([1,3), 6) -> [3,5)
    ///
    /// Fwd:
    ///  [-)
    /// 012345
    ///
    /// Rev:
    ///    [-)
    /// 012345
    ///
    pub fn reverse(&mut self, size: i64) {
        let istart = self.start;
        self.start = size - self.end;
        self.end = size - istart;
    }

    pub fn get_reverse_range(&self, size: i64) -> Self {
        Self {
            start: size - self.end,
            end: size - self.start,
        }
    }
}

impl fmt::Debug for IntRange {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{}-{})", self.start, self.end)
    }
}

/// Get the distance between 2 ranges
///
/// The notion of 'distance' between two ranges here means the gap between two ranges that don't intersect.
///
///    [---------)            [-----------)
///        R1    -------------     R2
///                R1/R2 dist
///
/// The distance is 0 if the ranges intersect or are adjacent
///
pub fn get_int_range_distance(ir1: &IntRange, ir2: &IntRange) -> usize {
    use std::cmp::max;
    max(max(ir2.start - ir1.end, ir1.start - ir2.end), 0) as usize
}

/// Get the direction and distance between 2 ranges
///
/// The direction is true if ir2 is 'ahead of' ir1, and false otherwise.
///
/// Return a distance of (true, 0) if the ranges intersect or are adjacent
pub fn get_int_range_dir_distance(ir1: &IntRange, ir2: &IntRange) -> (bool, usize) {
    use std::cmp::max;
    let d21 = ir2.start - ir1.end;
    let d12 = ir1.start - ir2.end;
    let (dir, dist) = if d12 > d21 {
        (d12 <= 0, d12)
    } else {
        (true, d21)
    };
    (dir, max(dist, 0) as usize)
}

pub fn get_recip_overlap(r1: &IntRange, r2: &IntRange) -> f64 {
    let min_span = 100;

    let r1_min_end = std::cmp::max(r1.end, r1.start + min_span);
    let r2_min_end = std::cmp::max(r2.end, r2.start + min_span);

    let olap = std::cmp::max(
        std::cmp::min(r1_min_end, r2_min_end) - std::cmp::max(r2.start, r1.start),
        0,
    );
    let span1 = r1_min_end - r1.start;
    let span2 = r2_min_end - r2.start;
    let span = std::cmp::max(span1, span2);

    olap as f64 / span as f64
}

#[allow(dead_code)]
pub fn get_overlap_range(r1: &IntRange, r2: &IntRange) -> Option<IntRange> {
    if !r1.intersect_range(r2) {
        return None;
    }
    Some(IntRange {
        start: std::cmp::max(r1.start, r2.start),
        end: std::cmp::min(r1.end, r2.end),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_int_range_distance() {
        let r1 = IntRange::from_pair(1, 4);
        let r2 = IntRange::from_pair(6, 8);
        let r3 = IntRange::from_pair(8, 9);
        let r4 = IntRange::from_pair(7, 9);

        assert_eq!(get_int_range_distance(&r1, &r2), 2);
        assert_eq!(get_int_range_distance(&r2, &r1), 2);
        assert_eq!(get_int_range_distance(&r2, &r3), 0);
        assert_eq!(get_int_range_distance(&r3, &r2), 0);
        assert_eq!(get_int_range_distance(&r2, &r4), 0);
        assert_eq!(get_int_range_distance(&r4, &r2), 0);
    }

    #[test]
    fn test_get_int_range_dir_distance() {
        let r1 = IntRange::from_pair(1, 4);
        let r2 = IntRange::from_pair(6, 8);
        let r3 = IntRange::from_pair(8, 9);
        let r4 = IntRange::from_pair(7, 9);

        assert_eq!(get_int_range_dir_distance(&r1, &r2), (true, 2));
        assert_eq!(get_int_range_dir_distance(&r2, &r1), (false, 2));
        assert_eq!(get_int_range_dir_distance(&r2, &r3), (true, 0));
        assert_eq!(get_int_range_dir_distance(&r3, &r2), (true, 0));
        assert_eq!(get_int_range_dir_distance(&r2, &r4), (true, 0));
        assert_eq!(get_int_range_dir_distance(&r4, &r2), (true, 0));
    }

    #[test]
    fn test_get_recip_overlap() {
        let r1 = IntRange::from_pair(1000, 2000);

        let r2 = IntRange::from_pair(0, 1000);
        let ro = get_recip_overlap(&r1, &r2);
        approx::assert_ulps_eq!(ro as f32, 0_f32, max_ulps = 4);

        let r2 = IntRange::from_pair(500, 1500);
        let ro = get_recip_overlap(&r1, &r2);
        approx::assert_ulps_eq!(ro as f32, 0.5_f32, max_ulps = 4);

        let r2 = IntRange::from_pair(1000, 2000);
        let ro = get_recip_overlap(&r1, &r2);
        approx::assert_ulps_eq!(ro as f32, 1.0_f32, max_ulps = 4);

        let r2 = IntRange::from_pair(1500, 2500);
        let ro = get_recip_overlap(&r1, &r2);
        approx::assert_ulps_eq!(ro as f32, 0.5_f32, max_ulps = 4);

        let r2 = IntRange::from_pair(1500, 2000);
        let ro = get_recip_overlap(&r1, &r2);
        approx::assert_ulps_eq!(ro as f32, 0.5_f32, max_ulps = 4);
    }

    #[test]
    fn test_reverse() {
        let size = 10_000;

        let r1 = IntRange::from_pair(1000, 2000);

        let mut r2 = r1.clone();

        r2.reverse(size);
        assert_ne!(r1, r2);
        assert_eq!(r2, r1.get_reverse_range(size));

        r2.reverse(size);
        assert_eq!(r1, r2);
    }
}
