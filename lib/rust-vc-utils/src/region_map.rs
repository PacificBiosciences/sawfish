//! A light wrapper on top of an interval tree for genomic region applications
//!

use bio::data_structures::interval_tree::{IntervalTree, IntervalTreeIterator};

/// A set of regions in a contiguous coordinate space (such as a chromosome) which can be efficiently queried
///
#[derive(Clone, Default)]
pub struct RegionMap {
    regions: IntervalTree<i64, u8>,
}

impl RegionMap {
    /// Return true if the start-end range intersects with any regions stored in this object
    ///
    pub fn intersect(&self, start: i64, end: i64) -> bool {
        //find any region that overlaps
        self.regions.find(start..end).next().is_some()
    }

    /// Return true if `pos` intersects with any regions stored in this object
    ///
    pub fn intersect_pos(&self, pos: i64) -> bool {
        self.intersect(pos, pos + 1)
    }

    pub fn find_overlaps(&self, start: i64, end: i64) -> IntervalTreeIterator<'_, i64, u8> {
        self.regions.find(start..end)
    }

    /// Add region
    ///
    /// Adds a region with the default value, regions are not collapsed
    ///
    pub fn add_region(&mut self, start: i64, end: i64) {
        self.regions.insert(start..end, Default::default());
    }

    /// Add region value
    ///
    /// Adds a value for a particular region, regions are not collapsed
    ///
    pub fn add_region_value(&mut self, start: i64, end: i64, value: u8) {
        self.regions.insert(start..end, value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intersect() {
        let mut regions = RegionMap::default();

        regions.add_region(100, 101);
        assert!(regions.intersect_pos(100));
        assert!(!regions.intersect_pos(99));
        assert!(!regions.intersect_pos(101));
    }

    #[test]
    fn test_find() {
        let mut regions = RegionMap::default();

        regions.add_region_value(100, 200, 0);
        regions.add_region_value(200, 250, 1);
        regions.add_region_value(250, 300, 2);
        let mut values = regions
            .find_overlaps(175, 275)
            .map(|x| *x.data())
            .collect::<Vec<_>>();
        values.sort();

        assert_eq!(values, vec![0, 1, 2]);
    }
}
