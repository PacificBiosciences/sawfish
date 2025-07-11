use std::collections::HashMap;

use bio::data_structures::interval_tree::{IntervalTree, IntervalTreeIterator};
use camino::Utf8Path;
use log::info;
use rust_vc_utils::ChromList;
use unwrap::unwrap;

use crate::genome_segment::parse_samtools_region_string;

/// A set of chromosome regions which can be efficiently queried
///
#[derive(Clone)]
pub struct ChromRegions {
    regions: IntervalTree<i64, u8>,
}

impl ChromRegions {
    pub fn new() -> Self {
        Self {
            regions: IntervalTree::new(),
        }
    }

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

    pub fn find_overlaps(&self, start: i64, end: i64) -> IntervalTreeIterator<i64, u8> {
        self.regions.find(start..end)
    }

    #[allow(dead_code)]
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

#[derive(Clone)]
pub struct GenomeRegions {
    pub chroms: HashMap<String, ChromRegions>,
    pub overlaps_allowed: bool,
}

impl GenomeRegions {
    /// Creates a new genome region lookup
    /// # Arguments
    /// * `overlaps_allowed` - if false, attempts to insert overlapping regions will panic
    pub fn new(overlaps_allowed: bool) -> Self {
        Self {
            chroms: HashMap::new(),
            overlaps_allowed,
        }
    }

    /// Create new object from bed file
    ///
    /// # Arguments
    ///
    /// * `label` - Used to error messages to describe what type of regions file this is
    /// * `overlaps_allowed` - if False, this will panic if any loaded regions overlap with each other
    ///
    pub fn from_bed(
        filename: &str,
        label: &str,
        overlaps_allowed: bool,
        payload_required: bool,
    ) -> Self {
        use rust_htslib::bgzf;
        use std::io::Read;

        info!("Reading {label} regions from file '{filename}'");

        let mut regions = GenomeRegions::new(overlaps_allowed);
        let mut reader = unwrap!(
            bgzf::Reader::from_path(filename),
            "Unable to open {label} regions file: '{filename}'"
        );

        let mut content = String::new();
        unwrap!(
            reader.read_to_string(&mut content),
            "Can't parse text from {label} regions file: '{filename}'"
        );

        for line in content.split('\n') {
            // The last line is expected to be empty
            if line.is_empty() {
                break;
            }

            let words = line.split('\t').collect::<Vec<_>>();

            assert!(words.len() >= 3);
            let chrom = words[0];
            let start = words[1].parse::<i64>().unwrap();
            let end = words[2].parse::<i64>().unwrap();

            //check if we have a value to load for the track
            if words.len() >= 5 {
                let value = words[4].parse::<u8>().unwrap();
                regions.add_region_value(chrom, start, end, value);
            } else {
                assert!(!payload_required);
                regions.add_region(chrom, start, end);
            }
        }

        regions
    }

    /// Create object from strings in 'samtools' region format (e.g. chr20:100-200)
    ///
    pub fn from_target_regions(
        chrom_list: &ChromList,
        target_regions: &Vec<String>,
        overlaps_allowed: bool,
    ) -> Self {
        let mut regions = Self::new(overlaps_allowed);
        for target_region in target_regions {
            let (chrom_index, start, end) = parse_samtools_region_string(chrom_list, target_region);

            let chrom_label = &chrom_list.data[chrom_index].label;
            regions.add_region(chrom_label, start, end);
        }
        regions
    }

    pub fn is_empty(&self) -> bool {
        self.chroms.is_empty()
    }

    /// This will add a region with the default value
    /// # Arguments
    /// * `chrom` - the contig string
    /// * `start` - the start coordinate (included)
    /// * `end` - the end coordinates (excluded)
    pub fn add_region(&mut self, chrom: &str, start: i64, end: i64) {
        self.add_region_value(chrom, start, end, Default::default());
    }

    /// This will add a region with the default value
    /// # Arguments
    /// * `chrom` - the contig string
    /// * `start` - the start coordinate (included)
    /// * `end` - the end coordinates (excluded)
    /// * `value` - the value associated with the region
    pub fn add_region_value(&mut self, chrom: &str, start: i64, end: i64, value: u8) {
        let chrom_regions = self
            .chroms
            .entry(chrom.to_owned())
            .or_insert_with(ChromRegions::new);

        if !self.overlaps_allowed && chrom_regions.intersect(start, end) {
            //TODO: may want to handle this problem as an error instead of a panic, so that e.g.
            // filename can get injected into error msg at higher level
            panic!("Overlaps are not allowed but were detected: {chrom} {start} {end}");
        }

        //check if we have a value to load for the track
        chrom_regions.add_region_value(start, end, value);
    }

    /// This will return an iterator over any overlapping regions. None if the contig does not exist.
    /// # Arguments
    /// * `chrom` - the contig string
    /// * `start` - the start coordinate (included)
    /// * `end` - the end coordinates (excluded)
    #[allow(unused)]
    pub fn find_overlaps(
        &self,
        chrom: &str,
        start: i64,
        end: i64,
    ) -> Option<IntervalTreeIterator<i64, u8>> {
        self.chroms
            .get(chrom)
            .map(|chrom_region| chrom_region.find_overlaps(start, end))
    }
}

pub fn write_genome_regions_to_bed(
    label: &str,
    filename: &Utf8Path,
    chrom_list: &ChromList,
    genome_regions: &[ChromRegions],
) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    info!("Writing {label} genome regions to bed file: '{filename}'");

    let f = unwrap!(
        File::create(filename),
        "Unable to create {label} genome regions bed file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    for (chrom_index, chrom_info) in chrom_list.data.iter().enumerate() {
        let chrom_regions = &genome_regions[chrom_index];
        for entry in chrom_regions.regions.find(0..chrom_info.length as i64) {
            let interval = entry.interval();
            writeln!(
                f,
                "{}\t{}\t{}",
                chrom_info.label, interval.start, interval.end
            )
            .unwrap();
        }
    }
}

/// This reads bed files into a vector of ChromRegions, indexed by chromosome id of the given
/// chrom_list
///
pub fn read_genome_regions_from_bed(
    label: &str,
    filename: &Utf8Path,
    chrom_list: &ChromList,
    quiet: bool,
    payload_required: bool,
) -> Vec<ChromRegions> {
    use rust_htslib::bgzf;
    use std::io::Read;

    if !quiet {
        info!("Reading {label} genome regions from bed file: '{filename}'");
    }

    let mut reader = unwrap!(
        bgzf::Reader::from_path(filename),
        "Unable to open {label} genome regions file: '{filename}'"
    );

    let mut content = String::new();
    unwrap!(
        reader.read_to_string(&mut content),
        "Can't parse text from {label} genome regions file: '{filename}'"
    );

    let mut chrom_regions = vec![ChromRegions::new(); chrom_list.data.len()];

    let mut last_chrom_index = 0;
    let mut last_chrom = "";
    for line in content.split('\n') {
        // The last line is expected to be empty
        if line.is_empty() {
            break;
        }

        let words = line.split('\t').collect::<Vec<_>>();

        assert!(words.len() >= 3);
        let chrom = words[0];
        let start = words[1].parse::<i64>().unwrap();
        let end = words[2].parse::<i64>().unwrap();

        if chrom != last_chrom {
            last_chrom = chrom;
            last_chrom_index = *unwrap!(
                chrom_list.label_to_index.get(chrom),
                "{label} genome regions include unknown chromosome name` `{chrom}` in file: '{filename}'"
            );
        }
        let regions = &mut chrom_regions[last_chrom_index];

        //check if we have a value to load for the track
        if words.len() >= 5 {
            let value = words[4].parse::<u8>().unwrap();
            regions.add_region_value(start, end, value);
        } else {
            assert!(!payload_required);
            regions.add_region(start, end);
        }
    }
    chrom_regions
}

/// This is a variant of GenomeRegions that allows fast lookup of chromosomes by index
#[derive(Clone)]
pub struct GenomeRegionsByChromIndex {
    pub chroms: Vec<ChromRegions>,
}

impl GenomeRegionsByChromIndex {
    /// Creates a chrom index version of an existing GenomeRegions object
    ///
    /// Any chromosomes in genome_regions not present in chrom_list will be dropped.
    ///
    #[allow(dead_code)]
    pub fn new(chrom_list: &ChromList, genome_regions: &GenomeRegions) -> Self {
        let chrom_count = chrom_list.data.len();
        let mut chroms = vec![ChromRegions::new(); chrom_count];
        for (chrom, chrom_regions) in genome_regions.chroms.iter() {
            if let Some(&chrom_index) = chrom_list.label_to_index.get(chrom) {
                chroms[chrom_index] = chrom_regions.clone();
            }
        }
        Self { chroms }
    }

    /// Return an iterator over any overlapping regions
    ///
    /// # Arguments
    /// * `chrom_index` - the contig index
    /// * `start` - the start coordinate (included)
    /// * `end` - the end coordinates (excluded)
    pub fn find_overlaps(
        &self,
        chrom_index: usize,
        start: i64,
        end: i64,
    ) -> IntervalTreeIterator<i64, u8> {
        self.chroms[chrom_index].find_overlaps(start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intersect() {
        let mut regions = ChromRegions::new();

        regions.add_region(100, 101);
        assert!(regions.intersect_pos(100));
        assert!(!regions.intersect_pos(99));
        assert!(!regions.intersect_pos(101));
    }

    /*
    #[test]
    fn test_find() {
        let mut regions = ChromRegions::new();

        regions.add_region_value(100, 200, 0);
        regions.add_region_value(200, 250, 1);
        regions.add_region_value(250, 300, 2);
        for x in regions.find_overlaps(175, 275) {
            let i = x.interval().clone();
            let s = i.end - i.start;
            eprintln!("{x:?}");
            eprintln!("{}", s);
        }
        assert!(false);
    }
    */

    #[test]
    fn test_overlap_nopanic() {
        let mut genome_regions = GenomeRegions::new(true);
        genome_regions.add_region("chr1", 10, 20);
        genome_regions.add_region("chr1", 19, 30);

        assert_eq!(
            genome_regions.find_overlaps("chr1", 9, 10).unwrap().count(),
            0
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 10, 11)
                .unwrap()
                .count(),
            1
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 18, 19)
                .unwrap()
                .count(),
            1
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 19, 20)
                .unwrap()
                .count(),
            2
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 21, 22)
                .unwrap()
                .count(),
            1
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 29, 30)
                .unwrap()
                .count(),
            1
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 30, 31)
                .unwrap()
                .count(),
            0
        );
        assert_eq!(
            genome_regions
                .find_overlaps("chr1", 0, 100)
                .unwrap()
                .count(),
            2
        );

        //let's check a non-existent chromosome also to make sure we get "None"
        assert!(genome_regions.find_overlaps("chr2", 0, 100).is_none());
    }

    #[test]
    fn test_payload() {
        let mut genome_regions = GenomeRegions::new(true);
        genome_regions.add_region_value("chr1", 10, 20, 12);
        genome_regions.add_region_value("chr1", 19, 30, 10);

        let values = genome_regions
            .find_overlaps("chr1", 10, 30)
            .unwrap()
            .map(|x| x.data())
            .copied()
            .collect::<Vec<_>>();
        assert_eq!(values[0], 12u8);
        assert_eq!(values[1], 10u8);
    }

    #[test]
    #[should_panic]
    /// Test scenarios that don't allow overlaps
    fn test_overlap_panic() {
        let mut genome_regions = GenomeRegions::new(false);
        genome_regions.add_region("chr1", 10, 20);
        genome_regions.add_region("chr1", 19, 30);
    }
}
