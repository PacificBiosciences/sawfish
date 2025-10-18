use std::collections::BTreeMap;

use camino::Utf8Path;
use log::info;
use rust_htslib::bam;
use rust_vc_utils::{ArraySegmenter, ChromList, MeanTracker, RegionMap, bigwig_utils};
use serde::{Deserialize, Serialize};
use unwrap::unwrap;

use crate::bam_utils::get_alignment_ref_segments;
use crate::filenames::DEPTH_BINS_MESSAGEPACK_FILENAME;
use crate::genome_regions::GenomeRegions;

/// Return the number of complete bins of size `bin_size` in `total_size`
///
/// Any incomplete bin at the end of the chromosome is excluded
///
pub fn get_complete_bin_count(total_size: u64, bin_size: u32) -> usize {
    (total_size / bin_size as u64) as usize
}

/// Return the zero-indexed bin number of position `pos` given bins of size `bin_size`
///
pub fn get_bin_index(pos: u64, bin_size: u32) -> usize {
    (pos / bin_size as u64) as usize
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum DepthBin {
    Excluded,
    Depth(f64),
}

pub type ChromDepthBins = Vec<DepthBin>;

/// All depth bin and related summary info from read scanning
#[derive(Clone, Default)]
pub struct AllChromDepthBinInfo {
    pub depth_bins: ChromDepthBins,
    pub read_length: MeanTracker,
}

/// Manages all temporary data structures required to produce depth tracks for a single chromosome
///
#[derive(Clone, Default)]
pub struct DepthBinsBuilder {
    /// Intermediate structure used to build standrd depth bins, which are used for segmetnation and
    /// most intuitive for visualization outputs for the user.
    ///
    /// Key is chrom position, value is depth delta
    ///
    depth_edges: BTreeMap<u64, i32>,

    read_length: MeanTracker,
}

impl DepthBinsBuilder {
    /// Merge data from another DepthBinsBuilder into this one
    pub fn merge(&mut self, other: &Self) {
        for (pos, delta) in other.depth_edges.iter() {
            *self.depth_edges.entry(*pos).or_insert(0) += delta;
        }
    }

    /// Process depth bin information from one bam record
    ///
    /// # Arguments:
    /// * min_del_size - The smallest alignment CIGAR deletion size that will be accounted for while building the depth
    ///   track (split reads are always accounted for)
    ///
    pub fn process_bam_record(&mut self, min_del_size: u32, record: &bam::Record) {
        // Update depth edges
        {
            let ref_segs = get_alignment_ref_segments(record.pos(), &record.cigar(), min_del_size);
            for ref_seg in ref_segs {
                assert!(ref_seg.start <= ref_seg.end);

                let start = std::cmp::max(ref_seg.start, 0) as u64;
                let end = std::cmp::max(ref_seg.end, 0) as u64;

                *self.depth_edges.entry(start).or_insert(0) += 1;
                *self.depth_edges.entry(end).or_insert(0) -= 1;
            }
        }

        // Update read_length tracker
        //
        // Ensure that each split read is only counted once
        if !record.is_supplementary() {
            self.read_length.insert(record.seq_len() as f64);
        }
    }

    /// Convert processed bam record data into depth bins
    ///
    pub fn get_depth_bins(&self, bin_size: u32, chrom_size: u64) -> AllChromDepthBinInfo {
        /// The fraction of positions in the bin containing `pos` that are not less than `pos`
        ///
        fn get_bin_fraction(pos: u64, bin_size: u32) -> f64 {
            let bin_size = bin_size as u64;
            (bin_size - (pos % bin_size)) as f64 / bin_size as f64
        }

        let bin_count = get_complete_bin_count(chrom_size, bin_size);
        let mut depth_bins = vec![0f64; bin_count];

        // Index of the lowest bin number that hasn't been processed yet:
        let mut min_unprocessed_bin = 0usize;
        let mut total_depth = 0;

        for (&pos, &delta) in self.depth_edges.iter() {
            let pos_bin = get_bin_index(pos, bin_size);

            assert!((pos_bin + 1) >= min_unprocessed_bin);

            // Equality is allowed because bin_count does not include any incomplete bin at the end of the chromosome.
            // For instance if bin_size is 1000, a chromosome 1999 bases long will have a bin_count of 1.
            //
            assert!(pos_bin <= bin_count);

            // Add previous total_depth to all unprocessed bins not greater than pos_bin:
            //
            // Note this is done even when total_depth is zero so that min_unprocessed_bin
            // is updated.
            //
            while (min_unprocessed_bin <= pos_bin) && (min_unprocessed_bin < bin_count) {
                depth_bins[min_unprocessed_bin] += total_depth as f64;
                min_unprocessed_bin += 1;
            }

            if pos_bin < bin_count {
                // Add the fractional component of delta to depth pos_bin:
                depth_bins[pos_bin] += delta as f64 * get_bin_fraction(pos, bin_size);
            }
            total_depth += delta;
            assert!(total_depth >= 0);
        }
        assert_eq!(total_depth, 0);

        let depth_bins = depth_bins.into_iter().map(DepthBin::Depth).collect();

        AllChromDepthBinInfo {
            depth_bins,
            read_length: self.read_length.clone(),
        }
    }

    /// Convert processed bam record data into depth bins
    ///
    pub fn get_regions_over_max_depth(&self, chrom_size: u64, max_depth: u32) -> RegionMap {
        assert!(max_depth > 0);

        let mut chrom_region_map = RegionMap::default();

        let mut total_depth = 0;
        let mut start_pos = None;
        for (pos, delta) in self.depth_edges.iter() {
            total_depth += delta;
            let is_max_depth = total_depth > max_depth as i32;
            if start_pos.is_some() {
                if !is_max_depth {
                    chrom_region_map.add_region(start_pos.unwrap(), *pos as i64);
                    start_pos = None;
                }
            } else if is_max_depth {
                start_pos = Some(*pos as i64);
            }
        }
        if let Some(start_pos) = start_pos {
            chrom_region_map.add_region(start_pos, chrom_size as i64);
        }

        assert_eq!(total_depth, 0);

        chrom_region_map
    }
}

#[derive(Deserialize, Serialize)]
pub struct GenomeDepthBins {
    pub bin_size: u32,
    pub depth_bins: Vec<ChromDepthBins>,
    pub mean_read_length: f64,
}

/// Mark excluded regions in depth bins track
///
pub fn apply_excluded_regions(
    chrom_list: &ChromList,
    excluded_regions: Option<&GenomeRegions>,
    bin_size: u32,
    depth_bins: &mut [ChromDepthBins],
) {
    if let Some(excluded_regions) = excluded_regions {
        for (chrom_index, chrom_depth_bins) in depth_bins.iter_mut().enumerate() {
            let chrom_label = &chrom_list.data[chrom_index].label;
            let chrom_excluded_regions = excluded_regions.chroms.get(chrom_label);

            if let Some(chrom_excluded_regions) = chrom_excluded_regions {
                for (bin_index, depth_bin) in chrom_depth_bins.iter_mut().enumerate() {
                    let start = bin_index as i64 * bin_size as i64;
                    let end = start + bin_size as i64;
                    if chrom_excluded_regions.intersect(start, end) {
                        *depth_bin = DepthBin::Excluded;
                    }
                }
            }
        }
    }
}

/// Write depth bin data from one sample to a wiggle file
///
#[allow(dead_code)]
pub fn write_depth_wig_file(
    filename: &str,
    genome_depth_bins: &GenomeDepthBins,
    chrom_list: &ChromList,
) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    info!("Writing wiggle depth track to file: '{filename}'");

    let f = unwrap!(
        File::create(filename),
        "Unable to create wiggle depth track file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    for (chrom_index, chrom_depth_bins) in genome_depth_bins.depth_bins.iter().enumerate() {
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        for depth_bin_segment_range in
            ArraySegmenter::new(chrom_depth_bins, |v| matches!(v, DepthBin::Excluded))
        {
            writeln!(
                f,
                "fixedStep chrom={} start={} step={bin_size} span={bin_size}",
                chrom_name,
                (depth_bin_segment_range.start as u32 * genome_depth_bins.bin_size) + 1,
                bin_size = genome_depth_bins.bin_size,
            )
            .unwrap();

            for bin_index in depth_bin_segment_range.clone() {
                match chrom_depth_bins[bin_index] {
                    DepthBin::Depth(depth) => {
                        writeln!(f, "{depth:.2}").unwrap();
                    }
                    DepthBin::Excluded => {
                        panic!("Unexpected depth bin state");
                    }
                }
            }
        }
    }
}

/// Write genome_depth_bins from one sample to a bigwig depth track file
///
pub fn write_depth_bigwig_file(
    filename: &Utf8Path,
    bin_size: u32,
    depth_bins: &[ChromDepthBins],
    chrom_list: &ChromList,
    label: &str,
) {
    info!("Writing {label} track to bigwig file: '{filename}'");

    let mut bigwig_writer = bigwig_utils::get_new_writer(filename.as_str(), chrom_list);

    for (chrom_index, chrom_depth_bins) in depth_bins.iter().enumerate() {
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        for depth_bin_segment_range in
            ArraySegmenter::new(chrom_depth_bins, |v| matches!(v, DepthBin::Excluded))
        {
            let mut depths = chrom_depth_bins[depth_bin_segment_range.clone()]
                .iter()
                .map(|v| match v {
                    DepthBin::Depth(d) => *d as f32,
                    _ => panic!("Unexpected depth bin state while writing {label} track"),
                })
                .collect::<Vec<_>>();

            bigwig_writer
                .add_interval_span_steps(
                    chrom_name,
                    depth_bin_segment_range.start as u32 * bin_size,
                    bin_size,
                    bin_size,
                    &mut depths,
                )
                .unwrap();
        }
    }
}

pub fn serialize_genome_depth_bins(discover_dir: &Utf8Path, genome_depth_bins: &GenomeDepthBins) {
    let mut buf = Vec::new();
    genome_depth_bins
        .serialize(&mut rmp_serde::Serializer::new(&mut buf))
        .unwrap();

    let filename = discover_dir.join(DEPTH_BINS_MESSAGEPACK_FILENAME);

    info!("Writing depth bins to binary file: '{filename}'");

    unwrap!(
        std::fs::write(&filename, buf.as_slice()),
        "Unable to open and write genome depth bins binary file: '{filename}'"
    );
}

pub fn deserialize_genome_depth_bins(discover_dir: &Utf8Path) -> GenomeDepthBins {
    let filename = discover_dir.join(DEPTH_BINS_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read genome depth bins binary file: '{filename}'"
    );
    rmp_serde::from_slice(&buf).unwrap()
}
