use std::collections::BTreeMap;
use std::path::Path;

use log::info;
use rust_htslib::bam;
use rust_vc_utils::{bigwig_utils, get_alignment_end, ArraySegmenter, ChromList};
use serde::{Deserialize, Serialize};
use unwrap::unwrap;

use crate::discover::DEPTH_BINS_MESSAGEPACK_FILENAME;
use crate::genome_regions::ChromRegions;

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

#[derive(Clone, Deserialize, Serialize)]
pub enum DepthBin {
    Excluded,
    Depth(f64),
}

impl DepthBin {
    pub fn increment(&mut self, value: f64) {
        if let DepthBin::Depth(depth) = self {
            *depth += value;
        }
    }
}

pub type ChromDepthBins = Vec<DepthBin>;

/// Manages intermediate data structures required to produce binned depth bins from bam
/// records of a single chromosome
///
pub struct DepthBinsBuilder {
    depth_edges: BTreeMap<u64, i32>,
}

impl DepthBinsBuilder {
    pub fn new() -> Self {
        Self {
            depth_edges: BTreeMap::new(),
        }
    }

    /// Merge data from another DepthBinsBuilder into this one
    pub fn merge(&mut self, other: &Self) {
        for (pos, delta) in other.depth_edges.iter() {
            *self.depth_edges.entry(*pos).or_insert(0) += delta;
        }
    }

    pub fn process_bam_record(&mut self, record: &bam::Record) {
        let begin = record.pos();
        let end = get_alignment_end(record);
        assert!(begin <= end);

        let begin = std::cmp::max(begin, 0) as u64;
        let end = std::cmp::max(end, 0) as u64;

        *self.depth_edges.entry(begin).or_insert(0) += 1;
        *self.depth_edges.entry(end).or_insert(0) -= 1;
    }

    /// Convert processed bam record data into depth bins
    ///
    pub fn get_depth_bins(
        &self,
        bin_size: u32,
        chrom_size: u64,
        excluded_regions: Option<&ChromRegions>,
    ) -> ChromDepthBins {
        /// The fraction of positions in the bin containing `pos` that are not less than `pos`
        ///
        fn get_bin_fraction(pos: u64, bin_size: u32) -> f64 {
            let bin_size = bin_size as u64;
            (bin_size - (pos % bin_size)) as f64 / bin_size as f64
        }

        let bin_count = get_complete_bin_count(chrom_size, bin_size);
        let mut depth_bins = vec![DepthBin::Depth(0f64); bin_count];

        // Mark which bins are excluded
        for (bin_index, depth_bin) in depth_bins.iter_mut().enumerate() {
            let start = bin_index as i64 * bin_size as i64;
            let end = start + bin_size as i64;
            if let Some(r) = excluded_regions {
                if r.intersect(start, end) {
                    *depth_bin = DepthBin::Excluded;
                }
            }
        }

        let mut min_unprocessed_bin = 0usize;
        let mut total_depth = 0;
        for (pos, delta) in self.depth_edges.iter() {
            let pos_bin = get_bin_index(*pos, bin_size);
            assert!((pos_bin + 1) >= min_unprocessed_bin);

            // Equality is allowed because any incomplete bin at the end of the chromosome is not
            // included in bin_count:
            assert!(pos_bin <= bin_count);

            // Add previous total_depth to all unprocessed bins not greater than pos_bin:
            //
            // Note this is done even when total_depth is zero so that min_unprocessed_bin
            // is updated.
            //
            while (min_unprocessed_bin <= pos_bin) && (min_unprocessed_bin < bin_count) {
                depth_bins[min_unprocessed_bin].increment(total_depth as f64);
                min_unprocessed_bin += 1;
            }

            // Add the fractional component of delta to pos_bin:
            if pos_bin < bin_count {
                depth_bins[pos_bin].increment(*delta as f64 * get_bin_fraction(*pos, bin_size));
            }
            total_depth += delta;
        }
        assert_eq!(total_depth, 0);

        depth_bins
    }

    /// Convert processed bam record data into depth bins
    ///
    pub fn get_regions_over_max_depth(&self, chrom_size: u64, max_depth: u32) -> ChromRegions {
        assert!(max_depth > 0);

        let mut chrom_regions = ChromRegions::new();

        let mut total_depth = 0;
        let mut start_pos = None;
        for (pos, delta) in self.depth_edges.iter() {
            total_depth += delta;
            let is_max_depth = total_depth > max_depth as i32;
            if start_pos.is_some() {
                if !is_max_depth {
                    chrom_regions.add_region(start_pos.unwrap(), *pos as i64);
                    start_pos = None;
                }
            } else if is_max_depth {
                start_pos = Some(*pos as i64);
            }
        }
        if let Some(start_pos) = start_pos {
            chrom_regions.add_region(start_pos, chrom_size as i64);
        }

        assert_eq!(total_depth, 0);

        chrom_regions
    }
}

/// Write depth bin data from one sample to a wiggle file
///
#[derive(Deserialize, Serialize)]
pub struct GenomeDepthBins {
    pub bin_size: u32,
    pub chroms: Vec<ChromDepthBins>,
}

#[allow(dead_code)]
pub fn write_depth_wig_file(
    filename: &str,
    genome_depth_bins: &GenomeDepthBins,
    chrom_list: &ChromList,
) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    info!("Writing wiggle depth track to file: '{}'", filename);

    let f = unwrap!(
        File::create(filename),
        "Unable to create wiggle depth track file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    for (chrom_index, chrom_depth_bins) in genome_depth_bins.chroms.iter().enumerate() {
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
                        writeln!(f, "{:.2}", depth).unwrap();
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
    filename: &Path,
    genome_depth_bins: &GenomeDepthBins,
    chrom_list: &ChromList,
    label: &str,
) {
    info!(
        "Writing {} track to bigwig file: '{}'",
        label,
        filename.display()
    );

    let mut bigwig_writer = bigwig_utils::get_new_writer(filename.to_str().unwrap(), chrom_list);

    for (chrom_index, chrom_depth_bins) in genome_depth_bins.chroms.iter().enumerate() {
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        for depth_bin_segment_range in
            ArraySegmenter::new(chrom_depth_bins, |v| matches!(v, DepthBin::Excluded))
        {
            let mut depths = chrom_depth_bins[depth_bin_segment_range.clone()]
                .iter()
                .map(|v| match v {
                    DepthBin::Depth(d) => *d as f32,
                    _ => panic!("Unexpected depth bin state"),
                })
                .collect::<Vec<_>>();

            bigwig_writer
                .add_interval_span_steps(
                    chrom_name,
                    depth_bin_segment_range.start as u32 * genome_depth_bins.bin_size,
                    genome_depth_bins.bin_size,
                    genome_depth_bins.bin_size,
                    &mut depths,
                )
                .unwrap();
        }
    }
}

pub fn serialize_genome_depth_bins(discover_dir: &Path, genome_depth_bins: &GenomeDepthBins) {
    let mut buf = Vec::new();
    genome_depth_bins
        .serialize(&mut rmp_serde::Serializer::new(&mut buf))
        .unwrap();

    let filename = discover_dir.join(DEPTH_BINS_MESSAGEPACK_FILENAME);

    info!(
        "Writing depth bins to binary file: '{}'",
        filename.display()
    );

    unwrap!(
        std::fs::write(&filename, buf.as_slice()),
        "Unable to open and write genome depth bins binary file: '{}'",
        filename.display()
    );
}

pub fn deserialize_genome_depth_bins(discover_dir: &Path) -> GenomeDepthBins {
    let filename = discover_dir.join(DEPTH_BINS_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read genome depth bins binary file: '{}'",
        filename.display()
    );
    rmp_serde::from_slice(&buf).unwrap()
}
