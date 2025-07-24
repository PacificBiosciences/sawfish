use std::sync::mpsc::channel;

use camino::Utf8Path;
use log::info;
use rust_vc_utils::{ArraySegmenter, ChromList, bigwig_utils, genome_ref};
use serde::{Deserialize, Serialize};
use unwrap::unwrap;

use crate::bam_scanner::SampleAlignmentScanResult;
use crate::depth_bins::{self, ChromDepthBins, DepthBin, GenomeDepthBins};
use crate::discover::{
    GC_BIAS_CORRECTED_DEPTH_BINS_BIGWIG_FILENAME, GENOME_GC_LEVELS_MESSAGEPACK_FILENAME,
    SAMPLE_GC_BIAS_MESSAGEPACK_FILENAME,
};

/// GC content information for a region
///
/// This is stored at the count level so that we have options to handle low counts downstream,
/// either with filtration or CIs
///
#[derive(Clone, Default)]
pub struct GCContent {
    pub gc: f64,
    pub at: f64,
}

impl GCContent {
    pub fn total(&self) -> f64 {
        self.gc + self.at
    }

    pub fn gc_frac(&self) -> Option<f64> {
        if self.total() < 500f64 {
            None
        } else {
            Some(self.gc / self.total())
        }
    }

    pub fn merge(&mut self, other: &Self) {
        self.gc += other.gc;
        self.at += other.at;
    }

    pub fn update_from_base(&mut self, base: &u8) {
        match *base {
            b'G' | b'C' => {
                self.gc += 1.0;
            }
            b'A' | b'T' => {
                self.at += 1.0;
            }
            _ => {}
        }
    }
}

pub type ChromGCBins = Vec<GCContent>;

#[derive(Clone)]
pub struct GenomeGCBins {
    /// `bin_size` should be the same as that used for the depth bins -- there should be a 1:1
    /// mapping from depth to gc content bins
    pub bin_size: u32,

    /// GC content for the entire genome
    ///
    /// This is used as a fallback for bins with low counts
    pub genome_bin: GCContent,

    /// Vector over chromosomes
    pub chroms: Vec<ChromGCBins>,
}

impl GenomeGCBins {
    pub fn new(bin_size: u32, chrom_count: usize) -> Self {
        Self {
            bin_size,
            genome_bin: GCContent::default(),
            chroms: vec![Vec::new(); chrom_count],
        }
    }
}

///
/// Return a 2-tuple of
/// (1) The GCContent data for all depth bins in the chromosome
/// (2) GCContent data for the entire chromosome
fn get_chrom_depth_bin_gc_content(
    chrom_ref: &[u8],
    depth_bin_size: u32,
    gc_genome_window_size: u32,
) -> (Vec<GCContent>, GCContent) {
    let chrom_size = chrom_ref.len();
    let bin_count = depth_bins::get_complete_bin_count(chrom_size as u64, depth_bin_size) as i64;

    // Local struct used to enumerate and sort all boundary points
    struct ChromGCSeg {
        pos: usize,
        bin_index: i64,
        is_end: bool,
    }

    // Get segmentation points from all bins in the chromosome and sort them
    //
    // Ultimately we want [[start, bin_id, is_start], [start, ...], ...] segmenting the full chrom
    let gc_segs = {
        let mut x = Vec::new();
        let leading_half_window = gc_genome_window_size as i64 / 2;
        let trailing_half_window = (gc_genome_window_size as i64) - leading_half_window;
        for bin_index in 0..bin_count {
            let center_pos = (((2 * bin_index) + 1) * depth_bin_size as i64) / 2;
            let start_pos = std::cmp::max(center_pos - leading_half_window, 0) as usize;
            let end_pos =
                std::cmp::min(center_pos + trailing_half_window, chrom_size as i64) as usize;
            x.push(ChromGCSeg {
                pos: start_pos,
                bin_index,
                is_end: false,
            });
            x.push(ChromGCSeg {
                pos: end_pos,
                bin_index,
                is_end: true,
            });
        }

        x.sort_by_key(|x| x.pos);
        x
    };

    // Compress any repeated segmentation points, and get start/end segment entries for each depth
    // bin:
    //
    let mut depth_bin_seg_range = vec![[0, 0]; bin_count as usize];
    let mut final_gc_segs = Vec::new();
    let mut last_start = 0;
    for gc_seg in gc_segs.into_iter() {
        if !final_gc_segs.is_empty() && last_start == gc_seg.pos {
            // This is a repeated start entry
        } else {
            final_gc_segs.push(gc_seg.pos);
            last_start = gc_seg.pos;
        }
        depth_bin_seg_range[gc_seg.bin_index as usize][gc_seg.is_end as usize] =
            final_gc_segs.len() - 1;
    }

    assert!(final_gc_segs.is_empty() || final_gc_segs[0] == 0);

    // Walk through chromosome sequence and produce a GCContent entry for each partition
    let fgcs_size = final_gc_segs.len();
    let mut fgcs_index_head = 0;
    let mut chrom_segment_gc_content = Vec::new();
    let mut chrom_gc_content = GCContent::default();
    for (pos, base) in chrom_ref.iter().enumerate() {
        if fgcs_index_head < fgcs_size {
            if final_gc_segs[fgcs_index_head] == pos {
                // Advance to next gcbin
                chrom_segment_gc_content.push(GCContent::default());
                fgcs_index_head += 1;
            }
            chrom_segment_gc_content[fgcs_index_head - 1].update_from_base(base);
        }
        chrom_gc_content.update_from_base(base);
    }

    // Sum segments for each depth bin
    let mut chrom_depth_bin_gc_content = Vec::new();
    for seg_range in depth_bin_seg_range {
        let mut gc_content = GCContent::default();
        for segment_gc_bin in chrom_segment_gc_content[seg_range[0]..seg_range[1]].iter() {
            gc_content.merge(segment_gc_bin);
        }
        chrom_depth_bin_gc_content.push(gc_content);
    }

    (chrom_depth_bin_gc_content, chrom_gc_content)
}

/// Return average gc content associated with each depth bin
///
/// GC content of each depth bin is found for a window size extending from the center of the bin. The window size may be
/// much larger than the size of the bin itself.
///
/// # Arguments
///
/// * `gc_genome_window_size` - Size of the window centered on each depth bin, from which that bin's GC-fraction is
///   taken. This window may be substantially different than the size of the bin itself (typically larger), so that we
///   consider the effect of the larger scale GC-content on the reads overlapping the depth bin itself.
///
pub fn get_depth_bin_gc_content(
    genome_ref: &genome_ref::GenomeRef,
    chrom_list: &ChromList,
    depth_bin_size: u32,
    gc_genome_window_size: u32,
    thread_count: usize,
) -> GenomeGCBins {
    info!("Getting depth bin GC content");

    // Check that chrom_list labels exist in the reference, and have the expected length
    for chrom_info in chrom_list.data.iter() {
        let chrom_ref = genome_ref.chroms.get(chrom_info.label.as_str());
        let chrom_ref = unwrap!(
            chrom_ref,
            "Can't find alignment file chromosome '{}' in the fasta reference genome input.",
            chrom_info.label
        );
        if chrom_ref.len() != chrom_info.length as usize {
            panic!(
                "Length of chromosome '{}' differs between alignment file header '{}' and fasta reference genome input '{}'",
                chrom_info.label,
                chrom_ref.len(),
                chrom_info.length
            );
        }
    }

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .unwrap();

    let (tx, rx) = channel();
    worker_pool.scope(move |scope| {
        for (chrom_index, chrom_info) in chrom_list.data.iter().enumerate() {
            let tx = tx.clone();
            let chrom_ref = genome_ref.chroms.get(chrom_info.label.as_str()).unwrap();
            scope.spawn(move |_| {
                let result = get_chrom_depth_bin_gc_content(
                    chrom_ref,
                    depth_bin_size,
                    gc_genome_window_size,
                );

                tx.send((chrom_index, result)).unwrap();
            });
        }
    });

    let chrom_count = chrom_list.data.len();
    let mut depth_bin_gc_content = GenomeGCBins::new(depth_bin_size, chrom_count);

    for (chrom_index, (chrom_depth_bin_gc_content, chrom_gc_content)) in rx {
        depth_bin_gc_content.chroms[chrom_index] = chrom_depth_bin_gc_content;
        depth_bin_gc_content.genome_bin.merge(&chrom_gc_content);
    }

    depth_bin_gc_content
}

/// Write gc fraction per depth bin to a bigwig file
///
pub fn write_gc_fraction_per_bin_bigwig_file(
    filename: &Utf8Path,
    depth_bin_gc_content: &GenomeGCBins,
    chrom_list: &ChromList,
) {
    info!("Writing gc fraction per bin to bigwig file: '{filename}'");

    let mut bigwig_writer = bigwig_utils::get_new_writer(filename.as_str(), chrom_list);

    for (chrom_index, chrom_gc_bins) in depth_bin_gc_content.chroms.iter().enumerate() {
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        for gc_bin_segment_range in ArraySegmenter::new(chrom_gc_bins, |v| v.gc_frac().is_none()) {
            let mut chrom_gc_bin_segment = chrom_gc_bins[gc_bin_segment_range.clone()]
                .iter()
                .map(|v| match v.gc_frac() {
                    Some(gc) => gc as f32,
                    None => panic!("Unexpected GC fraction"),
                })
                .collect::<Vec<_>>();

            bigwig_writer
                .add_interval_span_steps(
                    chrom_name,
                    gc_bin_segment_range.start as u32 * depth_bin_gc_content.bin_size,
                    depth_bin_gc_content.bin_size,
                    depth_bin_gc_content.bin_size,
                    &mut chrom_gc_bin_segment,
                )
                .unwrap();
        }
    }
}

#[derive(Clone)]
struct GCLevelBuilder {
    depth: f64,
    bins: u64,
}

impl GCLevelBuilder {
    fn new() -> Self {
        Self {
            depth: 0.0,
            bins: 0,
        }
    }
}

#[derive(Deserialize, Serialize)]
pub struct SampleGCBiasCorrectionData {
    pub gc_depth_reduction: Vec<f64>,
}

impl SampleGCBiasCorrectionData {
    fn new(gc_level_count: usize) -> Self {
        Self {
            gc_depth_reduction: vec![0.0; gc_level_count],
        }
    }

    #[allow(unused)]
    pub fn to_gc_depth_correction(&self) -> Vec<f64> {
        self.gc_depth_reduction
            .iter()
            .map(|x| 1.0 / x)
            .collect::<Vec<_>>()
    }
}

pub type GenomeGCLevels = Vec<Vec<usize>>;

pub struct GCBiasCorrectionData {
    /// GC level assignment vector for each chromosome, there is one gc level for each depth bin.
    ///
    /// This result depends on reference genome but not sample data
    pub genome_gc_levels: GenomeGCLevels,

    /// The gc-depth table used for GC bias correction in the sample
    pub sample_gc_bias_data: SampleGCBiasCorrectionData,
}

/// Convert GC fraction into a binned frequency level
///
fn gc_fraction_to_level(gc_fraction: f64, gc_level_count: usize) -> usize {
    std::cmp::max(
        0,
        std::cmp::min(
            gc_level_count - 1,
            (gc_fraction * gc_level_count as f64) as usize,
        ),
    )
}

/// Convert the GC and AT counts in each GC-bin into discrete GC-levels, over the whole genome
///
/// These values will be used for GC correction in all samples
///
fn get_genome_gc_levels(gc_bins: &GenomeGCBins, gc_level_count: usize) -> GenomeGCLevels {
    let genome_gc_frac = gc_bins.genome_bin.gc_frac().unwrap();

    let mut genome_gc_levels = Vec::new();
    for chrom_gc_bins in gc_bins.chroms.iter() {
        let mut chrom_gc_levels = Vec::new();
        for gc_bin in chrom_gc_bins.iter() {
            let gc_level =
                gc_fraction_to_level(gc_bin.gc_frac().unwrap_or(genome_gc_frac), gc_level_count);
            chrom_gc_levels.push(gc_level);
        }
        genome_gc_levels.push(chrom_gc_levels);
    }
    genome_gc_levels
}

/// Estimate GC-coverage relationship for sample
///
fn get_sample_gc_correction(
    chrom_list: &ChromList,
    genome_gc_levels: &GenomeGCLevels,
    sample_scan_result: &SampleAlignmentScanResult,
    coverage_est_regex: &str,
    gc_level_count: usize,
) -> SampleGCBiasCorrectionData {
    use regex::Regex;

    let chrom_include_regex = Regex::new(coverage_est_regex).unwrap();

    let mut gc_levels = vec![GCLevelBuilder::new(); gc_level_count];

    // Walk through depth and gc level bins together:
    for (chrom_index, (chrom_depth_bins, chrom_gc_levels)) in sample_scan_result
        .genome_depth_bins
        .depth_bins
        .iter()
        .zip(genome_gc_levels.iter())
        .enumerate()
    {
        if !chrom_include_regex.is_match(chrom_list.data[chrom_index].label.as_str()) {
            continue;
        }
        for (depth_bin, gc_level) in chrom_depth_bins.iter().zip(chrom_gc_levels.iter()) {
            if let depth_bins::DepthBin::Depth(depth) = depth_bin {
                gc_levels[*gc_level].depth += depth;
                gc_levels[*gc_level].bins += 1;
            }
        }
    }

    let debug = false;
    if debug {
        eprintln!("get_sample_gc_correction: intermediate debug output of depth totals:");
        eprintln!("gc_index\tgc_bins\tgc_depth\tgc_bin_depth");
        for (gc_index, gc_level) in gc_levels.iter().enumerate() {
            eprintln!(
                "{}\t{}\t{}\t{}",
                gc_index,
                gc_level.bins,
                gc_level.depth,
                gc_level.depth / gc_level.bins as f64
            );
        }
    }

    // Now convert depth totals to the depth correction factors that will be used downstream

    let mut sample_data = SampleGCBiasCorrectionData::new(gc_level_count);

    let min_base_coverage = 2_000_000;
    let min_bin_count = min_base_coverage / sample_scan_result.genome_depth_bins.bin_size as u64;

    // First find the max bin, or if no bins qualify just disable gc correction by setting all
    // values to 1.0
    //
    let mut max_depth = 0.0;
    let mut max_gc_index = 0;
    for (gc_index, gc_level) in gc_levels.iter().enumerate() {
        if gc_level.bins < min_bin_count {
            continue;
        }
        let depth = gc_level.depth / gc_level.bins as f64;
        sample_data.gc_depth_reduction[gc_index] = depth;
        if depth > max_depth {
            max_depth = depth;
            max_gc_index = gc_index;
        }
    }
    let max_depth = max_depth;
    let max_gc_index = max_gc_index;

    if max_depth > 0.0 {
        // Normalize by max depth:
        for depth in sample_data.gc_depth_reduction.iter_mut() {
            *depth /= max_depth;
        }

        // Fill in any unsupported high-gc bins with nearest supported value
        let mut last = 1.0;
        let mut is_fill = false;
        for gc_index in max_gc_index + 1..gc_level_count {
            if !is_fill {
                if sample_data.gc_depth_reduction[gc_index] > 0.0 {
                    last = sample_data.gc_depth_reduction[gc_index];
                } else {
                    is_fill = true;
                }
            }
            if is_fill {
                sample_data.gc_depth_reduction[gc_index] = last;
            }
        }

        // Fill in any unsupported low-gc bins with nearest supported value
        let mut last = 1.0;
        let mut is_fill = false;
        for gc_index in (0..max_gc_index).rev() {
            if !is_fill {
                if sample_data.gc_depth_reduction[gc_index] > 0.0 {
                    last = sample_data.gc_depth_reduction[gc_index];
                } else {
                    is_fill = true;
                }
            }
            if is_fill {
                sample_data.gc_depth_reduction[gc_index] = last;
            }
        }
    } else {
        // Disable gc correction:
        for depth in sample_data.gc_depth_reduction.iter_mut() {
            *depth = 1.0;
        }
    }

    sample_data
}

/// Estimate GC-coverage relationship for each sample and use this to compute GC correction factors
///
/// # Arguments
///
/// * `depth_bin_gc_content` - Count of GC and AT reference bases in a window centered on each depth bin
///
/// * `coverage_est_regex` - Regex used to select chromosomes for mean haploid coverage estimation. All selected
///   chromosomes are assumed diploid.
///
/// * `gc_level_count` - Number of equal divisions of the GC frequency range to use in the GC correction process
///
pub fn get_gc_correction(
    chrom_list: &ChromList,
    depth_bin_gc_content: &GenomeGCBins,
    sample_scan_result: &SampleAlignmentScanResult,
    coverage_est_regex: &str,
    gc_level_count: usize,
) -> GCBiasCorrectionData {
    info!("Getting GC depth correction levels");

    let genome_gc_levels = get_genome_gc_levels(depth_bin_gc_content, gc_level_count);
    let sample_gc_bias_data = get_sample_gc_correction(
        chrom_list,
        &genome_gc_levels,
        sample_scan_result,
        coverage_est_regex,
        gc_level_count,
    );
    GCBiasCorrectionData {
        genome_gc_levels,
        sample_gc_bias_data,
    }
}

/// Write a tsv file of gc depth reduction factors for one sample
///
/// This version of the GC correction table is intended to be simple to work with in methods
/// debugging scenarios.
///
fn write_gc_depth_reduction_debug_output(output_dir: &Utf8Path, gc_depth_reduction: &[f64]) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let filename = output_dir.join("gc_correction_table.tsv");

    info!("Writing gc correction table to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create gc correction table file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    let gc_levels = gc_depth_reduction.len();
    let gc_inc = 1.0 / gc_levels as f64;
    writeln!(f, "index\tmin_gc\tmax_gc\tgc_depth_factor").unwrap();
    for (index, gcr) in gc_depth_reduction.iter().enumerate() {
        writeln!(
            f,
            "{}\t{:.3}\t{:.3}\t{:.3}",
            index,
            index as f64 * gc_inc,
            (index + 1) as f64 * gc_inc,
            gcr
        )
        .unwrap();
    }
}

/// Given sample depth and gc correction information, return the gc-scaled depth
///
/// Note that scaled depth is not used by the segmentation model, this is only provided as a
/// debugging aid for humans.
///
fn gc_scale_sample_depth(
    depth_bins: &[ChromDepthBins],
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
    genome_gc_levels: &GenomeGCLevels,
) -> Vec<ChromDepthBins> {
    // This application requires the reciprocal of the depth reduction factor used for the
    // emission prob correction, so go ahead and setup the reciprocal value array here:
    let gc_depth_correction = sample_gc_bias_data
        .gc_depth_reduction
        .iter()
        .map(|x| 1.0 / x)
        .collect::<Vec<_>>();

    let mut gc_scaled_depth_bins = Vec::new();
    for (chrom_index, chrom_depth_bins) in depth_bins.iter().enumerate() {
        let chrom_gc_levels = &genome_gc_levels[chrom_index];

        let bin_count = chrom_depth_bins.len();
        let mut gc_scaled_chrom_depth_bins = Vec::with_capacity(bin_count);

        for bin_index in 0..bin_count {
            let depth = &chrom_depth_bins[bin_index];
            let gc_level = chrom_gc_levels[bin_index];

            gc_scaled_chrom_depth_bins.push(match depth {
                DepthBin::Depth(x) => DepthBin::Depth(x * gc_depth_correction[gc_level]),
                DepthBin::Excluded => DepthBin::Excluded,
            });
        }

        gc_scaled_depth_bins.push(gc_scaled_chrom_depth_bins);
    }

    gc_scaled_depth_bins
}

/// Scale depth by gc correction factors and write these scaled depths out to bigwig format
///
pub fn write_gc_scaled_depth_track_files(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
    genome_gc_levels: &GenomeGCLevels,
    genome_depth_bins: &GenomeDepthBins,
) {
    let gc_scaled_depth_filename = output_dir.join(GC_BIAS_CORRECTED_DEPTH_BINS_BIGWIG_FILENAME);

    let gc_scaled_depth_bins = gc_scale_sample_depth(
        &genome_depth_bins.depth_bins,
        sample_gc_bias_data,
        genome_gc_levels,
    );

    depth_bins::write_depth_bigwig_file(
        &gc_scaled_depth_filename,
        genome_depth_bins.bin_size,
        &gc_scaled_depth_bins,
        chrom_list,
        "gc scaled depth",
    );
}

/// This is a quick debug dump of gc reduction factors
///
/// To save time, result is written to depth bins even though these are not depths
///
fn get_gc_reduction_factor_track(
    depth_bins: &[ChromDepthBins],
    sample_gc_bias_data: &SampleGCBiasCorrectionData,
    gc_levels: &GenomeGCLevels,
) -> Vec<ChromDepthBins> {
    let mut rc_reduction_factor = Vec::new();

    for (chrom_index, chrom_depth_bins) in depth_bins.iter().enumerate() {
        let chrom_gc_levels = &gc_levels[chrom_index];

        let bin_count = chrom_depth_bins.len();
        let mut gc_scaled_chrom_depth_bins = Vec::with_capacity(bin_count);

        for bin_index in 0..bin_count {
            let depth = &chrom_depth_bins[bin_index];
            let gc_level = chrom_gc_levels[bin_index];

            gc_scaled_chrom_depth_bins.push(match depth {
                DepthBin::Depth(_) => {
                    DepthBin::Depth(sample_gc_bias_data.gc_depth_reduction[gc_level])
                }
                DepthBin::Excluded => DepthBin::Excluded,
            });
        }

        rc_reduction_factor.push(gc_scaled_chrom_depth_bins);
    }

    rc_reduction_factor
}

/// Write out gc reduction factors for each depth bin, over all samples
///
fn write_gc_reduction_track_files(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    gc_bias_data: &GCBiasCorrectionData,
    sample_scan_result: &SampleAlignmentScanResult,
) {
    let sample_gc_bias_data = &gc_bias_data.sample_gc_bias_data;
    let gc_scaled_depth_filename = output_dir.join("gc_reduction_factor.bw");

    let gc_reduction_factor = get_gc_reduction_factor_track(
        &sample_scan_result.genome_depth_bins.depth_bins,
        sample_gc_bias_data,
        &gc_bias_data.genome_gc_levels,
    );

    depth_bins::write_depth_bigwig_file(
        &gc_scaled_depth_filename,
        sample_scan_result.genome_depth_bins.bin_size,
        &gc_reduction_factor,
        chrom_list,
        "gc reduction factor",
    );
}

/// Write several debug output tracks/files related to the gc corrections process
///
/// Debug outputs include:
/// 1. A GC-fraction bigwig track (derived from reference so shared by all samples)
/// 2. GC-corrected depth bigwig track, one per sample
/// 3. GC-correction factor bigwig track, one per sample
/// 4. Tsv file of GC-correction factors as a function of gc-content, one per sample
///
pub fn write_gc_correction_debug_output(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    depth_bin_gc_content: &GenomeGCBins,
    sample_scan_result: &SampleAlignmentScanResult,
    gc_bias_data: &GCBiasCorrectionData,
) {
    let gc_fraction_filename = output_dir.join("gc_frac.bw");
    write_gc_fraction_per_bin_bigwig_file(&gc_fraction_filename, depth_bin_gc_content, chrom_list);

    write_gc_depth_reduction_debug_output(
        output_dir,
        &gc_bias_data.sample_gc_bias_data.gc_depth_reduction,
    );

    write_gc_scaled_depth_track_files(
        output_dir,
        chrom_list,
        &gc_bias_data.sample_gc_bias_data,
        &gc_bias_data.genome_gc_levels,
        &sample_scan_result.genome_depth_bins,
    );

    write_gc_reduction_track_files(output_dir, chrom_list, gc_bias_data, sample_scan_result);
}

/// Serialize the two parts of GC bias correction table into two separate files
pub fn serialize_gc_bias_correction_data(
    discover_dir: &Utf8Path,
    gc_bias_data: &GCBiasCorrectionData,
) {
    // serialize gc levels
    {
        let mut buf = Vec::new();
        gc_bias_data
            .genome_gc_levels
            .serialize(&mut rmp_serde::Serializer::new(&mut buf))
            .unwrap();

        let filename = discover_dir.join(GENOME_GC_LEVELS_MESSAGEPACK_FILENAME);
        info!("Writing genome gc levels to binary file: '{filename}'");

        unwrap!(
            std::fs::write(&filename, buf.as_slice()),
            "Unable to open and write genome gc levels to binary file: '{filename}'"
        );
    }

    // serialize sample gc bias data
    {
        let mut buf = Vec::new();
        gc_bias_data
            .sample_gc_bias_data
            .serialize(&mut rmp_serde::Serializer::new(&mut buf))
            .unwrap();

        let filename = discover_dir.join(SAMPLE_GC_BIAS_MESSAGEPACK_FILENAME);
        info!("Writing sample gc bias data to binary file: '{filename}'");

        unwrap!(
            std::fs::write(&filename, buf.as_slice()),
            "Unable to open and write sample gc bias data to binary file: '{filename}'"
        );
    }
}

pub fn deserialize_genome_gc_levels(discover_dir: &Utf8Path) -> GenomeGCLevels {
    let filename = discover_dir.join(GENOME_GC_LEVELS_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read genome gc levels binary file: '{filename}'"
    );
    unwrap!(
        rmp_serde::from_slice(&buf),
        "Unable to parse genome gc levels binary file: '{filename}'"
    )
}

pub fn deserialize_sample_gc_bias_data(discover_dir: &Utf8Path) -> SampleGCBiasCorrectionData {
    let filename = discover_dir.join(SAMPLE_GC_BIAS_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read sample gc bias binary file: '{filename}'"
    );
    rmp_serde::from_slice(&buf).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_vc_utils::GenomeRef;
    use std::collections::HashMap;

    #[test]
    fn test_get_depth_bin_gc_content() {
        let mut genome_ref = GenomeRef {
            chroms: HashMap::new(),
        };
        genome_ref.chroms.insert(
            "foo".to_string(),
            b"AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT".to_vec(),
        );
        genome_ref.chroms.insert(
            "bar".to_string(),
            b"CCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAA".to_vec(),
        );

        let mut chrom_list = ChromList::default();
        chrom_list.add_chrom("foo", 40);
        chrom_list.add_chrom("bar", 38);

        let depth_bin_size = 5;
        let result = get_depth_bin_gc_content(&genome_ref, &chrom_list, depth_bin_size, 9, 1);
        assert_eq!(result.bin_size, depth_bin_size);
        approx::assert_ulps_eq!(result.genome_bin.at, 38.0, max_ulps = 4);
        approx::assert_ulps_eq!(result.genome_bin.gc, 40.0, max_ulps = 4);
        assert_eq!(result.chroms.len(), 2);
        let test_bin = &result.chroms[0][0];
        approx::assert_ulps_eq!(test_bin.at, 7.0);
        approx::assert_ulps_eq!(test_bin.gc, 0.0);
        let test_bin = &result.chroms[0][2];
        approx::assert_ulps_eq!(test_bin.at, 2.0);
        approx::assert_ulps_eq!(test_bin.gc, 7.0);
        let test_bin = &result.chroms[0][7];
        approx::assert_ulps_eq!(test_bin.at, 7.0);
        approx::assert_ulps_eq!(test_bin.gc, 0.0);
    }

    #[test]
    fn test_gc_fraction_to_level() {
        assert_eq!(gc_fraction_to_level(-0.1, 2), 0);
        assert_eq!(gc_fraction_to_level(0.0, 2), 0);
        assert_eq!(gc_fraction_to_level(0.1, 2), 0);
        assert_eq!(gc_fraction_to_level(0.49, 2), 0);
        assert_eq!(gc_fraction_to_level(0.51, 2), 1);
        assert_eq!(gc_fraction_to_level(0.9, 2), 1);
        assert_eq!(gc_fraction_to_level(1.0, 2), 1);
        assert_eq!(gc_fraction_to_level(1.1, 2), 1);
    }
}
