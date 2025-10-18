use std::sync::mpsc::channel;

use log::{error, info};
use rust_htslib::bam::{self, Read};
use rust_vc_utils::{
    ChromList, GenomeRef, MeanTracker, ProgressReporter, RegionMap, filter_out_alignment_record,
    get_region_segments,
};
use unwrap::unwrap;

use crate::bam_utils;
use crate::cli;
use crate::cluster_breakpoints::{BreakBuilder, BreakObservations};
use crate::depth_bins::{AllChromDepthBinInfo, ChromDepthBins, DepthBinsBuilder, GenomeDepthBins};
use crate::genome_regions::GenomeRegions;
use crate::worker_thread_data::BamReaderWorkerThreadDataSet;

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segment(
    bam_reader: &mut bam::IndexedReader,
    scan_settings: &ScanSettings,
    chrom_list: &ChromList,
    chrom_index: usize,
    begin: i64,
    end: i64,
    chrom_sequence: &[u8],
    is_targeted_scan: bool,
    progress_reporter: &ProgressReporter,
) -> (BreakObservations, DepthBinsBuilder) {
    let mut depth_builder = DepthBinsBuilder::default();
    let mut break_builder = BreakBuilder::new(
        scan_settings.min_evidence_indel_size,
        scan_settings.min_sv_mapq,
        scan_settings.max_close_breakend_distance,
        begin,
        end,
        chrom_list,
        is_targeted_scan,
    );

    bam_reader
        .fetch(bam::FetchDefinition::Region(chrom_index as i32, begin, end))
        .unwrap();

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        if filter_out_alignment_record(&record) {
            continue;
        }

        if bam_utils::get_gap_compressed_identity(&record, chrom_sequence)
            < scan_settings.min_gap_compressed_identity
        {
            continue;
        }

        break_builder.process_bam_record(&record);

        // This filter keeps read evidence from being double-counted at region boundaries
        if record.pos() >= begin {
            depth_builder.process_bam_record(scan_settings.depth_bin_size, &record);
        }
    }

    let result = (break_builder.complete_processing(), depth_builder);

    let kb_completed = (end - begin) as u64 / 1000;
    progress_reporter.inc(kb_completed);

    result
}

/// Extend get_region_segments for use with intervals that start at non-zero positions by shifting
/// the results according to the given offset value.
pub fn get_region_segments_with_offset(
    offset: u64,
    size: u64,
    segment_size: u64,
) -> Vec<(u64, u64)> {
    get_region_segments(size, segment_size)
        .into_iter()
        .map(|(s, e)| (s + offset, e + offset))
        .collect()
}

/// Get the subsegments of the chromosome to distribute as individual worker thread jobs
///
/// Subsegments are either found from the full chromosome, or can be found for subsets of the
/// chromosome (typically in a debugging context)
///
fn get_chrom_worker_segments(
    chrom_size: u64,
    chrom_target_regions: Option<&RegionMap>,
    segment_size: u64,
) -> Vec<(u64, u64)> {
    match chrom_target_regions {
        Some(regions) => {
            let mut all_region_segments = Vec::new();
            for region in regions.find_overlaps(0, chrom_size as i64) {
                let interval = region.interval();
                let size = (interval.end - interval.start) as u64;
                let region_segments =
                    get_region_segments_with_offset(interval.start as u64, size, segment_size);
                all_region_segments.extend(region_segments);
            }
            all_region_segments
        }
        None => get_region_segments(chrom_size, segment_size),
    }
}

/// Scan the given chromosome for SV evidence
///
/// Delegate out segments of the chromosome to worker threads which scan each segment, then join
/// all results into chromosome scale data structures.
///
/// Returns a 3-tuple of:
/// 1. Breakpoint observations
/// 2. Binned depth data
/// 3. Regions exceeding max SV depth
///
#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segments(
    worker_thread_dataset: BamReaderWorkerThreadDataSet,
    scan_settings: &ScanSettings,
    chrom_list: &ChromList,
    chrom_index: usize,
    chrom_size: u64,
    chrom_sequence: &[u8],
    chrom_target_regions: Option<&RegionMap>,
    progress_reporter: &ProgressReporter,
) -> (BreakObservations, AllChromDepthBinInfo, DepthBinsBuilder) {
    let (tx, rx) = channel();

    let is_targeted_scan = chrom_target_regions.is_some();
    let chrom_segments =
        get_chrom_worker_segments(chrom_size, chrom_target_regions, scan_settings.segment_size);

    rayon::scope(move |scope| {
        for (begin, end) in chrom_segments {
            let worker_thread_dataset = worker_thread_dataset.clone();
            let tx = tx.clone();

            scope.spawn(move |_| {
                let worker_id = rayon::current_thread_index().unwrap();

                // Hardcode for single sample in discover mode:
                let sample_index = 0;
                let bam_reader =
                    &mut worker_thread_dataset[worker_id].lock().unwrap().bam_readers[sample_index];

                let result = scan_chromosome_segment(
                    bam_reader,
                    scan_settings,
                    chrom_list,
                    chrom_index,
                    begin as i64,
                    end as i64,
                    chrom_sequence,
                    is_targeted_scan,
                    progress_reporter,
                );
                tx.send(result).unwrap();
            });
        }
    });

    let (mut chrom_break_observations, mut depth_builder) = rx.iter().next().unwrap();
    for (mut other_break_observations, other_depth_builder) in rx {
        chrom_break_observations.merge(&mut other_break_observations);
        depth_builder.merge(&other_depth_builder);
    }

    let chrom_depth_bins = depth_builder.get_depth_bins(scan_settings.depth_bin_size, chrom_size);

    (chrom_break_observations, chrom_depth_bins, depth_builder)
}

pub struct SampleAlignmentScanResult {
    pub sample_name: String,

    pub genome_break_observations: BreakObservations,

    /// A vector with depth bins for each chromosome in the genome
    pub genome_depth_bins: GenomeDepthBins,

    pub genome_depth_bins_builder: Vec<DepthBinsBuilder>,
}

struct ScanSettings {
    depth_bin_size: u32,
    min_gap_compressed_identity: f64,
    min_evidence_indel_size: u32,
    min_sv_mapq: u32,
    max_close_breakend_distance: usize,

    /// This defines how large of a chromosome segment should be processed by a single thread
    segment_size: u64,
}

/// Read single alignment file and translate into raw SV evidence and depth data structures
///
/// # Arguments
///
/// * `target_regions` - Optional regions for targetted genome scanning.(supported for debug only)
///
/// Returns all scanned sample data
///
#[allow(clippy::too_many_arguments)]
pub fn scan_sample_bam_for_sv_evidence(
    shared_settings: &cli::SharedSettings,
    settings: &cli::DiscoverSettings,
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    sample_name: String,
    chrom_list: &ChromList,
    reference: &GenomeRef,
    target_regions: Option<&GenomeRegions>,
) -> SampleAlignmentScanResult {
    let scan_settings = &ScanSettings {
        depth_bin_size: settings.depth_bin_size,
        min_gap_compressed_identity: settings.min_gap_compressed_identity,
        min_evidence_indel_size: settings.get_min_evidence_indel_size(),
        max_close_breakend_distance: settings.max_close_breakend_distance,
        min_sv_mapq: settings.min_sv_mapq,
        segment_size: 20_000_000,
    };

    assert!(shared_settings.thread_count > 0);

    info!("Processing alignment file '{}'", settings.bam_filename);

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();

    let chrom_count = chrom_list.data.len();

    let genome_kb = chrom_list.data.iter().map(|x| x.length).sum::<u64>() / 1000;
    let progress_reporter =
        ProgressReporter::new(genome_kb, "Processed alignments on", "ref genome kb", false);

    let (tx, rx) = channel();

    let progress_reporter = &progress_reporter;
    let chrom_list_ref = &chrom_list;
    worker_pool.scope(move |scope| {
        for chrom_index in 0..chrom_count {
            let chrom_info = &chrom_list_ref.data[chrom_index];
            let chrom_size = chrom_info.length;
            let chrom_label = &chrom_info.label;
            let chrom_target_regions = match target_regions {
                Some(x) => x.chroms.get(chrom_label),
                None => None,
            };

            if target_regions.is_some() && chrom_target_regions.is_none() {
                continue;
            }

            let worker_thread_dataset = worker_thread_dataset.clone();
            let chrom_sequence: &[u8] = match reference.chroms.get(chrom_label) {
                Some(x) => x,
                None => {
                    error!(
                        "Chromosome \"{chrom_label}\" detected in BAM file, but not in reference file"
                    );
                    std::process::exit(exitcode::DATAERR);
                }
            };

            let tx = tx.clone();
            scope.spawn(move |_| {
                let results = scan_chromosome_segments(
                    worker_thread_dataset,
                    scan_settings,
                    chrom_list_ref,
                    chrom_index,
                    chrom_size,
                    chrom_sequence,
                    chrom_target_regions,
                    progress_reporter,
                );
                tx.send((chrom_index, results)).unwrap();
            });
        }
    });

    let mut genome_break_observations = BreakObservations::new();
    let mut depth_bins = vec![ChromDepthBins::new(); chrom_count];
    let mut read_length = MeanTracker::default();
    let mut genome_depth_bins_builder = vec![DepthBinsBuilder::default(); chrom_count];
    for (chrom_index, (mut break_observations, chrom_depth_bins, chrom_depth_bins_builder)) in rx {
        genome_break_observations.merge(&mut break_observations);
        depth_bins[chrom_index] = chrom_depth_bins.depth_bins;
        read_length.merge(&chrom_depth_bins.read_length);
        genome_depth_bins_builder[chrom_index] = chrom_depth_bins_builder;
    }

    let genome_depth_bins = GenomeDepthBins {
        bin_size: scan_settings.depth_bin_size,
        depth_bins,
        mean_read_length: read_length.mean(),
    };

    progress_reporter.clear();
    info!("Finished processing all alignments");

    SampleAlignmentScanResult {
        sample_name,
        genome_break_observations,
        genome_depth_bins,
        genome_depth_bins_builder,
    }
}
