use std::sync::mpsc::channel;

use log::info;
use rust_htslib::bam::{self, Read};
use rust_vc_utils::{
    ChromList, MeanTracker, ProgressReporter, RegionMap, filter_out_alignment_record,
    get_region_segments,
};
use unwrap::unwrap;

use crate::cli;
use crate::cluster_breakpoints::{BreakBuilder, BreakObservations};
use crate::depth_bins::{AllChromDepthBinInfo, ChromDepthBins, DepthBinsBuilder, GenomeDepthBins};
use crate::genome_ref_utils::{self, GenomeRef};
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
    reference: &GenomeRef,
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

        // The chromosome sequence is only loaded once it's needed, so that chromosomes without any
        // alignments are never loaded
        if genome_ref_utils::get_gap_compressed_identity(&record, reference, chrom_index)
            < scan_settings.min_gap_compressed_identity
        {
            continue;
        }

        // Decode the CIGAR once for all of the record's consumers below. The cache survives into
        // the next record read, so it has to be set again for every record that is used.
        record.cache_cigar();

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

/// Scan several chromosomes at once, so that small chromosomes are scanned in parallel
///
/// Returns the results of each chromosome, in the order of `chroms`
///
#[allow(clippy::too_many_arguments)]
fn scan_chromosome_batch(
    worker_thread_dataset: BamReaderWorkerThreadDataSet,
    scan_settings: &ScanSettings,
    chrom_list: &ChromList,
    chroms: &[(usize, Option<&RegionMap>)],
    reference: &GenomeRef,
    progress_reporter: &ProgressReporter,
) -> Vec<(BreakObservations, AllChromDepthBinInfo, DepthBinsBuilder)> {
    let (tx, rx) = channel();
    rayon::scope(|scope| {
        for (slot, &(chrom_index, chrom_target_regions)) in chroms.iter().enumerate() {
            let chrom_size = chrom_list.data[chrom_index].length;
            let is_targeted_scan = chrom_target_regions.is_some();
            let segments = get_chrom_worker_segments(
                chrom_size,
                chrom_target_regions,
                scan_settings.segment_size,
            );
            for (begin, end) in segments {
                let worker_thread_dataset = worker_thread_dataset.clone();
                let tx = tx.clone();
                scope.spawn(move |_| {
                    let worker_id = rayon::current_thread_index().unwrap();
                    let sample_index = 0;
                    let bam_reader =
                        &mut worker_thread_dataset[worker_id].lock().unwrap().bam_readers
                            [sample_index];
                    let result = scan_chromosome_segment(
                        bam_reader,
                        scan_settings,
                        chrom_list,
                        chrom_index,
                        begin as i64,
                        end as i64,
                        reference,
                        is_targeted_scan,
                        progress_reporter,
                    );
                    tx.send((slot, result)).unwrap();
                });
            }
        }
    });
    drop(tx);

    let mut per_chrom: Vec<Option<(BreakObservations, DepthBinsBuilder)>> =
        chroms.iter().map(|_| None).collect();
    for (slot, (mut obs, depth)) in rx {
        match per_chrom[slot].as_mut() {
            None => per_chrom[slot] = Some((obs, depth)),
            Some((chrom_obs, chrom_depth)) => {
                chrom_obs.merge(&mut obs);
                chrom_depth.merge(&depth);
            }
        }
    }
    chroms
        .iter()
        .zip(per_chrom)
        .map(|(&(chrom_index, _), x)| {
            let (obs, depth_builder) =
                x.unwrap_or_else(|| (BreakObservations::new(), DepthBinsBuilder::default()));
            let chrom_size = chrom_list.data[chrom_index].length;
            let depth_bins = depth_builder.get_depth_bins(scan_settings.depth_bin_size, chrom_size);
            (obs, depth_bins, depth_builder)
        })
        .collect()
}

/// Index of the first of the small contigs (decoys, unplaced) which follow the primary chromosomes
const FIRST_SMALL_CONTIG_INDEX: usize = 24;

type ScanBatch<'a> = Vec<(usize, Option<&'a RegionMap>)>;

/// Group the chromosomes to scan into batches, which are scanned one after another
///
/// The small contigs after the primary chromosomes are scanned first, all at once, so that
/// refinement of the clusters which involve them can start early, and the scan of single segment
/// contigs doesn't leave threads idle. Primary chromosomes follow one at a time, in order.
///
/// `target_of` returns None for a chromosome that isn't scanned, or its target regions, which are
/// None for a whole-chromosome scan.
///
fn get_scan_batches<'a>(
    chrom_count: usize,
    target_of: impl Fn(usize) -> Option<Option<&'a RegionMap>>,
) -> Vec<ScanBatch<'a>> {
    let mut batches: Vec<ScanBatch> = Vec::new();
    if chrom_count > FIRST_SMALL_CONTIG_INDEX {
        let small = (FIRST_SMALL_CONTIG_INDEX..chrom_count)
            .filter_map(|i| target_of(i).map(|t| (i, t)))
            .collect::<Vec<_>>();
        if !small.is_empty() {
            batches.push(small);
        }
    }
    for chrom_index in 0..chrom_count.min(FIRST_SMALL_CONTIG_INDEX) {
        if let Some(t) = target_of(chrom_index) {
            batches.push(vec![(chrom_index, t)]);
        }
    }
    batches
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
    scan_threads: usize,
    on_chrom_complete: &dyn Fn(usize, &BreakObservations),
) -> SampleAlignmentScanResult {
    let scan_settings = &ScanSettings {
        depth_bin_size: settings.depth_bin_size,
        min_gap_compressed_identity: settings.min_gap_compressed_identity,
        min_evidence_indel_size: settings.get_min_evidence_indel_size(),
        max_close_breakend_distance: settings.max_close_breakend_distance,
        min_sv_mapq: settings.min_sv_mapq,
        segment_size: 20_000_000,
    };

    assert!(scan_threads > 0 && scan_threads <= shared_settings.thread_count);

    info!("Processing alignment file '{}'", settings.bam_filename);

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(scan_threads)
        .build()
        .unwrap();

    let chrom_count = chrom_list.data.len();

    let genome_kb = chrom_list.data.iter().map(|x| x.length).sum::<u64>() / 1000;
    let progress_reporter =
        ProgressReporter::new(genome_kb, "Processed alignments on", "ref genome kb", false);

    // Chromosomes are scanned one at a time, in order, with the segments of each chromosome
    // scanned in parallel, so that the evidence from each chromosome is complete as early as
    // possible
    //
    let mut genome_break_observations = BreakObservations::new();
    let mut depth_bins = vec![ChromDepthBins::new(); chrom_count];
    let mut chrom_read_lengths = vec![MeanTracker::default(); chrom_count];
    let mut genome_depth_bins_builder = vec![DepthBinsBuilder::default(); chrom_count];

    let target_of = |chrom_index: usize| -> Option<Option<&RegionMap>> {
        let chrom_label = &chrom_list.data[chrom_index].label;
        let chrom_target_regions = match target_regions {
            Some(x) => x.chroms.get(chrom_label),
            None => None,
        };
        if target_regions.is_some() && chrom_target_regions.is_none() {
            None
        } else {
            Some(chrom_target_regions)
        }
    };
    let batches = get_scan_batches(chrom_count, target_of);

    for batch in batches {
        let results = worker_pool.install(|| {
            scan_chromosome_batch(
                worker_thread_dataset.clone(),
                scan_settings,
                chrom_list,
                &batch,
                reference,
                &progress_reporter,
            )
        });
        for (
            &(chrom_index, _),
            (mut break_observations, chrom_depth_bins, chrom_depth_bins_builder),
        ) in batch.iter().zip(results)
        {
            on_chrom_complete(chrom_index, &break_observations);
            genome_break_observations.merge(&mut break_observations);
            depth_bins[chrom_index] = chrom_depth_bins.depth_bins;
            chrom_read_lengths[chrom_index] = chrom_depth_bins.read_length;
            genome_depth_bins_builder[chrom_index] = chrom_depth_bins_builder;
        }
    }

    // Merged in chromosome order, so that the result doesn't depend on the scan order
    let mut read_length = MeanTracker::default();
    for x in chrom_read_lengths.iter() {
        read_length.merge(x);
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

#[cfg(test)]
mod tests {
    use super::*;

    fn chroms(batches: &[ScanBatch]) -> Vec<Vec<usize>> {
        batches
            .iter()
            .map(|b| b.iter().map(|x| x.0).collect())
            .collect()
    }

    #[test]
    fn test_get_scan_batches() {
        let all = |_: usize| Some(None);

        // Small contigs form the first batch, then primary chromosomes one at a time in order
        let batches = get_scan_batches(28, all);
        let mut expected = vec![vec![24, 25, 26, 27]];
        expected.extend((0..24).map(|i| vec![i]));
        assert_eq!(chroms(&batches), expected);

        // No small contigs
        assert_eq!(
            chroms(&get_scan_batches(FIRST_SMALL_CONTIG_INDEX, all)),
            (0..24).map(|i| vec![i]).collect::<Vec<_>>()
        );
        assert_eq!(
            chroms(&get_scan_batches(3, all)),
            vec![vec![0], vec![1], vec![2]]
        );
        assert!(get_scan_batches(0, all).is_empty());

        // One small contig
        assert_eq!(chroms(&get_scan_batches(25, all))[0], vec![24]);

        // Chromosomes without target regions are skipped, and a batch left empty is dropped
        let only_some = |i: usize| (i == 2 || i == 26).then_some(None);
        assert_eq!(
            chroms(&get_scan_batches(28, only_some)),
            vec![vec![26], vec![2]]
        );
        let only_primary = |i: usize| [0, 10, 20].contains(&i).then_some(None);
        assert_eq!(
            chroms(&get_scan_batches(28, only_primary)),
            vec![vec![0], vec![10], vec![20]]
        );

        // Every chromosome is in exactly one batch
        let mut all_chroms = chroms(&get_scan_batches(40, all))
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        all_chroms.sort();
        assert_eq!(all_chroms, (0..40).collect::<Vec<_>>());
    }
}
