use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::channel;
use std::time::{Duration, Instant};

use log::{Level, info, log_enabled};
use rust_vc_utils::{GenomeSegment, ProgressReporter};

use super::bam_prefetch::{ChunkProgress, get_prefetch_chunk_lead, prefetch_chunks};
use super::{SampleJointCallData, SharedJointCallData};
use crate::cli::{JointCallSettings, SharedSettings};
use crate::log_utils::debug_msg;
use crate::refine_sv::SVFilterType;
use crate::run_stats::ScoreStats;
use crate::score_sv::{
    MAX_WORKER_CACHE_BYTES, SampleReadCache, SampleScoreData, ScoreDebugSettings, ScoreSVSettings,
    get_search_segments, score_and_assess_refined_sv_group,
};
use crate::sv_group::SVGroup;
use crate::worker_thread_data::{BamReaderWorkerThreadDataSet, get_bam_reader_worker_thread_data};

/// Index of the SV group in the merged list, scoring time, and scored SV group
type ScoreWorkerReturnType = (usize, Duration, SVGroup);

/// Score a chunk of SV groups, which are adjacent in the genome so that they can share reads
/// through the worker's read cache
#[allow(clippy::too_many_arguments)]
fn score_sv_candidate_cluster_chunk(
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    worker_read_caches: &[Mutex<Vec<SampleReadCache>>],
    tx: std::sync::mpsc::Sender<ScoreWorkerReturnType>,
    score_settings: &ScoreSVSettings,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleScoreData],
    sv_groups: Vec<(usize, SVGroup)>,
    plan: &[GenomeSegment],
    progress_reporter: &ProgressReporter,
) {
    let worker_id = rayon::current_thread_index().unwrap();

    let bam_readers = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_readers;
    let bam_readers_ref = &mut bam_readers.iter_mut().collect::<Vec<_>>();
    let read_caches = &mut worker_read_caches[worker_id].lock().unwrap();

    // Plan the read cache windows from the regions of all groups in the chunk
    for read_cache in read_caches.iter_mut() {
        read_cache.set_plan(plan.to_vec());
    }

    for (group_index, sv_group) in sv_groups {
        score_sv_candidate_cluster(
            group_index,
            bam_readers_ref,
            read_caches,
            tx.clone(),
            score_settings,
            shared_data,
            all_sample_data,
            sv_group,
            progress_reporter,
        );
        SampleReadCache::trim_to_budget(read_caches, MAX_WORKER_CACHE_BYTES);
    }
}

#[allow(clippy::too_many_arguments)]
fn score_sv_candidate_cluster(
    group_index: usize,
    bam_readers_ref: &mut [&mut rust_htslib::bam::IndexedReader],
    read_caches: &mut [SampleReadCache],
    tx: std::sync::mpsc::Sender<ScoreWorkerReturnType>,
    score_settings: &ScoreSVSettings,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleScoreData],
    mut sv_group: SVGroup,
    progress_reporter: &ProgressReporter,
) {
    // Get total processing time for this cluster:
    let cluster_start_time = Instant::now();

    let debug_settings = ScoreDebugSettings {
        debug: false,
        debug_sample_index: None,
        debug_sv_index: None,
    };

    let any_debug = debug_settings.debug;
    if any_debug || log_enabled!(Level::Debug) {
        debug_msg!(
            any_debug,
            "Starting score for new SVGroup on worker {}",
            rayon::current_thread_index().unwrap()
        );
        for (region_index, group_region) in sv_group.group_regions.iter().enumerate() {
            debug_msg!(any_debug, "Group region {region_index}: {group_region:?}");
        }
        for (rsv_index, refined_sv) in sv_group.refined_svs.iter().enumerate() {
            debug_msg!(any_debug, "Refined SV {rsv_index}: {:?}", refined_sv.bp);
        }
    }

    score_and_assess_refined_sv_group(
        score_settings,
        &shared_data.reference,
        &shared_data.chrom_list,
        all_sample_data,
        bam_readers_ref,
        read_caches,
        &mut sv_group,
        &debug_settings,
    );

    let duration = cluster_start_time.elapsed();
    let result = (group_index, duration, sv_group);
    tx.send(result).unwrap();
    progress_reporter.inc(1);
}

/// Run joint scoring of all samples over multiple threads
pub(super) fn joint_genotype_all_samples(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    enable_phasing: bool,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    merged_sv_groups: Vec<SVGroup>,
) -> (Vec<SVGroup>, ScoreStats) {
    // Create a projection of all_sample_data into the corresponding scoring data structure
    let all_sample_scoring_data = &all_sample_data
        .iter()
        .map(|x| x.to_sample_score_data())
        .collect::<Vec<_>>();

    let bam_filenames = all_sample_data
        .iter()
        .map(|x| x.bam_filename.as_str())
        .collect::<Vec<_>>();

    let worker_thread_dataset = &get_bam_reader_worker_thread_data(
        shared_settings,
        &shared_data.ref_filename,
        &bam_filenames,
    );

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();

    info!("Starting SV candidate cluster genotyping");

    // Turn off the progress bar if we're logging debug info during genotyping:
    let force_periodic_updates = shared_settings.debug;
    let progress_reporter = ProgressReporter::new(
        merged_sv_groups.len() as u64,
        "Genotyped",
        "SV candidate clusters",
        force_periodic_updates,
    );
    let progress_reporter = &progress_reporter;

    let score_settings = &ScoreSVSettings::new(
        settings.min_sv_mapq,
        settings.min_gap_compressed_identity,
        settings.max_qscore,
        settings.report_supporting_reads,
        enable_phasing,
        settings.treat_single_copy_as_haploid,
    );

    let sample_count = all_sample_data.len();
    let worker_read_caches = &(0..shared_settings.thread_count)
        .map(|_| Mutex::new((0..sample_count).map(|_| SampleReadCache::new()).collect()))
        .collect::<Vec<_>>();

    // Score SV groups in chunks of groups adjacent in the genome, so that each worker can reuse
    // reads between the groups of a chunk. Results don't depend on the order in which groups are
    // scored.
    //
    // Chunk sizes shrink as the remaining groups do (guided scheduling), so that the last chunks
    // are small enough to balance the work between workers at the end.
    let mut merged_sv_groups = merged_sv_groups.into_iter().enumerate().collect::<Vec<_>>();
    merged_sv_groups.sort_by_key(|(_, x)| {
        let r = &x.group_regions[0];
        (r.chrom_index, r.range.start, r.range.end)
    });
    const MAX_CHUNK_SIZE: usize = 64;
    const MIN_CHUNK_SIZE: usize = 4;
    let mut chunks = Vec::new();
    let mut groups = merged_sv_groups.into_iter().peekable();
    let mut remaining = groups.len();
    while groups.peek().is_some() {
        let chunk_size =
            (remaining / (4 * shared_settings.thread_count)).clamp(MIN_CHUNK_SIZE, MAX_CHUNK_SIZE);
        let chunk = groups.by_ref().take(chunk_size).collect::<Vec<_>>();
        remaining -= chunk.len();
        chunks.push(chunk);
    }

    // Segments read by each chunk, for the read cache and read-ahead
    let chunk_segments = &chunks
        .iter()
        .map(|chunk| {
            chunk
                .iter()
                .flat_map(|(_, x)| get_search_segments(&shared_data.chrom_list, &x.group_regions))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    // Each worker takes the next chunk in genome order when it is free. (Spawning one task per chunk
    // would run the spawning worker's own tasks last-in first-out, so it would start at the end.)
    let chunk_count = chunks.len();
    let chunks = &chunks
        .into_iter()
        .map(|chunk| Mutex::new(Some(chunk)))
        .collect::<Vec<_>>();
    let next_chunk = &AtomicUsize::new(0);

    // The read-ahead stays within a fixed lead of the lowest chunk that isn't scored yet
    let progress = &ChunkProgress::new(chunk_count);
    let stop_prefetch = &AtomicBool::new(false);

    let (tx, rx) = channel();
    std::thread::scope(|thread_scope| {
        thread_scope.spawn(|| {
            prefetch_chunks(
                &bam_filenames,
                chunk_segments,
                progress,
                get_prefetch_chunk_lead(shared_settings.thread_count),
                stop_prefetch,
            )
        });
        worker_pool.scope(move |scope| {
            for _ in 0..shared_settings.thread_count {
                let tx = tx.clone();
                scope.spawn(move |_| {
                    loop {
                        let chunk_index = next_chunk.fetch_add(1, Ordering::Relaxed);
                        if chunk_index >= chunk_count {
                            break;
                        }
                        let chunk = chunks[chunk_index].lock().unwrap().take().unwrap();
                        score_sv_candidate_cluster_chunk(
                            worker_thread_dataset,
                            worker_read_caches,
                            tx.clone(),
                            score_settings,
                            shared_data,
                            all_sample_scoring_data,
                            chunk,
                            &chunk_segments[chunk_index],
                            progress_reporter,
                        );

                        progress.mark_finished(chunk_index);
                    }
                });
            }
        });
        stop_prefetch.store(true, Ordering::Relaxed);
    });

    let mut score_stats = ScoreStats::default();
    let mut scored_sv_groups = Vec::new();
    for (group_index, chunk_duration, chunk_scored_svs) in rx {
        score_stats.total_scoring_time_secs += chunk_duration.as_secs_f64();
        scored_sv_groups.push((group_index, chunk_scored_svs));
    }
    scored_sv_groups.sort_by_key(|(group_index, _)| *group_index);
    let scored_sv_groups = scored_sv_groups
        .into_iter()
        .map(|(_, x)| x)
        .collect::<Vec<_>>();

    // Get filtered SV stats
    for refined_sv in scored_sv_groups.iter().flat_map(|x| x.refined_svs.iter()) {
        match refined_sv.get_sv_filter_type() {
            SVFilterType::GTExcluded => score_stats.sv_gt_exclusion_filter += 1,
            SVFilterType::Replicated => score_stats.sv_replicate_haplotype_filter += 1,
            SVFilterType::NoFilter => (),
        }
    }

    progress_reporter.clear();
    info!("Finished genotyping all SV candidate clusters");

    (scored_sv_groups, score_stats)
}
