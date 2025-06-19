use std::sync::mpsc::channel;
use std::time::{Duration, Instant};

use log::{Level, info, log_enabled};
use rust_vc_utils::ProgressReporter;

use super::{SampleJointCallData, SharedJointCallData};
use crate::cli::{JointCallSettings, SharedSettings};
use crate::log_utils::debug_msg;
use crate::refine_sv::SVFilterType;
use crate::run_stats::ScoreStats;
use crate::score_sv::{
    SampleScoreData, ScoreDebugSettings, ScoreSVSettings, score_and_assess_refined_sv_group,
};
use crate::sv_group::SVGroup;
use crate::worker_thread_data::{BamReaderWorkerThreadDataSet, get_bam_reader_worker_thread_data};

type ScoreWorkerReturnType = (Duration, SVGroup);

fn score_sv_candidate_cluster_wrapper(
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    tx: std::sync::mpsc::Sender<ScoreWorkerReturnType>,
    score_settings: &ScoreSVSettings,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleScoreData],
    mut sv_group: SVGroup,
    progress_reporter: &ProgressReporter,
) {
    // Get total processing time for this cluster:
    let cluster_start_time = Instant::now();

    let worker_id = rayon::current_thread_index().unwrap();

    let bam_readers = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_readers;
    let bam_readers_ref = &mut bam_readers.iter_mut().collect::<Vec<_>>();

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
            worker_id
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
        &shared_data.genome_ref,
        &shared_data.chrom_list,
        all_sample_data,
        bam_readers_ref,
        &mut sv_group,
        &debug_settings,
    );

    let duration = cluster_start_time.elapsed();
    let result = (duration, sv_group);
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

    // Setup shared worker thread data structures:
    let all_discover_settings_ref = all_sample_data
        .iter()
        .map(|x| &x.discover_settings)
        .collect::<Vec<_>>();

    let worker_thread_dataset =
        &get_bam_reader_worker_thread_data(shared_settings, &all_discover_settings_ref);

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

    let (tx, rx) = channel();
    worker_pool.scope(move |scope| {
        for cluster_rsvs in merged_sv_groups {
            let tx = tx.clone();
            scope.spawn(move |_| {
                score_sv_candidate_cluster_wrapper(
                    worker_thread_dataset,
                    tx,
                    score_settings,
                    shared_data,
                    all_sample_scoring_data,
                    cluster_rsvs,
                    progress_reporter,
                );
            });
        }
    });

    let mut score_stats = ScoreStats::default();
    let mut scored_sv_groups = Vec::new();
    for (chunk_duration, chunk_scored_svs) in rx {
        score_stats.total_scoring_time_secs += chunk_duration.as_secs_f64();
        scored_sv_groups.push(chunk_scored_svs);
    }

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
