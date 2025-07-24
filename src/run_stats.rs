//! Track stats for the whole sawfish run
//!

use std::fs::{File, remove_file};

use camino::Utf8Path;
use log::info;
use serde::{Deserialize, Serialize};
use simple_error::{SimpleResult, try_with};
use unwrap::unwrap;

use crate::discover::RUN_STATS_FILENAME;

#[derive(Default, Deserialize, Serialize)]
pub struct ClusterStats {
    pub total_breakpoint_observation_count: usize,
    pub total_breakpoint_cluster_count: usize,
    pub single_region_candidate_count: usize,
    pub consolidated_single_region_candidate_count: usize,
    pub multi_region_candidate_count: usize,
    pub large_insertion_candidate_count: usize,
}

#[derive(Default, Deserialize, Serialize)]
pub struct RefineStats {
    pub candidate_vcf_output_record_count: usize,

    /// Count of candidate VCF records filtered out as duplicates
    ///
    /// Candidates are not typically deduped, this is usually only done at the final output stage
    pub candidate_vcf_duplicate_record_count: usize,

    /// These are duplications which have already been represented by an insertion
    pub redundant_duplications_filtered_out: usize,

    pub total_refinement_time_secs: f64,
    pub total_refinement_assembly_time_secs: f64,
    pub total_refinement_alignment_time_secs: f64,
    pub total_refinement_scoring_time_secs: f64,
}

#[derive(Default, Deserialize, Serialize)]
pub struct MultiSampleCategoryMergeStats {
    pub sv_groups_merged_as_duplicates: usize,
    pub sv_haplotypes_merged_as_duplicates: usize,
    pub sv_candidates_merged_as_duplicates: usize,
    pub sv_groups_after_merge: usize,
    pub sv_haplotypes_after_merge: usize,
    pub sv_candidates_after_merge: usize,
}

impl MultiSampleCategoryMergeStats {
    pub fn merge(&mut self, other: &Self) {
        self.sv_groups_merged_as_duplicates += other.sv_groups_merged_as_duplicates;
        self.sv_haplotypes_merged_as_duplicates += other.sv_haplotypes_merged_as_duplicates;
        self.sv_candidates_merged_as_duplicates += other.sv_candidates_merged_as_duplicates;
        self.sv_groups_after_merge += other.sv_groups_after_merge;
        self.sv_haplotypes_after_merge += other.sv_haplotypes_after_merge;
        self.sv_candidates_after_merge += other.sv_candidates_after_merge;
    }
}

#[derive(Default, Deserialize, Serialize)]
pub struct RunStep {
    pub name: String,
    pub version: String,
}

#[derive(Default, Deserialize, Serialize)]
pub struct MultiSampleMergeStats {
    pub multi_region: MultiSampleCategoryMergeStats,
    pub single_region: MultiSampleCategoryMergeStats,
}

#[derive(Default, Deserialize, Serialize)]
pub struct ScoreStats {
    pub sv_replicate_haplotype_filter: usize,
    pub sv_gt_exclusion_filter: usize,

    pub vcf_output_record_count: usize,
    pub vcf_duplicate_record_count: usize,

    pub total_scoring_time_secs: f64,
}

#[derive(Deserialize, Serialize)]
pub struct DiscoverRunStats {
    // Temporarily allow the run_step to be missing when parsing the run stats file:
    #[serde(default)]
    pub run_step: RunStep,

    pub sample_name: String,
    pub cluster_stats: ClusterStats,
    pub refine_stats: RefineStats,
}

#[derive(Deserialize, Serialize)]
pub struct JointCallRunStats {
    pub run_step: RunStep,
    pub merge_stats: MultiSampleMergeStats,
    pub score_stats: ScoreStats,
}

/// Write run_stats structure out in json format
pub fn write_discover_run_stats(discover_dir: &Utf8Path, run_stats: &DiscoverRunStats) {
    let filename = discover_dir.join(RUN_STATS_FILENAME);

    info!("Writing run statistics to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create run statistics json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &run_stats).unwrap();
}

pub fn read_discover_run_stats(discover_dir: &Utf8Path) -> SimpleResult<DiscoverRunStats> {
    use std::io::BufReader;

    let filename = discover_dir.join(RUN_STATS_FILENAME);
    let file = try_with!(
        File::open(&filename),
        "Unable to read discover-mode run statistics json file: '{filename}'"
    );

    let reader = BufReader::new(file);
    let run_stats = try_with!(
        serde_json::from_reader(reader),
        "Unable to parse discover-mode run statistics from json file: '{filename}'"
    );

    Ok(run_stats)
}

/// Write run_stats structure out in json format
pub fn write_joint_call_run_stats(output_dir: &Utf8Path, run_stats: &JointCallRunStats) {
    let filename = output_dir.join(RUN_STATS_FILENAME);

    info!("Writing run statistics to file: '{filename}'");

    let f = unwrap::unwrap!(
        std::fs::File::create(&filename),
        "Unable to create run statistics json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &run_stats).unwrap();
}

/// Delete run stats file
///
/// Delete a run stats file if one exists. This is typically done during a clobber run,
/// to prevent the old run stats file from being misinterpreted in the event of a crash.
///
pub fn delete_run_stats(output_dir: &Utf8Path) {
    let filename = output_dir.join(RUN_STATS_FILENAME);

    if filename.exists() {
        unwrap!(
            remove_file(&filename),
            "Can't remove original run statistics json file: '{filename}'"
        );
    }
}
