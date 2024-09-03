//! Track stats for the whole sawfish run
//!

use std::fs::File;
use std::path::Path;

use log::info;
use serde::{Deserialize, Serialize};
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
    pub sample_name: String,
    pub cluster_stats: ClusterStats,
    pub refine_stats: RefineStats,
}

#[derive(Deserialize, Serialize)]
pub struct JointCallRunStats {
    pub merge_stats: MultiSampleMergeStats,
    pub score_stats: ScoreStats,
}

/// Write run_stats structure out in json format
pub fn write_discover_run_stats(discover_dir: &Path, run_stats: &DiscoverRunStats) {
    let filename = discover_dir.join(RUN_STATS_FILENAME);

    info!("Writing run statistics to file: '{}'", filename.display());

    let f = unwrap!(
        File::create(&filename),
        "Unable to create run statistics json file: '{}'",
        filename.display()
    );

    serde_json::to_writer_pretty(&f, &run_stats).unwrap();
}

pub fn read_discover_run_stats(discover_dir: &Path) -> DiscoverRunStats {
    use std::io::BufReader;

    let filename = discover_dir.join(RUN_STATS_FILENAME);
    let file = unwrap!(
        File::open(&filename),
        "Unable to read discover-mode run stats json file: `{}`",
        filename.display()
    );
    let reader = BufReader::new(file);
    unwrap!(
        serde_json::from_reader(reader),
        "Unable to parse discover-mode settings from json file: `{}`",
        filename.display()
    )
}

/// Write run_stats structure out in json format
pub fn write_joint_call_run_stats(output_dir: &Path, run_stats: &JointCallRunStats) {
    let filename = output_dir.join("run_stats.json");

    info!("Writing run statistics to file: '{}'", filename.display());

    let f = unwrap::unwrap!(
        std::fs::File::create(&filename),
        "Unable to create run statistics json file: '{}'",
        filename.display()
    );

    serde_json::to_writer_pretty(&f, &run_stats).unwrap();
}
