mod find_inversions;
mod get_refined_svs;
mod joint_call_all_samples;
mod merge_haplotypes;
mod read_sample_data;
mod supporting_read_names;

use self::find_inversions::find_inversions;
use self::joint_call_all_samples::joint_genotype_all_samples;
use self::merge_haplotypes::merge_haplotypes;
pub use self::read_sample_data::SampleJointCallData;
use self::read_sample_data::{read_all_sample_data, SharedJointCallData};

use crate::cli::{JointCallSettings, SharedSettings};
use crate::large_variant_output::{write_indexed_sv_vcf_file, VcfSettings};
use crate::run_stats::{write_joint_call_run_stats, JointCallRunStats};

pub fn run_joint_call(shared_settings: &SharedSettings, settings: &JointCallSettings) {
    let (shared_data, all_sample_data) = read_all_sample_data(shared_settings, settings);

    let (merged_sv_groups, merge_stats) =
        merge_haplotypes(shared_settings, &shared_data.chrom_list, &all_sample_data);

    let enable_phasing = true;
    let (mut scored_svs, mut score_stats) = joint_genotype_all_samples(
        shared_settings,
        settings,
        enable_phasing,
        &shared_data,
        &all_sample_data,
        merged_sv_groups,
    );

    find_inversions(&mut scored_svs);

    let vcf_settings = VcfSettings::new(
        &shared_data.ref_filename,
        &settings.output_dir,
        settings.min_qual,
        settings.no_vcf_dedup,
        enable_phasing,
    );

    let sample_names = all_sample_data
        .iter()
        .map(|x| x.sample_name.as_str())
        .collect::<Vec<_>>();
    let vcf_stats = write_indexed_sv_vcf_file(
        shared_settings,
        &vcf_settings,
        &shared_data.genome_ref,
        &shared_data.chrom_list,
        &sample_names,
        &scored_svs,
        false,
    );

    score_stats.vcf_output_record_count = vcf_stats.output_record_count;
    score_stats.vcf_duplicate_record_count = vcf_stats.duplicate_record_count;

    write_joint_call_run_stats(
        &settings.output_dir,
        &JointCallRunStats {
            merge_stats,
            score_stats,
        },
    );

    if settings.report_supporting_reads {
        supporting_read_names::write_supporting_read_names(
            &settings.output_dir,
            &sample_names,
            &scored_svs,
        );
    }
}
