use log::info;
use rust_vc_utils::{get_genome_ref_from_fasta, ChromList};

use crate::bam_scanner::scan_sample_bam_for_sv_evidence;
use crate::cli;
use crate::cluster_breakpoints;
use crate::cnv_output;
use crate::copy_number_segmentation::{segment_sample_copy_number, write_copy_number_segment_file};
use crate::depth_bins;
use crate::gc_correction::*;
use crate::genome_regions::{write_genome_regions_to_bed, GenomeRegions};
use crate::large_variant_output::{write_indexed_sv_vcf_file, VcfSettings};
use crate::maf_utils;
use crate::refine_sv;
use crate::run_stats::{write_discover_run_stats, DiscoverRunStats};
use crate::worker_thread_data::get_bam_reader_worker_thread_data;

pub const ASSEMBLY_REGIONS_FILENAME: &str = "assembly.regions.bed";
pub const CANDIDATE_SV_FILENAME: &str = "candidate.sv.bcf";
pub const CONTIG_ALIGNMENT_FILENAME: &str = "contig.alignment.bam";
pub const DEPTH_BINS_BIGWIG_FILENAME: &str = "depth.bw";
pub const DEPTH_BINS_MESSAGEPACK_FILENAME: &str = "depth.mpack";
pub const EXPECTED_COPY_NUMBER_FILENAME: &str = "expected.copy.number.bed";
pub const GENOME_GC_LEVELS_MESSAGEPACK_FILENAME: &str = "genome.gclevels.mpack";
pub const MAX_SV_DEPTH_FILENAME: &str = "max.depth.bed";
pub const RUN_STATS_FILENAME: &str = "run.stats.json";
pub const SAMPLE_GC_BIAS_MESSAGEPACK_FILENAME: &str = "sample.gcbias.mpack";
pub const SETTINGS_FILENAME: &str = "discover.settings.json";

/// Settings shared between discover processes that are not intended to be runtime settable in the cli
pub struct StaticDiscoverSettings {
    /// Distance in read coordinate from the candidate SV breakends to extend read trimming for read
    /// input to assembly
    pub assembly_read_flank_size: usize,
}

impl StaticDiscoverSettings {
    fn new() -> Self {
        Self {
            assembly_read_flank_size: 300,
        }
    }
}

pub fn run_discover(shared_settings: &cli::SharedSettings, settings: &cli::DiscoverSettings) {
    cli::validate_discover_settings_data(settings);
    cli::write_discover_settings(&settings.output_dir, settings);

    let static_settings = StaticDiscoverSettings::new();
    let chrom_list = ChromList::from_bam_filename(&settings.bam_filename);
    let target_regions =
        GenomeRegions::from_target_regions(&chrom_list, &shared_settings.target_region_list, true);

    let genome_ref = get_genome_ref_from_fasta(&settings.ref_filename);

    let default_empty_filename = "".to_string();
    let excluded_regions = GenomeRegions::from_bed(
        settings
            .exclude_filename
            .as_ref()
            .unwrap_or(&default_empty_filename),
        "excluded",
        true,
        false,
    );
    let expected_copy_number_map = GenomeRegions::from_bed(
        settings
            .expected_copy_number_filename
            .as_ref()
            .unwrap_or(&default_empty_filename),
        "expected copy number",
        false,
        true,
    );

    // Setup shared worker thread data structures:
    let bam_reader_worker_thread_dataset =
        get_bam_reader_worker_thread_data(shared_settings, &[settings]);

    let mut sample_scan_result = scan_sample_bam_for_sv_evidence(
        shared_settings,
        settings,
        &bam_reader_worker_thread_dataset,
        &chrom_list,
        &excluded_regions,
        &genome_ref,
        &target_regions,
    );

    {
        // Write max_sv_depth to a BED file
        let max_sv_depth_bed = settings.output_dir.join(MAX_SV_DEPTH_FILENAME);
        write_genome_regions_to_bed(
            "max sv depth",
            &max_sv_depth_bed,
            &chrom_list,
            &sample_scan_result.genome_max_sv_depth_regions,
        );
    }

    let depth_bin_gc_content = get_depth_bin_gc_content(
        &genome_ref,
        &chrom_list,
        settings.depth_bin_size,
        settings.gc_genome_window_size,
        shared_settings.thread_count,
    );

    let gc_bias_data = get_gc_correction(
        &chrom_list,
        &depth_bin_gc_content,
        &sample_scan_result,
        settings.coverage_est_regex.as_str(),
        settings.gc_level_count,
    );

    let (breakpoint_clusters, cluster_stats) = cluster_breakpoints::process_breakpoint_clusters(
        settings,
        static_settings.assembly_read_flank_size,
        &chrom_list,
        &mut sample_scan_result.genome_break_observations,
    );

    // Start refining the SVs
    let (sv_candidate_groups, mut refine_stats) = {
        refine_sv::refine_sv_candidates(
            shared_settings,
            settings,
            &static_settings,
            &bam_reader_worker_thread_dataset,
            &genome_ref,
            &chrom_list,
            &breakpoint_clusters,
        )
    };

    let no_vcf_dedup = false;
    let enable_phasing = false;
    let vcf_settings = VcfSettings::new(
        &settings.ref_filename,
        &settings.output_dir,
        settings.min_qual,
        no_vcf_dedup,
        enable_phasing,
    );

    // Write SV candidate BCF output
    let vcf_stats = write_indexed_sv_vcf_file(
        shared_settings,
        &vcf_settings,
        &genome_ref,
        &chrom_list,
        &[],
        &sv_candidate_groups,
        true,
    );

    refine_stats.candidate_vcf_output_record_count = vcf_stats.output_record_count;
    refine_stats.candidate_vcf_duplicate_record_count = vcf_stats.duplicate_record_count;

    // Write various CNV data structure required for joint genotyping:
    {
        // Write expected copy number regions
        if let Some(input_filename) = &settings.expected_copy_number_filename {
            let output_filname = settings.output_dir.join(EXPECTED_COPY_NUMBER_FILENAME);
            std::fs::copy(input_filename, output_filname).unwrap();
        }

        // Write depth bins out to bigwig track (for reading in IGV)
        let bigwig_filename = &settings.output_dir.join(DEPTH_BINS_BIGWIG_FILENAME);
        depth_bins::write_depth_bigwig_file(
            bigwig_filename,
            &sample_scan_result.genome_depth_bins,
            &chrom_list,
            "depth",
        );

        // Write depth bins out to binary serialization format (for reading back in by sawfish joint
        // genotyper)
        depth_bins::serialize_genome_depth_bins(
            &settings.output_dir,
            &sample_scan_result.genome_depth_bins,
        );

        serialize_gc_bias_correction_data(&settings.output_dir, &gc_bias_data);
    }

    // Disable CNV output until we design multi-sample and SV integration behaviors
    let enable_cnv_output = false;

    // Skip most other CNV-processing in targeted SV mode for now
    let output_cnv_result = enable_cnv_output && target_regions.is_empty();

    if output_cnv_result {
        if settings.debug_gc_correction {
            write_gc_correction_debug_output(
                &settings.output_dir,
                &chrom_list,
                &depth_bin_gc_content,
                &sample_scan_result,
                &gc_bias_data,
            );
        }

        let maf_scan_result = maf_utils::scan_maf_file(
            settings
                .maf_filename
                .as_ref()
                .unwrap_or(&default_empty_filename),
            &excluded_regions,
        );

        if let Some(maf_scan_result) = maf_scan_result {
            maf_utils::write_maf_track_files(&settings.output_dir, &chrom_list, &maf_scan_result);
        }

        info!("Segmenting copy number");
        let copy_number_segment_result =
            segment_sample_copy_number(settings, &chrom_list, &sample_scan_result, &gc_bias_data);

        write_copy_number_segment_file(&settings.output_dir, &copy_number_segment_result);
        cnv_output::write_indexed_cnv_vcf_file(
            shared_settings,
            settings,
            &copy_number_segment_result,
            &expected_copy_number_map,
        );
    }

    write_discover_run_stats(
        &settings.output_dir,
        &DiscoverRunStats {
            sample_name: sample_scan_result.sample_name.clone(),
            cluster_stats,
            refine_stats,
        },
    );
}
