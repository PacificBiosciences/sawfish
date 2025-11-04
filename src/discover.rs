use rust_vc_utils::{ChromList, GenomeRef, get_genome_ref_from_fasta, get_sample_name};

use crate::bam_scanner::SampleAlignmentScanResult;
use crate::bam_scanner::scan_sample_bam_for_sv_evidence;
use crate::cli;
use crate::cluster_breakpoints;
use crate::copy_number_segmentation::{
    HaploidCoverage, get_sample_copy_number_segments, serialize_copy_number_segments,
    write_copy_number_segment_file,
};
use crate::depth_bins;
use crate::filenames::EXPECTED_COPY_NUMBER_BED_FILENAME;
use crate::gc_correction::*;
use crate::genome_regions::GenomeRegions;
use crate::globals::PROGRAM_VERSION;
use crate::large_variant_output::{VcfSettings, write_indexed_sv_vcf_file};
use crate::maf_utils::{scan_maf_file, serialize_maf_data};
use crate::refine_sv;
use crate::run_data::{DiscoverRunData, write_discover_run_data};
use crate::run_stats::{DiscoverRunStats, RunStep, delete_run_stats, write_discover_run_stats};
use crate::sv_scoring_exclusion::generate_and_write_sv_scoring_exclusion_regions;
use crate::worker_thread_data::get_bam_reader_worker_thread_data;

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

/// These steps of discover are only run with CNV is enabled
///
fn run_discover_cnv_steps(
    shared_settings: &cli::SharedSettings,
    settings: &cli::DiscoverSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_scan_result: &SampleAlignmentScanResult,
) -> HaploidCoverage {
    let depth_bin_gc_content = get_depth_bin_gc_content(
        genome_ref,
        chrom_list,
        settings.depth_bin_size,
        settings.gc_genome_window_size,
        shared_settings.thread_count,
    );

    let gc_bias_data = get_gc_correction(
        chrom_list,
        &depth_bin_gc_content,
        sample_scan_result,
        settings.coverage_est_regex.as_str(),
        settings.gc_level_count,
    );

    serialize_gc_bias_correction_data(&settings.output_dir, &gc_bias_data);

    if settings.debug_gc_correction {
        write_gc_correction_debug_output(
            &settings.output_dir,
            chrom_list,
            &depth_bin_gc_content,
            sample_scan_result,
            &gc_bias_data,
        );
    }

    let (haploid_coverage, copy_number_segments) =
        get_sample_copy_number_segments(settings, chrom_list, sample_scan_result, &gc_bias_data);

    write_copy_number_segment_file(&settings.output_dir, chrom_list, &copy_number_segments);

    serialize_copy_number_segments(&settings.output_dir, &copy_number_segments);

    haploid_coverage
}

pub fn run_discover(shared_settings: &cli::SharedSettings, settings: &cli::DiscoverSettings) {
    cli::validate_discover_settings_data(settings);

    // Now that we're committed to a run, remove any possible older run stats file that could be present in case this is a clobber run
    //
    // The run stats file is used as a marker of a successfully finished run, so removing it here allows run completion to be determined
    // from whether the new file is written at the end of this discover step.
    //
    delete_run_stats(&settings.output_dir);

    cli::write_discover_settings(&settings.output_dir, settings);

    let static_settings = StaticDiscoverSettings::new();
    let chrom_list = ChromList::from_bam_filename(&settings.bam_filename);
    let target_regions = if shared_settings.target_region_list.is_empty() {
        None
    } else {
        Some(GenomeRegions::from_target_regions(
            &chrom_list,
            &shared_settings.target_region_list,
            true,
        ))
    };

    let genome_ref = {
        let mut x = get_genome_ref_from_fasta(&settings.ref_filename);
        x.simplify_ambiguous_dna_bases();
        x
    };

    let cnv_excluded_regions = settings
        .cnv_excluded_regions_filename
        .as_ref()
        .map(|x| GenomeRegions::from_bed(x, "excluded", true, false));

    // Setup shared worker thread data structures:
    let bam_reader_worker_thread_dataset = get_bam_reader_worker_thread_data(
        shared_settings,
        &settings.ref_filename,
        &[&settings.bam_filename],
    );

    // Get sample name from bam header
    let sample_name = {
        use rust_htslib::bam::Read;
        // Temporarily access the first bam_reader to extract sample name from the header
        let bam_reader = &bam_reader_worker_thread_dataset[0]
            .lock()
            .unwrap()
            .bam_readers[0];
        get_sample_name(bam_reader.header(), "UnknownSampleName")
    };

    // Write expected copy number regions to discover output directory, these aren't used until the joint-call step
    if let Some(input_filename) = &settings.expected_copy_number_filename {
        let output_filname = settings.output_dir.join(EXPECTED_COPY_NUMBER_BED_FILENAME);
        std::fs::copy(input_filename, output_filname).unwrap();
    }

    if let Some(maf_filename) = &settings.maf_filename {
        // Right now the MAF data is read in from VCF and written to an intermediate format for output in the joint-call
        // step, so we can just close the loop right in this block.
        //
        // Consider dumping this whole thing onto its own thread if the time is annoying
        //

        // By default the sample name searched for in the VCF is that extracted from the BAM header. An alternative
        // sample name can be specified on the command-line:
        let maf_sample_name = if let Some(maf_sample_name) = &settings.maf_sample_name {
            maf_sample_name
        } else {
            &sample_name
        };
        let mut maf_data = scan_maf_file(maf_filename, &chrom_list, &[maf_sample_name]);

        // Change sample name back to the bam sample name if necessary:
        assert!(!maf_data.samples.is_empty());
        if maf_data.samples[0].sample_name != sample_name {
            maf_data.samples[0].sample_name = sample_name.clone();
        }

        serialize_maf_data(&settings.output_dir, &maf_data);
    }

    // Scan input bam for breakend and depth information
    let mut sample_scan_result = scan_sample_bam_for_sv_evidence(
        shared_settings,
        settings,
        &bam_reader_worker_thread_dataset,
        sample_name,
        &chrom_list,
        &genome_ref,
        target_regions.as_ref(),
    );

    // Apply CNV excluded regions to sample_scan_result depth_bins
    depth_bins::apply_excluded_regions(
        &chrom_list,
        cnv_excluded_regions.as_ref(),
        sample_scan_result.genome_depth_bins.bin_size,
        &mut sample_scan_result.genome_depth_bins.depth_bins,
    );

    // Write depth bins out to binary serialization format (for reading back in by sawfish joint genotyper)
    depth_bins::serialize_genome_depth_bins(
        &settings.output_dir,
        &sample_scan_result.genome_depth_bins,
    );

    let haploid_coverage = if !settings.disable_cnv {
        Some(run_discover_cnv_steps(
            shared_settings,
            settings,
            &genome_ref,
            &chrom_list,
            &sample_scan_result,
        ))
    } else {
        None
    };

    let max_sv_scoring_depth = generate_and_write_sv_scoring_exclusion_regions(
        settings,
        &chrom_list,
        haploid_coverage.as_ref(),
        sample_scan_result.genome_depth_bins_builder,
    );

    // Cluster breakpoint signatures
    let (breakpoint_clusters, cluster_stats) = cluster_breakpoints::process_breakpoint_clusters(
        settings,
        static_settings.assembly_read_flank_size,
        &chrom_list,
        &mut sample_scan_result.genome_break_observations,
    );

    // Refine the SVs
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

    let vcf_settings = {
        let no_vcf_dedup = false;
        let enable_phasing = false;
        let treat_single_copy_as_haploid = false;
        VcfSettings::new(
            &settings.ref_filename,
            &settings.output_dir,
            settings.min_qual,
            no_vcf_dedup,
            enable_phasing,
            treat_single_copy_as_haploid,
        )
    };

    // Write SV candidate BCF output
    let vcf_stats = write_indexed_sv_vcf_file(
        shared_settings,
        &vcf_settings,
        &genome_ref,
        &chrom_list,
        &[],
        &sv_candidate_groups,
        &[],
        true,
    );

    refine_stats.candidate_vcf_output_record_count = vcf_stats.output_record_count;
    refine_stats.candidate_vcf_duplicate_record_count = vcf_stats.duplicate_record_count;

    let run_step = RunStep {
        name: "discover".to_string(),
        version: PROGRAM_VERSION.to_string(),
    };

    write_discover_run_data(
        &settings.output_dir,
        &DiscoverRunData {
            run_step: run_step.clone(),
            sample_name: sample_scan_result.sample_name.clone(),
            max_sv_scoring_depth,
        },
    );

    // In addition to useful statistics this file acts as a marker for a successfully completed run, so it must be written last.
    write_discover_run_stats(
        &settings.output_dir,
        &DiscoverRunStats {
            run_step,
            sample_name: sample_scan_result.sample_name.clone(),
            cluster_stats,
            refine_stats,
        },
    );
}
