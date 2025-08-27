use rust_vc_utils::{ChromList, get_genome_ref_from_fasta, get_sample_name};

use crate::bam_scanner::scan_sample_bam_for_sv_evidence;
use crate::cli;
use crate::cluster_breakpoints;
use crate::copy_number_segmentation::{
    get_sample_copy_number_segments, serialize_copy_number_segments, write_copy_number_segment_file,
};
use crate::depth_bins;
use crate::gc_correction::*;
use crate::genome_regions::{GenomeRegions, write_genome_regions_to_bed};
use crate::globals::PROGRAM_VERSION;
use crate::large_variant_output::{VcfSettings, write_indexed_sv_vcf_file};
use crate::maf_utils::{scan_maf_file, serialize_maf_data};
use crate::refine_sv;
use crate::run_stats::{DiscoverRunStats, RunStep, delete_run_stats, write_discover_run_stats};
use crate::worker_thread_data::get_bam_reader_worker_thread_data;

pub const ASSEMBLY_REGIONS_FILENAME: &str = "assembly.regions.bed";
pub const CANDIDATE_SV_FILENAME: &str = "candidate.sv.bcf";
pub const CONTIG_ALIGNMENT_FILENAME: &str = "contig.alignment.bam";
pub const COPYNUM_SEGMENT_BEDGRAPH_FILENAME: &str = "copynum.bedgraph";
pub const COPYNUM_SEGMENT_MESSAGEPACK_FILENAME: &str = "copynum.mpack";
pub const DEPTH_BINS_BIGWIG_FILENAME: &str = "depth.bw";
pub const DEPTH_BINS_MESSAGEPACK_FILENAME: &str = "depth.mpack";
pub const EXPECTED_COPY_NUMBER_BED_FILENAME: &str = "expected.copy.number.bed";
pub const GENOME_GC_LEVELS_MESSAGEPACK_FILENAME: &str = "genome.gclevels.mpack";
pub const GC_BIAS_CORRECTED_DEPTH_BINS_BIGWIG_FILENAME: &str = "gc_bias_corrected_depth.bw";
pub const MAF_BIGWIG_FILENAME: &str = "maf.bw";
pub const MAF_MESSAGEPACK_FILENAME: &str = "maf.mpack";
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
        let maf_data = scan_maf_file(maf_filename, &chrom_list, &[&sample_name]);
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

    // Write max_sv_depth to a BED file
    {
        let max_sv_depth_bed = settings.output_dir.join(MAX_SV_DEPTH_FILENAME);
        write_genome_regions_to_bed(
            "max sv depth",
            &max_sv_depth_bed,
            &chrom_list,
            &sample_scan_result.genome_max_sv_depth_regions,
        );
    }

    if !settings.disable_cnv {
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

        serialize_gc_bias_correction_data(&settings.output_dir, &gc_bias_data);

        if settings.debug_gc_correction {
            write_gc_correction_debug_output(
                &settings.output_dir,
                &chrom_list,
                &depth_bin_gc_content,
                &sample_scan_result,
                &gc_bias_data,
            );
        }

        let copy_number_segments = get_sample_copy_number_segments(
            settings,
            &chrom_list,
            &sample_scan_result,
            &gc_bias_data,
        );

        write_copy_number_segment_file(&settings.output_dir, &chrom_list, &copy_number_segments);

        serialize_copy_number_segments(&settings.output_dir, &copy_number_segments);
    }

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

    // In addition to useful statistics this file acts as a marker for a successfully completed run, so it must be written last.
    write_discover_run_stats(
        &settings.output_dir,
        &DiscoverRunStats {
            run_step: RunStep {
                name: "discover".to_string(),
                version: PROGRAM_VERSION.to_string(),
            },
            sample_name: sample_scan_result.sample_name.clone(),
            cluster_stats,
            refine_stats,
        },
    );
}
