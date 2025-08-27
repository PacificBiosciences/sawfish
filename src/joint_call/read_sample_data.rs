use std::sync::mpsc::channel;

use camino::{Utf8Path, Utf8PathBuf};
use log::info;
use rust_vc_utils::{ChromList, GenomeRef, get_genome_ref_from_fasta};
use unwrap::unwrap;

use super::get_refined_svs::get_sample_sv_groups;

use crate::cli::{
    DiscoverSettings, InputSampleData, JointCallDerivedSettings, JointCallSettings, SharedSettings,
    read_discover_settings,
};
use crate::copy_number_segmentation::{SampleCopyNumberSegments, deserialize_copy_number_segments};
use crate::depth_bins::{GenomeDepthBins, deserialize_genome_depth_bins};
use crate::discover;
use crate::discover::EXPECTED_COPY_NUMBER_BED_FILENAME;
use crate::gc_correction::{
    GenomeGCLevels, SampleGCBiasCorrectionData, deserialize_genome_gc_levels,
    deserialize_sample_gc_bias_data,
};
use crate::genome_regions::{
    ChromRegions, GenomeRegions, GenomeRegionsByChromIndex, read_genome_regions_from_bed,
};
use crate::maf_utils::{MafData, deserialize_maf_data};
use crate::run_stats::read_discover_run_stats;
use crate::score_sv::SampleScoreData;
use crate::sv_group::SVGroup;

pub(super) struct SharedJointCallData {
    pub ref_filename: String,
    pub chrom_list: ChromList,
    pub genome_ref: GenomeRef,
    pub genome_gc_levels: GenomeGCLevels,
}

pub struct SampleJointCallData {
    /// This is the actual path given on the joint-call command-line, not the one in the discover_settings json file
    pub discover_dir: Utf8PathBuf,

    /// Sample's order within the joint-call run
    pub sample_index: usize,

    /// The bam filename to use for joint-call, this may be different than the path from discover_settings.
    ///
    /// The path in discover settings should not be directly used except possibly by the methods initializing this struct
    pub bam_filename: String,

    pub sample_name: String,
    pub discover_settings: DiscoverSettings,
    pub genome_max_sv_depth_regions: Vec<ChromRegions>,
    pub sv_groups: Vec<SVGroup>,

    /// Vector indexed on cluster index, containing the assembly regions for that cluster
    //pub assembly_regions: Vec<Vec<GenomeSegment>>,
    pub genome_depth_bins: GenomeDepthBins,
    pub sample_gc_bias_data: Option<SampleGCBiasCorrectionData>,
    pub copy_number_segments: SampleCopyNumberSegments,

    pub expected_copy_number_regions: Option<GenomeRegionsByChromIndex>,

    /// Minor allele frequency data for this sample
    pub maf_data: Option<MafData>,
}

impl SampleJointCallData {
    pub fn to_sample_score_data(&self) -> SampleScoreData<'_> {
        SampleScoreData {
            genome_max_sv_depth_regions: &self.genome_max_sv_depth_regions,
        }
    }
}

fn read_expected_copy_number_regions(
    chrom_list: &ChromList,
    discover_dir: &Utf8Path,
) -> Option<GenomeRegionsByChromIndex> {
    let filename = discover_dir.join(EXPECTED_COPY_NUMBER_BED_FILENAME);
    if filename.exists() {
        let chroms =
            read_genome_regions_from_bed("expected copy number", &filename, chrom_list, true, true);
        Some(GenomeRegionsByChromIndex { chroms })
    } else {
        None
    }
}

fn get_sample_joint_call_data(
    settings: &JointCallSettings,
    sample_index: usize,
    chrom_list: &ChromList,
    target_regions: &GenomeRegions,
    disable_small_indels: bool,
    treat_single_copy_as_haploid: bool,
    early_sample_data: EarlySampleData,
) -> SampleJointCallData {
    let debug = false;
    if debug {
        eprintln!(
            "Reading sample discovery input from '{}'",
            early_sample_data.discover_dir
        );
    }

    let EarlySampleData {
        discover_dir,
        discover_settings,
        bam_filename,
    } = early_sample_data;

    let disable_cnv = settings.disable_cnv || discover_settings.disable_cnv;

    let discover_run_stats = unwrap!(read_discover_run_stats(&discover_dir));

    let genome_max_sv_depth_regions = {
        let filename = discover_dir.join(discover::MAX_SV_DEPTH_FILENAME);
        read_genome_regions_from_bed("max sv depth", &filename, chrom_list, true, false)
    };

    let expected_copy_number_regions = read_expected_copy_number_regions(chrom_list, &discover_dir);

    let sv_groups = get_sample_sv_groups(
        disable_small_indels,
        treat_single_copy_as_haploid,
        &discover_dir,
        chrom_list,
        target_regions,
        expected_copy_number_regions.as_ref(),
    );

    let genome_depth_bins = deserialize_genome_depth_bins(&discover_dir);

    let sample_gc_bias_data = if !disable_cnv {
        Some(deserialize_sample_gc_bias_data(&discover_dir))
    } else {
        None
    };

    let copy_number_segments = if !disable_cnv {
        deserialize_copy_number_segments(&discover_dir)
    } else {
        SampleCopyNumberSegments::new(genome_depth_bins.bin_size, chrom_list)
    };

    let maf_data = discover_settings
        .maf_filename
        .as_ref()
        .map(|_| deserialize_maf_data(&discover_dir));

    SampleJointCallData {
        discover_dir,
        sample_index,
        bam_filename,
        sample_name: discover_run_stats.sample_name,
        discover_settings,
        genome_max_sv_depth_regions,
        sv_groups,
        genome_depth_bins,
        sample_gc_bias_data,
        copy_number_segments,
        expected_copy_number_regions,
        maf_data,
    }
}

fn get_all_sample_joint_call_data(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    chrom_list: &ChromList,
    target_regions: &GenomeRegions,
    all_early_sample_data: Vec<EarlySampleData>,
) -> Vec<SampleJointCallData> {
    // Read and process all sample-specific discovery data
    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();

    let (tx, rx) = channel();

    worker_pool.scope(move |scope| {
        for (sample_index, early_sample_data) in all_early_sample_data.into_iter().enumerate() {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let sample_data = get_sample_joint_call_data(
                    settings,
                    sample_index,
                    chrom_list,
                    target_regions,
                    shared_settings.disable_small_indels,
                    settings.treat_single_copy_as_haploid,
                    early_sample_data,
                );
                tx.send(sample_data).unwrap();
            });
        }
    });

    let mut all_sample_data = rx.into_iter().collect::<Vec<_>>();
    all_sample_data.sort_by_key(|x| x.sample_index);
    all_sample_data
}

struct EarlySampleData {
    discover_dir: Utf8PathBuf,
    discover_settings: DiscoverSettings,
    bam_filename: String,
}

/// Get per-sample information that we need to aquire before the main sample information loop
///
fn get_all_early_sample_data(all_input_sample_data: &[InputSampleData]) -> Vec<EarlySampleData> {
    let mut x = Vec::new();
    for input_sample_data in all_input_sample_data.iter() {
        let discover_settings = read_discover_settings(&input_sample_data.discover_dir);
        let bam_filename = match input_sample_data.bam_filename.clone() {
            Some(x) => x,
            None => discover_settings.bam_filename.clone(),
        }
        .clone();

        x.push(EarlySampleData {
            discover_dir: input_sample_data.discover_dir.clone(),
            discover_settings,
            bam_filename,
        });
    }
    x
}

/// Check that bam headers match in all input samples, and return the first chrom list as the shared reference copy to be used for all samples
///
fn get_multi_sample_chrom_list(all_early_sample_data: &[EarlySampleData]) -> ChromList {
    let mut chrom_list = None;
    for bam_filename in all_early_sample_data.iter().map(|x| &x.bam_filename) {
        let sample_chrom_list = ChromList::from_bam_filename(bam_filename);
        if let Some(chrom_list) = &chrom_list {
            assert_eq!(
                chrom_list, &sample_chrom_list,
                "Input sample bams do not have matching chromosome lists"
            );
        } else {
            chrom_list = Some(sample_chrom_list.clone());
        }
    }
    chrom_list.unwrap()
}

/// Get the referenece filename either from joint-call settings or input sample data
///
fn get_multi_sample_ref_filename(
    settings: &JointCallSettings,
    all_sample_data: &[SampleJointCallData],
) -> String {
    if let Some(x) = settings.ref_filename.clone() {
        x
    } else {
        let mut ref_filename = None;
        for discover_settings in all_sample_data.iter().map(|x| &x.discover_settings) {
            let sample_ref = &discover_settings.ref_filename;
            if let Some(ref_filename) = ref_filename.as_ref() {
                assert_eq!(ref_filename, sample_ref);
            } else {
                ref_filename = Some(sample_ref.clone());
            }
        }
        ref_filename.unwrap()
    }
}

pub(super) fn read_all_sample_data(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    derived_settings: &JointCallDerivedSettings,
) -> (SharedJointCallData, Vec<SampleJointCallData>) {
    let all_early_sample_data = get_all_early_sample_data(&derived_settings.all_input_sample_data);

    // Now that the information is collated in early_sample_data, this is the obviuus point to log all discover_dir and
    // bam_filename
    info!(
        "Initializing joint-call process for {} input samples",
        all_early_sample_data.len()
    );
    for (sample_index, early_sample_data) in all_early_sample_data.iter().enumerate() {
        info!(
            "Input sample {} discover_dir: '{}'",
            sample_index + 1,
            early_sample_data.discover_dir.as_str(),
        );
        info!(
            "Input sample {} alignment_file: '{}'",
            sample_index + 1,
            early_sample_data.bam_filename
        );
    }

    let chrom_list = get_multi_sample_chrom_list(&all_early_sample_data);

    // To keep target region logic as simple as possible, it is only used to filter SV groups during
    // read-in, and not stored with other settings.
    let target_regions =
        GenomeRegions::from_target_regions(&chrom_list, &shared_settings.target_region_list, true);

    let all_sample_data = get_all_sample_joint_call_data(
        shared_settings,
        settings,
        &chrom_list,
        &target_regions,
        all_early_sample_data,
    );

    let ref_filename = get_multi_sample_ref_filename(settings, &all_sample_data);

    let genome_ref = {
        let mut x = get_genome_ref_from_fasta(&ref_filename);
        x.simplify_ambiguous_dna_bases();
        x
    };

    // if CNV is enabled, read in genome GC levels from the first sample with CNV enabled, else return an empty vector
    let genome_gc_levels = if settings.disable_cnv {
        Vec::new()
    } else {
        match all_sample_data
            .iter()
            .find(|x| !x.discover_settings.disable_cnv)
        {
            Some(x) => deserialize_genome_gc_levels(&x.discover_dir),
            None => Vec::new(),
        }
    };

    let shared_data = SharedJointCallData {
        ref_filename,
        chrom_list,
        genome_ref,
        genome_gc_levels,
    };

    (shared_data, all_sample_data)
}
