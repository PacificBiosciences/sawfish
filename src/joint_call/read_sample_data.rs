use std::path::Path;
use std::sync::mpsc::channel;

use log::info;
use rust_vc_utils::{get_genome_ref_from_fasta, ChromList, GenomeRef};

use super::get_refined_svs::get_sample_sv_groups;

use crate::cli::{read_discover_settings, DiscoverSettings, JointCallSettings, SharedSettings};
use crate::depth_bins::{deserialize_genome_depth_bins, GenomeDepthBins};
use crate::discover;
use crate::discover::EXPECTED_COPY_NUMBER_FILENAME;
use crate::gc_correction::{
    deserialize_genome_gc_levels, deserialize_sample_gc_bias_data, GenomeGCLevels,
    SampleGCBiasCorrectionData,
};
use crate::genome_regions::{
    read_genome_regions_from_bed, ChromRegions, GenomeRegions, GenomeRegionsByChromIndex,
};
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
    pub sample_name: String,
    pub discover_settings: DiscoverSettings,
    pub genome_max_sv_depth_regions: Vec<ChromRegions>,
    pub sv_groups: Vec<SVGroup>,

    /// Vector indexed on cluster index, containing the assembly regions for that cluster
    //pub assembly_regions: Vec<Vec<GenomeSegment>>,
    pub genome_depth_bins: GenomeDepthBins,
    pub sample_gc_bias_data: SampleGCBiasCorrectionData,

    pub expected_copy_number_regions: Option<GenomeRegionsByChromIndex>,
}

impl SampleJointCallData {
    pub fn to_sample_score_data(&self) -> SampleScoreData {
        SampleScoreData {
            genome_max_sv_depth_regions: &self.genome_max_sv_depth_regions,
            genome_depth_bins: &self.genome_depth_bins,
            gc_depth_correction: self.sample_gc_bias_data.to_gc_depth_correction(),
        }
    }
}

fn read_expected_copy_number_regions(
    chrom_list: &ChromList,
    discover_dir: &Path,
) -> Option<GenomeRegionsByChromIndex> {
    let filename = discover_dir.join(EXPECTED_COPY_NUMBER_FILENAME);
    if filename.exists() {
        let chroms =
            read_genome_regions_from_bed("expected copy number", &filename, chrom_list, true, true);
        Some(GenomeRegionsByChromIndex { chroms })
    } else {
        None
    }
}

fn get_sample_joint_call_data(
    chrom_list: &ChromList,
    target_regions: &GenomeRegions,
    disable_small_indels: bool,
    discover_dir: &Path,
    discover_settings: DiscoverSettings,
) -> SampleJointCallData {
    let debug = false;
    if debug {
        eprintln!(
            "Reading sample discovery input from '{}'",
            discover_dir.display()
        );
    }

    let disovery_run_stats = read_discover_run_stats(discover_dir);

    let genome_max_sv_depth_regions = {
        let filename = discover_dir.join(discover::MAX_SV_DEPTH_FILENAME);
        read_genome_regions_from_bed("max sv depth", &filename, chrom_list, true, false)
    };

    let expected_copy_number_regions = read_expected_copy_number_regions(chrom_list, discover_dir);

    let sv_groups = get_sample_sv_groups(
        disable_small_indels,
        discover_dir,
        chrom_list,
        target_regions,
        expected_copy_number_regions.as_ref(),
    );

    let genome_depth_bins = deserialize_genome_depth_bins(discover_dir);
    let sample_gc_bias_data = deserialize_sample_gc_bias_data(discover_dir);

    SampleJointCallData {
        sample_name: disovery_run_stats.sample_name,
        discover_settings,
        genome_max_sv_depth_regions,
        sv_groups,
        genome_depth_bins,
        sample_gc_bias_data,
        expected_copy_number_regions,
    }
}

fn get_all_sample_joint_call_data(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    chrom_list: &ChromList,
    target_regions: &GenomeRegions,
    all_discover_settings: Vec<DiscoverSettings>,
) -> Vec<SampleJointCallData> {
    info!("Reading all sample discovery input");

    // Read and process all sample-specific discovery data
    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shared_settings.thread_count)
        .build()
        .unwrap();

    let (tx, rx) = channel();

    worker_pool.scope(move |scope| {
        for (sample_index, (discover_dir, discover_settings)) in settings
            .sample
            .iter()
            .zip(all_discover_settings.into_iter())
            .enumerate()
        {
            let tx = tx.clone();
            scope.spawn(move |_| {
                let sample_data = get_sample_joint_call_data(
                    chrom_list,
                    target_regions,
                    shared_settings.disable_small_indels,
                    discover_dir,
                    discover_settings,
                );
                tx.send((sample_index, sample_data)).unwrap();
            });
        }
    });

    let mut all_sample_data = rx.into_iter().collect::<Vec<_>>();
    all_sample_data.sort_by_key(|(index, _)| *index);
    all_sample_data.into_iter().map(|(_, x)| x).collect()
}

pub(super) fn read_all_sample_data(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
) -> (SharedJointCallData, Vec<SampleJointCallData>) {
    // Open the discover mode settings for all input samples:
    let all_discover_settings = settings
        .sample
        .iter()
        .map(|x| read_discover_settings(x))
        .collect::<Vec<_>>();

    // Check that bam headers match in all input samples, and store the first chrom list
    let mut chrom_list = None;
    for discover_settings in all_discover_settings.iter() {
        let sample_chrom_list = ChromList::from_bam_filename(&discover_settings.bam_filename);
        if let Some(chrom_list) = &chrom_list {
            assert_eq!(
                chrom_list, &sample_chrom_list,
                "Sample bams do not have matching chromosome lists"
            );
        } else {
            chrom_list = Some(sample_chrom_list.clone());
        }
    }
    let chrom_list = chrom_list.unwrap();

    // To keep target region logic as simple as possible, it is only used to filter SV groups during
    // read-in, and not stored with other settings.
    let target_regions =
        GenomeRegions::from_target_regions(&chrom_list, &shared_settings.target_region_list, true);

    let ref_filename = {
        let mut ref_filename = None;
        for discover_settings in all_discover_settings.iter() {
            let sample_ref = &discover_settings.ref_filename;
            if let Some(ref_filename) = ref_filename.as_ref() {
                assert_eq!(ref_filename, sample_ref);
            } else {
                ref_filename = Some(sample_ref.clone());
            }
        }
        ref_filename.unwrap()
    };

    // Get genome reference sequence
    let genome_ref = get_genome_ref_from_fasta(&ref_filename);

    let all_sample_data = get_all_sample_joint_call_data(
        shared_settings,
        settings,
        &chrom_list,
        &target_regions,
        all_discover_settings,
    );

    // Read in genome GC levels from the first sample, since we've already checked for a
    // reference match this can be shared across all samples
    let first_sample_discover_dir = settings.sample.first().unwrap();
    let genome_gc_levels = deserialize_genome_gc_levels(first_sample_discover_dir);

    let shared_data = SharedJointCallData {
        ref_filename,
        chrom_list,
        genome_ref,
        genome_gc_levels,
    };

    (shared_data, all_sample_data)
}
