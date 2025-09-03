use std::collections::BTreeMap;
use std::fs::File;

use camino::Utf8Path;
use indexmap::IndexMap;
use rust_vc_utils::ChromList;
use serde::Serialize;
use unwrap::unwrap;

use super::sample_output::get_sample_output_dir;
use super::{SampleJointCallData, SharedJointCallData};
use crate::copy_number_segmentation::{
    CopyNumberState, SampleCopyNumberSegmentationInput, SampleCopyNumberSegments,
};
use crate::filenames::COPYNUM_SUMMARY_FILENAME;

/// Records copy number details for any region(s) of the genome
#[derive(Default, Serialize)]
struct RegionCopyNumberInfo {
    total_copy_number_bases: u64,
    bases_per_copy_number: BTreeMap<u32, u64>,
    most_common_copy_number: u32,
}

impl RegionCopyNumberInfo {
    pub fn set_most_common_copy_number(&mut self) {
        self.most_common_copy_number = self
            .bases_per_copy_number
            .iter()
            .max_by_key(|(_, v)| **v)
            .map(|(x, _)| *x)
            .unwrap_or_default();
    }
}

#[derive(Serialize)]
struct SampleCopyNumberInfo {
    sample_name: String,
    gc_bias_corrected_haploid_coverage: f64,
    chromosomes: IndexMap<String, RegionCopyNumberInfo>,
}

fn write_sample_copy_number_info(
    sample_dir: &Utf8Path,
    chrom_list: &ChromList,
    sample_seg_input: &SampleCopyNumberSegmentationInput,
    sample_data: &SampleJointCallData,
    sample_cn_segments: &SampleCopyNumberSegments,
) {
    // First create the per-chromosome component of cnv info:
    let mut chromosomes = IndexMap::new();

    let chrom_count = chrom_list.data.len();
    let bin_size = sample_cn_segments.bin_size;

    for chrom_index in 0..chrom_count {
        let chrom_label = &chrom_list.data[chrom_index].label;
        let chrom_cn_segments = &sample_cn_segments.cn_segments[chrom_index];

        let mut chrom_value = RegionCopyNumberInfo::default();
        for cn_segment in chrom_cn_segments.iter() {
            if cn_segment.copy_number_info.state == CopyNumberState::Unknown {
                continue;
            }
            let cn_segment_pos_range = cn_segment.to_range(bin_size);
            let cn_segment_size = std::cmp::max(0, cn_segment_pos_range.size()) as u64;
            chrom_value.total_copy_number_bases += cn_segment_size;
            let copy_number = cn_segment.copy_number_info.as_u32_depth();
            let map_entry = chrom_value
                .bases_per_copy_number
                .entry(copy_number)
                .or_default();
            *map_entry += cn_segment_size;
        }

        if chrom_value.total_copy_number_bases == 0 {
            continue;
        }

        chrom_value.set_most_common_copy_number();

        chromosomes.insert(chrom_label.clone(), chrom_value);
    }

    // Assemble chromosome into the top-level copy number info object
    let info = SampleCopyNumberInfo {
        sample_name: sample_data.sample_name.clone(),
        gc_bias_corrected_haploid_coverage: sample_seg_input.gc_corrected_haploid_coverage,
        chromosomes,
    };

    let filename = sample_dir.join(COPYNUM_SUMMARY_FILENAME);
    //info!("Writing sample copy number info to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create sample copy number info json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &info).unwrap();
}

/// Summarize additional copy number details which don't naturally fit into any of the existing VCF or bed outputs
///
/// # Arguments
/// * `seg_input` - data on segmentation configuration for each sample with CNV enabled
/// * `cn_segments` - cn segmentation results for each sample with CNV enabled
///
pub fn write_copy_number_info(
    output_dir: &Utf8Path,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    seg_input: &[SampleCopyNumberSegmentationInput],
    cn_segments: &[SampleCopyNumberSegments],
) {
    assert_eq!(seg_input.len(), cn_segments.len());
    for (sample_seg_input, sample_cn_segments) in seg_input.iter().zip(cn_segments.iter()) {
        let sample_index = sample_seg_input.sample_index;
        let sample_data = &all_sample_data[sample_index];
        let sample_dir = get_sample_output_dir(output_dir, sample_index, &sample_data.sample_name);
        write_sample_copy_number_info(
            &sample_dir,
            &shared_data.chrom_list,
            sample_seg_input,
            sample_data,
            sample_cn_segments,
        );
    }
}
