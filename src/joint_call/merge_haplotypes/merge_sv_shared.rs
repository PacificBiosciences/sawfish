//! Shared haplotype merging code for both single-region and multi-region cases
//!

use std::collections::BTreeMap;

use crate::expected_ploidy::get_max_haplotype_count_for_regions;
use crate::joint_call::SampleJointCallData;
use crate::sv_group::SVGroup;

/// Fill in the expected copy-number related values for sv_group
///
pub(super) fn set_sv_group_expected_cn_info(
    treat_single_copy_as_haploid: bool,
    all_sample_data: &[SampleJointCallData],
    sv_groups: &mut [SVGroup],
) {
    for sv_group in sv_groups.iter_mut() {
        sv_group.sample_expected_cn_info.clear();
        for (sample_index, sample_data) in all_sample_data.iter().enumerate() {
            let (max_haplotype_count, expected_cn_info) = get_max_haplotype_count_for_regions(
                treat_single_copy_as_haploid,
                sample_data.expected_copy_number_regions.as_ref(),
                &sv_group.group_regions,
            );
            sv_group.sample_expected_cn_info.push(expected_cn_info);

            // Trim sample_haplotype_list to be no more than max_haplotype_count if necessary.
            //
            // Technically this shouldn't be needed here, because the genotyping routine should account for
            // haplotype count to, but seems good to be consistent about handling throught the full pipe.
            sv_group.sample_haplotype_list[sample_index].truncate(max_haplotype_count);
        }
    }
}

#[derive(Debug)]
pub(super) struct AnnoSVGroup<'a> {
    pub sv_group: &'a SVGroup,
    pub sample_index: usize,
}

/// Reformat single-sample SV group into multi-sample format
///
pub(super) fn clone_sv_group_as_multisample(
    sample_count: usize,
    input_sv_group_info: &AnnoSVGroup,
) -> SVGroup {
    let mut multisample_sv_group = input_sv_group_info.sv_group.clone();
    if sample_count != 1 {
        // The single sample sv_group input has all information set to a sample_index of 0.
        // The following steps move the sample index to its new value for all sample_index
        // annotations within SVGroup
        //
        let haplotype_source_sample_index = input_sv_group_info.sample_index;

        // 1. Move sample_haplotypes to new sample_index value
        {
            let sample_haplotypes = multisample_sv_group.sample_haplotype_list.pop().unwrap();
            multisample_sv_group.sample_haplotype_list = vec![Vec::new(); sample_count];
            multisample_sv_group.sample_haplotype_list[haplotype_source_sample_index] =
                sample_haplotypes;
        }

        // 2. Move SV ids to new sample_index value
        for rsv in multisample_sv_group.refined_svs.iter_mut() {
            rsv.id.sample_index = haplotype_source_sample_index;
        }

        // 3. Move group haplotypes to the new sample_index value
        for group_haplotype in multisample_sv_group.group_haplotypes.iter_mut() {
            group_haplotype.hap_id.sample_index = haplotype_source_sample_index;
        }

        // 4. Move neighbors to new sample_index value:
        for rsv in multisample_sv_group.refined_svs.iter_mut() {
            if let Some(x) = rsv.bp.breakend1_neighbor.as_mut() {
                x.sample_index = haplotype_source_sample_index;
            }
            if let Some(x) = rsv.bp.breakend2_neighbor.as_mut() {
                x.sample_index = haplotype_source_sample_index;
            }
        }
    }
    multisample_sv_group
}

#[derive(Default, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct MultiSampleClusterId {
    pub sample_index: usize,
    pub cluster_index: usize,
}

/// Structure to keep track of all multi-region merges:
#[derive(Default)]
pub(super) struct ClusterMergeMap {
    /// Any cluster that gets merged is mapped to its new sample and cluster id here:
    pub data: BTreeMap<MultiSampleClusterId, MultiSampleClusterId>,
}

impl ClusterMergeMap {
    pub fn extend(&mut self, other: Self) {
        self.data.extend(other.data);
    }
}
