//! Shared haplotype merging code for both single-region and multi-region cases
//!

use crate::sv_group::SVGroup;

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
    }
    multisample_sv_group
}
