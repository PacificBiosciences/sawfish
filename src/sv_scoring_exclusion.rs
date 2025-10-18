use rust_vc_utils::{ChromList, RegionMap};

use crate::cli;
use crate::copy_number_segmentation::HaploidCoverage;
use crate::depth_bins;
use crate::filenames::MAX_SV_DEPTH_FILENAME;
use crate::genome_regions::write_genome_regions_to_bed;

fn get_max_sv_scoring_depth(
    max_copy_number_for_sv_scoring: u32,
    haploid_coverage: Option<&HaploidCoverage>,
) -> u32 {
    // max_sv_scoring_depth is typically found from a multiple of the haploid_coverage, but the following limit is
    // applied to this value in all cases, and used when CNV is disabled and no haploid coverage estimate is available.
    //
    let max_sv_scoring_depth_limit = 1_000u32;

    if let Some(haploid_coverage) = haploid_coverage {
        let sample_max_depth =
            (haploid_coverage.gc_corrected_depth * max_copy_number_for_sv_scoring as f64) as u32;
        std::cmp::min(sample_max_depth, max_sv_scoring_depth_limit)
    } else {
        max_sv_scoring_depth_limit
    }
}

/// Return a 2-tuple of:
/// 1. The max sv scoring depth
/// 2. The regions of sv scoring exclusion, derived from this max depth
///
fn get_genome_sv_scoring_exclusion_regions(
    settings: &cli::DiscoverSettings,
    chrom_list: &ChromList,
    haploid_coverage: Option<&HaploidCoverage>,
    genome_depth_bins_builder: Vec<depth_bins::DepthBinsBuilder>,
) -> (u32, Vec<RegionMap>) {
    let max_sv_scoring_depth =
        get_max_sv_scoring_depth(settings.max_copy_number_for_sv_scoring, haploid_coverage);

    let regions = genome_depth_bins_builder
        .into_iter()
        .zip(chrom_list.data.iter())
        .map(|(chrom_depth_bins_builder, chrom_info)| {
            chrom_depth_bins_builder
                .get_regions_over_max_depth(chrom_info.length, max_sv_scoring_depth)
        })
        .collect();

    (max_sv_scoring_depth, regions)
}

/// Generate genome-wide regions where SV scoring will be excluded in downstream analysis,and write these
/// regions to a BED file for use in the joint-call step
///
/// Returns the max sv scoring depth
///
pub fn generate_and_write_sv_scoring_exclusion_regions(
    settings: &cli::DiscoverSettings,
    chrom_list: &ChromList,
    haploid_coverage: Option<&HaploidCoverage>,
    genome_depth_bins_builder: Vec<depth_bins::DepthBinsBuilder>,
) -> u32 {
    let (max_sv_scoring_depth, genome_sv_scoring_exclusion_regions) =
        get_genome_sv_scoring_exclusion_regions(
            settings,
            chrom_list,
            haploid_coverage,
            genome_depth_bins_builder,
        );

    let max_sv_depth_bed = settings.output_dir.join(MAX_SV_DEPTH_FILENAME);
    write_genome_regions_to_bed(
        "max sv depth",
        &max_sv_depth_bed,
        chrom_list,
        &genome_sv_scoring_exclusion_regions,
    );

    max_sv_scoring_depth
}
