use crate::genome_regions::GenomeRegionsByChromIndex;
use crate::genome_segment::GenomeSegment;

/// A simple copy number determination routine given a global expected copy number region map, and
/// a small set of focal regions (such as the breakpoints of an SV
///
/// A very simple method is used here, there must be a single copy number in all regions, or else
/// the default copy number of 2 is returned
///
fn get_expected_copy_number_for_regions(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_regions: &[GenomeSegment],
) -> u8 {
    let default_copy_number = 2;
    if expected_copy_number_regions.is_none() {
        return default_copy_number;
    }
    let expected_copy_number_regions = expected_copy_number_regions.unwrap();

    let mut all_target_cns = Vec::new();
    for target_region in target_regions {
        let mut target_cns = expected_copy_number_regions
            .find_overlaps(
                target_region.chrom_index,
                target_region.range.start,
                target_region.range.end,
            )
            .map(|x| *x.data())
            .collect::<Vec<_>>();
        if target_cns.is_empty() {
            return default_copy_number;
        }
        all_target_cns.append(&mut target_cns);
    }
    if all_target_cns.is_empty() {
        return default_copy_number;
    }
    let copy_number = all_target_cns[0];
    if all_target_cns.iter().all(|&x| x == copy_number) {
        copy_number
    } else {
        default_copy_number
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SVLocusPloidy {
    Haploid,
    Diploid,
}

/// Convert expected copy number to expected ploidy
///
pub fn get_expected_ploidy_for_regions(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_regions: &[GenomeSegment],
) -> SVLocusPloidy {
    let expected_copy_number =
        get_expected_copy_number_for_regions(expected_copy_number_regions, target_regions);
    if expected_copy_number < 2 {
        SVLocusPloidy::Haploid
    } else {
        SVLocusPloidy::Diploid
    }
}

/// Use sawfish 2 vs 1 region SV type rules to resolve the maximum allowed haplotypes given the
/// expected copy number regions for one sample.
///
/// Returns a 2-tuple of (max_haplotype_count, ploidy)
///
pub fn get_max_haplotype_count_for_regions(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_regions: &[GenomeSegment],
) -> (usize, SVLocusPloidy) {
    let ploidy = get_expected_ploidy_for_regions(expected_copy_number_regions, target_regions);
    let max_haplotype_count = if target_regions.len() == 1 {
        // Smaller SVs tend to occur in regions where there may be overlapping SVs/indels, for this
        // reason we assemble up to the expected sample ploidy in the region
        match ploidy {
            SVLocusPloidy::Haploid => 1,
            SVLocusPloidy::Diploid => 2,
        }
    } else {
        // The approach for large SVs is to just assemble the SV allele, and assume there is not a
        // second large SV allele with an overlapping signature, which is very unlikely. For this
        // reason we always limit the analysis to a single haplotype assembly
        1
    };

    (max_haplotype_count, ploidy)
}
