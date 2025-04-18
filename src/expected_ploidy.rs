use strum::EnumCount;

use crate::copy_number_segmentation::CopyNumberState;
use crate::genome_regions::GenomeRegionsByChromIndex;
use crate::genome_segment::GenomeSegment;

/// Provide a copy number for the given target region
///
/// If the target region overlaps multiple expected copy number regions, this routine will
/// provide the one covering the majority of the region.
///
/// This works well for typical large CNV needs.
///
pub fn get_majority_expected_copy_number_for_region(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_region: &GenomeSegment,
) -> u8 {
    let default_copy_number = 2;
    if expected_copy_number_regions.is_none() {
        return default_copy_number;
    }
    let expected_copy_number_regions = expected_copy_number_regions.unwrap();

    let cn_state_count = CopyNumberState::COUNT;
    let mut cn_span = vec![0; cn_state_count];
    for ecn_region in expected_copy_number_regions.find_overlaps(
        target_region.chrom_index,
        target_region.range.start,
        target_region.range.end,
    ) {
        let span = ecn_region.interval().end - ecn_region.interval().start;
        let ecn = *ecn_region.data() as usize;
        let ecn = std::cmp::min(ecn, CopyNumberState::High as usize);
        cn_span[ecn] += span;
    }

    let mut max_cn_index = 0;
    for cn_index in 1..cn_state_count {
        if cn_span[cn_index] > cn_span[max_cn_index] {
            max_cn_index = cn_index;
        }
    }

    if cn_span[max_cn_index] == 0 {
        default_copy_number
    } else {
        max_cn_index as u8
    }
}

/// A simple copy number determination routine given a global expected copy number region map, and
/// a small set of focal regions (such as the breakpoints of an SV)
///
/// A very simple method is used here, there must be a single copy number in all regions, or else
/// the default copy number of 2 is returned. This works well for small SVs
///
fn get_simple_expected_copy_number_for_regions(
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

/// Consolidate the combination of expected copy number and ploidy for the SV locus in one structure
///
/// For SV genotyping, we simplify the notion of 'ploidy' to be a relatively conservative mapping
/// from copy number. The only real use case we currently support for expected copy number is human
/// sex chromosomes, so the relatiionship can be enriched as needed if other applications arise.
///
#[derive(Clone, Copy, Debug)]
pub struct SVLocusExpectedCNInfo {
    pub expected_copy_number: u8,
}

impl SVLocusExpectedCNInfo {
    /// Return the ploidy to use for SV calling as a function of expected copy number
    ///
    /// # Arguments
    /// * `treat_single_copy_as_haploid` - It true, treat expected single-copy regions as haploid
    ///
    pub fn ploidy(&self, treat_single_copy_as_haploid: bool) -> SVLocusPloidy {
        if treat_single_copy_as_haploid && self.expected_copy_number == 1 {
            SVLocusPloidy::Haploid
        } else {
            SVLocusPloidy::Diploid
        }
    }
}

/// Get expected copy number and convert it to expected ploidy
///
pub fn get_expected_copy_number_info_for_regions(
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_regions: &[GenomeSegment],
) -> SVLocusExpectedCNInfo {
    let expected_copy_number =
        get_simple_expected_copy_number_for_regions(expected_copy_number_regions, target_regions);
    SVLocusExpectedCNInfo {
        expected_copy_number,
    }
}

/// Use sawfish 2 vs 1 region SV type rules to resolve the maximum allowed haplotypes given the
/// expected copy number regions for one sample.
///
/// Returns a 2-tuple of (max_haplotype_count, expected_copy_number_info)
///
/// # Arguments
/// * `treat_single_copy_as_haploid` - It true, treat expected single-copy regions as haploid
///
pub fn get_max_haplotype_count_for_regions(
    treat_single_copy_as_haploid: bool,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    target_regions: &[GenomeSegment],
) -> (usize, SVLocusExpectedCNInfo) {
    let expected_cn_info =
        get_expected_copy_number_info_for_regions(expected_copy_number_regions, target_regions);
    let max_haplotype_count = if target_regions.len() == 1 {
        // Smaller SVs tend to occur in regions where there may be overlapping SVs/indels, for this
        // reason we assemble up to the expected sample ploidy in the region
        match expected_cn_info.ploidy(treat_single_copy_as_haploid) {
            SVLocusPloidy::Haploid => 1,
            SVLocusPloidy::Diploid => 2,
        }
    } else {
        // The approach for large SVs is to just assemble the SV allele, and assume there is not a
        // second large SV allele with an overlapping signature, which is very unlikely. For this
        // reason we always limit the analysis to a single haplotype assembly
        1
    };

    (max_haplotype_count, expected_cn_info)
}
