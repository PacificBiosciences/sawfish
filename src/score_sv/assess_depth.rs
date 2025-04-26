use crate::breakpoint::{get_breakpoint_vcf_sv_type, Breakpoint, VcfSVType};

use crate::refine_sv::RefinedSV;

use super::{SampleScoreData, ScoreSVSettings};
use crate::depth_bins::DepthBin;
use crate::gc_correction::GenomeGCLevels;
use crate::score_sv::Genotype;

/// Return true if the deletion or duplication is not supported by depth in this sample
///
fn reject_deldup_depth_support_in_sample(
    genome_gc_levels: &GenomeGCLevels,
    sample_data: &SampleScoreData,
    bp: &Breakpoint,
    sv_type: VcfSVType,
) -> bool {
    // This stuff is constant over all samples, but no need to precompute:
    let chrom_index = bp.breakend1.segment.chrom_index;
    let chrom_gc_levels = &genome_gc_levels[chrom_index];

    // Start per-sample logic:
    let genome_depth_bins = &sample_data.genome_depth_bins;
    let chrom_depth_bins = &genome_depth_bins.chroms[chrom_index];
    let chrom_bin_count = chrom_depth_bins.len() as i64;

    let start_bin_index = bp.breakend1.segment.range.start / genome_depth_bins.bin_size as i64;
    let end_bin_index =
        bp.breakend2.as_ref().unwrap().segment.range.end / genome_depth_bins.bin_size as i64;

    let gc_depth_correction = &sample_data.gc_depth_correction;

    // Sum depth observations in the interior of the deletion or duplication
    let interior_depth = {
        let mut sum = 0.0;
        let mut count = 0;
        let mut excluded = 0;
        for bin_index in start_bin_index + 2..end_bin_index - 1 {
            let depth_bin = &chrom_depth_bins[bin_index as usize];
            match depth_bin {
                DepthBin::Depth(depth) => {
                    let gc_level = chrom_gc_levels[bin_index as usize];
                    sum += *depth * gc_depth_correction[gc_level];
                    count += 1;
                }
                DepthBin::Excluded => {
                    excluded += 1;
                }
            }
        }

        //eprintln!("Getting interior depth. count {count} excluded {excluded} sum {sum}");

        // Very simple quality heuristics:
        //
        if count < 10 || excluded > count {
            None
        } else {
            Some(sum / count as f64)
        }
    };

    if interior_depth.is_none() {
        return false;
    }

    // Sum depth observations on the flanks of the deletion or duplication
    // TODO should we just pull the haploid depth estimate instead?
    let flank_depth = {
        let flank_bin_count = 3;
        let left_flank_start = start_bin_index - (1 + flank_bin_count);
        let right_flank_end = end_bin_index + (1 + flank_bin_count);

        let mut sum = 0.0;
        let mut count = 0;
        let mut excluded = 0;
        for bin_index in left_flank_start..start_bin_index - 1 {
            if bin_index < 0 {
                excluded += 1;
                continue;
            }

            let bin = &chrom_depth_bins[bin_index as usize];
            match bin {
                DepthBin::Depth(depth) => {
                    let gc_level = chrom_gc_levels[bin_index as usize];
                    sum += *depth * gc_depth_correction[gc_level];
                    count += 1;
                }
                DepthBin::Excluded => {
                    excluded += 1;
                }
            }
        }

        for bin_index in end_bin_index + 1..right_flank_end {
            if bin_index >= chrom_bin_count {
                excluded += 1;
                continue;
            }

            let bin = &chrom_depth_bins[bin_index as usize];
            match bin {
                DepthBin::Depth(depth) => {
                    let gc_level = chrom_gc_levels[bin_index as usize];
                    sum += *depth * gc_depth_correction[gc_level];
                    count += 1;
                }
                DepthBin::Excluded => {
                    excluded += 1;
                }
            }
        }

        //eprintln!("Getting flank depth. count {count} excluded {excluded} sum {sum}");

        // Very simple quality heuristics:
        //
        if count < flank_bin_count || excluded > count {
            None
        } else {
            Some(sum / count as f64)
        }
    };

    // Now apply a simple cutoff for the ratio of interior to flank depth:
    if let (Some(f), Some(i)) = (flank_depth, interior_depth) {
        match sv_type {
            VcfSVType::Deletion => f > 0.0 && (i / f) > 0.8,
            VcfSVType::Duplication => f > 0.0 && (i / f) < 1.2,
            _ => false,
        }
    } else {
        false
    }
}

/// Return None if this sample doesn't get a vote (SV breakpoint not present), else return
/// true if the deletion or duplication is not supported by depth in this sample
///
fn optional_reject_deldup_depth_support_in_sample(
    genome_gc_levels: &GenomeGCLevels,
    all_sample_data: &[SampleScoreData],
    refined_sv: &RefinedSV,
    sv_type: VcfSVType,
    sample_index: usize,
) -> Option<bool> {
    // First determine if this sample is eligible for voting
    //
    // Require that breakpoint evidence is strong enough induce a non-REF genotype OR that that this
    // is sample nominating the candidate
    //
    let is_eligible = if refined_sv.id.sample_index == sample_index {
        true
    } else {
        let gt = &refined_sv.score.samples[sample_index].gt;
        if let Some(gt) = gt {
            *gt != Genotype::Ref
        } else {
            false
        }
    };

    // If eligible determine rejection status
    if is_eligible {
        Some(reject_deldup_depth_support_in_sample(
            genome_gc_levels,
            &all_sample_data[sample_index],
            &refined_sv.bp,
            sv_type,
        ))
    } else {
        None
    }
}

/// Return true if the deletion or duplication has valid depth information but fails to show depth support evidence
/// in the majority of samples where the genotype is not REF
///
fn reject_deldup_depth_support(
    genome_gc_levels: &GenomeGCLevels,
    all_sample_data: &[SampleScoreData],
    refined_sv: &RefinedSV,
    sv_type: VcfSVType,
) -> bool {
    let sample_count = all_sample_data.len();

    let mut total = 0;
    let mut reject = 0;
    for sample_index in 0..sample_count {
        let is_reject = optional_reject_deldup_depth_support_in_sample(
            genome_gc_levels,
            all_sample_data,
            refined_sv,
            sv_type,
            sample_index,
        );
        if let Some(is_reject) = is_reject {
            total += 1;
            if is_reject {
                reject += 1;
            }
        }
    }
    assert_ne!(total, 0);
    (reject as f64 / total as f64) >= 0.5
}

/// Return true for any SV subject to depth assessment
fn is_depth_assessed_variant_type(
    score_settings: &ScoreSVSettings,
    refined_sv: &RefinedSV,
) -> bool {
    let bp = &refined_sv.bp;
    let sv_type = get_breakpoint_vcf_sv_type(bp);
    match sv_type {
        VcfSVType::Deletion | VcfSVType::Duplication => {
            let be1 = &bp.breakend1;
            let be2 = bp.breakend2.as_ref().unwrap();
            let sv_size =
                std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) as usize;

            sv_size >= score_settings.min_sv_depth_assess_size
        }
        _ => false,
    }
}

/// Evaluate large deletions and duplications to see if these include a consistent depth signature
///
pub(super) fn assess_sv_depth_support(
    score_settings: &ScoreSVSettings,
    genome_gc_levels: &GenomeGCLevels,
    all_sample_data: &[SampleScoreData],
    refined_svs: &mut [RefinedSV],
) {
    let debug = false;

    for refined_sv in refined_svs.iter_mut().filter(|x| !x.filter_sv()) {
        if debug {
            eprintln!(
                "Starting depth assessment process for SV id: {:?}",
                refined_sv.id
            );
        }

        // Don't run assessment if the SV has been forced to BND format already
        if refined_sv.ext.force_breakpoint_representation {
            continue;
        }

        // Determine if the refined SV pattern is a candidate for depth assessment:
        if !is_depth_assessed_variant_type(score_settings, refined_sv) {
            continue;
        }

        let sv_type = get_breakpoint_vcf_sv_type(&refined_sv.bp);
        refined_sv.ext.force_breakpoint_representation = match sv_type {
            VcfSVType::Deletion | VcfSVType::Duplication => {
                reject_deldup_depth_support(genome_gc_levels, all_sample_data, refined_sv, sv_type)
            }
            _ => false,
        };
    }
}
