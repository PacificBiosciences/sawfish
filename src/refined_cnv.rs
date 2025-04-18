use std::cmp::Ordering;

use strum::EnumCount;

use crate::breakpoint::VcfSVType;
use crate::genome_segment::GenomeSegment;
use crate::prob_utils::error_prob_to_phred;
use crate::refine_sv::Genotype;

#[derive(Clone, Debug)]
/// Per-sample scoring information shared on all variant types (SVs and CNVs)
pub struct SharedSampleScoreInfo {
    /// Called copy number of sample at the variant locus
    pub copy_number: Option<u32>,

    /// Copy number quality score, only defined if copy_number is
    pub copy_number_qscore: i32,

    /// Most likely genotype (per sample)
    ///
    pub gt: Option<Genotype>,

    /// Genotype quality score, only defined if gt is
    pub gt_qscore: i32,

    /// Likelihood for each genotype expressed as a quality score (following the VCF PL format), only defined if gt is
    pub gt_lhood_qscore: Vec<i32>,
}

impl Default for SharedSampleScoreInfo {
    fn default() -> Self {
        Self {
            copy_number: None,
            copy_number_qscore: 0,
            gt: None,
            gt_qscore: 0,
            gt_lhood_qscore: vec![0; Genotype::COUNT],
        }
    }
}

/// All scoring information for a CNV used for
#[derive(Clone, Default)]
pub struct CNVSampleScoreInfo {
    /// Expected copy number of sample at the CNV locus
    pub expected_copy_number: u32,

    /// Posterior probability of the expected copy number across this segment
    ///
    /// This is an intermediate value used to compute tha QUAL score, it is stored
    /// at the sample level to enable the recalculation of QUAL when merging results
    /// across samples
    ///
    /// This is only defined if the copy_number is defined for the sample
    ///
    pub expected_copy_number_prob: Option<f64>,

    /// All shared per-sample scoring information for any variant type
    pub shared: SharedSampleScoreInfo,
}

#[derive(Default)]
pub struct CNVScoreInfo {
    /// Quality score for any non-expected depth level in any sample
    pub alt_score: Option<f32>,

    /// Scores for each sample
    pub samples: Vec<CNVSampleScoreInfo>,
}

impl CNVScoreInfo {
    pub fn update_alt_score(&mut self, max_qscore: u32) {
        self.alt_score = self.get_alt_score(max_qscore);
    }

    fn get_alt_score(&mut self, max_qscore: u32) -> Option<f32> {
        let expected_copy_number_probs = self
            .samples
            .iter()
            .filter_map(|x| x.expected_copy_number_prob)
            .collect::<Vec<_>>();
        if expected_copy_number_probs.is_empty() {
            None
        } else {
            let joint_expected_copy_number_prob = expected_copy_number_probs.iter().product();
            let qscore = (error_prob_to_phred(joint_expected_copy_number_prob) as f32)
                .min(max_qscore as f32)
                .max(0.0);
            Some(qscore)
        }
    }
}

/// Copy number variants including all detailed breakpoint analysis and scoring information
///
/// This format is only for CNVs that can't be explained as a single-breakpoint event and
/// merged into an SV record.
///
/// This is the final format for CNVs in preparation for VCF output
///
pub struct RefinedCNV {
    /// Segment defining the CNV
    ///
    /// We still have a fixed notion of segment uncertainty, so no need to record it here
    pub segment: GenomeSegment,

    /// Index of the sample CNV was originally segmented from. This is only used to help
    /// generate a unique VCF ID for CNV records
    pub source_sample_index: usize,

    pub score: CNVScoreInfo,
}

impl RefinedCNV {
    /// Assign a type to the CNV segment, for the purpose of VCF alt-allele selection
    ///
    /// In a single-sample context, we just look to see whether CN is below or above the expected
    /// CN to determine if the CNV types is "Deletion" or "Duplication".
    ///
    /// In a multi-sample context, we require that (1) at least one sample has a CN value not equal to
    /// expected CN, and (2) if the CN goes both below and above the expected CN in different samples,
    /// we change the CNV type from "Deletion" or "Duplication" to the more general "CopyNumberVariant" type.
    ///
    /// Return None if no CNV found
    ///
    pub fn get_cnv_type(&self) -> Option<VcfSVType> {
        use VcfSVType::*;

        let mut del = 0;
        let mut dup = 0;
        for sample in self.score.samples.iter() {
            if let Some(x) = sample.shared.copy_number {
                match x.cmp(&sample.expected_copy_number) {
                    Ordering::Less => {
                        del += 1;
                    }
                    Ordering::Greater => {
                        dup += 1;
                    }
                    _ => (),
                }
            }
        }
        let total = del + dup;
        if total == 0 {
            None
        } else if del == total {
            Some(Deletion)
        } else if dup == total {
            Some(Duplication)
        } else {
            Some(CopyNumberVariant)
        }
    }
}
