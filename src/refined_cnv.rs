use std::cmp::Ordering;

use rust_vc_utils::GenomeSegment;
use strum::EnumCount;

use crate::breakpoint::VcfSVType;
use crate::prob_utils::error_prob_to_phred;
use crate::refine_sv::Genotype;

#[derive(Clone, Debug)]
pub struct SharedSampleCopyNumInfo {
    /// Called copy number of sample at the variant locus
    pub copy_number: u32,

    /// Copy number quality score (CNQ)
    ///
    /// This quality expresses the probability that CN is incorrect, given a naive CN prior
    pub copy_number_qscore: i32,

    /// Posterior probability of the expected copy number across this segment
    ///
    /// This is an intermediate value used to compute tha QUAL score, it is stored at the sample level to enable the
    /// recalculation of QUAL when merging results across samples
    ///
    pub expected_copy_number_prob: f64,
}

#[derive(Clone, Debug)]
/// Per-sample scoring information shared on all variant types (SVs and CNVs)
pub struct SharedSampleScoreInfo {
    /// Copy number information for CNV or merged SV/CNV events
    pub copy_info: Option<SharedSampleCopyNumInfo>,

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
            copy_info: None,
            gt: None,
            gt_qscore: 0,
            gt_lhood_qscore: vec![0; Genotype::COUNT],
        }
    }
}

pub trait SharedVariantScoreInfo {
    fn get_shared_sample_iter(&self) -> impl Iterator<Item = &SharedSampleScoreInfo>;

    fn get_cnv_alt_score(&mut self, max_qscore: u32) -> Option<f32> {
        let expected_copy_number_probs = self
            .get_shared_sample_iter()
            .filter_map(|x| x.copy_info.as_ref().map(|y| y.expected_copy_number_prob))
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

/// All scoring information for a CNV used for
#[derive(Clone, Default)]
pub struct CNVSampleScoreInfo {
    /// Expected copy number of sample at the CNV locus
    pub expected_copy_number: u32,

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
    pub fn update_cnv_alt_score(&mut self, max_qscore: u32) {
        self.alt_score = self.get_cnv_alt_score(max_qscore);
    }
}

impl SharedVariantScoreInfo for CNVScoreInfo {
    fn get_shared_sample_iter(&self) -> impl Iterator<Item = &SharedSampleScoreInfo> {
        self.samples.iter().map(|x| &x.shared)
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
            if let Some(x) = &sample.shared.copy_info {
                match x.copy_number.cmp(&sample.expected_copy_number) {
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
