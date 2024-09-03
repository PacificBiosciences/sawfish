use std::collections::HashSet;
use std::hash::Hash;

use log::error;

use crate::expected_ploidy::SVLocusPloidy;
use crate::genome_segment::GenomeSegment;
use crate::int_range::IntRange;
use crate::refine_sv::RefinedSV;
use crate::simple_alignment::SimpleAlignment;

/// Index values are not required but useful to trace the origin of the haplotype assembly for
/// debugging:
#[derive(Clone)]
pub struct GroupHaplotypeId {
    pub sample_index: usize,
    pub cluster_index: usize,
    pub assembly_index: usize,
}

impl std::fmt::Debug for GroupHaplotypeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}",
            self.sample_index, self.cluster_index, self.assembly_index
        )
    }
}

/// Alt allele contig and alignment for one split alignment segment
///
/// For an SV with generalized breakpoints multiple contig alignments are needed, one for each split alignment segment. Note that
/// the contig_seq might be reverse complemented between the different split representations.
///
#[derive(Clone)]
pub struct ClusterAssemblyAlignment {
    pub contig_seq: Vec<u8>,

    #[allow(unused)]
    pub chrom_index: usize,

    pub contig_alignment: SimpleAlignment,

    pub high_quality_contig_range: IntRange,
}

#[derive(Clone, Default)]
pub struct ClusterAssembly {
    /// Contig alignment for each split alignment segment
    ///
    /// For single region SVs there will be only one aligment. For SVs with breakpoints, there should
    /// be (1 + breakpoint count) alignments.
    ///
    pub contig_alignments: Vec<ClusterAssemblyAlignment>,
    pub supporting_read_count: usize,
}

impl ClusterAssembly {
    pub fn single_region_refinement(&self) -> bool {
        assert!(!self.contig_alignments.is_empty());
        self.contig_alignments.len() == 1
    }
}

#[derive(Debug, Default)]
pub struct SVGroupTestStatus {
    pub are_group_regions_valid: bool,
    pub are_group_haplotypes_valid: bool,
    pub is_sample_max_haplotype_count_valid: bool,
    pub is_sample_haplotype_list_valid: bool,
    pub is_sv_haplotype_map_valid: bool,
    pub is_multi_region_sv_count_valid: bool,
}

impl SVGroupTestStatus {
    pub fn is_valid(&self) -> bool {
        self.are_group_regions_valid
            && self.are_group_haplotypes_valid
            && self.is_sample_max_haplotype_count_valid
            && self.is_sample_haplotype_list_valid
            && self.is_sv_haplotype_map_valid
            && self.is_multi_region_sv_count_valid
    }
}

#[derive(Clone)]
pub struct SVGroupHaplotype {
    pub hap_id: GroupHaplotypeId,

    pub contig_info: ClusterAssembly,
}

/// A set of RefinedSVs that will be scored as a group because of overlap or other short-range
/// interactions.
///
/// This includes single RefinedSVs as a degenerate case.
///
#[derive(Clone)]
pub struct SVGroup {
    /// Region(s) to search for supporting reads during scoring
    ///
    /// This is a multi-sample generalization of the assembly regions attached to each RefinedSV
    ///
    /// This must have length >= 1 in a valid group.
    ///
    pub group_regions: Vec<GenomeSegment>,

    /// Every assembly contig used to generate candidate SVs for the group should be included in
    /// the group_contigs list.
    ///
    /// The contig's order in this list defines its 'haplotype index'
    ///
    pub group_haplotypes: Vec<SVGroupHaplotype>,

    /// A vector with one entry for each sample covered by the SVGroup.
    ///
    /// Each sample entry lists zero to 2 (or ploidy limit) haplotype indexes associated with the
    /// sample
    ///
    pub sample_haplotype_list: Vec<Vec<usize>>,

    /// The maximum ploidy to use for each sample
    ///
    /// If there's any ambiguity in the ploidy count this will always be set to the higher value.
    ///
    pub sample_ploidy: Vec<SVLocusPloidy>,

    /// A vector with one entry for each Refined SV providing the index of the SV source haplotype
    pub sv_haplotype_map: Vec<usize>,

    /// x
    pub refined_svs: Vec<RefinedSV>,
}

impl SVGroup {
    pub fn is_single_region(&self) -> bool {
        self.group_regions.len() == 1
    }

    #[allow(clippy::field_reassign_with_default)]
    pub fn get_validity_test(&self) -> SVGroupTestStatus {
        let mut result = SVGroupTestStatus::default();

        // tmp
        result.are_group_regions_valid = !self.group_regions.is_empty();
        result.are_group_haplotypes_valid = {
            let expected = if self.is_single_region() { 1 } else { 2 };
            self.group_haplotypes
                .iter()
                .all(|x| x.contig_info.contig_alignments.len() == expected)
        };

        let haplotype_count = self.group_haplotypes.len();

        let sample_count = self.sample_haplotype_list.len();

        result.is_sample_max_haplotype_count_valid = self.sample_ploidy.len() == sample_count;

        result.is_sample_haplotype_list_valid = {
            let valid_haplotype_indexes = !self
                .sample_haplotype_list
                .iter()
                .flatten()
                .any(|&x| x >= haplotype_count);

            fn has_repeats<T>(iter: T) -> bool
            where
                T: IntoIterator,
                T::Item: Eq + Hash,
            {
                let mut uniq = HashSet::new();
                !iter.into_iter().all(move |x| uniq.insert(x))
            }

            let valid_sample_sets = !self
                .sample_haplotype_list
                .iter()
                .any(|x| has_repeats(x.iter()));

            let valid_hap_count = self
                .sample_haplotype_list
                .iter()
                .zip(self.sample_ploidy.iter())
                .all(|(list, ploidy)| {
                    let max = match ploidy {
                        SVLocusPloidy::Diploid => 2,
                        SVLocusPloidy::Haploid => 1,
                    };
                    list.len() <= max
                });

            valid_haplotype_indexes && valid_sample_sets && valid_hap_count
        };

        result.is_sv_haplotype_map_valid =
            !self.sv_haplotype_map.iter().any(|&x| x >= haplotype_count);

        result.is_multi_region_sv_count_valid =
            self.is_single_region() || self.refined_svs.len() == 1;

        result
    }

    pub fn assert_validity(&self) {
        let test_result = self.get_validity_test();
        if !test_result.is_valid() {
            error!("Invalid sv_group: {:?}\nreason: {:?}", self, test_result);
            for refined_sv in self.refined_svs.iter() {
                error!("RefinedSVID: {:?}", refined_sv.id);
            }
            panic!();
        }
    }

    /// SVGroup stores sv->haplotype map already, this method inverts it
    pub fn haplotype_to_sv_map(&self) -> Vec<Vec<usize>> {
        let haplotype_count = self.group_haplotypes.len();
        let mut x = vec![Vec::new(); haplotype_count];
        for (sv_index, &haplotype_index) in self.sv_haplotype_map.iter().enumerate() {
            x[haplotype_index].push(sv_index);
        }
        x
    }

    /// SVGroup stores sample->haplotype list already, this method inverts it
    #[allow(dead_code)]
    pub fn haplotype_to_sample_list(&self) -> Vec<Vec<usize>> {
        let haplotype_count = self.group_haplotypes.len();
        let mut x = vec![Vec::new(); haplotype_count];
        for (sample_index, haplotype_indexes) in self.sample_haplotype_list.iter().enumerate() {
            for &haplotype_index in haplotype_indexes.iter() {
                x[haplotype_index].push(sample_index);
            }
        }
        x
    }

    pub fn sample_to_sv_map(&self) -> Vec<Vec<usize>> {
        let haplotype_to_sv_map = self.haplotype_to_sv_map();

        let mut res = Vec::new();
        for haplotype_indexes in self.sample_haplotype_list.iter() {
            let mut sample_res = Vec::new();
            for &haplotype_index in haplotype_indexes {
                sample_res.extend(haplotype_to_sv_map[haplotype_index].iter());
            }
            res.push(sample_res);
        }
        res
    }
}

impl std::fmt::Debug for SVGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SVGroup:\n\tregions {:?}\n\thap count {}\n\tsample_hap_map: {:?}\n\tsv_hap_map: {:?}\n\tsv count {}",
               self.group_regions, self.group_haplotypes.len(), self.sample_haplotype_list, self.sv_haplotype_map, self.refined_svs.len())
    }
}
