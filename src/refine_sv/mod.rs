mod align;
mod assemble;
pub mod assembly_regions;
mod debug_output;
mod filter_dups;
mod trimmed_reads;

use std::fmt;
use std::sync::{
    Arc, Mutex,
    atomic::{self, AtomicBool},
    mpsc::channel,
};
use std::time::Instant;

use itertools::Itertools;
use log::info;
use rust_htslib::bam;
use rust_vc_utils::{ChromList, GenomeRef, ProgressReporter, downsample_vector};
use strum::EnumCount;

use self::assemble::{AssemblyResult, assemble_and_cluster_reads};
use self::assembly_regions::{
    get_assembly_regions_from_breakpoint_cluster, write_assembly_regions_to_bed,
};
use self::debug_output::*;
use self::filter_dups::filter_redundant_duplications;
use self::trimmed_reads::{
    get_sv_breakpoint_candidate_trimmed_reads, get_sv_candidate_region_trimmed_reads,
};

use crate::breakpoint::{Breakpoint, BreakpointCluster};
use crate::cli;
use crate::contig_output::{ContigAlignmentInfo, write_contig_alignments};
use crate::discover::StaticDiscoverSettings;
use crate::expected_ploidy::{SVLocusExpectedCNInfo, get_max_haplotype_count_for_regions};
use crate::genome_segment::{GenomeSegment, IntRange};
use crate::log_utils::debug_msg;
use crate::refine_sv::SVFilterType::{GTExcluded, NoFilter, Replicated};
use crate::refined_cnv::{SharedSampleScoreInfo, SharedVariantScoreInfo};
use crate::run_stats::RefineStats;
use crate::sv_group::SVGroup;
use crate::sv_id::{SVUniqueIdData, get_sv_id_label};
use crate::utils::print_fasta;
use crate::worker_thread_data::BamReaderWorkerThreadDataSet;

/// If the full consensus is usable, then return the range of the core assembly within it, and
/// None otherwise
#[allow(dead_code)]
fn is_full_consensus_usable(core: &[u8], prefix: &[u8], full: &[u8]) -> Option<IntRange> {
    if prefix.len() < core.len() {
        return None;
    }
    if full.len() < prefix.len() {
        return None;
    }
    let core_start = prefix.len() - core.len();
    let core_end = prefix.len();

    /*
    const MATCH_CHECK_SIZE: usize = 16;
    if core.len() < MATCH_CHECK_SIZE {
        return None;
    }
    if core[..MATCH_CHECK_SIZE] != prefix[core_start..core_start + MATCH_CHECK_SIZE] {
        return None;
    }

    if core[(core.len() - MATCH_CHECK_SIZE)..] != full[(core_end - MATCH_CHECK_SIZE)..core_end] {
        return None;
    }
     */

    Some(IntRange::from_pair(core_start as i64, core_end as i64))
}

#[derive(Clone, Copy, PartialEq, strum::Display, strum::EnumCount)]
pub enum AlleleType {
    Ref,
    Alt,
    Overlap,
}

impl fmt::Debug for AlleleType {
    // Set Debug trait to copy Display
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Display::fmt(self, f)
    }
}

#[derive(Clone, Debug, Default)]
pub struct AlleleCounts {
    pub any_strand: usize,
}

impl AlleleCounts {
    fn clear(&mut self) {
        self.any_strand = 0;
    }
}

#[derive(Clone, Debug, PartialEq, strum::EnumCount, strum::FromRepr)]
pub enum Genotype {
    Ref,
    Het,
    Hom,
}

/// Phasing information within the scope of one SVGroup
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SVPhaseStatus {
    /// No phase information available
    Unknown,

    /// SV is present (almost) only on the primary haplotype of the SVGroup
    Primary,

    /// SV is present (almost) only on the secondary haplotype of the SVGroup
    Secondary,
}

#[derive(Clone, Debug)]
pub struct SVSampleScoreInfo {
    /// This is true for any SV exceeding the maximum scoring depth, indicating that no scoring
    /// was attempted for this sample
    pub is_max_scoring_depth: bool,

    /// True if this SV is present on a haplotype associated with this sample
    ///
    /// This item is used as part of the scoring process but not part of SV reporting output
    ///
    pub on_sample_haplotype: bool,

    /// True if this SV is known to be represented on both haplotype assemblies for this sample,
    ///
    /// Note the haplotypes may still be different otherwise (either biologically or due to assembly
    /// noise), so we still keep both assemblies and evaluate both in the scoring process. Note also
    /// that the variant may still be a het for the case where one haplotype is split due to
    /// assembly noise.
    ///
    /// This item is used as part of the scoring process but not part of SV reporting output
    ///
    pub on_multiple_sample_haplotypes: bool,

    /// True if on_multiple_sample_haplotypes is true and this SV is not on the primary sample
    /// haplotype.
    ///
    /// If true this SV will be filtered from the scoring procedure for the given sample. It may
    /// be scored in a different sample.
    ///
    /// This item is used as part of the scoring process but not part of SV reporting output
    ///
    pub is_filtered_sv_replicate_within_sample: bool,

    /// True if overlapping genotypes in this sample suggest that this SV should be filtered out
    pub is_gt_excluded: bool,

    /// Allele counts for all alleles considered in the scoring process
    ///
    /// The allele is encoded by the vector index. Index 0 is the reference, index 1 is the SV
    /// allele in question, and other overlapping alleles are greater than 1.
    ///
    pub allele_depth: Vec<AlleleCounts>,

    /// Expected copy number info/ploidy of the SV locus
    pub expected_cn_info: SVLocusExpectedCNInfo,

    pub phase: SVPhaseStatus,

    /// Optional list of read names supporting the alt allele, used for debugging
    pub supporting_read_names: Vec<String>,

    /// All shared per-sample scoring information for any variant type
    pub shared: SharedSampleScoreInfo,
}

impl Default for SVSampleScoreInfo {
    fn default() -> Self {
        let mut x = Self {
            is_max_scoring_depth: false,
            on_sample_haplotype: false,
            on_multiple_sample_haplotypes: false,
            is_filtered_sv_replicate_within_sample: false,
            is_gt_excluded: false,
            allele_depth: Vec::new(),
            expected_cn_info: SVLocusExpectedCNInfo {
                expected_copy_number: 2,
            },
            phase: SVPhaseStatus::Unknown,
            supporting_read_names: Vec::new(),
            shared: SharedSampleScoreInfo::default(),
        };
        x.allele_depth
            .resize_with(AlleleType::COUNT, Default::default);
        x
    }
}

impl SVSampleScoreInfo {
    pub fn total_allele_depth(&self) -> usize {
        self.allele_depth.iter().map(|x| x.any_strand).sum()
    }

    /// Clears all allele counts and items associated with allele count, like supporting read names
    ///
    /// This should be called for certain cases where the allele count loop is run twice to resolve read to haplotype
    /// assignment in the first pass
    ///
    pub fn clear_allele_depth(&mut self) {
        for x in self.allele_depth.iter_mut() {
            x.clear();
        }
        self.supporting_read_names.clear();
    }
}

#[derive(Clone, Default)]
pub struct SVScoreInfo {
    /// Quality score for any non-reference genotype in any sample based on read-support likelihoods
    pub sv_alt_score: Option<f32>,

    /// Quality score for any non-default copy number in any sample based on depth-support likelihoods
    ///
    /// This is only defined if the SV has been merged with a CNV
    pub cnv_alt_score: Option<f32>,

    /// If true this SV should be filtered as a replicate.
    ///
    /// At the global level, this flag indicates that all SV replicates across all sample haplotypes
    /// have been identified, and one representative SV has been selected for output, excluding this
    /// one.
    ///
    pub is_filtered_sv_replicate: bool,

    /// Defined if any samples in the sv_group are phased
    pub phase_set: Option<i32>,

    /// Scores for each sample
    pub samples: Vec<SVSampleScoreInfo>,
}

impl SVScoreInfo {
    /// Quality score for any non-reference genotype in any sample or any non-default copy number in any sample
    pub fn alt_score(&self) -> Option<f32> {
        match (self.sv_alt_score, self.cnv_alt_score) {
            (Some(sv), Some(cnv)) => Some(sv.max(cnv)),
            (Some(sv), None) => Some(sv),
            (None, Some(cnv)) => Some(cnv),
            (None, None) => None,
        }
    }

    pub fn update_cnv_alt_score(&mut self, max_qscore: u32) {
        self.cnv_alt_score = self.get_cnv_alt_score(max_qscore);
    }
}

impl SharedVariantScoreInfo for SVScoreInfo {
    fn get_shared_sample_iter(&self) -> impl Iterator<Item = &SharedSampleScoreInfo> {
        self.samples.iter().map(|x| &x.shared)
    }
}

#[derive(Clone)]
pub struct AnnotatedOverlappingHaplotype {
    pub assembly_index: usize,

    /// Number of assembly reads supporting this haplotype
    pub supporting_read_count: usize,

    /// Number of SV candidates generated from this haplotype
    pub sv_candidate_count: usize,
}

#[derive(PartialEq)]
pub enum SVFilterType {
    NoFilter,
    Replicated,
    GTExcluded,
}

/// Extended items for RefinedSV, they're set into a separate struct so they can be defaulted
#[derive(Clone, Default)]
pub struct RefinedSVExt {
    /// If true, do not represent this SV as a higher order event, and assess/report it as a
    /// raw breakpoint
    ///
    /// This will probably end up being shared across all samples
    ///
    /// This item is only used for output formatting after refinement, it is conceptually a form
    /// of scoring as well, given that it's used to reclassify SV types in the final output.
    ///
    pub force_breakpoint_representation: bool,

    /// If defined, this SV is part of a group of breakends represented as an inversion, either the
    /// original breakend records (which will be marked as filtered in the output) or the combined
    /// inversion record.
    ///
    pub inversion_id: Option<usize>,

    /// If true this SV represents a VCF inversion record, this is a special type of refined SV which
    /// is only valid for the purpose of VCF output.
    ///
    pub is_inversion: bool,

    /// True for inversions where the GT of the two breakends dissagree in the majority of cases
    pub is_conflicting_gt: bool,

    /// True for SVs with CNV support
    pub is_cnv_match: bool,
}

/// Structural variants including all detailed breakpoint analysis and scoring information
///
/// This is the final format for SVs in preperation for VCF output
#[derive(Clone)]
pub struct RefinedSV {
    pub id: SVUniqueIdData,
    pub bp: Breakpoint,

    /// Homology sequence, from the perspective of breakend1 in fully left-shifted position
    pub breakend1_homology_seq: Vec<u8>,

    /// Mark true to indicate single as opposed to dual region refinement
    pub single_region_refinement: bool,

    /// Position in the contig sequence immediately before breakend1 of the RefinedSV
    pub contig_pos_before_breakend1: usize,

    /// Assembly regions for the RefinedSV
    ///
    /// This should be shared by all RefinedSV entries corresponding to the same cluster
    ///
    pub assembly_regions: Vec<GenomeSegment>,

    pub ext: RefinedSVExt,

    /// All supporting read and quality score information
    ///
    /// For the case of an inverted breakpoint, forward and reverse strand counts in SVScoreInfo
    /// refer to breakpoint1
    ///
    /// This item is only used for scoring after refinement
    pub score: SVScoreInfo,
}

impl RefinedSV {
    fn filter_gt_excluded_sv(&self) -> bool {
        let mut filter = false;
        for sample in self.score.samples.iter() {
            if sample.on_sample_haplotype {
                if sample.is_gt_excluded {
                    filter = true;
                } else {
                    // A single sample where the SV is present but not marked as a gt_excluded will
                    // deactivate the filter for all samples:
                    return false;
                }
            }
        }
        filter
    }

    pub fn get_sv_filter_type(&self) -> SVFilterType {
        if self.score.is_filtered_sv_replicate {
            Replicated
        } else if self.filter_gt_excluded_sv() {
            GTExcluded
        } else {
            NoFilter
        }
    }

    pub fn filter_sv(&self) -> bool {
        self.get_sv_filter_type() != NoFilter
    }
}

pub fn get_rsv_id_label(rsv: &RefinedSV) -> String {
    if rsv.ext.is_inversion {
        let pname = env!("CARGO_PKG_NAME");
        assert!(rsv.ext.inversion_id.is_some());
        let id = rsv.ext.inversion_id.unwrap();
        format!("{pname}:INV{id}")
    } else {
        get_sv_id_label(&rsv.id)
    }
}

#[allow(clippy::too_many_arguments)]
fn refine_sv_candidate_breakpoint(
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    bpc: &BreakpointCluster,
    cluster_index: usize,
    worker_thread_index: usize,
    debug: bool,
) -> (
    Vec<ContigAlignmentInfo>,
    Option<SVGroup>,
    ClusterDurationInfo,
) {
    // Get the reference regions from which to take reads for breakpoint assembly
    //
    let assembly_segments = get_assembly_regions_from_breakpoint_cluster(bpc);

    debug_msg!(
        debug,
        "Cluster {cluster_index}: Starting refine_sv_candidate_breakpoint on worker thread {worker_thread_index}. Regions: {}",
        assembly_segments
            .iter()
            .map(|x| x.to_region_str(chrom_list))
            .join(", ")
    );

    // Don't use expected copy number regions at this point (here we change expected copy number regions to None)
    //
    // Punting on this and preserving up to 2 haplotypes per locus allows sex to be changed downstream during joint genotyping.
    //
    let treat_single_copy_as_haploid = false;
    let (max_assembly_count, expected_cn_info) =
        get_max_haplotype_count_for_regions(treat_single_copy_as_haploid, None, &assembly_segments);

    let trimmed_reads = if assembly_segments.len() == 1 {
        get_sv_candidate_region_trimmed_reads(
            refine_settings,
            reference,
            chrom_list,
            bam_reader,
            &assembly_segments[0],
            bpc,
        )
    } else {
        // If we're handling the nb allele as 2 assemblies, then in theory we will do 2 region scans.
        //
        // For now the only things we'll detect are:
        // 1. Very large deletions encoded in the CIGAR alignment string
        // 2. Other types of breakends encoded in the SA split read pattern.
        //
        // For now, just scan the first region to get the reads with an alignment pattern consistent
        // with the candidate alt allele, and return trimmed reads only for those cases around the
        // breakend in question. The trimmed reads should be fwd or rev oriented consistent with the
        // scanned region (region 1)
        //
        // -----
        //
        // Future workflow for this case:
        // 1. Scan 2 regions
        // 2. Scan all reads in regions
        // 3. Cluster these and anticipate they should form into 3 signal groups:
        //      Alt, non-alt1, non-alt2
        // 4. Also get soft-clipped read evidence
        //
        get_sv_breakpoint_candidate_trimmed_reads(
            refine_settings,
            reference,
            chrom_list,
            bam_reader,
            &assembly_segments,
            &bpc.breakpoint,
        )
    };

    debug_msg!(
        debug,
        "Cluster {cluster_index}: refine_sv_candidate_breakpoint: Trimmed read count {}",
        trimmed_reads.len()
    );
    if debug {
        eprintln!("trimmed read info:");
        for tr in trimmed_reads.iter() {
            eprintln!(
                "{}: rr: {:?} lics: {:?}",
                tr.qname, tr.read_range, tr.large_insertion_clip_state
            );
        }
    }

    let is_large_insertion_inferred_from_clipping = bpc.large_insertion_info.is_some();
    let downsample_target = if is_large_insertion_inferred_from_clipping {
        refine_settings.max_trimmed_reads_at_large_insertion
    } else {
        refine_settings.max_trimmed_reads
    };

    // When more than downsample_target reads are found, deterministically subsample:
    let trimmed_reads = downsample_vector(trimmed_reads, downsample_target);

    let mut cluster_duration_info = ClusterDurationInfo::default();

    let assemblies = {
        let start_time = Instant::now();
        let is_multi_region = assembly_segments.len() > 1;
        let rval = assemble_and_cluster_reads(
            refine_settings,
            cluster_index,
            max_assembly_count,
            is_large_insertion_inferred_from_clipping,
            trimmed_reads,
            is_multi_region,
        );
        cluster_duration_info.assembly = start_time.elapsed();
        rval
    };

    debug_msg!(
        debug,
        "Cluster {cluster_index}: refine_sv_candidate_breakpoint: Assemblies count {}",
        assemblies.len()
    );
    if debug {
        let print_contigs = false;
        if print_contigs {
            let seqs = assemblies
                .iter()
                .flat_map(|x| &x.contigs)
                .map(|x| x.seq.as_slice())
                .collect::<Vec<_>>();
            print_fasta(250, &seqs);
        }
    }

    let is_large_insertion = bpc.large_insertion_info.is_some();

    let (contig_alignments, sv_group) = {
        let start_time = Instant::now();
        let rval = align::align_assemblies(
            refine_settings,
            reference,
            chrom_list,
            &assembly_segments,
            &assemblies,
            expected_cn_info,
            bpc,
            cluster_index,
            is_large_insertion,
        );
        cluster_duration_info.alignment = start_time.elapsed();
        rval
    };

    debug_msg!(
        debug,
        "Cluster {cluster_index}: refine_sv_candidate_breakpoint: Alignment count {}",
        contig_alignments.len()
    );

    (contig_alignments, sv_group, cluster_duration_info)
}

pub struct RefineSVSettings {
    /// Copied from cli settings
    min_indel_size: u32,

    /// Copied from cli settings
    min_evidence_indel_size: u32,

    /// Copied from cli settings
    min_sv_mapq: u32,

    /// Copied from cli settings
    min_gap_compressed_identity: f64,

    /// Distance in read coordinate from the candidate SV breakends to extend read trimming for read
    /// input to assembly
    assembly_read_flank_size: usize,

    /// Minimum soft-clip length on the read edge for use in large insertion assembly
    min_large_insertion_soft_clip_len: usize,

    /// Distance in ref coordinates to extend the reference template for realigning assemblies to
    /// the SV locus
    hap_alignment_ref_flank_size: i64,

    /// Maximum amount of additional non-assembled read prefix or suffix which is appended to the
    /// end of the contig
    max_contig_fake_flank_size: i64,

    /// For breakpoints (ie. 2 assembly region loci), retry fake flank and alignment at a larger fake flank size
    /// if first alignment attempt fails.
    breakpoint_retry_max_contig_fake_flank_size: i64,

    /// When reads are selected for use as fake flanks, require that they align at least this close
    /// to the corresponding edge of the assembly contig.
    max_ref_clip_for_fake_flank_selection: usize,

    /// Exact match sequence length required on the left and right flanks of a candidate SV
    /// assembly contig when it is aligned back to a a model of the candidate alt haplotype.
    ///
    min_assembly_edge_anchor: u32,

    /// Min number of supporting reads in an SV haplotype assembly
    min_assembly_read_support: usize,

    /// During region assembly, only track up to this many candidate assemblies:
    max_candidate_assemblies: usize,

    /// During region assembly, only track up to this many candidate assemblies if this is a large
    /// insertion locus:
    max_candidate_assemblies_at_large_insertion: usize,

    /// The number of trimmed reads at a single assembly region are subsampled down to this value
    max_trimmed_reads: usize,

    /// The number of trimmed reads at a single large insertion region are subsampled down to this value
    max_trimmed_reads_at_large_insertion: usize,

    two_region_alt_haplotype_spacer_size: usize,

    /// In a two region alignment, this is the required alignment quality for the left or right side
    /// anchor sequences
    anchor_min_gap_compressed_identity: f64,

    /// Filter overlapping or very close SVs
    ///
    /// This option should not be on in a typical production run, but is useful for GIAB v0.6
    /// compatibility.
    ///
    reduce_overlapping_sv_alleles: bool,

    /// If close overlapping SV filtration is on, this is the min distance threshold from previously accepted SVs
    /// allowed for an SV to be accepted.
    ///
    min_close_sv_range: i64,

    /// If true, collapse insertions and SNPs near the breakpoint into the SV
    collapse_indels_at_breakpoint: bool,

    /// Minimum normalized alignment score for the alignment of a trimmed read to the POA assembly
    /// graph, for the read to be accepted and clustered into that graph
    min_norm_score_for_read_cluster: f64,

    /// Minimum normalized alignment score for the alignment of a trimmed read to the POA assembly
    /// graph, for the read to be accepted and clustered into that graph
    ///
    /// This score is only used if the first cluster process fails
    min_norm_score_for_read_cluster_2nd_pass: f64,
}

impl RefineSVSettings {
    pub fn new(
        min_indel_size: u32,
        min_evidence_indel_size: u32,
        min_sv_mapq: u32,
        min_gap_compressed_identity: f64,
        reduce_overlapping_sv_alleles: bool,
        assembly_read_flank_size: usize,
    ) -> Self {
        Self {
            min_indel_size,
            min_evidence_indel_size,
            min_sv_mapq,
            min_gap_compressed_identity,
            assembly_read_flank_size,
            min_large_insertion_soft_clip_len: 500,
            hap_alignment_ref_flank_size: 3_500,
            max_contig_fake_flank_size: 1_500,
            breakpoint_retry_max_contig_fake_flank_size: 2_500,
            max_ref_clip_for_fake_flank_selection: 200,
            min_assembly_edge_anchor: 50,
            min_assembly_read_support: 2,
            max_candidate_assemblies: 8,
            max_candidate_assemblies_at_large_insertion: 2,
            max_trimmed_reads: 100,
            max_trimmed_reads_at_large_insertion: 50,
            two_region_alt_haplotype_spacer_size: 100,
            anchor_min_gap_compressed_identity: 0.9,
            reduce_overlapping_sv_alleles,
            min_close_sv_range: 100,
            collapse_indels_at_breakpoint: false,
            min_norm_score_for_read_cluster: 0.96,
            min_norm_score_for_read_cluster_2nd_pass: 0.91,
        }
    }
}

type RefineWorkerReturnType = (ClusterDebugInfo, Vec<ContigAlignmentInfo>, Option<SVGroup>);

/// Complete process to run refinement on a single breakpoint cluster, including various debug
/// overheads and channeling the result back to parent thread
///
#[allow(clippy::too_many_arguments)]
fn refine_sv_candidate_breakpoint_wrapper(
    debug_worker_status_report: bool,
    worker_debug_info: &AllWorkerDebugInfo,
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    tx: std::sync::mpsc::Sender<RefineWorkerReturnType>,
    refine_settings: &RefineSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    cluster_index: usize,
    bpc: &BreakpointCluster,
    progress_reporter: Option<&ProgressReporter>,
    debug: bool,
) {
    // Get total processing time for this cluster:
    let cluster_start_time = Instant::now();

    let worker_thread_index = rayon::current_thread_index().unwrap();
    if debug_worker_status_report {
        // Update worker thread status to reflect the start of this job
        let worker_id_debug_info = &mut *worker_debug_info[worker_thread_index].lock().unwrap();
        *worker_id_debug_info = Some(WorkerDebugInfo {
            cluster_start_time,
            bp: bpc.breakpoint.clone(),
            car: bpc.consolidated_assembly_segment.clone(),
        });
    }

    // Hard-code one sample in discover mode:
    let sample_index = 0;
    let bam_reader = &mut worker_thread_dataset[worker_thread_index]
        .lock()
        .unwrap()
        .bam_readers[sample_index];

    let (contig_alignments, sv_group, mut duration) = refine_sv_candidate_breakpoint(
        refine_settings,
        reference,
        chrom_list,
        bam_reader,
        bpc,
        cluster_index,
        worker_thread_index,
        debug,
    );

    duration.total = cluster_start_time.elapsed();
    let cdi = ClusterDebugInfo {
        duration,
        bp: bpc.breakpoint.clone(),
    };
    let result = (cdi, contig_alignments, sv_group);
    tx.send(result).unwrap();
    if debug_worker_status_report {
        // Update worker thread status to reflect the end of this job:
        let worker_id_debug_info = &mut *worker_debug_info[worker_thread_index].lock().unwrap();
        *worker_id_debug_info = None;
    } else if let Some(progress_reporter) = progress_reporter {
        progress_reporter.inc(1);
    }
}

fn get_cluster_job_order(clusters: &[BreakpointCluster]) -> Vec<(usize, &BreakpointCluster)> {
    let mut job_order = Vec::new();
    for (cluster_index, bpc) in clusters.iter().enumerate() {
        // Skip everything except soft-clip clusters that root a large insertion pairing:
        if bpc.breakpoint.breakend2.is_some() || bpc.large_insertion_info.is_none() {
            continue;
        }

        job_order.push((cluster_index, bpc));
    }

    // Start standard candidate refinement jobs
    for (cluster_index, bpc) in clusters.iter().enumerate() {
        // Skip:
        // (1) indel clusters which have been joined into a consolidated cluster
        // (2) all soft-clip clusters
        if bpc.consolidated_cluster_index.is_some() || bpc.breakpoint.breakend2.is_none() {
            continue;
        }
        job_order.push((cluster_index, bpc));
    }

    job_order
}

/// Refine all SV candidates
///
/// SV candidate refinement is distributed over a threadpool
///
#[allow(clippy::too_many_arguments)]
pub fn refine_sv_candidates(
    shared_settings: &cli::SharedSettings,
    settings: &cli::DiscoverSettings,
    static_settings: &StaticDiscoverSettings,
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    clusters: &[BreakpointCluster],
) -> (Vec<SVGroup>, RefineStats) {
    info!("Refining SV Candidates");

    // When true a worker thread status report is reported at regular time intervals during
    // refinement.
    //
    // To reduce output noise the refinement progress bar is disabled when this status report is
    // enabled.
    let debug_worker_status_report = false;

    let refine_settings = &RefineSVSettings::new(
        settings.min_indel_size,
        settings.get_min_evidence_indel_size(),
        settings.min_sv_mapq,
        settings.min_gap_compressed_identity,
        settings.reduce_overlapping_sv_alleles,
        static_settings.assembly_read_flank_size,
    );

    let mut refine_stats = RefineStats::default();

    // Set the cluster processing order so that large insertions are handled first, to improve
    // thread utilization (b/c large insertion jobs are slow)
    //
    let cluster_job_order = get_cluster_job_order(clusters);

    let disable_small_indels = settings.fast_cnv_mode || shared_settings.disable_small_indels;

    // Optionally apply debug filter
    let cluster_job_order = cluster_job_order
        .into_iter()
        .filter(|(_, bpc)| (!disable_small_indels || (bpc.assembly_segments.len() > 1)))
        .collect::<Vec<_>>();

    write_assembly_regions_to_bed(&settings.output_dir, chrom_list, clusters);

    let refine_threads = shared_settings.thread_count;

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(refine_threads)
        .build()
        .unwrap();

    let worker_debug_info: AllWorkerDebugInfo =
        Arc::new((0..refine_threads).map(|_| Mutex::new(None)).collect());

    let run_worker_status_reporter = Arc::new(AtomicBool::new(true));
    if debug_worker_status_report {
        let run_worker_status_reporter = run_worker_status_reporter.clone();
        let worker_debug_info = worker_debug_info.clone();
        std::thread::spawn(|| {
            worker_status_reporter(run_worker_status_reporter, worker_debug_info);
        });
    }
    let worker_debug_info = &worker_debug_info;

    info!("Starting refinement of breakpoint evidence clusters ");
    let progress_reporter = if !debug_worker_status_report {
        // Turn off the progress bar if we're logging debug info during refinement:
        let force_periodic_updates = shared_settings.debug;
        Some(ProgressReporter::new(
            cluster_job_order.len() as u64,
            "Refined",
            "breakpoint evidence clusters",
            force_periodic_updates,
        ))
    } else {
        None
    };
    let progress_reporter = progress_reporter.as_ref();

    let (tx, rx) = channel();
    worker_pool.scope(move |scope| {
        for (cluster_index, bpc) in cluster_job_order {
            let mut debug = false;
            if let Some(target_cluster_index) = settings.target_cluster_index {
                if cluster_index != target_cluster_index {
                    continue;
                }
                debug = true;
            }
            let tx = tx.clone();
            scope.spawn(move |_| {
                refine_sv_candidate_breakpoint_wrapper(
                    debug_worker_status_report,
                    worker_debug_info,
                    worker_thread_dataset,
                    tx,
                    refine_settings,
                    reference,
                    chrom_list,
                    cluster_index,
                    bpc,
                    progress_reporter,
                    debug,
                );
            });
        }
    });

    let mut debug_info = Vec::new();
    let mut contig_alignments = Vec::new();
    let mut sv_groups = Vec::new();
    for (chunk_debug_info, chunk_contig_alignments, sv_group) in rx {
        debug_info.push(chunk_debug_info);
        contig_alignments.extend(chunk_contig_alignments);
        sv_groups.push(sv_group);
    }

    // To clear the progress bar, this needs to be called before any other log entries are written:
    if let Some(progress_reporter) = progress_reporter {
        progress_reporter.clear();
    }
    info!("Finished refining all breakpoint evidence clusters");

    // Filter sv group options down to just existing sv groups:
    let mut sv_groups = sv_groups.into_iter().flatten().collect::<Vec<_>>();

    // Add an extra step to identify redundant SVs to filter out
    refine_stats.redundant_duplications_filtered_out =
        filter_redundant_duplications(chrom_list, &mut sv_groups);

    if debug_worker_status_report {
        // Turn off the worker status reporter if it's running:
        run_worker_status_reporter.store(false, atomic::Ordering::Relaxed);
    }

    // Write out a debug bam file of contig alignments to the reference
    write_contig_alignments(
        &settings.output_dir,
        shared_settings.thread_count,
        chrom_list,
        "all",
        contig_alignments,
    );

    // Add total refinement times to refine_stats:
    for cdi in debug_info.iter() {
        refine_stats.total_refinement_time_secs += cdi.duration.total.as_secs_f64();
        refine_stats.total_refinement_assembly_time_secs += cdi.duration.assembly.as_secs_f64();
        refine_stats.total_refinement_alignment_time_secs += cdi.duration.alignment.as_secs_f64();
        refine_stats.total_refinement_scoring_time_secs += cdi.duration.scoring.as_secs_f64();
    }

    // Write out debug information for the N most time-consuming clusters:
    write_debug_cluster_info(&settings.output_dir, debug_info);

    (sv_groups, refine_stats)
}
