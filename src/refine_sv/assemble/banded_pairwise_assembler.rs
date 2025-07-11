use std::collections::BTreeMap;

use rust_htslib::bam::record::Cigar;
use rust_vc_utils::cigar::{
    get_cigar_ref_offset, get_complete_read_clip_positions, update_ref_and_hard_clipped_read_pos,
};
use unwrap::unwrap;

use super::add_flanks::{
    ConsensusSequence, ConsensusSequenceInfo, add_flanks, get_consensus_alignment_info,
};
use super::{AssemblyResult, RefineSVSettings};
use crate::bio_align_utils::{AlignmentWeights, OverlapPairwiseAligner, PairwiseAligner};
use crate::refine_sv::trimmed_reads::TrimmedReadInfo;
use crate::simple_alignment::SimpleAlignment;
use crate::utils::{drop_true, print_fasta};

/// Return left soft clip size, right soft clip size, and read size
///
fn get_alignment_soft_clips(cigar: &[Cigar]) -> (u32, u32, u32) {
    let (left, right, size) = get_complete_read_clip_positions(cigar);
    assert!(size >= right);
    (left as u32, (size - right) as u32, size as u32)
}

/// Metrics used to determine the quality of a sequence to backbone alingment
#[derive(Debug, Default)]
struct BackboneAlignmentQualityMetrics {
    /// Alignment score
    score: i32,

    /// Number of read bases used in the alignment (subtracting any clipping on the left and right sides)
    overlap: u32,
}

#[derive(Debug, Default)]
struct BackboneAlignmentInfo {
    /// Index of the backbone being aligned to
    index: usize,

    metrics: BackboneAlignmentQualityMetrics,

    /// Number of read bases clipped from the left side of the alignment
    _left_clip: u32,

    alignment: Option<SimpleAlignment>,
}

// Scheme:
// 1. Build insertion backbone(s):
// Transition through reads in the same order we input to POA
// Define the backbone as equal to the first read
// Align the next read to the backbone. If the read extends the backbone, then redefine the backbone with this extension.
// Repeat for all reads to form a fully extended backbone.
// Rejection scoring follows the same scheme as for POA, so we might be building multiple backbones.
//
// 2.Pairwise POA from each backbone
// For each backbone and its associated read set, pairwise align all reads in the set to the backbone, and take a majority vote on each edit to determine the consensus sequence.

/// Standard alignment process applied to trimmed reads against all backbones
///
/// Selects the best scoring backbone. More in-depth analysis of the alignment
/// quality for backbone join/split decisions are left to downstream steps
///
fn get_best_backbone_alignments_for_trimmed_read(
    alignment_weights: &AlignmentWeights,
    candidate_backbones: &[ConsensusSequence],
    read: &[u8],
) -> (BackboneAlignmentInfo, BackboneAlignmentInfo) {
    let debug = false;

    let mut aligner = OverlapPairwiseAligner::new_aligner(alignment_weights);

    let mut rank1_backbone = BackboneAlignmentInfo::default();
    let mut rank2_backbone = BackboneAlignmentInfo::default();

    for (candidate_backbone_index, candidate_backbone) in candidate_backbones.iter().enumerate() {
        let score = aligner.align(candidate_backbone.seq.as_slice(), read);
        let alignment = aligner.get_last_alignment().unwrap();
        let (left_clip, right_clip, read_size) = get_alignment_soft_clips(&alignment.cigar);
        let overlap = read_size - (left_clip + right_clip);

        if debug {
            eprintln!("backbone index {candidate_backbone_index} score {score} overlap {overlap}");
        }
        if overlap > 0
            && (rank1_backbone.alignment.is_none() || score > rank1_backbone.metrics.score)
        {
            let metrics = BackboneAlignmentQualityMetrics { score, overlap };

            rank2_backbone = rank1_backbone;

            rank1_backbone = BackboneAlignmentInfo {
                index: candidate_backbone_index,
                metrics,
                _left_clip: left_clip,
                alignment: Some(alignment),
            };
        }
    }
    (rank1_backbone, rank2_backbone)
}

fn get_norm_score(align_metrics: &BackboneAlignmentQualityMetrics) -> f64 {
    (align_metrics.score as f64) / (align_metrics.overlap as f64)
}

/// Determine if the trimmed read alignment to the given backbone is high enough quality to join the
/// the read into that backbone or if a new backbone should be started from the read instead.
///
/// Selection criteria:
/// (1) minimum overlap size
/// (2) minimum score/overlap size
///
/// Return true if the backbone alignment should be rejected
///
fn reject_backbone_alignment(
    min_norm_score: f64,
    align_metrics: &BackboneAlignmentQualityMetrics,
) -> bool {
    if align_metrics.overlap == 0 {
        true
    } else {
        const MIN_OVERLAP: u32 = 100;
        align_metrics.overlap < MIN_OVERLAP || get_norm_score(align_metrics) < min_norm_score
    }
}

fn abs_diff<U: std::cmp::PartialOrd + std::ops::Sub<Output = U>>(x1: U, x2: U) -> U {
    if x1 < x2 { x2 - x1 } else { x1 - x2 }
}

/// Determine if 2 backbone alignments are at a very similar quality level
///
/// This is a simple heuristic filter so the similarity definition is just a score diff cutoff.
///
fn are_metrics_close(
    match_score: i32,
    align_metrics1: &BackboneAlignmentQualityMetrics,
    align_metrics2: &BackboneAlignmentQualityMetrics,
) -> bool {
    if align_metrics1.overlap == 0 || align_metrics2.overlap == 0 {
        false
    } else {
        let max_similar_count_diff = 3;

        // This means max score difference as measured in units of the match scoring value:
        let max_similar_match_diff = 3;
        abs_diff(align_metrics1.overlap, align_metrics2.overlap) <= max_similar_count_diff
            && abs_diff(align_metrics1.score, align_metrics2.score)
                <= (max_similar_match_diff * match_score)
    }
}

/// Process the next trimmed read into the candidate backbone set
///
fn add_trimmed_read_to_candidate_backbones(
    max_candidates: usize,
    min_norm_score: f64,
    alignment_weights: &AlignmentWeights,
    trimmed_read_index: usize,
    trimmed_read_info: &TrimmedReadInfo,
    candidate_backbones: &mut Vec<ConsensusSequence>,
    assembly_debug: bool,
) {
    let extra_debug = false;

    if assembly_debug {
        eprintln!(
            "trimmed read index/size/qname: {trimmed_read_index}/{}/{}",
            &trimmed_read_info.trimmed_read.len(),
            &trimmed_read_info.qname,
        );
    }

    let (rank1_alignment, rank2_alignment) = get_best_backbone_alignments_for_trimmed_read(
        alignment_weights,
        candidate_backbones,
        trimmed_read_info.trimmed_read.as_slice(),
    );

    if extra_debug {
        eprintln!("rank1_alignment: {rank1_alignment:?}");
        eprintln!("rank2_alignment: {rank2_alignment:?}");
    }

    let reject_rank1_alignment =
        reject_backbone_alignment(min_norm_score, &rank1_alignment.metrics);

    if !reject_rank1_alignment {
        // Test if rank1 and rank2 alignments are very close:
        if are_metrics_close(
            alignment_weights.match_,
            &rank1_alignment.metrics,
            &rank2_alignment.metrics,
        ) {
            if assembly_debug {
                eprintln!("Skipping read due to similar alignment metrics to multiple backbones");
            }
        } else {
            if assembly_debug {
                eprintln!(
                    "Adding read to candidate backbone {}",
                    rank1_alignment.index
                );
            }
            add_read_to_backbone(
                trimmed_read_index,
                &rank1_alignment.alignment.unwrap(),
                &trimmed_read_info.trimmed_read,
                &mut candidate_backbones[rank1_alignment.index],
            );
        }
    } else if candidate_backbones.len() < max_candidates {
        // make new cluster
        if assembly_debug {
            eprintln!(
                "Starting new candidate assembly {}",
                candidate_backbones.len()
            );
        }
        candidate_backbones.push(new_backbone_from_read(
            trimmed_read_index,
            &trimmed_read_info.trimmed_read,
        ));
    } else {
        // filter read out
        if assembly_debug {
            eprintln!(
                "Skipping new candidate assembly {}",
                candidate_backbones.len()
            );
        }
    }
}

/// Create a new backbone from a single starting read
fn new_backbone_from_read(read_index: usize, read: &[u8]) -> ConsensusSequence {
    ConsensusSequence {
        source_reads: vec![read_index],
        seq: read.to_vec(),
    }
}

fn add_read_to_backbone(
    read_index: usize,
    alignment: &SimpleAlignment,
    read: &[u8],
    consensus_seq: &mut ConsensusSequence,
) {
    consensus_seq.source_reads.push(read_index);

    // Extend left or right side of backbone if read alignment contains large overlap on one side
    //
    // TODO: Beyond the feasibility stage we'll want a better way to reject unexpected read end patterns
    //
    let (left_clip, right_clip, _) = get_alignment_soft_clips(&alignment.cigar);

    let min_clip_size = 2;

    if left_clip < min_clip_size && right_clip < min_clip_size {
        return;
    }

    let max_ref_clip = 10;
    if left_clip > right_clip {
        // Check that left_ref_clip is near zero as expected
        if alignment.ref_offset > max_ref_clip {
            return;
        }
        consensus_seq
            .seq
            .splice(0..0, read[..left_clip as usize].iter().cloned());
    } else {
        // Check that right_ref_clip is near zero as expected
        let right_ref_offset = consensus_seq.seq.len() as i64
            - (alignment.ref_offset + get_cigar_ref_offset(&alignment.cigar));
        if right_ref_offset > max_ref_clip {
            return;
        }
        consensus_seq
            .seq
            .extend_from_slice(&read[(read.len() - right_clip as usize)..]);
    }
}

fn build_insertion_backbones(
    max_candidates: usize,
    min_norm_score: f64,
    alignment_weights: &AlignmentWeights,
    trimmed_reads: &[TrimmedReadInfo],
    assembly_debug: bool,
) -> Vec<ConsensusSequence> {
    let mut candidate_backbones = Vec::new();

    for (trimmed_read_index, trimmed_read_info) in trimmed_reads.iter().enumerate() {
        add_trimmed_read_to_candidate_backbones(
            max_candidates,
            min_norm_score,
            alignment_weights,
            trimmed_read_index,
            trimmed_read_info,
            &mut candidate_backbones,
            assembly_debug,
        );
    }

    let extra_debug = false;
    if extra_debug {
        for (candidate_backbone_index, candidate_backbone) in candidate_backbones.iter().enumerate()
        {
            let count = candidate_backbone.source_reads.len();
            let length = candidate_backbone.seq.len();
            eprintln!(
                "Candidate backbone {candidate_backbone_index} length: {length} read count: {count} reads: {:?}",
                candidate_backbone.source_reads
            );
        }
    }

    // Sort top expected backbones by sequence count and report these as final results:
    candidate_backbones.sort_by_key(|k| std::cmp::Reverse(k.source_reads.len()));

    candidate_backbones
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
struct IndelEdit {
    del_len: usize,
    ins_seq: Vec<u8>,
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
struct IndelEditTracker {
    ref_pos: i64,
    indel_edit: IndelEdit,
}

impl IndelEditTracker {
    fn new(ref_pos: i64) -> Self {
        Self {
            ref_pos,
            indel_edit: IndelEdit {
                del_len: 0,
                ins_seq: Vec::new(),
            },
        }
    }
}

#[derive(Clone, Eq, Ord, PartialEq, PartialOrd)]
enum ConsensusIndel {
    NoIndel,
    Indel(IndelEdit),
}

/// Record the consensus base at pos and the consensus indel between the previous and current positions
#[derive(Clone, Default)]
struct ConsensusInfo {
    base: BTreeMap<u8, f32>,
    indel: BTreeMap<ConsensusIndel, f32>,
}

fn update_indel_entry(
    indel_edit_tracker: Option<&IndelEditTracker>,
    ref_pos: i64,
    read_quality: f32,
    consensus_edit: &mut [ConsensusInfo],
) {
    if let Some(indel_edit_tracker) = indel_edit_tracker {
        let cpos = &mut consensus_edit[indel_edit_tracker.ref_pos as usize];
        let indel_entry = cpos
            .indel
            .entry(ConsensusIndel::Indel(indel_edit_tracker.indel_edit.clone()))
            .or_default();
        *indel_entry += read_quality;
    } else {
        let cpos = &mut consensus_edit[ref_pos as usize];
        let indel_entry = cpos.indel.entry(ConsensusIndel::NoIndel).or_default();
        *indel_entry += read_quality;
    }
}

/// read_quality is used here to break ties on even counts
///
fn update_consensus_edit_from_read_alignment(
    alignment: &SimpleAlignment,
    read_seq: &[u8],
    read_quality: f32,
    consensus_edit: &mut [ConsensusInfo],
) {
    use Cigar::*;
    let mut ref_pos = alignment.ref_offset;
    let mut read_pos = 0;
    let mut indel_edit_tracker = None;
    for c in alignment.cigar.iter() {
        match c {
            Match(len) | Diff(len) | Equal(len) => {
                update_indel_entry(
                    indel_edit_tracker.as_ref(),
                    ref_pos,
                    read_quality,
                    consensus_edit,
                );
                indel_edit_tracker = None;
                for offset in 0..(*len as usize) {
                    let read_offset = read_pos + offset;
                    let allele = read_seq[read_offset];

                    let ref_offset = (ref_pos as usize) + offset;
                    let cpos = &mut consensus_edit[ref_offset];
                    let base_entry = cpos.base.entry(allele).or_default();
                    *base_entry += read_quality;
                    if offset > 0 {
                        let indel_entry = cpos.indel.entry(ConsensusIndel::NoIndel).or_default();
                        *indel_entry += read_quality;
                    }
                }
            }
            Del(len) => {
                if read_pos > 0 && ref_pos > 0 {
                    let indel_edit = &mut indel_edit_tracker
                        .get_or_insert(IndelEditTracker::new(ref_pos))
                        .indel_edit;
                    indel_edit.del_len += *len as usize;
                }
            }
            Ins(len) => {
                if read_pos > 0 && ref_pos > 0 {
                    let indel_edit = &mut indel_edit_tracker
                        .get_or_insert(IndelEditTracker::new(ref_pos))
                        .indel_edit;
                    let ins_seq = &read_seq[read_pos..read_pos + (*len as usize)];
                    indel_edit.ins_seq.extend(ins_seq);
                }
            }
            _ => (),
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
}

fn polish_backbone(
    alignment_weights: &AlignmentWeights,
    trimmed_reads: &[TrimmedReadInfo],
    consensus_seq_info: &mut ConsensusSequenceInfo,
) {
    let candidate_backbone = &mut consensus_seq_info.consensus_seq;
    let source_read_alignment_info = &mut consensus_seq_info.source_read_alignment_info;

    let mut consensus_edit = vec![ConsensusInfo::default(); candidate_backbone.seq.len()];

    let mut aligner = PairwiseAligner::new_glocal_aligner(alignment_weights);
    source_read_alignment_info.clear();
    for trimmed_read_index in candidate_backbone.source_reads.iter().copied() {
        let trimmed_read_info = &trimmed_reads[trimmed_read_index];
        let _ = aligner.align(&candidate_backbone.seq, &trimmed_read_info.trimmed_read);
        let alignment = unwrap!(aligner.get_last_alignment());

        update_consensus_edit_from_read_alignment(
            &alignment,
            &trimmed_read_info.trimmed_read,
            trimmed_read_info.read_quality,
            &mut consensus_edit,
        );

        source_read_alignment_info.push(get_consensus_alignment_info(
            candidate_backbone.seq.len(),
            &alignment,
        ));
    }

    let extra_debug = false;
    if extra_debug {
        eprintln!("Backbone before:");
        print_fasta(250, &[&candidate_backbone.seq]);
    }

    // Apply consensus edits to backbone
    //
    // Indels apply between the previous and present base so are searched first
    //
    candidate_backbone.seq.clear();
    let mut edit_index = 0;
    while edit_index < consensus_edit.len() {
        let consensus_info = &consensus_edit[edit_index];

        let max_indel = consensus_info
            .indel
            .iter()
            .max_by(|&(_, &a), &(_, &b)| a.partial_cmp(&b).expect("Tried to compare a NaN"));
        if let Some((ConsensusIndel::Indel(indel_edit), _)) = max_indel {
            candidate_backbone.seq.extend(&indel_edit.ins_seq);
            if indel_edit.del_len > 0 {
                edit_index += indel_edit.del_len;
                continue;
            }
        }
        let max_base = consensus_info
            .base
            .iter()
            .max_by(|&(_, &a), &(_, &b)| a.partial_cmp(&b).expect("Tried to compare a NaN"));
        if let Some((&base, _)) = max_base {
            candidate_backbone.seq.push(base);
        }

        edit_index += 1;
    }

    if extra_debug {
        eprintln!("Backbone after:");
        print_fasta(250, &[&candidate_backbone.seq]);
    }
}

fn polish_backbones(
    alignment_weights: &AlignmentWeights,
    trimmed_reads: &[TrimmedReadInfo],
    consensus_seq_infos: &mut [ConsensusSequenceInfo],
) {
    for consensus_seq_info in consensus_seq_infos.iter_mut() {
        polish_backbone(alignment_weights, trimmed_reads, consensus_seq_info);
    }
}

fn merge_consensus_sequences(
    alignment_weights: &AlignmentWeights,
    min_norm_score: f64,
    consensus_seq_infos: &mut Vec<ConsensusSequenceInfo>,
) {
    let debug = false;

    let mut aligner = OverlapPairwiseAligner::new_aligner(alignment_weights);
    let mut remove_backbone = vec![false; consensus_seq_infos.len()];
    for check_backbone_index in (1..consensus_seq_infos.len()).rev() {
        let check_source_read_count = consensus_seq_infos[check_backbone_index]
            .consensus_seq
            .source_reads
            .len();
        if check_source_read_count >= 6 {
            // Method is not intended to merge well-supported contigs:
            continue;
        }

        for anchor_backbone_index in 0..check_backbone_index {
            let anchor_source_read_count = consensus_seq_infos[anchor_backbone_index]
                .consensus_seq
                .source_reads
                .len();
            if check_source_read_count >= anchor_source_read_count {
                continue;
            }

            // If either contig has already been polished at consensus 3 or greater, then our expectation for the min_norm_score
            // increases as most errors will be polished out:
            let mut min_norm_score = min_norm_score;
            if anchor_source_read_count > 2 {
                min_norm_score += 0.01;
            }
            if check_source_read_count > 2 {
                min_norm_score += 0.01;
            }

            // Align "check" backbone (that we suspect could be a false split) back to 'anchor' backbone, to see if it only split
            // because it was unpolished
            //
            let score = aligner.align(
                &consensus_seq_infos[anchor_backbone_index].consensus_seq.seq,
                &consensus_seq_infos[check_backbone_index].consensus_seq.seq,
            );
            let alignment = unwrap!(aligner.get_last_alignment());
            let (left_clip, right_clip, read_size) = get_alignment_soft_clips(&alignment.cigar);
            let overlap = read_size - (left_clip + right_clip);
            let align_metrics = BackboneAlignmentQualityMetrics { score, overlap };

            let is_reject = reject_backbone_alignment(min_norm_score, &align_metrics);

            if debug {
                eprintln!(
                    "Check merge of backbone {check_backbone_index} into {anchor_backbone_index}"
                );
                eprintln!(
                    "Score/Overlap/NormScore/IsReject {score}/{overlap}/{}/{is_reject}",
                    get_norm_score(&align_metrics)
                );
                eprintln!("Alignment: {:?}", &alignment);
                aligner.print_pairwise_alignment(
                    &consensus_seq_infos[anchor_backbone_index].consensus_seq.seq,
                    &consensus_seq_infos[check_backbone_index].consensus_seq.seq,
                );
            }

            if is_reject {
                continue;
            }

            remove_backbone[check_backbone_index] = true;
        }
    }

    drop_true(consensus_seq_infos, &remove_backbone);
}

/// This is a general assembly method provided as an alternative to spoa. It is based entirely on sparse
/// pairwise alignments.
///
/// The principal benefit over spoa is on memory usage for very large insertions, but this is also generally
/// faster.
///
/// Just like the spoa assembler (at this time at least) this method suffers from not having true clustering
/// of the top two alleles. It is progressive and greedy.
///
/// It is almost ready to be the assembler for general cases at this point, the performance is just slightly lower
/// than spoa. For now we use it for large insertions only.
///
/// Future ideas:
/// - An actual clustering routine
/// - Sort candidates by supported read count before checking for read membership
///
pub fn assemble_and_cluster_reads(
    refine_settings: &RefineSVSettings,
    cluster_index: usize,
    max_candidate_assemblies: usize,
    max_assembly_count: usize,
    trimmed_reads: Vec<TrimmedReadInfo>,
    is_multi_region: bool,
    assembly_debug: bool,
) -> Vec<AssemblyResult> {
    if assembly_debug {
        eprintln!("starting new backbone assembly for cluster index {cluster_index}");
    }

    // Weights for pairwise alignments used to first build and then polish backbone alignments
    let alignment_weights = AlignmentWeights {
        match_: 1,
        mismatch: -3,
        gap_open: 0,
        gap_extend: -1,
    };

    let min_norm_score =
        refine_settings.min_norm_score_for_read_cluster * alignment_weights.match_ as f64;
    let candidate_backbones = build_insertion_backbones(
        max_candidate_assemblies,
        min_norm_score,
        &alignment_weights,
        &trimmed_reads,
        assembly_debug,
    );

    // If the first clustering doesn't produce any contigs, make a second pass at clustering with
    // a more inclusive clustering criteria. This can help identify SVs in regions with an increased
    // level of local sequencing noise.
    //
    let candidate_backbones = if !candidate_backbones.is_empty()
        && candidate_backbones[0].source_reads.len() < refine_settings.min_assembly_read_support
    {
        let min_norm_score = refine_settings.min_norm_score_for_read_cluster_2nd_pass
            * alignment_weights.match_ as f64;
        build_insertion_backbones(
            max_candidate_assemblies,
            min_norm_score,
            &alignment_weights,
            &trimmed_reads,
            assembly_debug,
        )
    } else {
        candidate_backbones
    };

    // Go ahead and filter backbones by ploidy and read support before starting polish step
    let candidate_backbones = candidate_backbones
        .into_iter()
        .take(max_assembly_count)
        .filter(|x| x.source_reads.len() >= refine_settings.min_assembly_read_support)
        .collect::<Vec<_>>();

    // Reformat into structure designed for adding flanks later
    let mut consensus_seq_infos = {
        let mut x = Vec::new();
        for backbone in candidate_backbones.into_iter() {
            x.push(ConsensusSequenceInfo {
                consensus_seq: backbone,
                source_read_alignment_info: Vec::new(),
            });
        }
        x
    };

    polish_backbones(&alignment_weights, &trimmed_reads, &mut consensus_seq_infos);

    merge_consensus_sequences(&alignment_weights, min_norm_score, &mut consensus_seq_infos);

    add_flanks(
        refine_settings,
        cluster_index,
        max_assembly_count,
        &trimmed_reads,
        is_multi_region,
        consensus_seq_infos,
    )
}
