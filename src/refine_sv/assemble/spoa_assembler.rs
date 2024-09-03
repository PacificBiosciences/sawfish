use std::ffi::CStr;

use spoa::{get_alignment_clip_size, get_alignment_overlap_size, AlignmentEngine};

use super::add_flanks::{
    add_flanks, get_consensus_alignment_info, ConsensusSequence, ConsensusSequenceInfo,
    TrimmedReadConsensusAlignmentInfo,
};
use super::{AssemblyResult, RefineSVSettings};
use crate::refine_sv::trimmed_reads::TrimmedReadInfo;
use crate::spoa_utils::{
    get_alignment_to_graph_consensus, get_pairwise_alignment, get_query_alignment_from_msa,
};

#[derive(Default)]
struct AssemblyAlignmentInfo {
    /// Index of the assembly contig being aligned to
    index: usize,

    /// Alignment score
    score: i32,

    /// Number of read bases used in the alignment (subtracting any clipping on the left and right sides)
    overlap: u32,

    /// Number of read bases clipped from the left side of the alignment
    left_clip: u32,

    alignment: Option<spoa::Alignment>,
}

/// Standard alignment process applied to trimmed reads against all contigs
///
/// Selects the best scoring alignment. More in-depth analysis of the alignment
/// quality for cluster join/split decisions are left to downstream steps
///
fn get_best_assembly_alignment_for_trimmed_read(
    alignment_engine: &mut AlignmentEngine,
    cread: &CStr,
    candidate_assemblies: &[CandidateContigAssembly],
) -> AssemblyAlignmentInfo {
    let mut best_assembly = AssemblyAlignmentInfo::default();
    for (candidate_assembly_index, candidate_assembly) in candidate_assemblies.iter().enumerate() {
        let (alignment, score) = alignment_engine.align(cread, &candidate_assembly.graph);
        let overlap = get_alignment_overlap_size(&alignment);
        let clip = get_alignment_clip_size(&alignment);
        if clip >= 0 && (best_assembly.alignment.is_none() || score > best_assembly.score) {
            best_assembly = AssemblyAlignmentInfo {
                index: candidate_assembly_index,
                score,
                overlap,
                left_clip: clip as u32,
                alignment: Some(alignment),
            };
        }
    }
    best_assembly
}

fn get_norm_score(align_info: &AssemblyAlignmentInfo) -> f64 {
    (align_info.score as f64) / (align_info.overlap as f64)
}

/// Determine if the trimmed read alignment to the given assembly is high enough quality to join the
/// the read into that assembly or if a new assembly should be started from the read instead.
///
/// Selection criteria:
/// (1) minimum overlap size
/// (2) minimum score/overlap size
///
/// Return true if the assembly should be rejected
///
fn reject_assembly_alignment(min_norm_score: f64, align_info: &AssemblyAlignmentInfo) -> bool {
    if align_info.overlap == 0 {
        true
    } else {
        const MIN_OVERLAP: u32 = 100;
        align_info.overlap < MIN_OVERLAP || get_norm_score(align_info) < min_norm_score
    }
}

/// Determine if the read alignment to the POA graph should be accepted.
///
/// This is related to reject_assembly_alignment above, but adds extra checks that are helpful for
/// local alignments.
///
fn reject_assembly_backup_alignment(
    min_norm_score: f64,
    read_size: u32,
    align_info: &AssemblyAlignmentInfo,
) -> bool {
    let read_end_tolerance = 3;
    // Check that read is aligned end to end:
    if align_info.left_clip > read_end_tolerance {
        return true;
    }
    let left_side = align_info.left_clip + align_info.overlap;
    assert!(read_size >= left_side);
    let right_clip = read_size - left_side;
    if right_clip > read_end_tolerance {
        return true;
    }

    reject_assembly_alignment(min_norm_score, align_info)
}

/// Print out reports on a given read assembly alignment
fn debug_assembly_alignment_quality(
    min_norm_score: f64,
    alignment_engine: &mut AlignmentEngine,
    alignment_engine2: &mut AlignmentEngine,
    candidate_assemblies: &mut [CandidateContigAssembly],
    trimmed_read_info: &TrimmedReadInfo,
    cread: &CStr,
    align_info: &AssemblyAlignmentInfo,
) {
    use crate::spoa_utils::{print_fasta, print_msa};
    eprintln!(
        "Best candidate assembly overlap: {} score: {} Norm score: {} is_reject: {} qname: {} read_range: {:?}",
        align_info.overlap, align_info.score, get_norm_score(align_info), reject_assembly_alignment(min_norm_score, align_info), trimmed_read_info.qname, trimmed_read_info.read_range,
    );
    {
        let (score, query_start_clip, overlap, msa) = get_alignment_to_graph_consensus(
            alignment_engine,
            &mut candidate_assemblies[align_info.index].graph,
            cread,
        );
        eprintln!("Overlap alignment to consensus score: {score} query_start_clip: {query_start_clip} overlap: {overlap}\nmsa:\n");
        print_msa(&msa);
    }
    {
        let (score, query_start_clip, overlap, msa) = get_alignment_to_graph_consensus(
            alignment_engine2,
            &mut candidate_assemblies[align_info.index].graph,
            cread,
        );
        eprintln!("SW alignment to consensus score: {score} query_start_clip: {query_start_clip} overlap: {overlap}\nmsa:\n");
        print_msa(&msa);
    }

    let reads = vec![cread.to_owned()];
    eprintln!("Full Read:");
    print_fasta(&reads);
}

/// Process the next trimmed read into the candidate assembly set
///
#[allow(clippy::too_many_arguments)]
fn add_trimmed_read_to_candidate_assembles(
    max_candidate_assemblies: usize,
    min_norm_score: f64,
    ov_alignment_engine: &mut AlignmentEngine,
    sw_alignment_engine: &mut AlignmentEngine,
    candidate_assemblies: &mut Vec<CandidateContigAssembly>,
    trimmed_read_index: usize,
    trimmed_read_info: &TrimmedReadInfo,
    assembly_debug: bool,
) {
    if assembly_debug {
        eprintln!(
            "trimmed read index/size/qname: {trimmed_read_index}/{}/{}",
            &trimmed_read_info.trimmed_read.len(),
            &trimmed_read_info.qname,
        );
    }

    let cread = unsafe { CStr::from_bytes_with_nul_unchecked(&trimmed_read_info.trimmed_read) };

    let (best_assembly_alignment, reject_best_assembly_alignment) = {
        // Try the primary alignment strategy, and switch to a fallback strategy if that gives a
        // poor result:
        //
        let alignment = get_best_assembly_alignment_for_trimmed_read(
            ov_alignment_engine,
            cread,
            candidate_assemblies,
        );

        let reject_alignment = reject_assembly_alignment(min_norm_score, &alignment);

        if !reject_alignment {
            (alignment, reject_alignment)
        } else {
            if assembly_debug {
                eprintln!("Rejecting OV alignment and going to SW");
            }
            let alignment = get_best_assembly_alignment_for_trimmed_read(
                sw_alignment_engine,
                cread,
                candidate_assemblies,
            );

            let read_size = trimmed_read_info.trimmed_read_len() as u32;
            let reject_alignment =
                reject_assembly_backup_alignment(min_norm_score, read_size, &alignment);
            (alignment, reject_alignment)
        }
    };

    if assembly_debug && best_assembly_alignment.alignment.is_some() {
        debug_assembly_alignment_quality(
            min_norm_score,
            ov_alignment_engine,
            sw_alignment_engine,
            candidate_assemblies,
            trimmed_read_info,
            cread,
            &best_assembly_alignment,
        );
    }

    if !reject_best_assembly_alignment {
        if assembly_debug {
            eprintln!(
                "Adding read to candidate assembly {}",
                best_assembly_alignment.index
            );
        }
        candidate_assemblies[best_assembly_alignment.index].push_read(
            trimmed_read_index,
            &best_assembly_alignment.alignment.unwrap(),
            cread,
        );
    } else if candidate_assemblies.len() < max_candidate_assemblies {
        // make new cluster
        if assembly_debug {
            eprintln!(
                "Starting new candidate assembly {}",
                candidate_assemblies.len()
            );
        }
        candidate_assemblies.push(CandidateContigAssembly::new());
        let candidate_assembly = candidate_assemblies.last_mut().unwrap();
        candidate_assembly.push_read(trimmed_read_index, &spoa::Alignment::new(), cread);
    } else {
        // filter read out
        if assembly_debug {
            eprintln!(
                "Skipping new candidate assembly {}",
                candidate_assemblies.len()
            );
        }
    }
}

fn cluster_and_sort_candidate_assemblies(
    max_candidate_assemblies: usize,
    min_norm_score: f64,
    trimmed_reads: &[TrimmedReadInfo],
    ov_alignment_engine: &mut AlignmentEngine,
    sw_alignment_engine: &mut AlignmentEngine,
    assembly_debug: bool,
) -> Vec<CandidateContigAssembly> {
    let mut candidate_assemblies: Vec<CandidateContigAssembly> = Vec::new();

    if assembly_debug {
        eprintln!("Starting assembly round with min_norm_score {min_norm_score}");
    }

    for (trimmed_read_index, trimmed_read_info) in trimmed_reads.iter().enumerate() {
        add_trimmed_read_to_candidate_assembles(
            max_candidate_assemblies,
            min_norm_score,
            ov_alignment_engine,
            sw_alignment_engine,
            &mut candidate_assemblies,
            trimmed_read_index,
            trimmed_read_info,
            assembly_debug,
        );
    }

    let extra_debug = false;
    if extra_debug {
        for (candidate_assembly_index, candidate_assembly) in
            candidate_assemblies.iter_mut().enumerate()
        {
            let count = candidate_assembly.graph.get_sequence_count();
            let con_len = candidate_assembly.graph.consensus().into_bytes().len();
            eprintln!("Candidate assembly {candidate_assembly_index} read count {count} cons length {con_len}");
        }
    }

    // Get top expected clusters by sequence count and convert these into final assembly results:
    candidate_assemblies.sort_by_key(|k| std::cmp::Reverse(k.graph.get_sequence_count()));

    candidate_assemblies
}

/// Align trimmed read back to the core_consensus sequence, with SW this time, since the contig
/// should be a superset of all trimmed reads at this point.
///
/// Return some basic stats on the edges of this alignment
///
fn get_trimmed_read_consensus_alignment_info(
    sw_alignment_engine: &mut AlignmentEngine,
    trimmed_read_info: &TrimmedReadInfo,
    consensus: &CStr,
    contig_size: usize,
) -> TrimmedReadConsensusAlignmentInfo {
    let trimmed_cread =
        unsafe { CStr::from_bytes_with_nul_unchecked(&trimmed_read_info.trimmed_read) };
    let (_score, query_start_clip, _overlap, msa) =
        get_pairwise_alignment(sw_alignment_engine, consensus, trimmed_cread);
    let read_len = (trimmed_read_info.trimmed_read.len() - 1) as u32;
    let alignment = get_query_alignment_from_msa(query_start_clip, read_len, &msa);

    get_consensus_alignment_info(contig_size, &alignment)
}

/// Align all trimmed reads back to the core_consensus, with SW this time, since the contig
/// should be a superset of all trimmed reads at this point.
///
/// Return some basic stats on the edges of these alignments
///
fn get_consensus_alignment_info_for_all_trimmed_reads(
    sw_alignment_engine: &mut AlignmentEngine,
    trimmed_reads: &[TrimmedReadInfo],
    source_reads: &[usize],
    consensus: &CStr,
    consensus_size: usize,
) -> Vec<TrimmedReadConsensusAlignmentInfo> {
    let mut source_read_alignment_info = Vec::new();
    for read_index in source_reads.iter() {
        let trimmed_read_info = &trimmed_reads[*read_index];
        let alignment_info = get_trimmed_read_consensus_alignment_info(
            sw_alignment_engine,
            trimmed_read_info,
            consensus,
            consensus_size,
        );
        source_read_alignment_info.push(alignment_info);
    }
    source_read_alignment_info
}

/// Used to represent each candidate contig during the assembly process
struct CandidateContigAssembly {
    source_reads: Vec<usize>,
    graph: spoa::Graph,
}

impl CandidateContigAssembly {
    fn new() -> Self {
        Self {
            source_reads: Vec::new(),
            graph: spoa::Graph::new(),
        }
    }

    fn push_read(&mut self, read_index: usize, alignment: &spoa::Alignment, cread: &CStr) {
        self.source_reads.push(read_index);
        self.graph.add_alignment_noqual(alignment, cread);
    }
}

struct PoaAlignmentWeights {
    match_: i8,
    mismatch: i8,
    gap_open: i8,
    gap_extend: i8,
}

/// Assemble trimmed reads and return the top N assembly sequences using the spoa library
///
/// This strategy tests alignments one at a time on a set of POAs which are each supposed to represent one haplotype.
/// The candidate haplotypes are created each time a read is not sufficiently similar to one of the existing haplotypes.
///
pub fn assemble_and_cluster_reads(
    refine_settings: &RefineSVSettings,
    cluster_index: usize,
    max_candidate_assemblies: usize,
    max_assembly_count: usize,
    mut trimmed_reads: Vec<TrimmedReadInfo>,
    is_multi_region: bool,
    assembly_debug: bool,
) -> Vec<AssemblyResult> {
    // nul-terminate all reads for use in spoa
    for r in trimmed_reads.iter_mut() {
        r.add_trimmed_read_nul_term();
    }

    // Note that in spoa, gap_extend is only applied to gaps larger than 1. If gap_extend is
    // set lower than gap_open, gap_extend will be ignored and it will switch to linear mode using
    // the gap open score.
    //
    let contig_asm_weights = PoaAlignmentWeights {
        match_: 1,
        mismatch: -3,
        gap_open: -1,
        gap_extend: -1,
    };
    let mut ov_alignment_engine = AlignmentEngine::new(
        spoa::AlignmentType::kOV,
        contig_asm_weights.match_,
        contig_asm_weights.mismatch,
        contig_asm_weights.gap_open,
        contig_asm_weights.gap_extend,
        contig_asm_weights.gap_open,
        contig_asm_weights.gap_extend,
    );
    let mut sw_alignment_engine = AlignmentEngine::new(
        spoa::AlignmentType::kSW,
        contig_asm_weights.match_,
        contig_asm_weights.mismatch,
        contig_asm_weights.gap_open,
        contig_asm_weights.gap_extend,
        contig_asm_weights.gap_open,
        contig_asm_weights.gap_extend,
    );

    // sort trimmed reads from longest to shortest:
    //
    // TODO - this is debatable -- the spoa ovelap mode is really designed to go from left to right
    // Try taking this out
    //trimmed_reads.sort_by_key(|x| std::cmp::Reverse(x.trimmed_read.len()));

    let min_norm_score =
        refine_settings.min_norm_score_for_read_cluster * contig_asm_weights.match_ as f64;
    let candidate_assemblies = cluster_and_sort_candidate_assemblies(
        max_candidate_assemblies,
        min_norm_score,
        &trimmed_reads,
        &mut ov_alignment_engine,
        &mut sw_alignment_engine,
        assembly_debug,
    );

    // If the first clustering doesn't produce any contigs, make a second pass at clustering with
    // a more inclusive clustering criteria. This can help identify SVs in regions with an increased
    // level of local sequencing noise.
    //
    let candidate_assemblies = if !candidate_assemblies.is_empty()
        && candidate_assemblies[0].graph.get_sequence_count()
            < refine_settings.min_assembly_read_support as u32
    {
        let min_norm_score = refine_settings.min_norm_score_for_read_cluster_2nd_pass
            * contig_asm_weights.match_ as f64;
        cluster_and_sort_candidate_assemblies(
            max_candidate_assemblies,
            min_norm_score,
            &trimmed_reads,
            &mut ov_alignment_engine,
            &mut sw_alignment_engine,
            assembly_debug,
        )
    } else {
        candidate_assemblies
    };

    // Convert candidate assembly graph to the consensus sequence format used in all downstream steps:
    //
    let mut consensus_seq_infos = Vec::new();
    for mut candidate_assembly in candidate_assemblies.into_iter() {
        // Assert that graph is consistent with source read count:
        let read_count = candidate_assembly.graph.get_sequence_count() as usize;
        assert_eq!(read_count, candidate_assembly.source_reads.len());

        let core_consensus_cstring = candidate_assembly.graph.consensus();
        let core_consensus = core_consensus_cstring.as_bytes().to_vec();

        let source_read_alignment_info = get_consensus_alignment_info_for_all_trimmed_reads(
            &mut sw_alignment_engine,
            &trimmed_reads,
            &candidate_assembly.source_reads,
            &core_consensus_cstring,
            core_consensus.len(),
        );

        let consensus_seq = ConsensusSequence {
            source_reads: candidate_assembly.source_reads.clone(),
            seq: core_consensus,
        };
        let consensus_seq_info = ConsensusSequenceInfo {
            consensus_seq,
            source_read_alignment_info,
        };
        consensus_seq_infos.push(consensus_seq_info);
    }

    add_flanks(
        refine_settings,
        cluster_index,
        max_assembly_count,
        &trimmed_reads,
        is_multi_region,
        consensus_seq_infos,
    )
}
