mod read_cache;

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::default::Default;
use std::ops::Range;

use rust_htslib::bam;
use rust_vc_utils::{ChromList, GenomeSegment, IntRange, RegionMap};
use strum::EnumCount;

use self::read_cache::ScoringRead;
pub use self::read_cache::{MAX_WORKER_CACHE_BYTES, SampleReadCache};
use crate::bam_utils::{TargetMatchType, get_alignment_closest_to_target_ref_pos};
use crate::breakpoint::{Breakend, BreakendDirection, Breakpoint};
use crate::expected_ploidy::{SVLocusExpectedCNInfo, SVLocusPloidy};
use crate::genome_ref_utils::{GenomeRef, get_ref_segment_seq};
use crate::prob_utils::{error_prob_to_phred, ln_error_prob_to_qphred, normalize_ln_distro};
use crate::refine_sv::{
    AlleleCounts, AlleleType, Genotype, RefinedSV, SVPhaseStatus, SVSampleScoreInfo, SVScoreInfo,
};
use crate::sv_group::{ClusterAssemblyAlignment, GroupHaplotypeId, SVGroup, SVGroupHaplotype};
use crate::wfa2_utils::{
    BoundedScore, ConsensusScoreAligner, PairwiseAligner, print_pairwise_alignment,
};

/// Get the reference allele breakend flank sequence corresponding to the genotype alignment region of the target SV breakend
///
/// # Arguments
/// * `is_breakend1` - Use only for debug messages
///
/// Return the breakend flank sequence, or None if the flank region is too small (ie. trying to sample from the chromosome edge)
///
fn convert_breakend_to_ref_flank(
    score_settings: &ScoreSVSettings,
    breakend: &Breakend,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    is_breakend1: bool,
    debug: bool,
) -> Option<Vec<u8>> {
    let mut flank_segment = breakend.segment.clone();
    match breakend.dir {
        BreakendDirection::LeftAnchor => {
            flank_segment.range.start = flank_segment.range.end;
            flank_segment.asymmetric_expand_by(
                chrom_list,
                0,
                score_settings.read_support_flank_size,
            );
        }
        BreakendDirection::RightAnchor => {
            flank_segment.range.end = flank_segment.range.start + 1;
            flank_segment.asymmetric_expand_by(
                chrom_list,
                score_settings.read_support_flank_size,
                0,
            );
        }
    }

    if debug {
        let bel = if is_breakend1 { '1' } else { '2' };
        eprintln!("breakend{bel}: {breakend:?}");
        eprintln!("Ref breakend{bel} flank segment: {flank_segment:?}");
    }

    if flank_segment.range.start < 0
        || flank_segment.range.end <= flank_segment.range.start
        || flank_segment.range.size() < score_settings.min_read_support_flank_size
    {
        None
    } else {
        Some(get_ref_segment_seq(reference, &flank_segment))
    }
}

/// Get the reference allele breakpoint flank sequences corresponding to the genotype alignment regions of the target SV breakpoint
///
/// Returns the flank sequence for breakend1 and breakend2
///
fn get_ref_flanks(
    score_settings: &ScoreSVSettings,
    bp: &Breakpoint,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    debug: bool,
) -> [Option<Vec<u8>>; 2] {
    let breakend1_flank = convert_breakend_to_ref_flank(
        score_settings,
        &bp.breakend1,
        reference,
        chrom_list,
        true,
        debug,
    );

    let breakend2 = bp.breakend2.as_ref().unwrap();
    let breakend2_flank = convert_breakend_to_ref_flank(
        score_settings,
        breakend2,
        reference,
        chrom_list,
        false,
        debug,
    );

    [breakend1_flank, breakend2_flank]
}

/// Get the reference allele flank sequences from the flanks returned by get_ref_flanks
fn get_ref_flank_seqs(ref_flanks: &[Option<Vec<u8>>; 2]) -> BreakpointAlleleFlankSeqs<'_> {
    // We don't need any alternate registrations for reference sequence, so leave these as None
    //
    BreakpointAlleleFlankSeqs {
        breakend_flanks: [
            [ref_flanks[0].as_deref(), None],
            [ref_flanks[1].as_deref(), None],
        ],
    }
}

/// Get haplotype contig segment sequence corresponding to the breakend flank
///
/// In this case we have a contig and contig alignment for the haplotype of interest, so the logic
/// should be similar to extracting flank sequences from reads
///
/// # Arguments
///
/// * `assembly_segment` - This is optionally passed in to attempt to get the flank sequence by
///   registering the alignment point at the edge of the assembly region instead of the edge of the
///   candidate SV. This can be important to rescue some complex cases where the candidate SV
///   alignment is very different than the mapped read alignment.
///
/// * `hap_id`` - Haplotype id is only used to improve debugging messages
///
/// * `is_breakend1` - Only used to improve debug messages
///
fn convert_breakend_to_haplotype_contig_flank_seq<'a>(
    score_settings: &ScoreSVSettings,
    breakend: &Breakend,
    assembly_segment: Option<&GenomeSegment>,
    hap_id: &GroupHaplotypeId,
    breakend_contig_alignment: &'a ClusterAssemblyAlignment,
    is_breakend1: bool,
    debug: bool,
) -> Option<&'a [u8]> {
    let homology_len = breakend.get_homology_len();
    let contig_len = breakend_contig_alignment.contig_seq.len();

    // Check if flank range falls off either end of the contig or falls below min size, and if so
    // return None
    //
    let get_flank_range = |start: i64, end: i64| {
        if start < 0
            || end >= (contig_len as i64)
            || (end - start) < score_settings.min_read_support_flank_size
        {
            None
        } else {
            Some((start as usize)..(end as usize))
        }
    };

    let contig_alignment = &breakend_contig_alignment.contig_alignment;

    let flank_range = match breakend.dir {
        BreakendDirection::LeftAnchor => {
            let ref_target_pos = breakend.segment.range.start;
            let ref_to_contig_anchor_pos = {
                match assembly_segment {
                    Some(region) => region.range.start,
                    None => ref_target_pos,
                }
            };
            let read_to_ref_map = get_alignment_closest_to_target_ref_pos(
                contig_alignment.ref_offset,
                &contig_alignment.cigar,
                ref_to_contig_anchor_pos,
                TargetMatchType::TargetOrLower,
            )?;

            // Make the simple mapping adjustment if the closest ref pos is offset from the target
            let shift = ref_target_pos - read_to_ref_map.ref_pos;
            let shifted_read_pos = read_to_ref_map.read_pos as i64 + shift;

            // Add one to move to the position immediately after the breakend, then account for
            // exact homology len to get the intended flank region start position
            let flank_start = shifted_read_pos + 1 + homology_len;

            get_flank_range(
                flank_start,
                flank_start + score_settings.read_support_flank_size,
            )
        }
        BreakendDirection::RightAnchor => {
            let ref_target_pos = breakend.segment.range.end;
            let ref_to_contig_anchor_pos = {
                match assembly_segment {
                    Some(region) => region.range.end,
                    None => ref_target_pos,
                }
            };
            let read_to_ref_map = get_alignment_closest_to_target_ref_pos(
                contig_alignment.ref_offset,
                &contig_alignment.cigar,
                ref_to_contig_anchor_pos,
                TargetMatchType::TargetOrHigher,
            )?;

            // Make the simple mapping adjustment if the closest ref pos is offset from the target
            let shift = ref_target_pos - read_to_ref_map.ref_pos;
            let shifted_read_pos = read_to_ref_map.read_pos as i64 + shift;

            // Account for exact homology len to get the intended flank region end position
            let flank_end = shifted_read_pos - homology_len;

            get_flank_range(
                flank_end - score_settings.read_support_flank_size,
                flank_end,
            )
        }
    };
    if let Some(flank_range) = flank_range {
        if debug {
            let bel = if is_breakend1 { '1' } else { '2' };
            eprintln!("Haplotype {hap_id:?} breakend{bel} flank range: {flank_range:?}");
        }
        Some(&breakend_contig_alignment.contig_seq[flank_range])
    } else {
        None
    }
}

/// Get contig subsequence corresponding to both breakends of the given alt allele breakpoint
///
/// # Arguments
///
/// * `assembly_segment` - Reference coordinate bounds of interesting read evidence that lead to this
///   region being nominated for SV calling. This is used to pull back to the left or right end of the
///   entire SV affected region to more accurately extract the correct segment.
///
fn get_contig_flank_seqs<'a>(
    score_settings: &ScoreSVSettings,
    bp: &Breakpoint,
    assembly_segment: &GenomeSegment,
    haplotype_info: &'a SVGroupHaplotype,
    debug: bool,
) -> BreakpointAlleleFlankSeqs<'a> {
    let breakend1_contig_alignment = &haplotype_info.contig_info.contig_alignments[0];
    let contig_size = breakend1_contig_alignment.contig_seq.len();
    if debug {
        eprintln!("Contig size: {contig_size}");
    }

    let breakend1_flank = convert_breakend_to_haplotype_contig_flank_seq(
        score_settings,
        &bp.breakend1,
        None,
        &haplotype_info.hap_id,
        breakend1_contig_alignment,
        true,
        debug,
    );
    let breakend1_ext_flank = convert_breakend_to_haplotype_contig_flank_seq(
        score_settings,
        &bp.breakend1,
        Some(assembly_segment),
        &haplotype_info.hap_id,
        breakend1_contig_alignment,
        true,
        debug,
    );

    let breakend2_contig_alignment = {
        if haplotype_info.contig_info.single_region_refinement() {
            breakend1_contig_alignment
        } else {
            &haplotype_info.contig_info.contig_alignments[1]
        }
    };

    let breakend2 = bp.breakend2.as_ref().unwrap();
    let breakend2_flank = convert_breakend_to_haplotype_contig_flank_seq(
        score_settings,
        breakend2,
        None,
        &haplotype_info.hap_id,
        breakend2_contig_alignment,
        false,
        debug,
    );
    let breakend2_ext_flank = convert_breakend_to_haplotype_contig_flank_seq(
        score_settings,
        breakend2,
        Some(assembly_segment),
        &haplotype_info.hap_id,
        breakend2_contig_alignment,
        false,
        debug,
    );

    BreakpointAlleleFlankSeqs {
        breakend_flanks: [
            [breakend1_flank, breakend1_ext_flank],
            [breakend2_flank, breakend2_ext_flank],
        ],
    }
}

#[derive(Clone, Debug)]
struct ReadFlankInfo {
    pub range: Range<usize>,

    /// How much of the flank length is truncated due to the start of the read
    start_truncation: usize,

    /// How much of the flank length is truncated due to the end of the read
    end_truncation: usize,
}

/// Get read segment corresponding to breakend flank
///
/// * assembly_segment - This is optionally passed in to attempt to get the read flank by registering
///   the alignment point at the edge of the assembly region instead of the edge of the candidate SV.
///   This can be important to rescue some complex cases where the candidate SV alignment is very
///   different than the mapped read alignment.
///
fn get_read_breakend_flank_range(
    score_settings: &ScoreSVSettings,
    read: &ScoringRead,
    breakend: &Breakend,
    assembly_segment: Option<&GenomeSegment>,
) -> Option<ReadFlankInfo> {
    // TODO: Could we use the read to contig mapping when this simple process doesn't work?
    let ref_to_read_anchor_offset = 10;
    let homology_len = breakend.get_homology_len();

    // Check if flank range falls below the minimum size, and if so return None
    //
    let get_flank_range = |start: i64, end: i64, start_truncation: usize, end_truncation: usize| {
        if (end - start) < score_settings.min_read_support_flank_size {
            None
        } else {
            Some(ReadFlankInfo {
                range: start as usize..end as usize,
                start_truncation,
                end_truncation,
            })
        }
    };

    match breakend.dir {
        BreakendDirection::LeftAnchor => {
            let ref_target_pos = breakend.segment.range.start;
            let ref_to_read_anchor_pos = {
                match assembly_segment {
                    Some(region) => {
                        if region.range.start >= (ref_target_pos - ref_to_read_anchor_offset) {
                            return None;
                        }
                        region.range.start
                    }
                    None => ref_target_pos,
                }
            };
            let ref_to_read_offset_anchor_pos = ref_to_read_anchor_pos - ref_to_read_anchor_offset;

            let read_to_ref_map = get_alignment_closest_to_target_ref_pos(
                read.pos,
                &read.cigar,
                ref_to_read_offset_anchor_pos,
                TargetMatchType::TargetOrLower,
            )?;

            // Make the simple mapping adjustment if the closest ref pos is offset from the target
            let shift = ref_target_pos - read_to_ref_map.ref_pos;
            let shifted_read_pos = read_to_ref_map.read_pos as i64 + shift;

            if shifted_read_pos + 1 >= score_settings.read_support_anchor_min_edge_distance {
                // Add 1 to move to the position immediately after the breakend, then add exact
                // homology len to the read position to get the intended flank region start position
                let flank_start = shifted_read_pos + 1 + homology_len;

                let read_len = read.seq_len() as i64;
                let full_flank_end = flank_start + score_settings.read_support_flank_size;
                let (flank_end, end_truncation) = if full_flank_end > read_len {
                    (read_len, (full_flank_end - read_len) as usize)
                } else {
                    (full_flank_end, 0)
                };
                get_flank_range(flank_start, flank_end, 0, end_truncation)
            } else {
                None
            }
        }
        BreakendDirection::RightAnchor => {
            let ref_target_pos = breakend.segment.range.end;
            let ref_to_read_anchor_pos = {
                match assembly_segment {
                    Some(region) => {
                        if region.range.end <= (ref_target_pos + ref_to_read_anchor_offset) {
                            return None;
                        }
                        region.range.end
                    }
                    None => ref_target_pos,
                }
            };

            let ref_to_read_offset_anchor_pos = ref_to_read_anchor_pos + ref_to_read_anchor_offset;

            let read_to_ref_map = get_alignment_closest_to_target_ref_pos(
                read.pos,
                &read.cigar,
                ref_to_read_offset_anchor_pos,
                TargetMatchType::TargetOrHigher,
            )?;

            // Make the simple mapping adjustment if the closest ref pos is offset from the target
            let shift = ref_target_pos - read_to_ref_map.ref_pos;
            let shifted_read_pos = read_to_ref_map.read_pos as i64 + shift;

            let read_len = read.seq_len() as i64;
            if (read_len - shifted_read_pos) >= score_settings.read_support_anchor_min_edge_distance
            {
                // Account for exact homology len to get the intended flank region end position
                let flank_end = shifted_read_pos - homology_len;
                let full_flank_start = flank_end - score_settings.read_support_flank_size;
                let (flank_start, start_truncation) = if full_flank_start < 0 {
                    (0, (-full_flank_start) as usize)
                } else {
                    (full_flank_start, 0)
                };
                get_flank_range(flank_start, flank_end, start_truncation, 0)
            } else {
                None
            }
        }
    }
}

/// Range of flanking sequences around one breakend of a given allele breakpoint
///
/// The range format assumes a single source sequence, like a read or assembly contig
///
/// Sometimes multiple ranges are generated for one read, reflecting that we aren't sure exactly
/// which one is an accurate flank extraction
///
type BreakendAlleleFlankRange = Vec<Option<ReadFlankInfo>>;

/// Ranges of flanking sequences around a breakpoint for a given allele
///
/// The range format assumes a single source sequence, like a read or assembly contig
///
/// The 2 ranges correspond to the two breakends of a single breakpoint
///
#[derive(Debug)]
struct BreakpointAlleleFlankRanges {
    breakend_flanks: [BreakendAlleleFlankRange; BreakendFlankType::COUNT],
}

impl BreakpointAlleleFlankRanges {
    fn is_empty(&self) -> bool {
        self.breakend_flanks.iter().flatten().flatten().count() == 0
    }
}

type FlankSeq<'a> = Option<&'a [u8]>;

/// One or more flanking sequences around one breakend of a given allele breakpoint
///
/// Sometimes multiple sequences are generated with different registration methods, reflecting that
/// we aren't sure exactly which one is the most accurate flank extraction
///
type BreakendAlleleFlankSeq<'a> = [FlankSeq<'a>; BreakendFlankRegistrationType::COUNT];

/// Flanking sequences around a breakpoint for a given allele
///
/// The breakend_flanks array here has size 2, and represents breakend 1 and 2 of the breakpoint
///
struct BreakpointAlleleFlankSeqs<'a> {
    breakend_flanks: [BreakendAlleleFlankSeq<'a>; BreakendFlankType::COUNT],
}

struct OverlappingHaplotypeFlankSeqs<'a> {
    haplotype_index: usize,
    flank_seqs: BreakpointAlleleFlankSeqs<'a>,
}

struct BreakendMatchSplitReadInfo {
    seq_order_read_split_segments: Vec<GenomeSegment>,
    primary_split_segment_index: usize,
}

fn get_segment_breakend_distance(segment: &GenomeSegment, breakend: &Breakend) -> Option<i64> {
    if segment.chrom_index != breakend.segment.chrom_index {
        None
    } else {
        let dist = match breakend.dir {
            BreakendDirection::LeftAnchor => {
                (segment.range.end - breakend.segment.range.start).abs()
            }
            BreakendDirection::RightAnchor => {
                (segment.range.start - breakend.segment.range.start).abs()
            }
        };
        Some(dist)
    }
}

/// Return true if the split read segment at `segment_index` has a distance to the breakend that is not higher than that of any other split read segment
///
fn is_segment_breakend_match(
    segment_index: usize,
    breakend: &Breakend,
    split_read_info: &BreakendMatchSplitReadInfo,
) -> bool {
    let segment_breakend_dist = split_read_info
        .seq_order_read_split_segments
        .iter()
        .map(|x| get_segment_breakend_distance(x, breakend))
        .collect::<Vec<_>>();
    if let Some(segment_dist) = segment_breakend_dist[segment_index] {
        segment_breakend_dist
            .iter()
            .flatten()
            .all(|&x| x >= segment_dist)
    } else {
        false
    }
}

/// Check if BAM read segment from a split read alignment matches the given breakend
///
fn does_read_segment_match_breakend(
    read: &ScoringRead,
    target_breakend: &Breakend,
    mate_breakend: &Breakend,
    split_read_info: &BreakendMatchSplitReadInfo,
) -> bool {
    // Method:
    // For primary alignment, check that same-chrom end (left-anchored target_breakend) or pos (left-anchored target_breakend) of this segment is closer than any other segment.
    //
    // Find the adjacent split read segment corresponding to the other side of the breakend, check that chrom id matches, and make the similar 'closest position' check on that breakend.

    let target_segment_index = split_read_info.primary_split_segment_index;
    let is_primary_segment_match =
        is_segment_breakend_match(target_segment_index, target_breakend, split_read_info);
    if !is_primary_segment_match {
        return false;
    }

    // Find the adjacent split read segment corresponding to the other side of this breakend
    let is_left_anchor = target_breakend.dir == BreakendDirection::LeftAnchor;
    let is_plus_one_index = is_left_anchor ^ read.is_reverse;
    let mate_segment_index = if is_plus_one_index {
        // One segment downstream on the read
        if target_segment_index + 1 < split_read_info.seq_order_read_split_segments.len() {
            target_segment_index + 1
        } else {
            return false;
        }
    } else {
        // One segment upsteam on the read
        if target_segment_index > 0 {
            target_segment_index - 1
        } else {
            return false;
        }
    };

    is_segment_breakend_match(mate_segment_index, mate_breakend, split_read_info)
}

///
/// # Arguments
///
/// * `split_read_info` - split alignment segments for the bam record in read order, this is used to check the match of this split segment to the
///   breakend. If this vector is empty no check is needed.
///
fn get_read_flanks_for_breakend(
    score_settings: &ScoreSVSettings,
    read: &ScoringRead,
    target_breakend: &Breakend,
    mate_breakend: &Breakend,
    assembly_segment: &GenomeSegment,
    split_read_info: Option<&BreakendMatchSplitReadInfo>,
) -> BreakendAlleleFlankRange {
    if split_read_info.is_none()
        || does_read_segment_match_breakend(
            read,
            target_breakend,
            mate_breakend,
            split_read_info.unwrap(),
        )
    {
        let breakend_flank =
            get_read_breakend_flank_range(score_settings, read, target_breakend, None);
        let breakend_ext_flank = get_read_breakend_flank_range(
            score_settings,
            read,
            target_breakend,
            Some(assembly_segment),
        );
        vec![breakend_flank, breakend_ext_flank]
    } else {
        vec![None, None]
    }
}

/// Get read segment corresponding to both breakends of the given alt allele breakpoint
///
/// # Arguments
///
/// * `assembly_segment` - Bounds of interesting read evidence that lead to this region being
///   nominated for SV calling. This is used to pull back to the left or right end of the entire SV
///   affected region to more accurately extract the correct segment.
///
/// * `single_region_refinement` - True for single-region (indel-like) scoring regions
///
fn get_read_flanks_for_breakpoint(
    score_settings: &ScoreSVSettings,
    read: &ScoringRead,
    bp: &Breakpoint,
    assembly_segment: &GenomeSegment,
    single_region_refinement: bool,
    debug: bool,
) -> BreakpointAlleleFlankRanges {
    // Confirm that the flanks are extracted from the correct split read segment, except for
    // single-region (indel-like) SVs
    let split_read_info = if single_region_refinement {
        None
    } else {
        read.split_read_info.as_ref()
    };

    let breakend1 = &bp.breakend1;
    let breakend2 = bp.breakend2.as_ref().unwrap();

    let breakend1_flanks = get_read_flanks_for_breakend(
        score_settings,
        read,
        breakend1,
        breakend2,
        assembly_segment,
        split_read_info,
    );
    let breakend2_flanks = get_read_flanks_for_breakend(
        score_settings,
        read,
        breakend2,
        breakend1,
        assembly_segment,
        split_read_info,
    );

    if debug {
        eprintln!("Evaluating support from read {}", read.qname);
        eprintln!("read_breakend1_flank_range: {:?}", breakend1_flanks[0]);
        eprintln!("read_breakend1_ext_flank_range: {:?}", breakend1_flanks[1]);
        eprintln!("read_breakend2_flank_range: {:?}", breakend2_flanks[0]);
        eprintln!("read_breakend2_ext_flank_range: {:?}", breakend2_flanks[1]);
    }

    BreakpointAlleleFlankRanges {
        breakend_flanks: [breakend1_flanks, breakend2_flanks],
    }
}

/// Hasher for the flank alignment cache, faster than the default for long byte keys
#[derive(Default)]
struct FlankHasher(u64);

impl std::hash::Hasher for FlankHasher {
    fn write(&mut self, bytes: &[u8]) {
        const K: u64 = 0x517c_c1b7_2722_0a95;
        let (chunks, remainder) = bytes.as_chunks::<8>();
        for &c in chunks {
            let x = u64::from_le_bytes(c);
            self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(K);
        }
        for &b in remainder {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(K);
        }
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

type FlankHashMap<K, V> = HashMap<K, V, std::hash::BuildHasherDefault<FlankHasher>>;

/// Scores read flank alignments to allele flanks, reusing the score of any alignment already
/// computed for the same SV group
///
/// The same alignment comes up repeatedly: both registrations of a read flank are often the same
/// read segment, aligned to the same reference flank, and SVs of a group often share a breakend.
/// Scores are cached by the sequences aligned, so they are exactly those of a new alignment.
///
struct FlankAlignmentScorer {
    /// Taken from CONSENSUS_SCORE_ALIGNER, and returned to it on drop
    aligner: Option<ConsensusScoreAligner>,
    /// Current read flank, aligned to each allele flank
    text: Vec<u8>,
    text_id: usize,
    /// True if the aligner's text is the current read flank
    is_aligner_text_set: bool,
    text_ids: FlankHashMap<Vec<u8>, usize>,
    pattern_ids: FlankHashMap<Vec<u8>, usize>,
    /// Exact scores, and upper bounds of scores not computed exactly
    scores: FlankHashMap<(usize, usize), BoundedScore>,
    /// Bytes of sequence in the cache
    cache_bytes: usize,
}

thread_local! {
    /// Each thread's aligner, reused between SV groups to avoid setting up WFA2's memory each time
    static CONSENSUS_SCORE_ALIGNER: std::cell::RefCell<Option<ConsensusScoreAligner>> =
        const { std::cell::RefCell::new(None) };
}

impl Drop for FlankAlignmentScorer {
    fn drop(&mut self) {
        CONSENSUS_SCORE_ALIGNER.with(|x| *x.borrow_mut() = self.aligner.take());
    }
}

impl FlankAlignmentScorer {
    /// The cache is cleared when it holds more sequence than this
    const MAX_CACHE_BYTES: usize = 16 * 1024 * 1024;

    fn new() -> Self {
        let aligner = CONSENSUS_SCORE_ALIGNER
            .with(|x| x.borrow_mut().take())
            .unwrap_or_else(ConsensusScoreAligner::new);
        Self {
            aligner: Some(aligner),
            text: Vec::new(),
            text_id: 0,
            is_aligner_text_set: false,
            text_ids: FlankHashMap::default(),
            pattern_ids: FlankHashMap::default(),
            scores: FlankHashMap::default(),
            cache_bytes: 0,
        }
    }

    fn intern(
        ids: &mut FlankHashMap<Vec<u8>, usize>,
        cache_bytes: &mut usize,
        seq: &[u8],
    ) -> usize {
        if let Some(&id) = ids.get(seq) {
            id
        } else {
            let id = ids.len();
            ids.insert(seq.to_vec(), id);
            *cache_bytes += seq.len();
            id
        }
    }

    fn set_text(&mut self, text: &[u8]) {
        if self.cache_bytes > Self::MAX_CACHE_BYTES {
            self.text_ids.clear();
            self.pattern_ids.clear();
            self.scores.clear();
            self.cache_bytes = 0;
        }
        self.text.clear();
        self.text.extend_from_slice(text);
        self.text_id = Self::intern(&mut self.text_ids, &mut self.cache_bytes, text);
        self.is_aligner_text_set = false;
    }

    fn text(&self) -> &[u8] {
        &self.text
    }

    fn get_aligner(&mut self) -> &mut ConsensusScoreAligner {
        let aligner = self.aligner.as_mut().unwrap();
        if !self.is_aligner_text_set {
            aligner.set_text(&self.text);
            self.is_aligner_text_set = true;
        }
        aligner
    }

    /// Alignment score of pattern to the current read flank
    fn align_score(&mut self, pattern: &[u8]) -> i32 {
        let pattern_id = Self::intern(&mut self.pattern_ids, &mut self.cache_bytes, pattern);
        let key = (self.text_id, pattern_id);
        let score = match self.scores.get(&key) {
            Some(&BoundedScore::Exact(x)) => return x,
            Some(&BoundedScore::AtMost(_)) => self.get_aligner().dp_align_score(pattern),
            None => self.get_aligner().align_score(pattern),
        };
        self.scores.insert(key, BoundedScore::Exact(score));
        score
    }

    /// Alignment score of pattern to the current read flank, or an upper bound if computing it
    /// needs the DP matrix
    fn align_score_bounded(&mut self, pattern: &[u8]) -> BoundedScore {
        let pattern_id = Self::intern(&mut self.pattern_ids, &mut self.cache_bytes, pattern);
        let key = (self.text_id, pattern_id);
        if let Some(&score) = self.scores.get(&key) {
            return score;
        }
        let score = self.get_aligner().align_score_bounded(pattern);
        self.scores.insert(key, score);
        score
    }
}

/// # Arguments
/// * `label` - Only used for debug messages
///
/// * `breakend_index` - Only used for debug messagesq
///
fn get_read_flank_align_score(
    scorer: &mut FlankAlignmentScorer,
    seq: &[u8],
    label: &str,
    breakend_index: usize,
    debug: bool,
) -> f64 {
    let score = scorer.align_score(seq);
    let norm_score = score as f64 / scorer.text().len() as f64;

    if debug {
        eprintln!(
            "Read to {label} breakend{breakend_index} score: {norm_score}. Alignment (read above):"
        );
        let mut aligner = PairwiseAligner::new_consensus_aligner();
        aligner.set_text(scorer.text());
        aligner.align(seq);
        print_pairwise_alignment(&aligner, seq);
    }

    norm_score
}

/// # Arguments
/// * `label` - Only used for debug messages
///
/// * `breakend_index` - Only used for debug messages
///
#[allow(clippy::too_many_arguments)]
fn get_read_flank_align_score_option(
    scorer: &mut FlankAlignmentScorer,
    seq: Option<&[u8]>,
    label: &str,
    breakend_index: usize,
    read_start_truncation: usize,
    read_end_truncation: usize,
    min_read_support_flank_size: i64,
    debug: bool,
) -> Option<f64> {
    let seq = get_truncated_flank_seq(
        seq,
        read_start_truncation,
        read_end_truncation,
        min_read_support_flank_size,
    )?;
    Some(get_read_flank_align_score(
        scorer,
        seq,
        label,
        breakend_index,
        debug,
    ))
}

/// Get the part of an allele flank sequence aligned to a read flank with the given truncation,
/// or None if it is too short
fn get_truncated_flank_seq(
    seq: Option<&[u8]>,
    read_start_truncation: usize,
    read_end_truncation: usize,
    min_read_support_flank_size: i64,
) -> Option<&[u8]> {
    let seq = seq?;
    if (read_start_truncation + read_end_truncation + min_read_support_flank_size as usize)
        > seq.len()
    {
        return None;
    }
    Some(&seq[read_start_truncation..seq.len() - read_end_truncation])
}

/// Get the normalized alignment scores of the current read flank to the ref and alt allele flanks,
/// and the higher of the two scores that are exact
///
/// These two scores are only used to find the best allele for this read flank: they're compared
/// to each other and to the score of an overlapping haplotype's flank, the best score must be at
/// least zero, and ties for the best score are handled specially. So a score that is lower than
/// the other score can't change the result, provided it stays lower.
///
/// The scores of read flanks that don't match an allele need the DP matrix (see
/// ConsensusScoreAligner), which is slow. When such a score's upper bound is already lower than the
/// other allele's exact score, the upper bound is used as its value instead.
///
fn get_ref_alt_flank_scores(
    scorer: &mut FlankAlignmentScorer,
    ref_flank: Option<&[u8]>,
    alt_flank: Option<&[u8]>,
) -> (Option<f64>, Option<f64>, Option<i32>) {
    let ref_bound = ref_flank.map(|x| scorer.align_score_bounded(x));
    let alt_bound = alt_flank.map(|x| scorer.align_score_bounded(x));
    let exact = |x: Option<BoundedScore>| match x {
        Some(BoundedScore::Exact(x)) => Some(x),
        _ => None,
    };

    // A score known to be lower than the other's exact score, or else the exact score
    let mut resolve =
        |bound: Option<BoundedScore>, flank: Option<&[u8]>, other_exact: Option<i32>| match bound {
            Some(BoundedScore::Exact(x)) => Some((x, true)),
            Some(BoundedScore::AtMost(x)) if other_exact.is_some_and(|y| x < y) => Some((x, false)),
            Some(BoundedScore::AtMost(_)) => Some((scorer.align_score(flank.unwrap()), true)),
            None => None,
        };
    let ref_score = resolve(ref_bound, ref_flank, exact(alt_bound));
    let ref_exact = ref_score.and_then(|(x, is_exact)| is_exact.then_some(x));
    let alt_score = resolve(alt_bound, alt_flank, ref_exact);

    let alt_exact = alt_score.and_then(|(x, is_exact)| is_exact.then_some(x));
    let text_len = scorer.text().len() as f64;
    (
        ref_score.map(|(x, _)| x as f64 / text_len),
        alt_score.map(|(x, _)| x as f64 / text_len),
        ref_exact.max(alt_exact),
    )
}

/// Get the normalized alignment score of the current read flank to an overlapping haplotype's
/// allele flank, when this score isn't averaged (see `are_overlap_scores_averaged`)
///
/// The score is then only used to find the best allele for this read flank, together with the ref
/// and alt scores, as in `get_ref_alt_flank_scores`. So if it needs the DP matrix and its upper bound
/// is already lower than the ref or alt allele's exact score, the upper bound is used instead.
///
fn get_overlap_flank_score(
    scorer: &mut FlankAlignmentScorer,
    flank: &[u8],
    max_ref_alt_exact_score: Option<i32>,
) -> f64 {
    let score = match scorer.align_score_bounded(flank) {
        BoundedScore::Exact(x) => x,
        BoundedScore::AtMost(x) if max_ref_alt_exact_score.is_some_and(|y| x < y) => x,
        BoundedScore::AtMost(_) => scorer.align_score(flank),
    };
    score as f64 / scorer.text().len() as f64
}

/// True if the overlapping haplotype scores of this sample are averaged to select a guest
/// haplotype, see `select_overlapping_guest_haplotype_index`
///
fn is_overlapping_guest_haplotype_selected(
    treat_single_copy_as_haploid: bool,
    sv_group: &SVGroup,
    sample_index: usize,
) -> bool {
    let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];
    let max_haplotype_count =
        match sv_group.sample_expected_cn_info[sample_index].ploidy(treat_single_copy_as_haploid) {
            SVLocusPloidy::Diploid => 2,
            SVLocusPloidy::Haploid => 1,
        };
    let sv_group_haplotype_count = sv_group.group_haplotypes.len();
    sample_haplotype_list.len() == 1 && max_haplotype_count > 1 && sv_group_haplotype_count > 1
}

/// Scan one search segment to process read support for the SV allele
///
/// # Arguments
///
/// * `assembly_segment` - The original assembly segment associated with this SV
///
/// * `search_segment` - An expansion of target_segment, intended for setting the range of read scanning
///
/// * `are_overlap_scores_averaged` - True if overlapping haplotype scores must be exact, see
///   `get_overlap_flank_score`
///
/// * `scores` - All read support information indexed on read qname, all read scores are appended to this structure.
///
#[allow(clippy::too_many_arguments)]
fn get_read_support_scores_from_segment(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    read_cache: &mut SampleReadCache,
    target_sv: &RefinedSV,
    ref_flank_seqs: &BreakpointAlleleFlankSeqs,
    contig_flank_seqs: &BreakpointAlleleFlankSeqs,
    overlapping_haplotypes_flank_seqs: &[OverlappingHaplotypeFlankSeqs],
    assembly_segment: &GenomeSegment,
    search_segment: &GenomeSegment,
    are_overlap_scores_averaged: bool,
    scorer: &mut FlankAlignmentScorer,
    scores: &mut LocusSupportScores,
    debug: bool,
) {
    let standard_registration_index = BreakendFlankRegistrationType::Standard as usize;

    let single_region_refinement = target_sv.single_region_refinement;
    let bp = &target_sv.bp;

    let mut read_flank_buffer = Vec::new();
    read_cache.for_each_read(
        score_settings,
        reference,
        chrom_list,
        bam_reader,
        search_segment,
        |scoring_read| {
            let read_flanks = get_read_flanks_for_breakpoint(
                score_settings,
                scoring_read,
                bp,
                assembly_segment,
                single_region_refinement,
                debug,
            );
            if read_flanks.is_empty() {
                return;
            }

            let read_scores = scores.entry(scoring_read.qname.clone()).or_default();

            for breakend_index in 0..2 {
                // This filter is added for cases where the same qname is encountered multiple times as
                // both primary and supplementary alignments, as would be typical in a duplication
                // (N=2+), large non-dup insertion (N=2), inversion or other complex SV type.
                //
                // Only the reference allele is checked because the ref and alt alleles should always be
                // filled in as a pair and the intent of subsequent read updates is not to fill in a
                // 'better' version of the allele score if it's defined in both cases, but rather just
                // to fill in any scores that were missing.
                //
                if read_scores.ref_allele.breakend_flanks[breakend_index].is_any_condition_scored()
                {
                    continue;
                }

                for (registration_index, breakend_flank) in read_flanks.breakend_flanks
                    [breakend_index]
                    .iter()
                    .enumerate()
                    .filter_map(|(i, x)| x.clone().map(|x| (i, x)))
                {
                    scoring_read.get_seq(breakend_flank.range, &mut read_flank_buffer);
                    let read_breakend_flank = read_flank_buffer.as_slice();

                    scorer.set_text(read_breakend_flank);

                    // This activates more detailed scoring output than we normally want from standard
                    // scoring debug.
                    let debug_score_alignments = false;
                    let (ref_score, alt_score, max_ref_alt_exact_score) = if debug_score_alignments
                    {
                        let ref_score = get_read_flank_align_score_option(
                            scorer,
                            ref_flank_seqs.breakend_flanks[breakend_index]
                                [standard_registration_index],
                            "ref",
                            breakend_index,
                            breakend_flank.start_truncation,
                            breakend_flank.end_truncation,
                            score_settings.min_read_support_flank_size,
                            debug_score_alignments,
                        );
                        let alt_score = get_read_flank_align_score_option(
                            scorer,
                            contig_flank_seqs.breakend_flanks[breakend_index][registration_index],
                            "alt",
                            breakend_index,
                            breakend_flank.start_truncation,
                            breakend_flank.end_truncation,
                            score_settings.min_read_support_flank_size,
                            debug_score_alignments,
                        );
                        (ref_score, alt_score, None)
                    } else {
                        let truncated_flank_seq = |seq| {
                            get_truncated_flank_seq(
                                seq,
                                breakend_flank.start_truncation,
                                breakend_flank.end_truncation,
                                score_settings.min_read_support_flank_size,
                            )
                        };
                        get_ref_alt_flank_scores(
                            scorer,
                            truncated_flank_seq(
                                ref_flank_seqs.breakend_flanks[breakend_index]
                                    [standard_registration_index],
                            ),
                            truncated_flank_seq(
                                contig_flank_seqs.breakend_flanks[breakend_index]
                                    [registration_index],
                            ),
                        )
                    };
                    read_scores.ref_allele.breakend_flanks[breakend_index].registrations
                        [registration_index] = ref_score;
                    read_scores.alt_allele.breakend_flanks[breakend_index].registrations
                        [registration_index] = alt_score;

                    // Expand overlap score structure to match the size of overlap flanking sequences
                    //
                    for haplotype_index in (read_scores.overlapping_alleles.len()
                        ..overlapping_haplotypes_flank_seqs.len())
                        .map(|x| overlapping_haplotypes_flank_seqs[x].haplotype_index)
                    {
                        read_scores
                            .overlapping_alleles
                            .push(OverlappingHaplotypeSupportScores {
                                haplotype_index,
                                ..Default::default()
                            });
                    }

                    for (overlap_index, overlapping_haplotype_flank_seqs) in
                        overlapping_haplotypes_flank_seqs.iter().enumerate()
                    {
                        let overlap_flank_seq = overlapping_haplotype_flank_seqs
                            .flank_seqs
                            .breakend_flanks[breakend_index][registration_index];
                        let overlap_hap_score = if are_overlap_scores_averaged
                            || debug_score_alignments
                        {
                            get_read_flank_align_score_option(
                                scorer,
                                overlap_flank_seq,
                                "overlap_hap",
                                breakend_index,
                                breakend_flank.start_truncation,
                                breakend_flank.end_truncation,
                                score_settings.min_read_support_flank_size,
                                debug_score_alignments,
                            )
                        } else {
                            get_truncated_flank_seq(
                                overlap_flank_seq,
                                breakend_flank.start_truncation,
                                breakend_flank.end_truncation,
                                score_settings.min_read_support_flank_size,
                            )
                            .map(|x| get_overlap_flank_score(scorer, x, max_ref_alt_exact_score))
                        };
                        read_scores.overlapping_alleles[overlap_index]
                            .scores
                            .breakend_flanks[breakend_index]
                            .registrations[registration_index] = overlap_hap_score;
                    }
                }
            }
        },
    );
}

#[derive(Clone, Debug)]
struct FlankAlleleScore {
    score: Option<f64>,
    allele_type: AlleleType,
}

#[derive(Clone, Copy, Debug, strum::EnumCount)]
enum BreakendFlankType {
    Breakend1,
    Breakend2,
}

#[derive(Clone, Copy, Debug, strum::EnumCount)]
enum BreakendFlankRegistrationType {
    Standard,
    AsmRegionEdge,
}

#[derive(Clone, Debug)]
struct BestAlleleInfo {
    best_allele_score: FlankAlleleScore,

    /// Option is defined as the other allele type if the best score is a tie:
    tying_score_allele_type: Option<AlleleType>,

    /// Flank type is stored here for debugging
    #[allow(dead_code)]
    flank_type: BreakendFlankType,

    /// Flank registration type is stored here for debugging
    #[allow(dead_code)]
    flank_registration_type: BreakendFlankRegistrationType,
}

/// Sort best allele conditions so as to workaround common reference bias patterns.
///
/// Pick the best allele condition with the highest top score among all non-reference best allele
/// conditions.
///
/// Only pick the reference allele when all test conditions favor the reference, in which case the
/// top scoring reference allele condition is picked (for consistency even though this won't change
/// the outcome).
///
/// Sort is oriented from best to worst choice
///
fn sort_best_allele_score_conditions(a: &FlankAlleleScore, b: &FlankAlleleScore) -> Ordering {
    if a.allele_type != AlleleType::Ref && b.allele_type == AlleleType::Ref {
        Ordering::Less
    } else if b.allele_type != AlleleType::Ref && a.allele_type == AlleleType::Ref {
        Ordering::Greater
    } else {
        b.score.partial_cmp(&a.score).unwrap()
    }
}

fn sort_best_allele_conditions(a: &BestAlleleInfo, b: &BestAlleleInfo) -> Ordering {
    // currently the sorting only depends on the flank allele score field
    sort_best_allele_score_conditions(&a.best_allele_score, &b.best_allele_score)
}

/// Given the allele scores from a single flank test condition, determine which allele is best
/// under this condition.
///
fn get_best_allele_for_flank_condition(
    read_scores: &ReadSupportScores,
    overlapping_allele_index: Option<usize>,
    is_alt_on_top_haplotype_rank: bool,
    flank_type: BreakendFlankType,
    flank_registration_type: BreakendFlankRegistrationType,
    read_to_hap_allele_assignment: Option<AlleleType>,
    debug: bool,
) -> BestAlleleInfo {
    // Flank alignments must score at least this high to count for ANY allele
    let min_flank_score = 0.0;

    let mut allele_scores = Vec::new();
    allele_scores.push(FlankAlleleScore {
        score: read_scores.ref_allele.breakend_flanks[flank_type as usize].registrations
            [flank_registration_type as usize],
        allele_type: AlleleType::Ref,
    });
    allele_scores.push(FlankAlleleScore {
        score: read_scores.alt_allele.breakend_flanks[flank_type as usize].registrations
            [flank_registration_type as usize],
        allele_type: AlleleType::Alt,
    });
    if let Some(overlapping_allele_index) = overlapping_allele_index {
        allele_scores.push(FlankAlleleScore {
            score: read_scores.overlapping_alleles[overlapping_allele_index]
                .scores
                .breakend_flanks[flank_type as usize]
                .registrations[flank_registration_type as usize],
            allele_type: AlleleType::Overlap,
        });
    };

    if debug {
        eprintln!(
            "flank_type {flank_type:?} condition {flank_registration_type:?} allele_scores {allele_scores:?}"
        );
    }

    allele_scores.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    let mut tying_score_allele_type = if allele_scores[0].score == allele_scores[1].score {
        Some(allele_scores[1].allele_type)
    } else {
        None
    };

    // Special handling for the case of tying non-ref allele type scores
    let mut best_allele_index = 0;
    if tying_score_allele_type.is_some()
        && allele_scores[0].allele_type != AlleleType::Ref
        && allele_scores[1].allele_type != AlleleType::Ref
    {
        if let Some(hap_allele_type) = read_to_hap_allele_assignment {
            // First try to resolve the case with read_to_hap:
            if allele_scores[1].allele_type == hap_allele_type {
                best_allele_index = 1;
            }
        } else if !is_alt_on_top_haplotype_rank
            && allele_scores[0].allele_type == AlleleType::Alt
            && allele_scores[1].allele_type == AlleleType::Overlap
        {
            // Next resolve the case based on top/bottom haplotypes
            //
            // Define the conditions under which an ALT vs OVERLAP allele score tie is broken in favor of
            // OVERLAP:
            best_allele_index = 1;
        }
    }

    if best_allele_index == 1 {
        tying_score_allele_type = Some(allele_scores[0].allele_type);
    };

    // Get best allele score but filter for min_flank_score at this point.
    let best_allele_score = {
        let mut x = allele_scores[best_allele_index].clone();
        let filter_score = if let Some(score) = x.score {
            score < min_flank_score
        } else {
            false
        };
        if filter_score {
            x.score = None;
        }
        x
    };

    BestAlleleInfo {
        best_allele_score,
        tying_score_allele_type,
        flank_type,
        flank_registration_type,
    }
}

/// Check whether this read count is usable in the enumeration of allele depth for the
/// SV allele in question.
///
/// This is mainly a test of whether the allele assignment would conflict with a count
/// already assigned to another haplotype.
///
/// Note this has side effects on `read_to_hap``
///
fn check_if_read_allele_depth_count_is_usable(
    is_assigned_allele_score_tie: bool,
    on_multiple_sample_haplotypes: bool,
    assigned_allele: AlleleType,
    qname: &str,
    sv_haplotype_index: usize,
    overlapping_haplotype_index: Option<usize>,
    read_to_hap: &mut HashMap<String, usize>,
) -> bool {
    // The read_to_hap object is only used to record the haplotype of reads assigned to the
    // alt allele, and then requires any subsequent alt allele assignments to be to the same
    // haplotype.
    //
    // Haplotype assignment needs to be skipped when we can't clarify which haplotype the
    // read is on (such as reads supporting SV1 in the diagram below), because there can be
    // other SVs (SV2 in the diagram below) on one of the haplotypes, which are not be
    // shared on both assembly haplotypes.
    //
    //           -SV1-           -SV2-
    // HAP1: ABC-------ABCABCABCABCABCABC
    // HAP2: ABC-------ABCABCABC------ABC
    //
    let dont_assign_read_haplotype = is_assigned_allele_score_tie || on_multiple_sample_haplotypes;

    if dont_assign_read_haplotype {
        true
    } else if assigned_allele == AlleleType::Alt {
        let read_hap = read_to_hap
            .entry(qname.to_string())
            .or_insert(sv_haplotype_index);
        *read_hap == sv_haplotype_index
    } else if assigned_allele == AlleleType::Overlap {
        if let Some(overlapping_haplotype_index) = overlapping_haplotype_index {
            let read_hap = read_to_hap
                .entry(qname.to_string())
                .or_insert(overlapping_haplotype_index);
            *read_hap == overlapping_haplotype_index
        } else {
            false
        }
    } else {
        true
    }
}

/// Translate read support scores from locus into supporting read counts for a single allele in one
/// sample
///
/// Among other uses, the read count found here will be the VCF FORMAT/AD (allele depth) output for
/// each allele.
///
/// In many cases like large insertions, the ALT reads will perfectly support the REF allele on one
/// of the two breakends, so the policy for now is that the read supports ALT if the alt allele is
/// supported from either of the two tested flanks, but we should have a more model-based method to
/// translate all flank scores to final counts in future.
///
#[allow(clippy::too_many_arguments)]
fn set_allele_depth_counts_for_single_sample_and_sv_allele(
    score_settings: &ScoreSVSettings,
    locus_scores: &LocusSupportScores,
    debug: bool,
    sv_group: &mut SVGroup,
    sample_sv_info: &SampleSVInfo,
    sample_index: usize,
    read_to_hap: &mut HashMap<String, usize>,
) {
    let sv_index = sample_sv_info.sv_index;
    if debug {
        let target_sv = &sv_group.refined_svs[sv_index];
        let id = &target_sv.id;
        eprintln!("Getting AD counts for refined SV bp: {:?}", target_sv.bp);
        eprintln!(
            "SV ID: {}:{}:{}:{}",
            id.sample_index, id.cluster_index, id.assembly_index, id.alignment_index
        );
        eprintln!("Single region: {}", sv_group.is_single_region());
        eprintln!("Starting read_to_hap size: {}", read_to_hap.len());
        for (qname, hap) in read_to_hap.iter() {
            eprintln!("read_to_hap: hap: {hap} qname: {qname}");
        }
    }

    // The haplotype assembly index from which the target SV was proposed:
    let sv_haplotype_index = sv_group.sv_haplotype_map[sv_index];
    let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];

    let is_alt_on_top_haplotype_rank = if !sample_haplotype_list.is_empty() {
        sv_haplotype_index == sample_haplotype_list[0]
    } else {
        false
    };

    let overlapping_haplotype_index = {
        // Get full overlap list:
        let overlapping_haplotype_indexes =
            get_overlapping_haplotype_indexes(sv_group, sample_index, sample_sv_info);

        // Convert this into the relevant overlapping haplotype index for this case:
        if overlapping_haplotype_indexes.len() == 1 {
            overlapping_haplotype_indexes.first().cloned()
        } else {
            // TODO add more docs here about the reasoning for handling multiple guest overlaps this way?
            None
        }
    };

    let sample_score = &mut sv_group.refined_svs[sv_index].score.samples[sample_index];
    sample_score.clear_allele_depth();

    // Translate read scores into counts
    for (qname, read_scores) in locus_scores.iter() {
        if debug {
            eprintln!("Translating read scores to AD counts for: {qname}");
        }

        let overlapping_allele_index = {
            if let Some(overlapping_haplotype_index) = overlapping_haplotype_index {
                read_scores
                    .overlapping_alleles
                    .iter()
                    .enumerate()
                    .find(|(_, x)| x.haplotype_index == overlapping_haplotype_index)
                    .map(|(i, _)| i)
            } else {
                None
            }
        };

        // If there's an overlapping allele, we need to check whether this read has already been assinged as supporting the other allele
        let read_to_hap_allele_assignment = if overlapping_allele_index.is_some() {
            match read_to_hap.get(qname) {
                Some(hap) => {
                    if *hap == sv_haplotype_index {
                        Some(AlleleType::Alt)
                    } else {
                        // The qname has already been assigned to 'some other haplotype', so we mark
                        // it as an overlap count even thought we don't check that it is the same
                        // index as any current overlap haplotype... they all go into the same overlap
                        // count pool.
                        Some(AlleleType::Overlap)
                    }
                }
                None => None,
            }
        } else {
            None
        };

        // Gather details of the best allele found under different test conditions.
        //
        // Test conditions are different regions of alignment like breakend1 vs breakend2, or
        // different ways of attempting to register reads to the breakends.
        //
        let mut best_allele_tests = {
            let mut x = Vec::new();
            for flank_type in [BreakendFlankType::Breakend1, BreakendFlankType::Breakend2] {
                for flank_registration_type in [
                    BreakendFlankRegistrationType::Standard,
                    BreakendFlankRegistrationType::AsmRegionEdge,
                ] {
                    let best_allele = get_best_allele_for_flank_condition(
                        read_scores,
                        overlapping_allele_index,
                        is_alt_on_top_haplotype_rank,
                        flank_type,
                        flank_registration_type,
                        read_to_hap_allele_assignment,
                        debug,
                    );
                    x.push(best_allele);
                }
            }
            x
        };

        if debug {
            eprintln!("Best alleles under all test conditions:");
            for (best_allele_index, best_allele) in best_allele_tests.iter().enumerate() {
                eprintln!("best_allele{best_allele_index} {best_allele:?}");
            }
        }

        // Sort best allele findings under various test conditions so as to workaround common
        // reference bias patterns.
        //
        // Pick the best allele condition with the highest top score amont all non-reference best
        // allele choices.
        //
        // Only pick the reference allele when all test conditions favor the reference.
        //
        best_allele_tests.sort_by(sort_best_allele_conditions);

        if debug {
            eprintln!("Best alleles under all test conditions after sort:");
            for (best_allele_index, best_allele) in best_allele_tests.iter().enumerate() {
                eprintln!("best_allele{best_allele_index} {best_allele:?}");
            }
        }

        // After the above sort, the best_allele should be found at sorted condition index zero in a typical case
        //
        // The additional logic below handles score tie-breaking logic that will sometimes lead to picking an
        // index greater than zero.
        //
        let mut best_allele_condition_index = 0;
        if let Some(tying_allele_type) = best_allele_tests[0].tying_score_allele_type {
            // Special processing when the best condition includes a scoring tie.
            //
            // In this case we check the next best scoring conditions, to see if one of the two
            // tying alleles under the best condition becomes the untied winner. If this changes
            // the allele type, then change this to the assigned allele.
            //
            for (allele_index, allele_tests) in best_allele_tests.iter().enumerate().skip(1) {
                if allele_tests.tying_score_allele_type.is_none() {
                    if allele_tests.best_allele_score.score.is_some()
                        && allele_tests.best_allele_score.allele_type == tying_allele_type
                    {
                        best_allele_condition_index = allele_index;
                    }
                    break;
                }
            }
        }

        // Now that we have the best allele index, translate this into allele type and score status.
        //
        // From this point it is refered to as the 'assigned allele'
        //
        let (assigned_allele, is_assigned_allele_score_tie) = {
            let assigned = &best_allele_tests[best_allele_condition_index];
            (
                assigned
                    .best_allele_score
                    .score
                    .map(|_| assigned.best_allele_score.allele_type),
                assigned.tying_score_allele_type.is_some(),
            )
        };

        if let Some(assigned_allele) = assigned_allele {
            // The 'count allele' type represents what type will be used to represent this count in the output (most often the VCF AD field).
            // We may count an allele as something different than its actual assignment in some cases.
            //
            // Represent overlapping alleles as reference, unless the overlapping haplotype is
            // known to have a matching alt allele.
            //
            let count_allele = if assigned_allele == AlleleType::Overlap {
                if sample_score.on_multiple_sample_haplotypes {
                    AlleleType::Alt
                } else {
                    AlleleType::Overlap
                }
            } else {
                assigned_allele
            };

            let use_count = check_if_read_allele_depth_count_is_usable(
                is_assigned_allele_score_tie,
                sample_score.on_multiple_sample_haplotypes,
                assigned_allele,
                qname,
                sv_haplotype_index,
                overlapping_haplotype_index,
                read_to_hap,
            );

            if debug {
                eprintln!(
                    "AD Outcome: use_count: {use_count} count_allele: {count_allele} assigned_allele {assigned_allele} is_assigned_allele_score_tie {is_assigned_allele_score_tie}"
                );
            }

            if use_count {
                let counts = &mut sample_score.allele_depth[count_allele as usize];
                counts.any_strand += 1;

                /*
                if read_scores.is_reverse {
                    counts.rev_strand += 1;
                } else {
                    counts.fwd_strand += 1;
                }
                 */

                if score_settings.report_supporting_read_names && (count_allele == AlleleType::Alt)
                {
                    sample_score.supporting_read_names.push(qname.clone());
                }
            }
        } else {
            //
            if debug {
                eprintln!("AD Outcome: No assigned allele.");
            }
        }
    }
}

#[derive(Default, Clone)]
struct BreakendFlankScores {
    registrations: [Option<f64>; BreakendFlankRegistrationType::COUNT],
}

impl BreakendFlankScores {
    fn is_any_condition_scored(&self) -> bool {
        self.registrations[BreakendFlankRegistrationType::Standard as usize].is_some()
            || self.registrations[BreakendFlankRegistrationType::AsmRegionEdge as usize].is_some()
    }
}

/// Read support for a specific allele
#[derive(Default, Clone)]
struct ReadAlleleSupportScores {
    /// Array over each of the two breakend flanks
    breakend_flanks: [BreakendFlankScores; BreakendFlankType::COUNT],
}

#[derive(Default, Clone)]
struct OverlappingHaplotypeSupportScores {
    haplotype_index: usize,
    scores: ReadAlleleSupportScores,
}

/// Read support for all SV alleles at one locus
#[derive(Clone, Default)]
struct ReadSupportScores {
    ref_allele: ReadAlleleSupportScores,
    alt_allele: ReadAlleleSupportScores,
    overlapping_alleles: Vec<OverlappingHaplotypeSupportScores>,
}

/// Alleles represented by the quality score model only
///
#[derive(strum::EnumCount)]
pub enum QualityModelAlleles {
    Ref,
    Alt,
}

const ALLELE_GENOTYPE_PROBS: [[f64; QualityModelAlleles::COUNT]; Genotype::COUNT] =
    [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]];

/// Generate all genotype quality values for one SV allele in one sample
///
/// All quality scores are approximated from the allele counts rather than direct read alignments.
///
/// Returns the posterior prob for the reference allele, to facilitate QUAL estimation for all
/// samples
///
fn get_sv_sample_qualities_from_allele_counts(
    settings: &ScoreSVSettings,
    expected_cn_info: SVLocusExpectedCNInfo,
    sample_score: &mut SVSampleScoreInfo,
) -> f64 {
    sample_score.expected_cn_info = expected_cn_info;

    let is_haploid = sample_score
        .expected_cn_info
        .ploidy(settings.treat_single_copy_as_haploid)
        == SVLocusPloidy::Haploid;

    // Hardcode allele likelihood for each read matched to an allele
    let wrong_read_allele_lhood = 1e-5;

    // Get genotype likelihoods P(D|G)
    let mut gt_ln_lhoods = [0.0; Genotype::COUNT];
    for genotype_index in 0..Genotype::COUNT {
        if is_haploid && genotype_index == Genotype::Het as usize {
            gt_ln_lhoods[genotype_index] = f64::NEG_INFINITY;
            continue;
        }

        let allele_probs = ALLELE_GENOTYPE_PROBS[genotype_index];
        for (read_allele_index, read_allele_counts) in sample_score.allele_depth.iter().enumerate()
        {
            let mut read_allele_lhood = [1.0; 2];
            let wrong_read_allele_index = if read_allele_index == 1 {
                QualityModelAlleles::Ref
            } else {
                QualityModelAlleles::Alt
            } as usize;
            read_allele_lhood[wrong_read_allele_index] = wrong_read_allele_lhood;
            let read_lhood =
                read_allele_lhood[0] * allele_probs[0] + read_allele_lhood[1] * allele_probs[1];

            gt_ln_lhoods[genotype_index] += read_lhood.ln() * read_allele_counts.any_strand as f64;
        }
    }

    let mut gt_pprob = {
        let prior = if is_haploid {
            &settings.derived_score_parameters.haploid_sv_ln_prior
        } else {
            &settings.derived_score_parameters.diploid_sv_ln_prior
        };
        (0..Genotype::COUNT)
            .map(|x| gt_ln_lhoods[x] + prior[x])
            .collect::<Vec<_>>()
    };

    // iterator reversal here allows this to tie-break to '0/0' for zero AD counts
    let max_gt_lhood_index = gt_ln_lhoods
        .iter()
        .enumerate()
        .rev()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .unwrap();

    let _max_gt_pprob_index = normalize_ln_distro(&mut gt_pprob).unwrap();

    // Get phred ln_lhood value to use for PL outputs:
    for genotype_index in 0..Genotype::COUNT {
        sample_score.shared.gt_lhood_qscore[genotype_index] = std::cmp::min(
            ln_error_prob_to_qphred(
                gt_ln_lhoods[genotype_index] - gt_ln_lhoods[max_gt_lhood_index],
            ),
            settings.max_qscore as i32,
        );
    }

    {
        let total_allele_depth = sample_score
            .allele_depth
            .iter()
            .map(|x| x.any_strand)
            .sum::<usize>();

        // If the allele count is zero, GT remains None and will be output as './.' (or '.') in VCF
        if total_allele_depth != 0 {
            sample_score.shared.gt = Genotype::from_repr(max_gt_lhood_index);
        }
    }

    // Extract the second-highest genotype likelihood as the GQ value
    let gt_qscore = sample_score
        .shared
        .gt_lhood_qscore
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != max_gt_lhood_index)
        .min_by(|(_, a), (_, b)| a.cmp(b))
        .map(|(_, a)| a)
        .unwrap();

    sample_score.shared.gt_qscore = *gt_qscore;

    gt_pprob[Genotype::Ref as usize]
}

/// Generate diploid quality scores.
///
/// Use a simplified quality model built directly from the read support counts.
///
fn get_sv_qualities_from_allele_counts(
    settings: &ScoreSVSettings,
    sv_index: usize,
    sample_expected_cn_info: &[SVLocusExpectedCNInfo],
    sv_score: &mut SVScoreInfo,
    debug_settings: &ScoreDebugSettings,
) {
    let mut joint_ref_prob = 1.0;
    for (sample_index, sample_score) in sv_score.samples.iter_mut().enumerate() {
        let _debug = debug_settings.focused(sample_index, sv_index);
        let sample_ref_prob = get_sv_sample_qualities_from_allele_counts(
            settings,
            sample_expected_cn_info[sample_index],
            sample_score,
        );
        joint_ref_prob *= sample_ref_prob;
    }

    sv_score.sv_alt_score = Some(
        (error_prob_to_phred(joint_ref_prob) as f32)
            .min(settings.max_qscore as f32)
            .max(0.0),
    );
}

/// Get the overlapping haplotype indexes for this sample/SV combination
///
/// Eligible overlapping haplotypes are:
/// - For a sample with 2 native haplotypes, where one native is being processed as alt:
///   - Add one overlapping haplotype corresponding to case not being processed as the alt.
/// - For a sample with 1 native haplotype, where the native is being processed as alt:
///   - Add all other 'guest' haplotypes in the sv_groups as overlapping.
/// - For a sample with 1 native haplotype, where a guest is being processed as alt:
///   - Add the native haplotpye as overlapping.
///
fn get_overlapping_haplotype_indexes(
    sv_group: &SVGroup,
    sample_index: usize,
    sample_sv_info: &SampleSVInfo,
) -> Vec<usize> {
    let mut overlapping_haplotype_indexes = Vec::new();

    let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];
    let sv_haplotype_index = sv_group.sv_haplotype_map[sample_sv_info.sv_index];
    match sample_sv_info.relation {
        SampleSVRelation::Native => {
            if sample_haplotype_list.len() == 2 {
                overlapping_haplotype_indexes.extend(
                    sample_haplotype_list
                        .iter()
                        .cloned()
                        .filter(|&x| x != sv_haplotype_index),
                );
            } else if sample_haplotype_list.len() == 1 {
                overlapping_haplotype_indexes.extend(
                    (0..sv_group.group_haplotypes.len()).filter(|&x| x != sv_haplotype_index),
                );
            }
        }
        SampleSVRelation::Guest if sample_haplotype_list.len() == 1 => {
            overlapping_haplotype_indexes.push(sample_haplotype_list[0]);
        }
        _ => (),
    };
    overlapping_haplotype_indexes
}

/// Capture support information between different reads and SV alleles at a given locus for one
/// sample
type LocusSupportScores = BTreeMap<String, ReadSupportScores>;

/// Expand SV group regions into the larger search segments used to find candidate supporting reads
pub fn get_search_segments(
    chrom_list: &ChromList,
    group_regions: &[GenomeSegment],
) -> Vec<GenomeSegment> {
    let search_segment_expansion_size = 50;
    let mut search_segments = group_regions.to_vec();
    for search_segment in search_segments.iter_mut() {
        search_segment.expand_by(chrom_list, search_segment_expansion_size);
    }
    search_segments
}

/// Get read support scores for a single sample and a single SV allele, but considered in the context of
/// other overlapping SV haplotypes in the same sv_group
///
/// High-level scoring idea is to compare alignment of subsegments from each read to the ref and alt
/// allele(s), and attribute the read to the best alignment (or reject read if all alignments are
/// poor).
///
/// Most of the complexity here is from trying to define the segments to align. We try to define
/// the following 2 distinguishing regions:
/// - The right side of breakend1
/// - The left side of breakend2
///
/// For the first case, the region of interest is after breakend1+homology. We need to find the same
/// region in the reference, the contig and the read.
///
#[allow(clippy::too_many_arguments)]
fn get_locus_support_scores(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    genome_max_sv_depth_regions: &[RegionMap],
    bam_readers: &mut [&mut bam::IndexedReader],
    read_caches: &mut [SampleReadCache],
    sv_group: &mut SVGroup,
    sample_sv_info: &SampleSVInfo,
    sample_index: usize,
    scorer: &mut FlankAlignmentScorer,
    debug: bool,
) -> Option<LocusSupportScores> {
    let sv_index = sample_sv_info.sv_index;
    let target_segments = &sv_group.group_regions;
    let target_sv = &sv_group.refined_svs[sv_index];

    if debug {
        eprintln!(
            "Getting locus support scores for sample {sample_index} refined SV bp: {:?}",
            target_sv.bp
        );
        eprintln!(
            "SV ID: {}:{}:{}:{}",
            target_sv.id.sample_index,
            target_sv.id.cluster_index,
            target_sv.id.assembly_index,
            target_sv.id.alignment_index
        );
        eprintln!(
            "Single region refine: {}",
            target_sv.single_region_refinement
        );
    }

    let search_segments = get_search_segments(chrom_list, target_segments);

    // Check whether AD counting for this SV should be skipped due to high depth:
    for search_segment in search_segments.iter() {
        if is_segment_over_max_depth(search_segment, genome_max_sv_depth_regions) {
            let target_sv = &mut sv_group.refined_svs[sv_index];
            target_sv.score.samples[sample_index].is_max_scoring_depth = true;
            return None;
        }
    }

    // Drop mutable reference after this point:
    let sv_group: &SVGroup = sv_group;

    // Overlapping haplotype scores need to be exact if they're averaged
    let are_overlap_scores_averaged = sample_sv_info.relation == SampleSVRelation::Native
        && is_overlapping_guest_haplotype_selected(
            score_settings.treat_single_copy_as_haploid,
            sv_group,
            sample_index,
        );

    let ref_flanks = get_ref_flanks(score_settings, &target_sv.bp, reference, chrom_list, debug);
    let ref_flank_seqs = get_ref_flank_seqs(&ref_flanks);

    let mut scores = LocusSupportScores::new();

    let segment_count = target_segments.len();
    for segment_index in 0..segment_count {
        let search_segment = &search_segments[segment_index];
        let target_segment = &target_segments[segment_index];
        if debug {
            eprintln!("Starting read scan for search_segment: {search_segment:?}")
        }

        // Get the target SV contig breakpoint flanks
        let sv_contig_info = &sv_group.group_haplotypes[sv_group.sv_haplotype_map[sv_index]];
        let contig_flank_seqs = get_contig_flank_seqs(
            score_settings,
            &target_sv.bp,
            target_segment,
            sv_contig_info,
            debug,
        );

        // Get flanking sequences for all eligible haplotypes overlapping the alt haplotype
        //
        let overlapping_haplotypes_flank_seqs = {
            let mut x = Vec::new();
            for haplotype_index in
                get_overlapping_haplotype_indexes(sv_group, sample_index, sample_sv_info)
            {
                let flank_seqs = get_contig_flank_seqs(
                    score_settings,
                    &target_sv.bp,
                    target_segment,
                    &sv_group.group_haplotypes[haplotype_index],
                    debug,
                );
                x.push(OverlappingHaplotypeFlankSeqs {
                    haplotype_index,
                    flank_seqs,
                });
            }
            x
        };

        get_read_support_scores_from_segment(
            score_settings,
            reference,
            chrom_list,
            bam_readers[sample_index],
            &mut read_caches[sample_index],
            target_sv,
            &ref_flank_seqs,
            &contig_flank_seqs,
            &overlapping_haplotypes_flank_seqs,
            target_segment,
            search_segment,
            are_overlap_scores_averaged,
            scorer,
            &mut scores,
            debug,
        );
    }
    Some(scores)
}

fn is_segment_over_max_depth(
    segment: &GenomeSegment,
    genome_max_sv_depth_regions: &[RegionMap],
) -> bool {
    let chrom_max_sv_depth_regions = &genome_max_sv_depth_regions[segment.chrom_index];
    chrom_max_sv_depth_regions.intersect(segment.range.start, segment.range.end)
}

/// Identify and mark any matching SVs across multiple haplotypes in the same sample
///
/// The matching criteria used here are based on simple Breakpoint object equality
///
fn mark_within_sample_replicate_svs(sv_group: &mut SVGroup) {
    let sample_to_sv = sv_group.sample_to_sv_map();
    for (sample_index, sv_indexes) in sample_to_sv.iter().enumerate() {
        // Mark whether each SV is present on a haplotype associated with the sample
        for &sv_index in sv_indexes {
            let refined_sv = &mut sv_group.refined_svs[sv_index];
            let sample_score = &mut refined_sv.score.samples[sample_index];
            sample_score.on_sample_haplotype = true;
        }

        let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];
        if sample_haplotype_list.len() < 2 {
            // Filter out everything except multi-haplotype samples for the remaining analysis
            continue;
        }

        let mut bp_count = BTreeMap::new();
        let mut bp_found = BTreeSet::new();

        // Count up instances of SVs with a matching breakpoint in this sample:
        for &sv_index in sv_indexes {
            let bp = sv_group.refined_svs[sv_index].bp.clone();
            let bpe = bp_count.entry(bp).or_insert(0);
            *bpe += 1;
        }

        // Mark all SVs where the breakpoint count is greater than 1, add a separate mark to all
        // but the first so that the replicates can be filtered.
        for &sv_index in sv_indexes {
            let refined_sv = &mut sv_group.refined_svs[sv_index];
            let sample_score = &mut refined_sv.score.samples[sample_index];
            if let Some(&count) = bp_count.get(&refined_sv.bp)
                && count > 1
            {
                sample_score.on_multiple_sample_haplotypes = true;
                if !bp_found.insert(refined_sv.bp.clone()) {
                    sample_score.is_filtered_sv_replicate_within_sample = true;
                }
            }
        }
    }
}

/// Resolve duplicate SVs across all samples in the sv_group
///
fn resolve_cross_sample_replicate_svs(svs: &mut [RefinedSV]) {
    // Count matching breakpoint SV instances in this sv_group:
    let bp_count = {
        let mut x = BTreeMap::new();
        for bp in svs.iter().map(|x| x.bp.clone()) {
            *x.entry(bp).or_insert(0) += 1;
        }
        x
    };

    // Build a map specifying which SVs should be consolidated into others:
    let consolidate_from_to_sv_replicate = {
        let mut x = Vec::new();
        // Arbitrarily choose the first matching SV to be the representative of the matching set
        let mut first_bp_in_group = BTreeMap::new();
        for (sv_index, refined_sv) in svs.iter().enumerate() {
            if let Some(&count) = bp_count.get(&refined_sv.bp)
                && count > 1
            {
                match first_bp_in_group.get(&refined_sv.bp) {
                    Some(&first_sv_index) => {
                        x.push((sv_index, first_sv_index));
                    }
                    None => {
                        first_bp_in_group.insert(refined_sv.bp.clone(), sv_index);
                    }
                }
            }
        }
        x
    };

    // Use the from->to map to update sv_group
    for (from_index, to_index) in consolidate_from_to_sv_replicate {
        // Mark the from index as filtered:
        svs[from_index].score.is_filtered_sv_replicate = true;

        // Copy AD counts between the 'from' and 'to' sv as follows: If the 'from' SV is present in
        // the sample, and not duplicate filtered within the sample, then copy the counts.
        let sample_count = svs[from_index].score.samples.len();
        for sample_index in 0..sample_count {
            let move_counts = {
                let from_sample_score = &svs[from_index].score.samples[sample_index];
                from_sample_score.on_sample_haplotype
                    && !from_sample_score.is_filtered_sv_replicate_within_sample
            };

            if move_counts {
                svs[to_index].score.samples[sample_index].allele_depth =
                    svs[from_index].score.samples[sample_index]
                        .allele_depth
                        .clone();
            }
        }
    }
}

/// This method applies some additional post-hoc heuristics to fix haplotype inconsistencies
///
/// The longer term intent is that a more deeply haplotype-focused scoring system will eliminate the
/// need for this step.
///
fn refine_sample_overlapping_haplotype_scoring(
    score_settings: &ScoreSVSettings,
    sv_group: &mut SVGroup,
    haplotype_to_sv_map: &[Vec<usize>],
    sample_index: usize,
    debug_settings: &ScoreDebugSettings,
) {
    let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];

    // Skip any samples without overlapping haplotypes in this SV group
    if sample_haplotype_list.len() < 2 {
        return;
    }

    let top_haplotype_index = sample_haplotype_list[0];
    let overlap_haplotype_index = sample_haplotype_list[1];
    let overlap_sv_indexes = &haplotype_to_sv_map[overlap_haplotype_index];

    // Refine the scoring among the cluster SV group for consistency:
    // 1. If any 'top' haplotype SV is genotyped as 1/1 in a sample, then:
    //    - Put overlap AD counts back into alt for the 'top' haplotype SV in question
    //    - Zero out the support counts for the 'botton' haplotype in this sample
    //    - Rescore
    //    - Nominate the overlapping SVs on the bottom haplotype for filtration -- this will only
    //      occur if the nomination is consistent across all samples in which the SV occurs
    //
    // ...add more consistency operations as we encounter them

    // Filter for SV indexes on the top haplotype of this sample:
    let refined_sv_count = sv_group.refined_svs.len();
    for refined_sv_index in
        (0..refined_sv_count).filter(|&x| sv_group.sv_haplotype_map[x] == top_haplotype_index)
    {
        if sv_group.refined_svs[refined_sv_index].filter_sv() {
            continue;
        }

        let is_homozygous_gt = {
            let refined_sv = &sv_group.refined_svs[refined_sv_index];
            let sample_score = &refined_sv.score.samples[sample_index];
            sample_score.shared.gt == Some(Genotype::Hom)
        };

        if !is_homozygous_gt {
            continue;
        }

        // Fix overlap AD counts and qualities in this SV:
        {
            let refined_sv = &mut sv_group.refined_svs[refined_sv_index];
            let sample_score = &mut refined_sv.score.samples[sample_index];
            sample_score.allele_depth[AlleleType::Alt as usize].any_strand +=
                sample_score.allele_depth[AlleleType::Overlap as usize].any_strand;
            sample_score.allele_depth[AlleleType::Overlap as usize].any_strand = 0;

            get_sv_qualities_from_allele_counts(
                score_settings,
                refined_sv_index,
                &sv_group.sample_expected_cn_info,
                &mut refined_sv.score,
                debug_settings,
            );
        }

        // Search for any SVs on this sample's overlapping haplotype which intersect this SV, and
        // zero them out:
        //
        let base_sv_range = {
            let bp = &sv_group.refined_svs[refined_sv_index].bp;
            IntRange::from_pair(
                bp.breakend1.segment.range.start,
                bp.breakend2.as_ref().unwrap().segment.range.start + 1,
            )
        };

        for &overlap_sv_index in overlap_sv_indexes.iter() {
            let overlap_sv_range = {
                let bp = &sv_group.refined_svs[overlap_sv_index].bp;
                IntRange::from_pair(
                    bp.breakend1.segment.range.start,
                    bp.breakend2.as_ref().unwrap().segment.range.start + 1,
                )
            };
            if !base_sv_range.intersect_range(&overlap_sv_range) {
                continue;
            }

            // Zero-out the SV in this sample:
            let overlap_sv = &mut sv_group.refined_svs[overlap_sv_index];
            {
                let move_counts =
                    |adepth: &mut [AlleleCounts], from: AlleleType, to: AlleleType| {
                        adepth[to as usize].any_strand += adepth[from as usize].any_strand;
                        adepth[from as usize].any_strand = 0;
                    };

                let overlap_sample_score = &mut overlap_sv.score.samples[sample_index];
                let adepth = &mut overlap_sample_score.allele_depth;
                move_counts(adepth, AlleleType::Alt, AlleleType::Ref);
                move_counts(adepth, AlleleType::Overlap, AlleleType::Ref);
                overlap_sample_score.is_gt_excluded = true;
            }
            get_sv_qualities_from_allele_counts(
                score_settings,
                overlap_sv_index,
                &sv_group.sample_expected_cn_info,
                &mut overlap_sv.score,
                debug_settings,
            );
        }
    }
}

struct DerivedScoreParameters {
    haploid_sv_ln_prior: [f64; Genotype::COUNT],
    diploid_sv_ln_prior: [f64; Genotype::COUNT],
}

/// By convention the haploid prior uses the diploid genotype states but with zero prob for het
///
fn get_haploid_ln_prior(prior: f64) -> [f64; Genotype::COUNT] {
    use Genotype::*;
    let mut p = [0.0; Genotype::COUNT];
    p[Hom as usize] = prior;
    p[Ref as usize] = 1.0 - p[Hom as usize];
    for p in p.iter_mut() {
        *p = p.ln();
    }
    p
}

fn get_diploid_ln_prior(prior: f64) -> [f64; Genotype::COUNT] {
    use Genotype::*;
    let mut p = [0.0; Genotype::COUNT];
    p[Het as usize] = prior;
    p[Hom as usize] = prior / 2.0;
    p[Ref as usize] = 1.0 - (p[Het as usize] + p[Hom as usize]);
    for p in p.iter_mut() {
        *p = p.ln();
    }
    p
}

impl DerivedScoreParameters {
    fn new(sv_prior: f64) -> Self {
        Self {
            haploid_sv_ln_prior: get_haploid_ln_prior(sv_prior),
            diploid_sv_ln_prior: get_diploid_ln_prior(sv_prior),
        }
    }
}

pub struct ScoreSVSettings {
    /// Copied from cli settings
    min_sv_mapq: u32,

    /// Copied from cli settings
    min_gap_compressed_identity: f64,

    /// Length of sequence aligned between reads and candidate alleles during scoring
    ///
    read_support_flank_size: i64,

    /// Min read support length for scoring, must not be greater than read_support_flank_size
    min_read_support_flank_size: i64,

    /// Min allowed distance from the read edge to the anchoring position used for read support
    read_support_anchor_min_edge_distance: i64,

    max_qscore: u32,

    /// Additional scoring parameters precomputed from sv_prior
    derived_score_parameters: DerivedScoreParameters,

    /// parallels command-line debug option to enumerate all read-names supporting the SV allele for each VCF record output
    report_supporting_read_names: bool,

    /// If true, enable local phased genotype output
    enable_phasing: bool,

    treat_single_copy_as_haploid: bool,
}

impl ScoreSVSettings {
    pub fn new(
        min_sv_mapq: u32,
        min_gap_compressed_identity: f64,
        max_qscore: u32,
        report_supporting_read_names: bool,
        enable_phasing: bool,
        treat_single_copy_as_haploid: bool,
    ) -> Self {
        // Prior SV prob (theta)
        let sv_prior = 5e-4;
        let read_support_flank_size = 500;
        let min_read_support_flank_size = std::cmp::min(100, read_support_flank_size);

        Self {
            min_sv_mapq,
            min_gap_compressed_identity,
            read_support_flank_size,
            min_read_support_flank_size,
            read_support_anchor_min_edge_distance: 200,
            max_qscore,
            derived_score_parameters: DerivedScoreParameters::new(sv_prior),
            report_supporting_read_names,
            enable_phasing,
            treat_single_copy_as_haploid,
        }
    }
}

/// SV state terminology used to original SV/sample relationships during multi-sample analysis
///
#[derive(Clone, Copy, Debug, PartialEq)]
enum SampleSVRelation {
    /// A haplotype/SV already assigned to/discovered from at least one sample
    Native,

    /// A non-native haplotype/SV being trialed on another sample, when the sample has less than 2 native haplotypes
    Guest,

    /// A haplotype/SV is a non-native item from another sample, when the sample in question already has 2 native haplotypes.
    Blocked,
}

#[derive(Debug)]
struct SampleSVInfo {
    sv_index: usize,
    relation: SampleSVRelation,
}

/// Build an order for SV evaluation, which can be different for each sample
///
/// The sample-specific SV order is currently designed to ensure that native haplotypes are evaluated first
///
fn get_sample_sv_order(
    treat_single_copy_as_haploid: bool,
    sv_group: &SVGroup,
) -> Vec<Vec<SampleSVInfo>> {
    let sv_count = sv_group.refined_svs.len();
    let haplotype_to_sv_map = sv_group.haplotype_to_sv_map();
    let mut sample_sv_order = Vec::new();
    for (sample_index, native_haplotype_list) in sv_group.sample_haplotype_list.iter().enumerate() {
        let mut added_to_sv_order = vec![false; sv_count];
        let mut sv_order = Vec::new();
        {
            // Step 1: Add native SVs
            let relation = SampleSVRelation::Native;
            for &sv_index in native_haplotype_list
                .iter()
                .flat_map(|&hap_index| haplotype_to_sv_map[hap_index].iter())
            {
                if added_to_sv_order[sv_index] {
                    // TODO: Make this an assertion once sv_group is cleaned up to prevent it
                    continue;
                }
                sv_order.push(SampleSVInfo { sv_index, relation });
                added_to_sv_order[sv_index] = true;
            }
        }

        {
            // Step 2: Add all remaining SVs:
            let max_haplotype_count = match sv_group.sample_expected_cn_info[sample_index]
                .ploidy(treat_single_copy_as_haploid)
            {
                SVLocusPloidy::Diploid => 2,
                SVLocusPloidy::Haploid => 1,
            };
            let relation = if native_haplotype_list.len() >= max_haplotype_count {
                SampleSVRelation::Blocked
            } else {
                SampleSVRelation::Guest
            };
            for sv_index in (0..sv_count).filter(|&sv_index| !added_to_sv_order[sv_index]) {
                sv_order.push(SampleSVInfo { sv_index, relation });
            }
        }
        sample_sv_order.push(sv_order);
    }
    sample_sv_order
}

/// Get read alignment scores to each allele in each sample
///
/// # Arguments
/// * `sample_sv_order` - Specify the order in which SV alleles should be processed in each sample. This enables SV alleles 'native' to the sample to be processed first.
///
/// Returns:
/// - A vector over all samples of:
///   - A vector over all SVs of:
///     - All read alignment scores to the SV for that sample, or None
///
#[allow(clippy::too_many_arguments)]
fn get_sv_group_read_alignment_scores(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    all_sample_genome_max_sv_depth_regions: &[&[RegionMap]],
    bam_readers: &mut [&mut bam::IndexedReader],
    read_caches: &mut [SampleReadCache],
    sv_group: &mut SVGroup,
    sample_sv_order: &[Vec<SampleSVInfo>],
    debug_settings: &ScoreDebugSettings,
) -> Vec<Vec<Option<LocusSupportScores>>> {
    let mut all_sample_locus_scores = Vec::new();

    let sv_count = sv_group.refined_svs.len();
    let mut scorer = FlankAlignmentScorer::new();

    for (sample_index, genome_max_sv_depth_regions) in
        all_sample_genome_max_sv_depth_regions.iter().enumerate()
    {
        let mut all_sv_locus_scores_for_sample = vec![None; sv_count];
        for sample_sv_info in sample_sv_order[sample_index]
            .iter()
            .filter(|x| x.relation != SampleSVRelation::Blocked)
        {
            let sv_index = sample_sv_info.sv_index;
            if sv_group.refined_svs[sv_index].filter_sv() {
                continue;
            }
            let debug = debug_settings.focused(sample_index, sv_index);

            all_sv_locus_scores_for_sample[sv_index] = get_locus_support_scores(
                score_settings,
                reference,
                chrom_list,
                genome_max_sv_depth_regions,
                bam_readers,
                read_caches,
                sv_group,
                sample_sv_info,
                sample_index,
                &mut scorer,
                debug,
            );
        }
        all_sample_locus_scores.push(all_sv_locus_scores_for_sample);
    }
    all_sample_locus_scores
}

/// For samples with only one native assembled haplotype and a max haplotype count of two, check all
/// guest haplotype scores to determine which one (if any) will be assigned as the 'guest'
/// overlapping haplotype.
///
/// Still not clear the best way to do this so procedure is a bit of a first-pass hack.. avg
/// every score for the overlaps of one sv in the haplotype, and pick the highest avg
///
fn select_overlapping_guest_haplotype_index(
    treat_single_copy_as_haploid: bool,
    sv_group: &SVGroup,
    sample_sv_order: &[Vec<SampleSVInfo>],
    all_sv_locus_scores_for_sample: &[Option<LocusSupportScores>],
    sample_index: usize,
) -> Option<usize> {
    let sv_group_haplotype_count = sv_group.group_haplotypes.len();
    if is_overlapping_guest_haplotype_selected(treat_single_copy_as_haploid, sv_group, sample_index)
    {
        for sample_sv_info in sample_sv_order[sample_index]
            .iter()
            .filter(|x| x.relation == SampleSVRelation::Native)
        {
            let sv_index = sample_sv_info.sv_index;
            if sv_group.refined_svs[sv_index].filter_sv() {
                continue;
            }

            #[derive(Clone, Default)]
            struct AvgInfo {
                total: f64,
                count: i32,
            }

            // Not sure about a good way to determine best, so start out with a simple average over all reads, breakends and conditions
            let mut avg_per_hap = vec![AvgInfo::default(); sv_group_haplotype_count];
            for (_, read_scores) in all_sv_locus_scores_for_sample[sv_index].iter().flatten() {
                for overlapping_allele in read_scores.overlapping_alleles.iter() {
                    let avg_info = &mut avg_per_hap[overlapping_allele.haplotype_index];
                    for read_alignment_score in overlapping_allele
                        .scores
                        .breakend_flanks
                        .iter()
                        .flat_map(|x| x.registrations)
                        .flatten()
                    {
                        avg_info.total += read_alignment_score;
                        avg_info.count += 1;
                    }
                }
            }

            struct HapAvg {
                haplotype_index: usize,
                avg: i64,
                count: i32,
            }

            let mut avg_list = Vec::new();
            let native_haplotype_index = sv_group.sv_haplotype_map[sample_sv_info.sv_index];
            for haplotype_index in
                (0..sv_group_haplotype_count).filter(|&x| x != native_haplotype_index)
            {
                if avg_per_hap[haplotype_index].count > 1 {
                    let count = avg_per_hap[haplotype_index].count;
                    let favg = avg_per_hap[haplotype_index].total / count as f64;

                    // Switch avg to an integer to keep the selection deterministic
                    let avg = (favg * 1e9) as i64;
                    avg_list.push(HapAvg {
                        haplotype_index,
                        avg,
                        count,
                    });
                }
            }

            // Sort for the highest avg
            //
            // Highest avg should be unique most of the time, but add count
            // and hap index to keep it deterministic in case of a tie
            //
            avg_list.sort_by(|a, b| {
                b.avg
                    .cmp(&a.avg)
                    .then(b.count.cmp(&a.count))
                    .then(a.haplotype_index.cmp(&b.haplotype_index))
            });
            return avg_list.first().map(|x| x.haplotype_index);
        }
    }
    None
}

/// After scoring but before enumerating AD counts, we determine if additional, 'guest' haplotypes
/// can be assigned to any samples with less than max assembled haplotype count (ussually this is 2)
///
fn add_guest_haplotypes(
    treat_single_copy_as_haploid: bool,
    sv_group: &mut SVGroup,
    sample_sv_order: &[Vec<SampleSVInfo>],
    all_sample_locus_scores: &[Vec<Option<LocusSupportScores>>],
    debug_settings: &ScoreDebugSettings,
) {
    for (sample_index, all_sv_locus_scores_for_sample) in all_sample_locus_scores.iter().enumerate()
    {
        let debug = debug_settings.sample(sample_index);
        let overlapping_guest_haplotype_index = select_overlapping_guest_haplotype_index(
            treat_single_copy_as_haploid,
            sv_group,
            sample_sv_order,
            all_sv_locus_scores_for_sample,
            sample_index,
        );
        if let Some(haplotype_index) = overlapping_guest_haplotype_index {
            if debug {
                eprintln!(
                    "For sample {}: adding guest haplotype {} to sample_hap_list {:?}",
                    sample_index, haplotype_index, sv_group.sample_haplotype_list[sample_index]
                );
            }
            sv_group.sample_haplotype_list[sample_index].push(haplotype_index);
        }
    }
}

/// Convert read to SV haplotype alignment scores into supporting read counts, for all alleles in an sv_group, over all samples
///
fn set_allele_depth_counts_for_sv_group(
    score_settings: &ScoreSVSettings,
    sv_group: &mut SVGroup,
    sample_sv_order: &[Vec<SampleSVInfo>],
    all_sample_locus_scores: &[Vec<Option<LocusSupportScores>>],
    debug_settings: &ScoreDebugSettings,
) {
    for (sample_index, all_sv_locus_scores_for_sample) in all_sample_locus_scores.iter().enumerate()
    {
        // To keep the AD counts across the locus consistent, record once a qname has been counted
        // towards a certain haplotype, so that it cannot be counted towards another
        //
        // The key of the map is qname, the value of the map is the haplotype index from the sv
        // group. Nothing gets reserved for the reference haplotype in this structure.
        //
        let mut read_to_hap = HashMap::new();

        // For some overlapping haplotype cases, it's more accurate to build a complete read_to_hap
        // structure first, and then get the AD counts based off this haplotype assignment.
        //
        // To simplify things we just run the AD count function twice in these cases for now.
        //
        let ad_count_passes = {
            let sample_haplotype_list = &sv_group.sample_haplotype_list[sample_index];
            let is_multihaplotype_sample = sample_haplotype_list.len() > 1;
            if is_multihaplotype_sample { 2 } else { 1 }
        };
        let mut max_native_support_count = 0;
        for _ in 0..ad_count_passes {
            for sample_sv_info in sample_sv_order[sample_index]
                .iter()
                .filter(|x| x.relation != SampleSVRelation::Blocked)
            {
                let sv_index = sample_sv_info.sv_index;
                if sv_group.refined_svs[sv_index].filter_sv() {
                    continue;
                }

                if sv_group.refined_svs[sv_index].score.samples[sample_index]
                    .is_filtered_sv_replicate_within_sample
                {
                    continue;
                }

                if let Some(locus_scores) = &all_sv_locus_scores_for_sample[sv_index] {
                    let debug = debug_settings.focused(sample_index, sv_index);

                    set_allele_depth_counts_for_single_sample_and_sv_allele(
                        score_settings,
                        locus_scores,
                        debug,
                        sv_group,
                        sample_sv_info,
                        sample_index,
                        &mut read_to_hap,
                    );

                    if sample_sv_info.relation == SampleSVRelation::Native {
                        let sample_score =
                            &sv_group.refined_svs[sv_index].score.samples[sample_index];
                        max_native_support_count = std::cmp::max(
                            max_native_support_count,
                            sample_score.total_allele_depth(),
                        );
                    }
                }
            }
        }

        // Fill in REF counts for blocked alleles
        for sample_sv_info in sample_sv_order[sample_index]
            .iter()
            .filter(|x| x.relation == SampleSVRelation::Blocked)
        {
            let sv_index = sample_sv_info.sv_index;
            if sv_group.refined_svs[sv_index].filter_sv() {
                continue;
            }

            let sample_score = &mut sv_group.refined_svs[sv_index].score.samples[sample_index];
            sample_score.allele_depth[AlleleType::Ref as usize].any_strand =
                max_native_support_count;
        }
    }
}

pub struct ScoreDebugSettings {
    pub debug: bool,

    // If debug is true, this provides a way to narrow the scope of debug output to a specific
    // sample. Setting to None will print all samples.
    pub debug_sample_index: Option<usize>,

    // Same as above but for SVs in the sv_group
    pub debug_sv_index: Option<usize>,
}

impl ScoreDebugSettings {
    pub fn focused(&self, sample_index: usize, sv_index: usize) -> bool {
        self.debug
            && match self.debug_sample_index {
                Some(x) => x == sample_index,
                None => true,
            }
            && match self.debug_sv_index {
                Some(x) => x == sv_index,
                None => true,
            }
    }

    pub fn sample(&self, sample_index: usize) -> bool {
        self.debug
            && match self.debug_sample_index {
                Some(x) => x == sample_index,
                None => true,
            }
    }
}

/// Setup information for phased genotype output within individual SV Groups
fn set_local_phasing(sv_group: &mut SVGroup) {
    let sample_count = sv_group.sample_haplotype_list.len();

    // First determine, for each sample, if the sv_group contains more than one het:
    let phase_sample = {
        let mut x = Vec::new();
        for sample_index in 0..sample_count {
            let sample_het_count = sv_group
                .refined_svs
                .iter()
                .map(|x| x.score.samples[sample_index].shared.gt.as_ref())
                .filter(|x| *x == Some(&Genotype::Het))
                .count();
            x.push(sample_het_count > 1);
        }
        x
    };

    // Skip further steps if no samples are phased
    if !phase_sample.iter().any(|&x| x) {
        return;
    }

    // Determine what value to use for phase set, for now use the minimum sv_group pos
    let phase_set = sv_group
        .refined_svs
        .iter()
        .map(|x| x.bp.breakend1.segment.range.start + 1)
        .min()
        .unwrap() as i32;

    // Setup phased genotypes in samples with multiple hets
    for (sv_index, sv_score_info) in sv_group
        .refined_svs
        .iter_mut()
        .map(|x| &mut x.score)
        .enumerate()
    {
        // Setup the phase_set tag for this SV group
        sv_score_info.phase_set = Some(phase_set);

        // Get the haplotype index of this SV
        let sv_hap_index = sv_group.sv_haplotype_map[sv_index];
        for (sample_index, sample_score_info) in sv_score_info.samples.iter_mut().enumerate() {
            if !phase_sample[sample_index] {
                continue;
            }

            for (sample_hap_rank, hap_index) in sv_group.sample_haplotype_list[sample_index]
                .iter()
                .enumerate()
            {
                if sv_hap_index == *hap_index {
                    if sample_hap_rank == 0 {
                        sample_score_info.phase = SVPhaseStatus::Primary;
                    } else if sample_hap_rank == 1 {
                        sample_score_info.phase = SVPhaseStatus::Secondary;
                    } else {
                        panic!("Invalid sample_hap_rank {sample_hap_rank}");
                    }
                    break;
                }
            }
        }
    }
}

/// Score all refined_svs from a single group
///
/// Score each SV in the context of all other SVs in the group
///
#[allow(clippy::too_many_arguments)]
fn score_refined_sv_group(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    all_sample_genome_max_sv_depth_regions: &[&[RegionMap]],
    bam_readers: &mut [&mut bam::IndexedReader],
    read_caches: &mut [SampleReadCache],
    sv_group: &mut SVGroup,
    debug_settings: &ScoreDebugSettings,
) {
    let sample_count = sv_group.sample_haplotype_list.len();

    // Initialize a scoring structure for each sample:
    for refined_sv in sv_group.refined_svs.iter_mut() {
        refined_sv
            .score
            .samples
            .resize(sample_count, SVSampleScoreInfo::default());
    }

    let sample_sv_order =
        get_sample_sv_order(score_settings.treat_single_copy_as_haploid, sv_group);

    let all_sample_locus_scores = get_sv_group_read_alignment_scores(
        score_settings,
        reference,
        chrom_list,
        all_sample_genome_max_sv_depth_regions,
        bam_readers,
        read_caches,
        sv_group,
        &sample_sv_order,
        debug_settings,
    );

    if debug_settings.debug {
        eprintln!("sv_group before guest: {sv_group:?}");
    }

    add_guest_haplotypes(
        score_settings.treat_single_copy_as_haploid,
        sv_group,
        &sample_sv_order,
        &all_sample_locus_scores,
        debug_settings,
    );

    if debug_settings.debug {
        eprintln!("sv_group after guest: {sv_group:?}");
    }

    // Redo sample_sv_order after add_guest_haplotypes
    let sample_sv_order =
        get_sample_sv_order(score_settings.treat_single_copy_as_haploid, sv_group);

    // These marks are used by the AD counting routine
    mark_within_sample_replicate_svs(sv_group);

    set_allele_depth_counts_for_sv_group(
        score_settings,
        sv_group,
        &sample_sv_order,
        &all_sample_locus_scores,
        debug_settings,
    );

    resolve_cross_sample_replicate_svs(&mut sv_group.refined_svs);

    for (sv_index, refined_sv) in sv_group
        .refined_svs
        .iter_mut()
        .enumerate()
        .filter(|(_, x)| !x.filter_sv())
    {
        get_sv_qualities_from_allele_counts(
            score_settings,
            sv_index,
            &sv_group.sample_expected_cn_info,
            &mut refined_sv.score,
            debug_settings,
        );
    }

    let haplotype_to_sv_map = sv_group.haplotype_to_sv_map();
    for sample_index in 0..sample_count {
        refine_sample_overlapping_haplotype_scoring(
            score_settings,
            sv_group,
            &haplotype_to_sv_map,
            sample_index,
            debug_settings,
        );
    }

    // Setup phasing data
    //
    if score_settings.enable_phasing {
        set_local_phasing(sv_group);
    }
}

pub struct SampleScoreData<'a> {
    pub genome_sv_scoring_exclusion_regions: &'a [RegionMap],
}

/// Run breakpoint scoring on all refined_svs from a single group
///
/// Each SV candidate is evaluated in the context of all other SVs in the same group
///
#[allow(clippy::too_many_arguments)]
pub fn score_and_assess_refined_sv_group(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    all_sample_data: &[SampleScoreData],
    bam_readers: &mut [&mut bam::IndexedReader],
    read_caches: &mut [SampleReadCache],
    sv_group: &mut SVGroup,
    debug_settings: &ScoreDebugSettings,
) {
    sv_group.assert_validity(score_settings.treat_single_copy_as_haploid);

    let sample_count = all_sample_data.len();
    assert_ne!(sample_count, 0);
    assert_eq!(sample_count, bam_readers.len());
    assert_eq!(sample_count, read_caches.len());
    assert_eq!(sample_count, sv_group.sample_haplotype_list.len());
    assert!(!sv_group.refined_svs.is_empty());
    assert!(!sv_group.group_regions.is_empty());

    let all_sample_genome_max_sv_depth_regions = all_sample_data
        .iter()
        .map(|x| x.genome_sv_scoring_exclusion_regions)
        .collect::<Vec<_>>();

    score_refined_sv_group(
        score_settings,
        reference,
        chrom_list,
        &all_sample_genome_max_sv_depth_regions,
        bam_readers,
        read_caches,
        sv_group,
        debug_settings,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sort_best_allele_conditions() {
        // These values are ignored by the sort, so just setup placeholders for all cases:
        let tying_score_allele_type = None;
        let flank_type = BreakendFlankType::Breakend1;
        let flank_registration_type = BreakendFlankRegistrationType::Standard;

        let bnone = BestAlleleInfo {
            best_allele_score: FlankAlleleScore {
                score: None,
                allele_type: AlleleType::Ref,
            },
            tying_score_allele_type,
            flank_type,
            flank_registration_type,
        };

        let bsomeref = BestAlleleInfo {
            best_allele_score: FlankAlleleScore {
                score: Some(0.5),
                allele_type: AlleleType::Ref,
            },
            tying_score_allele_type,
            flank_type,
            flank_registration_type,
        };

        // Test that Some(x) is ranked before None:
        let mut v = [bnone.clone(), bsomeref.clone()];
        v.sort_by(sort_best_allele_conditions);
        assert!(v[0].best_allele_score.score.is_some());

        let bsomealt = BestAlleleInfo {
            best_allele_score: FlankAlleleScore {
                score: Some(0.3),
                allele_type: AlleleType::Alt,
            },
            tying_score_allele_type,
            flank_type,
            flank_registration_type,
        };

        // Test that Alt is ranked before Ref, even if it has a lower score
        let mut v = [bnone.clone(), bsomeref.clone(), bsomealt.clone()];
        v.sort_by(sort_best_allele_conditions);
        assert_eq!(v[0].best_allele_score.allele_type, AlleleType::Alt);

        let bsomealt2 = BestAlleleInfo {
            best_allele_score: FlankAlleleScore {
                score: Some(0.4),
                allele_type: AlleleType::Alt,
            },
            tying_score_allele_type,
            flank_type,
            flank_registration_type,
        };

        // Test that best Alt score is selected
        let mut v = [
            bnone.clone(),
            bsomeref.clone(),
            bsomealt.clone(),
            bsomealt2.clone(),
        ];
        v.sort_by(sort_best_allele_conditions);
        assert_eq!(v[0].best_allele_score.allele_type, AlleleType::Alt);
        assert!(v[0].best_allele_score.score.is_some());
        assert!(v[0].best_allele_score.score.unwrap() > 0.35);

        // Test that Alt is selected over higher scoring Ref (simpler version this time)
        let mut v = [bsomealt.clone(), bsomeref.clone()];
        v.sort_by(sort_best_allele_conditions);
        assert_eq!(v[0].best_allele_score.allele_type, AlleleType::Alt);
    }

    fn random_seq(rand: &mut impl FnMut(usize) -> usize, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rand(4)]).collect()
    }

    fn xorshift(seed: u64) -> impl FnMut(usize) -> usize {
        let mut state = seed;
        move |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        }
    }

    /// Flanks related to the read flank with a few edits, so that their scores are found exactly,
    /// or unrelated, so that they need the DP matrix
    fn random_allele_flank(rand: &mut impl FnMut(usize) -> usize, text: &[u8]) -> Vec<u8> {
        match rand(3) {
            0 => {
                let mut x = text.to_vec();
                for _ in 0..rand(6) {
                    let i = rand(x.len());
                    x[i] = b"ACGT"[rand(4)];
                }
                x
            }
            1 => text[rand(text.len() / 3)..].to_vec(),
            _ => {
                let len = 30 + rand(text.len());
                random_seq(rand, len)
            }
        }
    }

    /// Scores from the cache are those of a new alignment, whether the first request was for the
    /// exact score or a bound, and after the cache is cleared
    #[test]
    fn test_flank_alignment_scorer_cache() {
        let mut rand = xorshift(0x1357_9bdf_2468_ace0);
        let mut scorer = FlankAlignmentScorer::new();
        let mut fresh = ConsensusScoreAligner::new();
        for case in 0..30 {
            let text_len = 80 + rand(300);
            let text = random_seq(&mut rand, text_len);
            if case == 15 {
                // The next text clears the cache
                scorer.cache_bytes = FlankAlignmentScorer::MAX_CACHE_BYTES + 1;
            }
            scorer.set_text(&text);
            if case == 15 {
                assert_eq!(scorer.cache_bytes, text.len());
                assert!(scorer.scores.is_empty());
            }
            fresh.set_text(&text);
            let flanks = (0..4)
                .map(|_| random_allele_flank(&mut rand, &text))
                .collect::<Vec<_>>();
            for round in 0..2 {
                for flank in flanks.iter() {
                    let expected = fresh.align_score(flank);
                    // The first round asks for the bound first, and the second for the exact score
                    if round == 0 {
                        match scorer.align_score_bounded(flank) {
                            BoundedScore::Exact(x) => assert_eq!(x, expected),
                            BoundedScore::AtMost(x) => assert!(expected <= x),
                        }
                    }
                    assert_eq!(scorer.align_score(flank), expected, "case {case}");
                    assert!(matches!(
                        scorer.align_score_bounded(flank),
                        BoundedScore::Exact(x) if x == expected
                    ));
                }
            }
        }
    }

    /// The scores from get_ref_alt_flank_scores may be upper bounds, but only for a score whose bound
    /// is lower than the other's exact score, so the comparison of the two, and the best exact score, are
    /// those of the exact scores
    #[test]
    fn test_get_ref_alt_flank_scores_preserve_comparison() {
        let mut rand = xorshift(0x2468_ace0_1357_9bdf);
        let mut scorer = FlankAlignmentScorer::new();
        let mut exact_scorer = ConsensusScoreAligner::new();
        let mut bound_count = 0;
        for case in 0..60 {
            let text_len = 80 + rand(300);
            let text = random_seq(&mut rand, text_len);
            scorer.set_text(&text);
            exact_scorer.set_text(&text);
            let ref_flank = (rand(5) != 0).then(|| random_allele_flank(&mut rand, &text));
            let alt_flank = (rand(5) != 0).then(|| random_allele_flank(&mut rand, &text));

            let (ref_score, alt_score, best_exact) =
                get_ref_alt_flank_scores(&mut scorer, ref_flank.as_deref(), alt_flank.as_deref());

            let text_len = text.len() as f64;
            let mut exact =
                |x: &Option<Vec<u8>>| x.as_ref().map(|y| exact_scorer.dp_align_score(y));
            let ref_exact = exact(&ref_flank);
            let alt_exact = exact(&alt_flank);
            assert_eq!(ref_score.is_some(), ref_exact.is_some());
            assert_eq!(alt_score.is_some(), alt_exact.is_some());
            let to_score = |x: Option<f64>| x.map(|y| (y * text_len).round() as i32);
            let (ref_score, alt_score) = (to_score(ref_score), to_score(alt_score));

            // A returned score is the exact score, or an upper bound of it
            for (returned, expected) in [(ref_score, ref_exact), (alt_score, alt_exact)] {
                if let (Some(r), Some(e)) = (returned, expected) {
                    assert!(r >= e, "case {case}");
                    bound_count += usize::from(r > e);
                }
            }
            if let (Some(r), Some(a)) = (ref_score, alt_score) {
                assert_eq!(
                    r.cmp(&a),
                    ref_exact.cmp(&alt_exact),
                    "case {case}: {ref_score:?} {alt_score:?} vs {ref_exact:?} {alt_exact:?}"
                );
            }
            assert_eq!(best_exact, ref_exact.max(alt_exact), "case {case}");
        }
        // Some scores must have been replaced by bounds
        assert!(bound_count > 0);
    }
}
