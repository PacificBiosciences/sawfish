mod add_flanks;
mod banded_pairwise_assembler;
mod spoa_assembler;

use super::{RefineSVSettings, trimmed_reads::TrimmedReadInfo};
use crate::int_range::IntRange;

pub struct AssemblyResultContig {
    /// Region of the contig assembled into a consensus sequence by POA
    ///
    /// The contig may also include optional left and right flanks at lower consensus quality
    ///
    pub high_quality_range: IntRange,

    /// Full contig sequence
    pub seq: Vec<u8>,

    /// Fake flank size used to construct this contig (only used for debug at present)
    #[allow(unused)]
    pub fake_flank_size: i64,
}

pub struct AssemblyResult {
    /// Index of reads used to assemble the consensus sequence (from the input trimmed reads list)
    pub source_reads: Vec<usize>,

    /// Typical indel-like SVs should have only one contig, but for some cases the aligner will iterate
    /// to a larger contig if the initial contig alignment fails.
    ///
    pub contigs: Vec<AssemblyResultContig>,
}

impl AssemblyResult {
    pub fn supporting_read_count(&self) -> usize {
        self.source_reads.len()
    }
}

/// Assemble trimmed reads and return the top N assembly sequences
///
/// * cluster_index - Used for debug messages only
/// * max_assembly_count - Maximum number of assemblies to report out. Often this is set to
///   sample ploidy. For some cases it's limited to less than sample ploidy to reduce false positives.
/// * is_large_insertion_inferred_from_clipping - If true, this may be used to control the runtime of large assemblies
/// * is_multi_region - false if this is a single region indel assembly job
///
// TODO: runtime optimization for large insertions should be automated instead of a function arg
//
pub(super) fn assemble_and_cluster_reads(
    refine_settings: &RefineSVSettings,
    cluster_index: usize,
    max_assembly_count: usize,
    is_large_insertion_inferred_from_clipping: bool,
    trimmed_reads: Vec<TrimmedReadInfo>,
    is_multi_region: bool,
) -> Vec<AssemblyResult> {
    let assembly_debug = false;
    if assembly_debug {
        eprintln!("Starting trimmed read assembly process for cluster {cluster_index}");
    }

    if trimmed_reads.is_empty() {
        return Vec::new();
    }

    // This limits total candidate contig assemblies to control runtime in high depth regions
    let max_candidate_assemblies = if is_large_insertion_inferred_from_clipping {
        refine_settings.max_candidate_assemblies_at_large_insertion
    } else {
        refine_settings.max_candidate_assemblies
    };

    // Disable banded_pairwise assembler for all cases for now. It may use less memory from larger insertions but still
    // unclear if this is how we want to do it
    //
    let big_asm_mode = false;

    /*
    // Determine if we use alt assembly routine for larger reads:
    let min_big_asm_len = 8_000;
    let big_asm_mode = trimmed_reads
        .iter()
        .any(|x| x.trimmed_read_len() > min_big_asm_len);
    */

    if big_asm_mode {
        banded_pairwise_assembler::assemble_and_cluster_reads(
            refine_settings,
            cluster_index,
            max_candidate_assemblies,
            max_assembly_count,
            trimmed_reads,
            is_multi_region,
            assembly_debug,
        )
    } else {
        spoa_assembler::assemble_and_cluster_reads(
            refine_settings,
            cluster_index,
            max_candidate_assemblies,
            max_assembly_count,
            trimmed_reads,
            is_multi_region,
            assembly_debug,
        )
    }
}
