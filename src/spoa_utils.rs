//! Utilities to help interface with the spoa library via the spoa-rs crate
//!

use std::ffi::{CStr, CString};

use rust_htslib::bam::record::Cigar;
use spoa::{get_alignment_clip_size, get_alignment_overlap_size, AlignmentEngine, Graph};

use crate::simple_alignment::SimpleAlignment;

/// Create a new poa graph with only the given sequence
///
pub fn get_single_seq_graph(seq: &CStr) -> Graph {
    let mut new_graph = Graph::new();
    new_graph.add_alignment_noqual(&spoa::Alignment::new(), seq);
    new_graph
}

/// Align query to text sequence, and return the corresponding ascii MSA with additional alignment
/// details
///
/// Returns a 4-tuple of:
/// 1. alignment score
/// 2. left-side query soft-clip
/// 3. overlap size
/// 4. strings formatted as an ascii MSA
///
pub fn get_pairwise_alignment(
    engine: &mut AlignmentEngine,
    text: &CStr,
    query: &CStr,
) -> (i32, u32, u32, Vec<CString>) {
    let mut tmp_graph = get_single_seq_graph(text);

    let (aln, score) = engine.align(query, &tmp_graph);
    let query_clip = get_alignment_clip_size(&aln);
    let overlap = get_alignment_overlap_size(&aln);
    assert!(query_clip >= 0);
    tmp_graph.add_alignment_noqual(&aln, query);

    (
        score,
        query_clip as u32,
        overlap,
        tmp_graph.multiple_sequence_alignment(false),
    )
}

/// Align sequence to the consensus of the given graph, and return the corresponding ascii MSA
///
/// Returns a 4-tuple of:
/// 1. alignment score
/// 2. left-side query soft-clip
/// 3. overlap size
/// 4. strings formatted as an ascii MSA
///
pub fn get_alignment_to_graph_consensus(
    engine: &mut AlignmentEngine,
    graph: &mut Graph,
    query: &CStr,
) -> (i32, u32, u32, Vec<CString>) {
    let consensus = graph.consensus();
    get_pairwise_alignment(engine, &consensus, query)
}

/// Handles wrapping long alignments into multiple rows
///
#[allow(dead_code)]
pub fn print_msa(alignments: &[CString]) {
    const WIDTH: usize = 100;

    if alignments.is_empty() {
        return;
    }

    let bytes = alignments.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();

    let len = bytes.first().unwrap().len();
    for &seq in bytes.iter() {
        assert_eq!(seq.len(), len);
    }

    let rows = (len + WIDTH - 1) / WIDTH;

    let ruler = {
        let mut ruler = String::new();
        for i in 0..WIDTH {
            if (i + 1) % 10 == 0 {
                ruler.push('.');
            } else {
                ruler.push(' ');
            }
        }
        ruler
    };

    for row_index in 0..rows {
        let start = WIDTH * row_index;
        let end = std::cmp::min(start + WIDTH, len);
        eprintln!("{}", ruler);
        for &seq in bytes.iter() {
            eprintln!("{}", std::str::from_utf8(&seq[start..end]).unwrap());
        }
        // Print line to highlight mismatches:
        let mut mismatches = String::new();
        for pos in start..end {
            let mut is_mismatch = false;
            let mut cons = None;
            for &seq in bytes.iter() {
                if let Some(cons) = cons {
                    if seq[pos] != cons {
                        is_mismatch = true;
                        break;
                    }
                } else {
                    cons = Some(seq[pos]);
                }
            }
            mismatches.push(if is_mismatch { 'X' } else { ' ' });
        }
        eprintln!("{mismatches}");
        eprintln!();
    }
}

/// Handles wrapping long alignments into multiple rows
///
#[allow(dead_code)]
pub fn print_fasta(alignments: &[CString]) {
    const WIDTH: usize = 100;

    if alignments.is_empty() {
        return;
    }

    let bytes = alignments.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();

    for (seq_index, seq) in bytes.iter().enumerate() {
        let len = seq.len();
        let rows = (len + WIDTH - 1) / WIDTH;

        eprintln!("> {}", seq_index);
        for row_index in 0..rows {
            let start = WIDTH * row_index;
            let end = std::cmp::min(start + WIDTH, len);
            eprintln!("{}", std::str::from_utf8(&seq[start..end]).unwrap());
        }
    }
}

/// Use spoa to get a consensus of two reads using the given alignment engine
///
#[allow(dead_code)]
fn get_two_read_consensus(
    alignment_engine: &mut AlignmentEngine,
    read1: &CStr,
    read2: &CStr,
) -> CString {
    let mut tmp_graph = get_single_seq_graph(read1);

    let (alignment, _score) = alignment_engine.align(read2, &tmp_graph);
    tmp_graph.add_alignment_noqual(&alignment, read2);
    tmp_graph.consensus()
}

/// Produces a CIGAR-like alignment summary from a 2 sequence MSA
///
/// This is only intended as a stop-gap test measure before trying to get spoa to output this type
/// of result in a more efficient direct way
///
/// # Arguments
/// * query_start_clip - The query position (zero-indexed) where the query-to ref alignment starts
///   in the `alignments` strings
/// * query_size - Full query length
/// * alignments - 2 strings formatted to pretty-print an alignment between ref (1st) and query (2nd)
///   sequences
///
/// Returns a 2-tuple of (1) the reference sequence offset for the start of the alignment and (2)
/// the alignment cigar. The alignment cigar includes left and right soft clipping of the query
/// sequence so that the full query length is represented.
///
pub fn get_query_alignment_from_msa(
    query_start_clip: u32,
    query_size: u32,
    alignments: &[CString],
) -> SimpleAlignment {
    assert_eq!(alignments.len(), 2);

    let bytes = alignments.iter().map(|x| x.as_bytes()).collect::<Vec<_>>();

    let len = bytes.first().unwrap().len();
    for &seq in bytes.iter() {
        assert_eq!(seq.len(), len);
    }

    // Find the start and end of the alignment (any edge will not be marked as soft-slip)
    let mut ref_offset = 0;
    let mut first_match_index = None;
    for msa_index in 0..len {
        let r = bytes[0][msa_index];
        let q = bytes[1][msa_index];

        if r != b'-' && q != b'-' {
            first_match_index = Some(msa_index);
            break;
        }
        if r != b'-' {
            ref_offset += 1;
        }
    }

    let first_match_index = match first_match_index {
        Some(i) => i,
        None => {
            return SimpleAlignment::new();
        }
    };

    let mut last_match_index = 0;
    for msa_index in (0..len).rev() {
        let r = bytes[0][msa_index];
        let q = bytes[1][msa_index];

        if r != b'-' && q != b'-' {
            last_match_index = msa_index;
            break;
        }
    }

    let mut query_pos = query_start_clip;
    let mut state = Cigar::Match(0);
    let mut cigar = Vec::new();
    if query_start_clip > 0 {
        cigar.push(Cigar::SoftClip(query_start_clip));
    }
    for msa_index in first_match_index..=last_match_index {
        let r = bytes[0][msa_index];
        let q = bytes[1][msa_index];

        let new_state = if r != b'-' && q != b'-' {
            query_pos += 1;
            if r == q {
                Cigar::Equal(1)
            } else {
                Cigar::Diff(1)
            }
        } else if r != b'-' {
            Cigar::Del(1)
        } else {
            query_pos += 1;
            Cigar::Ins(1)
        };

        if std::mem::discriminant(&new_state) == std::mem::discriminant(&state) {
            use Cigar::*;
            if let Equal(ref mut n) | Diff(ref mut n) | Del(ref mut n) | Ins(ref mut n) = state {
                *n += 1;
            }
        } else {
            if state != Cigar::Match(0) {
                cigar.push(state)
            }
            state = new_state;
        }
    }
    if state != Cigar::Match(0) {
        cigar.push(state)
    }

    assert!(query_size >= query_pos);
    let query_end_clip = query_size - query_pos;
    if query_end_clip > 0 {
        cigar.push(Cigar::SoftClip(query_end_clip));
    }

    SimpleAlignment { ref_offset, cigar }
}

/// Hide details behind a standard-ish interface to make it easy to change aligners in future
///
/// Use spoa to align for convenience even though we don't need the graph
#[allow(dead_code)]
pub struct PairwiseAligner {
    engine: AlignmentEngine,
    graph: Graph,
}

#[allow(dead_code)]
impl PairwiseAligner {
    pub fn new_large_indel_aligner() -> Self {
        // For this case we need a convex alignment that allows very long indels.
        // We also lower the mismatch cost compared to assembly to allow for SNPs
        Self {
            engine: AlignmentEngine::new(spoa::AlignmentType::kSW, 1, -2, -2, -1, -30, 0),
            graph: Graph::new(),
        }
    }

    /// can we delete get_spoa_graph_for_assembly_region now?
    pub fn set_text(&mut self, text: &[u8]) {
        // Convert text to a CString for use in spoa
        let mut t = text.to_vec();
        t.push(0u8);
        let t = unsafe { CStr::from_bytes_with_nul_unchecked(&t) };

        let (alignment, _score) = self.engine.align(t, &self.graph);
        self.graph.add_alignment_noqual(&alignment, t);
    }

    pub fn align_pattern(&mut self, pattern: &[u8]) -> Option<SimpleAlignment> {
        // Convert pattern to a CString for use in spoa
        let mut p = pattern.to_vec();
        p.push(0u8);
        let p = unsafe { CStr::from_bytes_with_nul_unchecked(&p) };

        let (_score, query_start_clip, _query_overlap, msa) =
            get_alignment_to_graph_consensus(&mut self.engine, &mut self.graph, p);
        Some(get_query_alignment_from_msa(
            query_start_clip,
            pattern.len() as u32,
            &msa,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_cigar_from_msa() {
        let test_align = [
            CString::from_vec_with_nul(b"TATTTC-CCCGTT-\0".to_vec()).unwrap(),
            CString::from_vec_with_nul(b"-A--TGTCC-GTTA\0".to_vec()).unwrap(),
        ];

        let qa = get_query_alignment_from_msa(0, 10, &test_align);

        assert_eq!(qa.ref_offset, 1);

        use Cigar::*;
        let expected_cigar = vec![
            Equal(1),
            Del(2),
            Equal(1),
            Diff(1),
            Ins(1),
            Equal(2),
            Del(1),
            Equal(3),
            SoftClip(1),
        ];
        assert_eq!(qa.cigar, expected_cigar);
    }
}
