use spoa_sys::ffi;
use std::ffi::{CStr, CString};

pub use ffi::AlignmentType;

/// An opaque struct containing an alignment of a sequence to a `Graph`
pub struct Alignment(cxx::UniquePtr<ffi::Alignment>);

impl Alignment {
    /// Construct a new, empty `Alignment`. This is useful for adding the first sequence to a graph
    pub fn new() -> Self {
        Alignment(ffi::create_alignment())
    }

    /// True if the alignment has no aligned positions
    pub fn is_empty(&self) -> bool {
        ffi::alignment_is_empty(self.0.as_ref().unwrap())
    }
}
/// An opaque struct for aligning sequences to a `Graph`
pub struct AlignmentEngine(cxx::UniquePtr<ffi::AlignmentEngine>);

impl AlignmentEngine {
    /// Construct a new `AlignmentEngine`, given an `AlignmentType`, a match score (`m`), a
    /// mismatch score (`n`), a gap open score (`g`), a gap extension score (`e`), a second gap
    /// open score (`q`), and a second gap extension score (`c`).
    pub fn new(typ: AlignmentType, m: i8, n: i8, g: i8, e: i8, q: i8, c: i8) -> Self {
        AlignmentEngine(ffi::create_alignment_engine(typ, m, n, g, e, q, c))
    }

    /// Align a `sequence` to a `graph`, returning an `Alignment`
    pub fn align(&mut self, sequence: &CStr, graph: &Graph) -> (Alignment, i32) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        let mut score = 0;
        let aln = Alignment(unsafe {
            ffi::align(
                self.0.pin_mut(),
                sequence.as_ptr(),
                sequence_len,
                graph.0.as_ref().unwrap(),
                &mut score,
            )
        });
        (aln, score)
    }

    /// Align a `sequence` to a `graph`, with the same result as `align`, computed more efficiently
    /// if the alignment score is at least `min_score`
    ///
    /// Returns the alignment, its score, and true if the alignment is non-empty with a score of at
    /// least `min_score`. Without `full_fallback`, the alignment and score are only valid if the
    /// alignment reached `min_score`, or if the alignment is non-empty, which means that pruning
    /// was abandoned and the full alignment was computed.
    pub fn align_pruned(
        &mut self,
        sequence: &CStr,
        graph: &Graph,
        min_score: i32,
        full_fallback: bool,
    ) -> (Alignment, i32, bool) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        let mut score = 0;
        let mut reached = false;
        let aln = Alignment(unsafe {
            ffi::align_pruned(
                self.0.pin_mut(),
                sequence.as_ptr(),
                sequence_len,
                graph.0.as_ref().unwrap(),
                min_score,
                full_fallback,
                &mut score,
                &mut reached,
            )
        });
        (aln, score, reached)
    }

    /// Upper bound of the score of any alignment of a sequence of length `sequence_len` to `graph`
    pub fn max_possible_score(&mut self, sequence_len: usize, graph: &Graph) -> i64 {
        ffi::max_possible_score(
            self.0.pin_mut(),
            u32::try_from(sequence_len).unwrap(),
            graph.0.as_ref().unwrap(),
        )
    }

    /// Align a `sequence` to a `graph`, with the same result as `align`, computed more efficiently
    /// if the alignment score is close to one of `expected_fracs` of the maximum possible score
    pub fn align_pruned_auto(
        &mut self,
        sequence: &CStr,
        graph: &Graph,
        expected_fracs: &[f64],
    ) -> (Alignment, i32) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        let mut score = 0;
        let aln = Alignment(unsafe {
            ffi::align_pruned_auto(
                self.0.pin_mut(),
                sequence.as_ptr(),
                sequence_len,
                graph.0.as_ref().unwrap(),
                expected_fracs.as_ptr(),
                expected_fracs.len() as u32,
                &mut score,
            )
        });
        (aln, score)
    }
}

/// An opaque struct for the partial order alignment graph
pub struct Graph(cxx::UniquePtr<ffi::Graph>);

impl Graph {
    /// Construct a new, empty `Graph`.
    pub fn new() -> Self {
        Graph(ffi::create_graph())
    }

    /// Clear graph
    pub fn clear(&mut self) {
        unsafe {
            ffi::clear(self.0.pin_mut());
        }
    }

    /// Add an `Alignment` to a `Graph`, along with its `sequence` and per-base `quality` scores.
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &CStr, quality: &CStr) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        let quality_len = u32::try_from(quality.to_bytes().len()).unwrap();
        assert!(sequence_len == quality_len);
        unsafe {
            ffi::add_alignment(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr(),
                sequence_len,
                quality.as_ptr(),
                quality_len,
            )
        }
    }

    /// Add an `Alignment` to a `Graph`, along with its `sequence`
    ///
    /// The `alignment` should be derived from calling `AlignmentEngine::align` with the
    /// `sequence`.
    pub fn add_alignment_noqual(&mut self, alignment: &Alignment, sequence: &CStr) {
        let sequence_len = u32::try_from(sequence.to_bytes().len()).unwrap();
        unsafe {
            ffi::add_alignment_noqual(
                self.0.pin_mut(),
                alignment.0.as_ref().unwrap(),
                sequence.as_ptr(),
                sequence_len,
            )
        }
    }

    /// Generate a consenus sequence from the partial order `Graph`.
    pub fn consensus(&mut self) -> CString {
        let mut buf = Vec::from(
            ffi::generate_consensus(self.0.pin_mut())
                .as_ref()
                .unwrap()
                .as_bytes(),
        );
        buf.push(0);
        CString::from_vec_with_nul(buf).unwrap()
    }

    /// Generate a multiple sequence alignment for all sequences added to the `Graph`.
    ///
    /// If `include_consensus` is provided, the consensus sequence is provided, also aligned, at
    /// the end.
    pub fn multiple_sequence_alignment(&mut self, include_consensus: bool) -> Vec<CString> {
        let msa = ffi::generate_multiple_sequence_alignment(self.0.pin_mut(), include_consensus);
        let mut result = vec![];
        for aln in msa.as_ref().unwrap().iter() {
            let mut buf = Vec::from(aln.as_bytes());
            buf.push(0);
            result.push(CString::from_vec_with_nul(buf).unwrap());
        }
        result
    }

    pub fn get_sequence_count(&self) -> u32 {
        ffi::get_sequence_count(self.0.as_ref().unwrap())
    }
}

impl Default for Graph {
    fn default() -> Self {
        Graph::new()
    }
}

/// Identify the number of read bases used in the alignment (subtracting any clipping on the left and right sides)
///
pub fn get_alignment_overlap_size(alignment: &Alignment) -> u32 {
    ffi::get_alignment_overlap_size(alignment.0.as_ref().unwrap())
}

/// Return the read base index at the start of the alignment, or -1 if no aligned base can be found
///
pub fn get_alignment_clip_size(alignment: &Alignment) -> i32 {
    ffi::get_alignment_clip_size(alignment.0.as_ref().unwrap())
}

#[cfg(test)]
mod fastq;

#[cfg(test)]
mod tests {
    use super::*;

    // thanks to pjedge for this small
    const SMALL_SEQS: &[&CStr] = unsafe {
        &[
            CStr::from_bytes_with_nul_unchecked(b"ATTGCCCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCCCGAT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AACGCCCGTC\0"),
            CStr::from_bytes_with_nul_unchecked(b"AGTGCTCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"AATGCTCGTT\0"),
        ]
    };

    const PROG_SEQS: &[&CStr] = unsafe {
        &[
            CStr::from_bytes_with_nul_unchecked(b"ATTGCCCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"GCCGTTTGG\0"),
            CStr::from_bytes_with_nul_unchecked(b"CGATGGCGT\0"),
            CStr::from_bytes_with_nul_unchecked(b"CGTCGGCGTT\0"),
            CStr::from_bytes_with_nul_unchecked(b"TGGCCTTACC\0"),
            CStr::from_bytes_with_nul_unchecked(b"GCCATACCC\0"),
        ]
    };

    #[test]
    fn consensus_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let qual = {
                let mut qual = vec![34u8; seq.to_bytes().len()];
                qual.push(0);
                CString::from_vec_with_nul(qual).unwrap()
            };
            let (aln, _score) = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, &qual);
        }

        let consensus = graph.consensus();
        assert!(consensus.as_bytes() == b"AATGCCCGTT");
    }

    #[test]
    fn consensus_small_noqual() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let (aln, _score) = eng.align(seq, &graph);
            graph.add_alignment_noqual(&aln, seq);
        }

        let consensus = graph.consensus();
        assert!(consensus.as_bytes() == b"AATGCCCGTT");
    }

    #[test]
    fn test_get_sequence_count() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let qual = {
                let mut qual = vec![34u8; seq.to_bytes().len()];
                qual.push(0);
                CString::from_vec_with_nul(qual).unwrap()
            };
            let (aln, _score) = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, &qual);
        }

        assert_eq!(graph.get_sequence_count() as usize, SMALL_SEQS.len());
    }

    #[test]
    fn test_alignment_overlap_size() {
        let mut eng = AlignmentEngine::new(AlignmentType::kOV, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        let mut vals = Vec::new();
        for seq in PROG_SEQS {
            let (aln, _score) = eng.align(seq, &graph);
            let val = get_alignment_overlap_size(&aln);
            vals.push(val);
            graph.add_alignment_noqual(&aln, seq);
        }
        assert_eq!(vals, vec![0, 6, 6, 9, 7, 8])
    }

    #[test]
    fn test_alignment_clip_size() {
        let mut eng = AlignmentEngine::new(AlignmentType::kSW, 1, -2, -3, -1, -3, -1);
        let mut graph = Graph::new();

        let mut vals = Vec::new();
        for seq in PROG_SEQS.iter().rev() {
            let (aln, _score) = eng.align(seq, &graph);
            let val = get_alignment_clip_size(&aln);
            vals.push(val);
            graph.add_alignment_noqual(&aln, seq);
        }
        /*
                let msa = graph.multiple_sequence_alignment(false);
                for line in msa.iter() {
                    println!("{:?}", line);
                }
        */
        assert_eq!(vals, vec![-1, 2, 4, 3, 2, 5])
    }

    #[test]
    fn consensus_spoa() {
        // test borrowed from spoa itself, obviously
        let file = std::fs::File::open("spoa-sys/spoa/test/data/sample.fastq.gz").unwrap();
        let gz = flate2::read::GzDecoder::new(file);
        let reader = fastq::FastqReader::new(std::io::BufReader::new(gz));

        let mut eng = AlignmentEngine::new(AlignmentType::kSW, 5, -4, -8, -6, -8, -6);
        let mut graph = Graph::new();

        for record in reader {
            let record = record.unwrap();
            let (aln, _score) = eng.align(&record.seq, &graph);
            graph.add_alignment(&aln, &record.seq, &record.qual);
        }

        let consensus = graph.consensus();
        let expected = b"AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAATGCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGTTGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGCTTTACACGCGCAACCAAGGATTTCGG";
        assert!(consensus.as_bytes() == expected);
    }

    #[test]
    fn test_msa_small() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
        let mut graph = Graph::new();

        for seq in SMALL_SEQS {
            let qual = {
                let mut qual = vec![34u8; seq.to_bytes().len()];
                qual.push(0);
                CString::from_vec_with_nul(qual).unwrap()
            };
            let (aln, _score) = eng.align(seq, &graph);
            graph.add_alignment(&aln, seq, &qual);
        }

        let msa = graph.multiple_sequence_alignment(true);
        assert!(msa.len() == 7);
        assert!(msa[0].as_bytes() == b"ATTGCC-CGTT");
        assert!(msa[1].as_bytes() == b"AATG-C-CGTT");
        assert!(msa[2].as_bytes() == b"AATGCC-CGAT");
        assert!(msa[3].as_bytes() == b"AACGCC-CGTC");
        assert!(msa[4].as_bytes() == b"AGTG-CTCGTT");
        assert!(msa[5].as_bytes() == b"AATG-CTCGTT");
        assert!(msa[6].as_bytes() == b"AATGCC-CGTT");
    }

    /// Verify how gap scoring works
    #[test]
    fn test_gap_scoring() {
        let mut eng = AlignmentEngine::new(AlignmentType::kNW, 1, -4, -3, -2, -3, -2);
        let mut graph = Graph::new();

        let seq1 = unsafe { CStr::from_bytes_with_nul_unchecked(b"ATTGCCCGTT\0") };
        let seq2 = unsafe { CStr::from_bytes_with_nul_unchecked(b"ATTGCCGTT\0") };

        let (aln, _score) = eng.align(seq1, &graph);
        graph.add_alignment_noqual(&aln, seq1);
        let (_aln, score) = eng.align(seq2, &graph);

        assert_eq!(score, 6);
    }
}
