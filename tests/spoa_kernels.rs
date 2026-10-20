//! Check spoa's low-memory linear alignment and its pruned alignments against the full-matrix
//! alignment, on random graphs and around the lengths where scores switch from 16 to 32 bits
//!
//! `SPOA_VERIFY_LOW_MEMORY_LINEAR` makes spoa compute the full-matrix alignment for every low-memory
//! linear alignment and abort if the alignment or score differ. The pruned alignments are compared
//! with the unpruned alignment here, by their scores and by the graphs built from them.
//!
//! The AVX2 kernels are used where the CPU has them. `spoa_kernels_without_avx2` repeats the tests
//! with `SPOA_DISABLE_AVX2` set, for the kernels that are used otherwise (SSE or NEON, which is
//! also what runs on arm64 machines).

use std::ffi::CString;
use std::sync::Once;

use spoa::{Alignment, AlignmentEngine, AlignmentType, Graph};

fn init() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        // Safe because this runs before any alignment, and every test calls it first
        unsafe { std::env::set_var("SPOA_VERIFY_LOW_MEMORY_LINEAR", "1") };
    });
}

struct Rng(u64);

impl Rng {
    fn next(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }

    fn seq(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[self.next(4)]).collect()
    }

    /// Copy of seq with SNVs and small indels (at per_thousand per 1000 bases), and occasional
    /// deletions of up to 400 bases, which give the graph edges between distant nodes
    fn mutate(&mut self, seq: &[u8], per_thousand: usize, long_deletions: bool) -> Vec<u8> {
        let mut x = Vec::with_capacity(seq.len());
        let mut i = 0;
        while i < seq.len() {
            if self.next(1000) < per_thousand {
                match self.next(5) {
                    0 | 1 => x.push(b"ACGT"[self.next(4)]),
                    2 => (0..1 + self.next(5)).for_each(|_| x.push(b"ACGT"[self.next(4)])),
                    3 => i += self.next(5),
                    _ if long_deletions => i += 20 + self.next(380),
                    _ => i += self.next(5),
                }
            }
            if i < seq.len() {
                x.push(seq[i]);
            }
            i += 1;
        }
        if x.is_empty() {
            x.push(b'A');
        }
        x
    }
}

/// Weights as (match, mismatch, gap), used as a linear gap model
const WEIGHT_SETS: [(i8, i8, i8); 3] = [
    // The weights used for contig assembly
    (1, -3, -1),
    (2, -5, -3),
    (5, -4, -8),
];

fn engines() -> Vec<(String, AlignmentEngine)> {
    let mut x = Vec::new();
    for (typ, label) in [(AlignmentType::kSW, "SW"), (AlignmentType::kOV, "OV")] {
        for (m, n, g) in WEIGHT_SETS {
            x.push((
                format!("{label} {m} {n} {g}"),
                AlignmentEngine::new(typ, m, n, g, g, g, g),
            ));
        }
    }
    x
}

/// Graphs built from the same sequences, by different alignment methods that must agree
struct Graphs {
    full: Graph,
    pruned: Graph,
    auto: Graph,
}

impl Graphs {
    fn new() -> Self {
        Self {
            full: Graph::new(),
            pruned: Graph::new(),
            auto: Graph::new(),
        }
    }

    /// Align seq to the graphs, using the full alignment for one and the pruned alignments, with
    /// minimum scores below, at and above the optimum, for the others
    fn add(&mut self, engine: &mut AlignmentEngine, seq: &[u8], step: usize, label: &str) {
        let cseq = CString::new(seq).unwrap();
        if step == 0 {
            for graph in [&mut self.full, &mut self.pruned, &mut self.auto] {
                graph.add_alignment_noqual(&Alignment::new(), &cseq);
            }
            return;
        }

        let (full_alignment, full_score) = engine.align(&cseq, &self.full);
        let full_found = !full_alignment.is_empty();
        let optimum = if full_found { full_score } else { 0 };

        // Minimum scores at and around the optimum, and values that don't prune at all
        let min_scores = [i32::MIN / 4, -7, 0, 1, optimum - 1, optimum, optimum + 1];
        let min_score = min_scores[step % min_scores.len()];

        // With a full fallback the result is always that of the full alignment
        let (alignment, score, reached) = engine.align_pruned(&cseq, &self.pruned, min_score, true);
        assert_eq!(
            reached,
            !alignment.is_empty() && score >= min_score,
            "{label} step {step} min score {min_score}"
        );
        if full_found {
            assert_eq!(score, full_score, "{label} step {step} min {min_score}");
        } else {
            assert!(alignment.is_empty(), "{label} step {step}");
        }
        self.pruned.add_alignment_noqual(&alignment, &cseq);

        // Without it, the result is only valid if the minimum score was reached, or if the
        // pruning was abandoned
        let (unpruned, score, reached) = engine.align_pruned(&cseq, &self.pruned, min_score, false);
        if reached || !unpruned.is_empty() {
            let (_, check_score) = engine.align(&cseq, &self.pruned);
            assert_eq!(score, check_score, "{label} step {step} min {min_score}");
        }

        let expected_fracs: &[&[f64]] = &[&[0.5], &[0.1, 0.9], &[0.0, 0.3, 0.6, 0.99]];
        let (alignment, score) =
            engine.align_pruned_auto(&cseq, &self.auto, expected_fracs[step % 3]);
        if full_found {
            assert_eq!(score, full_score, "{label} step {step} auto");
        } else {
            assert!(alignment.is_empty(), "{label} step {step} auto");
        }
        self.auto.add_alignment_noqual(&alignment, &cseq);
        self.full.add_alignment_noqual(&full_alignment, &cseq);
    }

    /// The graphs built by every method must be the same
    fn check_same(&mut self, label: &str) {
        let consensus = self.full.consensus();
        assert!(consensus == self.pruned.consensus(), "{label}");
        assert!(consensus == self.auto.consensus(), "{label}");
        let msa = self.full.multiple_sequence_alignment(true);
        assert!(
            msa == self.pruned.multiple_sequence_alignment(true),
            "{label}"
        );
        assert!(
            msa == self.auto.multiple_sequence_alignment(true),
            "{label}"
        );
    }
}

/// Graphs of related and unrelated sequences, with long-range edges and nodes with several
/// predecessors
#[test]
fn spoa_random_graphs() {
    init();
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for (label, mut engine) in engines() {
        for graph_index in 0..10 {
            let len = [40, 150, 400, 900, 2500][graph_index % 5] + rng.next(60);
            let base = rng.seq(len);
            let mut graphs = Graphs::new();
            for step in 0..9 {
                let seq = match step % 4 {
                    // Related with small variants, and with long deletions
                    0 | 1 => rng.mutate(&base, 30, false),
                    2 => rng.mutate(&base, 60, true),
                    // Unrelated sequence, or a short one, which may not align at all
                    _ if rng.next(2) == 0 => rng.seq(len / 2 + 1),
                    _ => rng.mutate(&base[len / 4..], 20, true),
                };
                graphs.add(&mut engine, &seq, step, &label);
            }
            graphs.check_same(&format!("{label} graph {graph_index}"));
        }
    }
}

/// Sequences around the length where scores no longer fit in 16 bits, which depends on the
/// weights, aligned to a small graph so that the full-matrix check stays small
#[test]
fn spoa_int16_boundary() {
    init();
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let base = rng.seq(800);
    // Lengths either side of the one at which scores need 32 bits (28,930 for the assembly
    // weights, from IsInt16ScoreRange in sisd_alignment_engine.cpp)
    let boundary_lengths = [[28_900, 28_960], [9_550, 9_650], [3_550, 3_650]];
    for (typ, label) in [(AlignmentType::kSW, "SW"), (AlignmentType::kOV, "OV")] {
        for ((m, n, g), lengths) in WEIGHT_SETS.into_iter().zip(boundary_lengths) {
            let label = format!("{label} {m} {n} {g}");
            let mut engine = AlignmentEngine::new(typ, m, n, g, g, g, g);
            let mut graphs = Graphs::new();
            graphs.add(&mut engine, &base, 0, &label);
            let variant = rng.mutate(&base, 40, true);
            graphs.add(&mut engine, &variant, 1, &label);
            for (step, len) in lengths.into_iter().enumerate() {
                // Copies of the graph's sequence with variants, up to the length wanted
                let mut seq = Vec::new();
                while seq.len() < len {
                    seq.extend(rng.mutate(&base, 30, false));
                }
                seq.truncate(len);
                graphs.add(&mut engine, &seq, 2 + step, &label);
            }
            graphs.check_same(&format!("{label} boundary"));
        }
    }
}

/// Repeat the tests with the non-AVX2 kernels
#[test]
#[cfg(target_arch = "x86_64")]
fn spoa_kernels_without_avx2() {
    let status = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["spoa_", "--skip", "without_avx2"])
        .env("SPOA_DISABLE_AVX2", "1")
        .status()
        .unwrap();
    assert!(status.success());
}
