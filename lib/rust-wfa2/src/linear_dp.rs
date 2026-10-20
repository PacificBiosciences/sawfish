//! Exact score of WFA2 gap-linear ends-free alignment with a match score, by dynamic programming
//!
//! This returns the same score as WFA2-lib for one aligner configuration: gap-linear penalties with a
//! negative match score, ends-free alignment and no heuristic. Only the score is computed, not the
//! CIGAR.
//!
//! WFA2's work grows with the square of the alignment penalty, so it is slow for dissimilar
//! sequences, where it visits each DP cell many times. This instead fills the whole DP matrix once,
//! by anti-diagonals so that the inner loop vectorizes.
//!
//! # Why the DP gives WFA2's result
//!
//! With match score M < 0, mismatch X and indel G, WFA2 aligns with penalties 0 for a match,
//! X' = 2X - 2M for a mismatch and G' = 2G - M for an indel. An ends-free begin at text position
//! h (cell (0,h)) or pattern position v (cell (v,0)) enters the wavefronts at penalty -M*h or -M*v.
//! WFA2 computes wavefronts of increasing penalty s and, after extending the matches of each one,
//! checks its diagonals k = h - v from low to high. It stops at the first diagonal whose furthest
//! point is an end cell: h = tlen with at most pattern_end_free pattern bases left, or v = plen with
//! at most text_end_free text bases left. The reported score is (-M*(v+h) - s)/2 for that end cell
//! (v,h).
//!
//! Let D(v,h) be the lowest penalty of a path to cell (v,h) within the DP matrix. On each diagonal
//! the end condition holds from one point onwards, so a diagonal's furthest point at penalty s is an
//! end cell exactly when an end cell on that diagonal can be reached at penalty s. WFA2's wavefronts
//! can also hold points outside the matrix (h > tlen or v > plen), because it only trims these at
//! the wavefront ends. Such a point can never return into the matrix, and its path can be replaced
//! by one along the matrix edge (indels instead of mismatches, no indels beyond the edge) that ends
//! in an end cell with a lower penalty. So points outside the matrix never meet the end condition
//! before the best end cell inside it. WFA2 therefore stops at s = the lowest D over all end cells,
//! at the end cell on the lowest diagonal among those with that D, which is what is computed here.
//!
//! `SAWFISH_VERIFY_FAST_WFA` in sawfish, and the tests below, compare this against WFA2.
//!

/// Penalties above this use the i32 version of the DP
const I16_MAX_PENALTY: i64 = 30_000;

/// Reusable buffers for the gap-linear ends-free DP
#[derive(Default)]
pub struct LinearEndsFreeDp {
    buf16: [Vec<i16>; 3],
    buf32: [Vec<i32>; 3],
    text_rev: Vec<u8>,
    /// Penalty of each cell in the last pattern row, (plen,h), indexed by h
    last_row: Vec<i32>,
    /// Penalty of each cell in the last text column, (v,tlen), indexed by v
    last_col: Vec<i32>,
}

/// Penalties after WFA2's transformation for a negative match score
#[derive(Clone, Copy)]
struct Penalties {
    /// Penalty per ends-free begin position (-M)
    begin: i32,
    mismatch: i32,
    indel: i32,
}

trait Cell: Copy + Ord {
    fn from_i32(x: i32) -> Self;
    fn to_i32(self) -> i32;
    fn add(self, x: Self) -> Self;
}

impl Cell for i16 {
    #[inline(always)]
    fn from_i32(x: i32) -> Self {
        x as i16
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        self as i32
    }
    #[inline(always)]
    fn add(self, x: Self) -> Self {
        self.wrapping_add(x)
    }
}

impl Cell for i32 {
    #[inline(always)]
    fn from_i32(x: i32) -> Self {
        x
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        self
    }
    #[inline(always)]
    fn add(self, x: Self) -> Self {
        self.wrapping_add(x)
    }
}

/// Fill the DP matrix by anti-diagonals d = v + h, recording the last row and column
///
/// Cell (v,h) of anti-diagonal d is stored at index v of its buffer. The text is reversed so that
/// the text base compared with pattern base v-1 on anti-diagonal d, text[d-v-1], is
/// text_rev[tlen-d+v], which moves forward with v like the other inputs.
///
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn fill<T: Cell>(
    bufs: &mut [Vec<T>; 3],
    pattern: &[u8],
    text_rev: &[u8],
    pen: Penalties,
    pattern_begin_free: usize,
    text_begin_free: usize,
    last_row: &mut [i32],
    last_col: &mut [i32],
) {
    let plen = pattern.len();
    let tlen = text_rev.len();
    for b in bufs.iter_mut() {
        b.clear();
        b.resize(plen + 1, T::from_i32(0));
    }
    let [b0, b1, b2] = bufs;
    // a: anti-diagonal d-2, b: d-1, c: d
    let (mut a, mut b, mut c) = (b0, b1, b2);
    let mismatch = T::from_i32(pen.mismatch);
    let indel = T::from_i32(pen.indel);

    for d in 0..=(plen + tlen) {
        // Interior cells, with v >= 1 and h >= 1
        let lo = std::cmp::max(1, d.saturating_sub(tlen));
        let hi = std::cmp::min(plen, d.saturating_sub(1));
        if lo <= hi {
            let n = hi + 1 - lo;
            let t0 = tlen + lo - d;
            let a_diag = &a[lo - 1..lo - 1 + n];
            let b_up = &b[lo - 1..lo - 1 + n];
            let b_left = &b[lo..lo + n];
            let p = &pattern[lo - 1..lo - 1 + n];
            let t = &text_rev[t0..t0 + n];
            let out = &mut c[lo..lo + n];
            for i in 0..n {
                let diag = a_diag[i].add(if p[i] == t[i] {
                    T::from_i32(0)
                } else {
                    mismatch
                });
                let gap = std::cmp::min(b_up[i], b_left[i]).add(indel);
                out[i] = std::cmp::min(diag, gap);
            }
        }

        // Edge cells (0,h=d) and (v=d,0)
        if d == 0 {
            c[0] = T::from_i32(0);
        } else {
            if d <= tlen {
                let gap = b[0].to_i32() + pen.indel;
                c[0] = T::from_i32(if d <= text_begin_free {
                    std::cmp::min(gap, pen.begin * d as i32)
                } else {
                    gap
                });
            }
            if d <= plen {
                let gap = b[d - 1].to_i32() + pen.indel;
                c[d] = T::from_i32(if d <= pattern_begin_free {
                    std::cmp::min(gap, pen.begin * d as i32)
                } else {
                    gap
                });
            }
        }

        if d >= plen && d - plen <= tlen {
            last_row[d - plen] = c[plen].to_i32();
        }
        if d >= tlen && d - tlen <= plen {
            last_col[d - tlen] = c[d - tlen].to_i32();
        }

        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut b, &mut c);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn fill_i16_avx2(
    bufs: &mut [Vec<i16>; 3],
    pattern: &[u8],
    text_rev: &[u8],
    pen: Penalties,
    pattern_begin_free: usize,
    text_begin_free: usize,
    last_row: &mut [i32],
    last_col: &mut [i32],
) {
    fill(
        bufs,
        pattern,
        text_rev,
        pen,
        pattern_begin_free,
        text_begin_free,
        last_row,
        last_col,
    )
}

#[cfg(target_arch = "x86_64")]
fn has_avx2() -> bool {
    static AVX2: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *AVX2.get_or_init(|| std::is_x86_feature_detected!("avx2"))
}

impl LinearEndsFreeDp {
    pub fn new() -> Self {
        Self::default()
    }

    /// Score of WFA2's gap-linear ends-free alignment of pattern to text
    ///
    /// The arguments are the same as for a WFA2 gap-linear aligner created with these match,
    /// mismatch and indel scores, followed by `align_ends_free`, and the result is the same as that
    /// aligner's `score()`. The match score must be negative and the sequences non-empty.
    ///
    #[allow(clippy::too_many_arguments)]
    pub fn score(
        &mut self,
        match_: i32,
        mismatch: i32,
        indel: i32,
        pattern: &[u8],
        text: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) -> i32 {
        assert!(match_ < 0 && mismatch > 0 && indel > 0);
        assert!(!pattern.is_empty() && !text.is_empty());
        let pen = Penalties {
            begin: -match_,
            mismatch: 2 * mismatch - 2 * match_,
            indel: 2 * indel - match_,
        };
        let plen = pattern.len();
        let tlen = text.len();
        let to_len = |x: i32| x.max(0) as usize;
        let (pbf, pef) = (to_len(pattern_begin_free), to_len(pattern_end_free));
        let (tbf, tef) = (to_len(text_begin_free), to_len(text_end_free));

        self.text_rev.clear();
        self.text_rev.extend(text.iter().rev());
        self.last_row.clear();
        self.last_row.resize(tlen + 1, 0);
        self.last_col.clear();
        self.last_col.resize(plen + 1, 0);

        // Every cell's penalty is at most the begin penalty of the larger begin-free length plus
        // the largest penalty for each base of both sequences
        let max_penalty = pen.begin as i64 * std::cmp::max(pbf, tbf) as i64
            + std::cmp::max(pen.mismatch, pen.indel) as i64 * (plen + tlen) as i64;
        if max_penalty <= I16_MAX_PENALTY {
            #[cfg(target_arch = "x86_64")]
            if has_avx2() {
                unsafe {
                    fill_i16_avx2(
                        &mut self.buf16,
                        pattern,
                        &self.text_rev,
                        pen,
                        pbf,
                        tbf,
                        &mut self.last_row,
                        &mut self.last_col,
                    )
                };
                return self.best_end_score(pen.begin, plen, tlen, pef, tef);
            }
            fill(
                &mut self.buf16,
                pattern,
                &self.text_rev,
                pen,
                pbf,
                tbf,
                &mut self.last_row,
                &mut self.last_col,
            );
        } else {
            fill(
                &mut self.buf32,
                pattern,
                &self.text_rev,
                pen,
                pbf,
                tbf,
                &mut self.last_row,
                &mut self.last_col,
            );
        }
        self.best_end_score(pen.begin, plen, tlen, pef, tef)
    }

    /// Find the end cell with the lowest penalty, then the lowest diagonal, and return its score
    fn best_end_score(&self, begin: i32, plen: usize, tlen: usize, pef: usize, tef: usize) -> i32 {
        // (penalty, diagonal, v + h) of the best end cell so far
        let mut best = (i32::MAX, i64::MAX, 0i64);
        let mut check = |penalty: i32, v: usize, h: usize| {
            let k = h as i64 - v as i64;
            if (penalty, k) < (best.0, best.1) {
                best = (penalty, k, (v + h) as i64);
            }
        };
        // Text is aligned, with at most pattern_end_free pattern bases left
        for v in plen.saturating_sub(pef)..=plen {
            check(self.last_col[v], v, tlen);
        }
        // Pattern is aligned, with at most text_end_free text bases left
        for h in tlen.saturating_sub(tef)..=tlen {
            check(self.last_row[h], plen, h);
        }
        ((begin as i64 * best.2 - best.0 as i64) / 2) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aligner::{
        AlignmentScope, AlignmentStatus, Heuristic, MemoryModel, WFAlignerGapLinear,
    };

    fn wfa2_score(penalties: (i32, i32, i32), pattern: &[u8], text: &[u8], clips: [i32; 4]) -> i32 {
        let mut aligner = WFAlignerGapLinear::new_with_match(
            penalties.0,
            penalties.1,
            penalties.2,
            AlignmentScope::Alignment,
            MemoryModel::MemoryHigh,
        );
        aligner.set_heuristic(Heuristic::None);
        let status = aligner.align_ends_free(pattern, text, clips[0], clips[1], clips[2], clips[3]);
        assert_eq!(status, AlignmentStatus::StatusAlgCompleted);
        aligner.score()
    }

    /// Compare against WFA2 on random sequence pairs, from near-identical to unrelated
    ///
    /// This passed for 200,000 pairs.
    #[test]
    fn test_linear_dp_matches_wfa2() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let bases = b"ACGTN";
        let mut dp = LinearEndsFreeDp::new();
        for iter in 0..5000 {
            let penalties = if iter % 3 == 0 {
                (-1, 3, 2)
            } else {
                (-2, 3, 4)
            };
            let tlen = 1 + rand(if iter % 10 == 0 { 600 } else { 120 });
            let text = (0..tlen).map(|_| bases[rand(4)]).collect::<Vec<_>>();
            // Mutate the text into the pattern at a random rate, so that alignments range from
            // near-identical to unrelated
            let rate = rand(101);
            let mut pattern = Vec::new();
            for &b in text.iter() {
                if rand(100) < rate {
                    match rand(4) {
                        0 => pattern.push(bases[rand(5)]),
                        1 => {}
                        2 => {
                            pattern.push(b);
                            pattern.push(bases[rand(4)]);
                        }
                        _ => pattern.push(bases[rand(4)]),
                    }
                } else {
                    pattern.push(b);
                }
            }
            if pattern.is_empty() {
                pattern.push(b'A');
            }
            // sawfish always uses equal clips of 40% of the shorter sequence. With unequal clips,
            // WFA2's result can depend on the aligner's earlier use, so those aren't compared.
            let c = (std::cmp::min(pattern.len(), text.len()) as f64 * 0.4) as i32;
            let clips = [c, c, c, c];
            let expected = wfa2_score(penalties, &pattern, &text, clips);
            let score = dp.score(
                penalties.0,
                penalties.1,
                penalties.2,
                &pattern,
                &text,
                clips[0],
                clips[1],
                clips[2],
                clips[3],
            );
            assert_eq!(
                score,
                expected,
                "iter {iter} clips {clips:?} pattern {} text {}",
                String::from_utf8_lossy(&pattern),
                String::from_utf8_lossy(&text)
            );
        }
    }
}

#[cfg(test)]
mod bench {
    use super::*;
    use crate::aligner::{
        AlignmentScope, AlignmentStatus, Heuristic, MemoryModel, WFAlignerGapLinear,
    };

    /// Timing of the DP and of score-only WFA2 with step limits, by sequence divergence
    ///
    /// Run with: cargo test --release bench_linear_dp -- --ignored --nocapture
    #[test]
    #[ignore]
    fn bench_linear_dp() {
        let mut state = 0x1234_5678_9abc_def1u64;
        let mut rand = |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let bases = b"ACGT";
        let n = 200;
        let mut dp = LinearEndsFreeDp::new();
        for rate in [0, 1, 2, 5, 10, 20, 30, 50, 100] {
            let mut pairs = Vec::new();
            for _ in 0..n {
                let text = (0..500).map(|_| bases[rand(4)]).collect::<Vec<_>>();
                let mut pattern = Vec::new();
                for &b in text.iter() {
                    if rand(100) < rate {
                        match rand(3) {
                            0 => {}
                            1 => {
                                pattern.push(b);
                                pattern.push(bases[rand(4)]);
                            }
                            _ => pattern.push(bases[rand(4)]),
                        }
                    } else {
                        pattern.push(b);
                    }
                }
                pairs.push((pattern, text));
            }
            let clip = |p: &[u8], t: &[u8]| (std::cmp::min(p.len(), t.len()) as f64 * 0.4) as i32;
            let t0 = std::time::Instant::now();
            let mut sum = 0i64;
            for (p, t) in pairs.iter() {
                let c = clip(p, t);
                sum += dp.score(-1, 3, 2, p, t, c, c, c, c) as i64;
            }
            let dp_us = t0.elapsed().as_secs_f64() * 1e6 / n as f64;
            let mut line = format!("rate {rate:3} dp {dp_us:7.1}us sum {sum}");
            for cap in [i32::MAX, 200, 120, 80, 50] {
                let mut a = WFAlignerGapLinear::new_with_match(
                    -1,
                    3,
                    2,
                    AlignmentScope::Alignment,
                    MemoryModel::MemoryHigh,
                );
                a.set_heuristic(Heuristic::None);
                if cap != i32::MAX {
                    a.set_max_alignment_steps(cap);
                }
                let t0 = std::time::Instant::now();
                let mut done = 0;
                for (p, t) in pairs.iter() {
                    let c = clip(p, t);
                    let st = a.align_ends_free(p, t, c, c, c, c);
                    if st == AlignmentStatus::StatusAlgCompleted {
                        done += 1;
                    }
                }
                line += &format!(
                    " | wfa2 cap {} {:6.1}us done {done}",
                    if cap == i32::MAX { 0 } else { cap },
                    t0.elapsed().as_secs_f64() * 1e6 / n as f64
                );
            }
            eprintln!("{line}");
        }
    }
}
