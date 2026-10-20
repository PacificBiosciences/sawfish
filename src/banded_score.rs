//! Score-only version of rust-bio's banded aligner (bio 2.3.0, `alignment::pairwise::banded`)
//!
//! rust-bio's banded aligner computes its DP only within a band around k-mer matches, but always
//! allocates and initializes a full (m+1) x (n+1) traceback matrix. For the multi-kilobase contigs
//! compared when merging haplotypes across samples, this takes hundreds of MB to GB per alignment
//! and most of the merge time. Haplotype merging only uses the alignment score.
//!
//! This file copies the band construction unchanged, and the DP of `Aligner::compute_alignment`
//! with every traceback statement removed. The traceback is only written by the DP, never read
//! for a score, so the score is the same as `Aligner::custom_with_matches(..).score`.
//!
//! Original code: Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks, licensed
//! under the MIT license (http://opensource.org/licenses/MIT).
//!

// The copied code is kept as close to the original as possible
#![allow(clippy::manual_repeat_n, clippy::explicit_counter_loop)]

use std::cmp::{Ordering, max, min};
use std::iter::repeat;
use std::ops::Range;

use bio::alignment::pairwise::{MIN_SCORE, MatchFunc, Scoring};
use bio::alignment::sparse;
use bio::utils::TextSlice;

const MAX_CELLS: usize = 5_000_000;
const DEFAULT_MATCH_SCORE: i32 = 2;

/// Score of rust-bio's banded aligner, `banded::Aligner::with_scoring(scoring, k, w)`, for
/// `custom_with_matches(x, y, matches)`
#[allow(non_snake_case)]
pub struct BandedScoreAligner<F: MatchFunc> {
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    Sn: Vec<i32>,
    scoring: Scoring<F>,
    band: Band,
    k: usize,
    w: usize,
}

impl<F: MatchFunc> BandedScoreAligner<F> {
    pub fn with_scoring(scoring: Scoring<F>, k: usize, w: usize) -> Self {
        Self {
            S: [Vec::new(), Vec::new()],
            I: [Vec::new(), Vec::new()],
            D: [Vec::new(), Vec::new()],
            Sn: Vec::new(),
            scoring,
            band: Band::default(),
            k,
            w,
        }
    }

    pub fn custom_with_matches(
        &mut self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        matches: &[(u32, u32)],
    ) -> i32 {
        self.band = Band::create_with_matches(x, y, self.k, self.w, &self.scoring, matches);
        self.compute_score(x, y)
    }

    fn compute_score(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> i32 {
        if self.band.num_cells() > MAX_CELLS {
            return MIN_SCORE;
        }

        let (m, n) = (x.len(), y.len());

        for k in 0..2 {
            self.I[k].clear();
            self.D[k].clear();
            self.S[k].clear();
            self.D[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.I[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.S[k].extend(repeat(MIN_SCORE).take(m + 1));
        }
        self.Sn.clear();
        self.Sn.extend(repeat(MIN_SCORE).take(m + 1));

        {
            // Handle j = 0
            let curr = 0;
            let i_start = self.band.ranges[0].start;
            let i_end = self.band.ranges[0].end;
            if i_start == 0 {
                self.S[curr][0] = 0;
            }

            for i in max(1, i_start)..i_end {
                if i == 1 {
                    self.I[curr][i] = self.scoring.gap_open + self.scoring.gap_extend;
                } else {
                    // Insert all i characters
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score =
                        self.scoring.xclip_prefix + self.scoring.gap_open + self.scoring.gap_extend; // Clip then insert
                    if i_score > c_score {
                        self.I[curr][i] = i_score;
                    } else {
                        self.I[curr][i] = c_score;
                    }
                }

                if self.I[curr][i] > self.S[curr][i] {
                    self.S[curr][i] = self.I[curr][i];
                }

                if self.scoring.xclip_prefix > self.S[curr][i] {
                    self.S[curr][i] = self.scoring.xclip_prefix;
                }

                // Track the score if we do a suffix clip (x) after this character
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                }
            }

            for i in i_end..min(m + 1, self.band.ranges[min(n, 1)].end) {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
            }

            if i_end < (m + 1) {
                self.S[curr][m] = MIN_SCORE;
            }
            // Track the score if we do clip (y) from origin
            if self.scoring.yclip_prefix > self.scoring.yclip_suffix {
                self.Sn[0] = self.scoring.yclip_prefix;
            } else {
                self.Sn[0] = self.scoring.yclip_suffix;
            }
        }

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            let i_start = self.band.ranges[j].start;
            let i_end = self.band.ranges[j].end;

            if i_start == 0 {
                // Handle i = 0
                self.I[curr][0] = MIN_SCORE;

                if j == 1 {
                    self.D[curr][0] = self.scoring.gap_open + self.scoring.gap_extend;
                } else {
                    // Delete all j characters
                    let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
                    let c_score =
                        self.scoring.yclip_prefix + self.scoring.gap_open + self.scoring.gap_extend;
                    if d_score > c_score {
                        self.D[curr][0] = d_score;
                    } else {
                        self.D[curr][0] = c_score;
                    }
                }

                if self.D[curr][0] > self.scoring.yclip_prefix {
                    self.S[curr][0] = self.D[curr][0];
                } else {
                    self.S[curr][0] = self.scoring.yclip_prefix;
                }

                // Track the score if we do suffix clip (y) from here
                if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
                    self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
                }
            }

            for i in i_start.saturating_sub(1)..i_start {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
                self.D[curr][i] = MIN_SCORE;
            }
            self.S[curr][m] = MIN_SCORE;

            let q = y[j - 1];
            let xclip_score = self.scoring.xclip_prefix
                + max(
                    if j == n {
                        max(self.scoring.yclip_prefix, self.Sn[0])
                    } else {
                        self.scoring.yclip_prefix
                    },
                    self.scoring.gap_open + self.scoring.gap_extend * (j as i32),
                );

            for i in max(1, i_start)..i_end {
                let p = x[i - 1];

                let m_score = self.S[prev][i - 1] + self.scoring.match_fn.score(p, q);

                let i_score = self.I[curr][i - 1] + self.scoring.gap_extend;
                let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
                let mut best_i_score = if i_score > s_score { i_score } else { s_score };
                if j == n {
                    let clip_score =
                        self.Sn[i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
                    if clip_score > best_i_score {
                        best_i_score = clip_score;
                    }
                }

                let d_score = self.D[prev][i] + self.scoring.gap_extend;
                let s_score = self.S[prev][i] + self.scoring.gap_open + self.scoring.gap_extend;
                let best_d_score = if d_score > s_score { d_score } else { s_score };

                if i != m {
                    self.S[curr][i] = MIN_SCORE;
                }
                let mut best_s_score = self.S[curr][i];

                if m_score > best_s_score {
                    best_s_score = m_score;
                }

                if best_i_score > best_s_score {
                    best_s_score = best_i_score;
                }

                if best_d_score > best_s_score {
                    best_s_score = best_d_score;
                }

                if xclip_score > best_s_score {
                    best_s_score = xclip_score;
                }

                let yclip_score = self.scoring.yclip_prefix
                    + self.scoring.gap_open
                    + self.scoring.gap_extend * (i as i32);
                if yclip_score > best_s_score {
                    best_s_score = yclip_score;
                }

                self.S[curr][i] = best_s_score;
                self.I[curr][i] = best_i_score;
                self.D[curr][i] = best_d_score;

                // Track the score if we do suffix clip (x) from here
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                }

                // Track the score if we do suffix clip (y) from here
                if self.S[curr][i] + self.scoring.yclip_suffix > self.Sn[i] {
                    self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
                }
            }

            // Suffix clip (y) from i = m and reset Sn[m] if required
            if self.S[curr][m] + self.scoring.yclip_suffix > self.Sn[m] {
                self.Sn[m] = self.S[curr][m] + self.scoring.yclip_suffix;
            }
            if i_end < (m + 1) {
                self.S[curr][m] = MIN_SCORE;
            }

            for i in i_end..min(m + 1, self.band.ranges[min(n, j + 1)].end) {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
                self.D[curr][i] = MIN_SCORE;
            }
        }

        // Handle suffix clipping in the j=n case
        for i in 0..=m {
            let j = n;
            let curr = j % 2;
            // These entries are not set in the loop above and could contain leftover
            // values from previous columns. Reset them to MIN_SCORE
            if i != m && (i < self.band.ranges[j].start || i > self.band.ranges[j].end) {
                self.S[curr][i] = MIN_SCORE;
            }
            if self.Sn[i] > self.S[curr][i] {
                self.S[curr][i] = self.Sn[i];
            }
            if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last column of I as this could also change
        for i in max(1, self.band.ranges[n].start)..self.band.ranges[n].end {
            let j = n;
            let curr = j % 2;
            let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if s_score > self.I[curr][i] {
                self.I[curr][i] = s_score;
            }
            if s_score > self.S[curr][i] {
                self.S[curr][i] = s_score;
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                }
            }
        }

        if n > 0 {
            let j = n;
            let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
            let mut best_score = max(d_score, self.scoring.yclip_prefix);
            if self.scoring.yclip_suffix > best_score {
                best_score = self.scoring.yclip_suffix;
            }
            if (self.scoring.xclip_suffix + best_score) > self.S[n % 2][m] {
                self.S[n % 2][m] = self.scoring.xclip_suffix + best_score;
            }
        }

        if m > 0 {
            let i = m;
            let c_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
            let mut best_score = max(c_score, self.scoring.xclip_prefix);
            if self.scoring.xclip_suffix > best_score {
                best_score = self.scoring.xclip_suffix;
            }
            if (self.scoring.yclip_suffix + best_score) > self.S[n % 2][m] {
                self.S[n % 2][m] = self.scoring.yclip_suffix + best_score;
            }
        }

        self.S[n % 2][m]
    }
}

trait MatchPair {
    fn continues(&self, p: Option<(u32, u32)>) -> bool;
}

impl MatchPair for (u32, u32) {
    fn continues(&self, p: Option<(u32, u32)>) -> bool {
        match p {
            Some(_p) => self.0 == _p.0 + 1 && self.1 == _p.1 + 1,
            None => false,
        }
    }
}

#[derive(Default, Clone, Debug)]
struct Band {
    rows: usize,
    cols: usize,
    ranges: Vec<Range<usize>>,
}

impl Band {
    // Create new Band instance with given size
    //
    // # Arguments
    //
    // * `m` - the expected size of x
    // * `n` - the expected size of y
    //
    fn new(m: usize, n: usize) -> Self {
        Band {
            rows: m + 1,
            cols: n + 1,
            ranges: vec![m + 1..0; n + 1],
        }
    }

    // Add cells around a kmer of length 'k', starting at 'start', which are within a
    // distance of 'w' in x or y directions to the band.
    fn add_kmer(&mut self, start: (u32, u32), k: usize, w: usize) {
        let (r, c) = (start.0 as usize, start.1 as usize);
        // println!("{} {} {}", r, k, self.rows);
        debug_assert!(r + k <= self.rows);
        debug_assert!(c + k <= self.cols);

        if k == 0 {
            return;
        }

        let i = r.saturating_sub(w);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
        }

        let mut i = r.saturating_sub(w);
        for j in min(c + w, self.cols)..min(c + k + w, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
            i += 1;
        }

        let mut i = r + w + k;
        let mut j = (c + k - 1).saturating_sub(w);
        loop {
            if j <= c.saturating_sub(w) {
                break;
            }
            j -= 1;
            i -= 1;
            self.ranges[j].end = max(self.ranges[j].end, min(i, self.rows));
        }

        let i = min(r + w + k, self.rows);
        for j in (c + k - 1).saturating_sub(w)..min(c + k + w, self.cols) {
            self.ranges[j].end = max(self.ranges[j].end, i);
        }
    }

    // Add cells around a specific position to the band. An cell which is within 'w' distance
    // in x or y directions are added
    fn add_entry(&mut self, pos: (u32, u32), w: usize) {
        let (r, c) = (pos.0 as usize, pos.1 as usize);

        let istart = r.saturating_sub(w);
        let iend = min(r + w + 1, self.rows);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, istart);
            self.ranges[j].end = max(self.ranges[j].end, iend);
        }
    }

    // Each gap generates a line from the start to end.
    fn add_gap(&mut self, start: (u32, u32), end: (u32, u32), w: usize) {
        let nrows = end.0 - start.0;
        let ncols = end.1 - start.1;
        if nrows > ncols {
            for r in start.0..end.0 {
                let c = start.1 + (end.1 - start.1) * (r - start.0) / (end.0 - start.0);
                self.add_entry((r, c), w);
            }
        } else {
            for c in start.1..end.1 {
                let r = start.0 + (end.0 - start.0) * (c - start.1) / (end.1 - start.1);
                self.add_entry((r, c), w);
            }
        }
    }

    // The band needs to start either at (0,0) or at a point that is zero score from (0,0).
    // This naturally sets the start positions correctly for global, semiglobal and local
    // modes. Similarly the band has to either end at (m,n) or at a point from which there is
    // a zero score path to (m,n).
    //
    // At the minimum, irrespective of the score (0,0)->start or end->(m,n), we extend the band
    // diagonally for a length "lazy_extend"(2k) or when it hits the corner, whichever happens first
    //
    // start - the index of the first matching kmer in LCSk++
    // end - the index of the last matching kmer in LCSk++
    //
    fn set_boundaries<F: MatchFunc>(
        &mut self,
        start: (u32, u32),
        end: (u32, u32),
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
    ) {
        let lazy_extend: usize = 2 * k;

        // -------------- START --------------
        // Nothing to do if the start is already at (0,0)
        let (r, c) = (start.0 as usize, start.1 as usize);
        if !(r == 0usize && c == 0usize) {
            let mut score_to_start = if r > 0 { scoring.xclip_prefix } else { 0i32 };
            score_to_start += if c > 0 { scoring.yclip_prefix } else { 0i32 };

            if score_to_start == 0 {
                // Just do a "lazy_extend"
                // First diagonally
                let d = min(lazy_extend, min(r, c));
                self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);

                // If we hit one of the edges before completing lazy_extend
                self.add_gap(
                    (
                        r.saturating_sub(lazy_extend) as u32,
                        c.saturating_sub(lazy_extend) as u32,
                    ),
                    ((r - d) as u32, (c - d) as u32),
                    w,
                );
            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let diagonal_score = match r.cmp(&c) {
                    // We will hit (r-c, 0)
                    Ordering::Greater => scoring.xclip_prefix,
                    // We will hit (0, c-r)
                    Ordering::Less => scoring.yclip_prefix,
                    Ordering::Equal => 0,
                };

                if diagonal_score == 0 {
                    let d = min(r, c);
                    self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let start = (
                        r.saturating_sub(lazy_extend) as u32,
                        c.saturating_sub(lazy_extend) as u32,
                    );
                    let end = ((r - d) as u32, (c - d) as u32);
                    if (start.0 <= end.0) && (start.1 <= end.1) {
                        self.add_gap(start, end, w);
                    }
                } else {
                    // Band to origin
                    self.add_gap((0u32, 0u32), start, w);
                }
            }
        }

        // -------------- END --------------
        // Nothing to do if the last kmer ends at (m, n)
        let (r, c) = (end.0 as usize + k, end.1 as usize + k);
        debug_assert!(r <= self.rows);
        debug_assert!(c <= self.cols);
        if !(r == self.rows && c == self.cols) {
            let mut score_from_end = if r == self.rows {
                0
            } else {
                scoring.xclip_suffix
            };
            score_from_end += if c == self.cols {
                0
            } else {
                scoring.yclip_suffix
            };

            if score_from_end == 0 {
                // Just a lazy_extend
                let d = min(lazy_extend, min(self.rows - r, self.cols - c));
                self.add_kmer((r as u32, c as u32), d, w);

                let r1 = min(self.rows, r + d) - 1;
                let c1 = min(self.cols, c + d) - 1;
                let r2 = min(self.rows, r + lazy_extend);
                let c2 = min(self.cols, c + lazy_extend);
                if (r1 <= r2) && (c1 <= c2) {
                    self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                }
            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let dr = self.rows - r;
                let dc = self.cols - c;
                let diagonal_score = match dr.cmp(&dc) {
                    // We will hit (r+dc, self.cols)
                    Ordering::Greater => scoring.xclip_suffix,
                    // We will hit (self.rows, c+dr)
                    Ordering::Less => scoring.yclip_suffix,
                    // We will hit the corner
                    Ordering::Equal => 0,
                };

                if diagonal_score == 0 {
                    let d = min(dr, dc);
                    self.add_kmer((r as u32, c as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let r1 = min(self.rows, r + d) - 1;
                    let c1 = min(self.cols, c + d) - 1;
                    let r2 = min(self.rows, r + lazy_extend);
                    let c2 = min(self.cols, c + lazy_extend);
                    if (r1 <= r2) && (c1 <= c2) {
                        self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                    }
                } else {
                    // Band to lower right corner
                    let rows = self.rows as u32;
                    let cols = self.cols as u32;
                    self.add_gap((r as u32, c as u32), (rows, cols), w);
                }
            }
        }
    }

    fn create_with_matches<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
        matches: &[(u32, u32)],
    ) -> Band {
        if matches.is_empty() {
            let mut band = Band::new(x.len(), y.len());
            band.full_matrix();
            return band;
        }

        let match_score = match scoring.match_scores {
            Some((m, _)) => m,
            None => DEFAULT_MATCH_SCORE,
        };

        let res = sparse::sdpkpp(
            matches,
            k,
            match_score as u32,
            scoring.gap_open,
            scoring.gap_extend,
        );
        Band::create_from_match_path(x, y, k, w, scoring, &res.path, matches)
    }

    fn create_from_match_path<F: MatchFunc>(
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        k: usize,
        w: usize,
        scoring: &Scoring<F>,
        path: &[usize],
        matches: &[(u32, u32)],
    ) -> Band {
        let mut band = Band::new(x.len(), y.len());

        if matches.is_empty() {
            band.full_matrix();
            return band;
        }

        let ps = path[0];
        let pe = path[path.len() - 1];

        // Set the boundaries
        band.set_boundaries(matches[ps], matches[pe], k, w, scoring);
        let mut prev: Option<(u32, u32)> = None;

        for &idx in path {
            let curr = matches[idx];
            if curr.continues(prev) {
                let p = prev.unwrap();
                band.add_entry((p.0 + k as u32, p.1 + k as u32), w);
            } else {
                if let Some(p) = prev {
                    band.add_gap((p.0 + (k - 1) as u32, p.1 + (k - 1) as u32), curr, w)
                }
                band.add_kmer(curr, k, w);
            }
            prev = Some(curr);
        }
        band
    }

    fn full_matrix(&mut self) {
        self.ranges.clear();
        self.ranges.resize(self.cols, 0..self.rows);
    }

    fn num_cells(&self) -> usize {
        let mut banded_cells = 0;
        for j in 0..self.ranges.len() {
            banded_cells += self.ranges[j].end.saturating_sub(self.ranges[j].start);
        }
        banded_cells
    }
}
