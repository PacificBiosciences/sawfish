//! Exact low-memory reimplementation of WFA2 gap-affine-2p ends-free alignment with a match score
//!
//! This computes exactly the same score and CIGAR as WFA2-lib for one aligner configuration: gap-affine
//! 2-piece penalties with a negative match score, ends-free alignment, no heuristic, and either the
//! MemoryHigh or MemoryMed memory model (these two backtrace the same wavefronts in different ways).
//!
//! WFA2 keeps every wavefront for its backtrace (MemoryHigh), or a backtrace buffer with every
//! wavefront in a ring of the last max_score_scope scores (MemoryMed). For long gaps this is a large
//! amount of memory, written once and read back from DRAM, which limits how many alignments can run
//! on each CPU at the same time. Here, only the input wavefronts needed for the next score are kept,
//! and each cell of each score stores just the backtrace choice it would make and its number of
//! matches, so that the working set mostly stays in cache.
//!
//! The wavefront computation follows WFA2 step by step, including its rules for wavefront limits,
//! trimming of out-of-bounds cells at the wavefront ends, out-of-bounds M cells, ends-free
//! initialization at each score, and the order of the termination check, since these all affect which
//! alignment is returned among those with the optimal score.
//!

use std::collections::HashMap;

// Backtrace types, with the same values as WFA2's backtrace_type, so that ties between candidate
// predecessors are broken in the same way
const BT_NONE: u8 = 0;
const BT_I1_OPEN: u8 = 1;
const BT_I1_EXT: u8 = 2;
const BT_I2_OPEN: u8 = 3;
const BT_I2_EXT: u8 = 4;
const BT_D1_OPEN: u8 = 5;
const BT_D1_EXT: u8 = 6;
const BT_D2_OPEN: u8 = 7;
const BT_D2_EXT: u8 = 8;
const BT_M: u8 = 9;

// Flags for the gap components in the cell backtrace byte, set when the gap is extended rather
// than opened
const I1_EXT_FLAG: u8 = 1 << 4;
const I2_EXT_FLAG: u8 = 1 << 5;
const D1_EXT_FLAG: u8 = 1 << 6;
const D2_EXT_FLAG: u8 = 1 << 7;

/// Wavefront offset type
///
/// WFA2 uses i32 offsets, with WAVEFRONT_OFFSET_NULL = i32::MIN/2. Offsets are only ever compared
/// or incremented by one per score, so i16 offsets give exactly the same results as long as no
/// value can overflow. That holds when the sequences are short enough and the score stays below
/// a limit, see I16_MAX_SEQ_LEN and i16_max_score. Using i16 halves the size of the wavefronts,
/// which are mostly limited by memory bandwidth.
///
trait Offset: Copy + Ord + std::fmt::Debug {
    const NULL: Self;
    fn from_i32(x: i32) -> Self;
    fn to_i32(self) -> i32;
    /// True if offset on diagonal k is out of bounds (h > tlen or v > plen), with WFA2's unsigned
    /// comparison, which also treats NULL-like (very negative) offsets as out of bounds
    fn is_oob(self, k: i32, plen: i32, tlen: i32) -> bool;
}

impl Offset for i32 {
    const NULL: Self = i32::MIN / 2;
    #[inline(always)]
    fn from_i32(x: i32) -> Self {
        x
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        self
    }
    #[inline(always)]
    fn is_oob(self, k: i32, plen: i32, tlen: i32) -> bool {
        (self as u32) > tlen as u32 || (self.wrapping_sub(k) as u32) > plen as u32
    }
}

impl Offset for i16 {
    const NULL: Self = i16::MIN / 2;
    #[inline(always)]
    fn from_i32(x: i32) -> Self {
        x as i16
    }
    #[inline(always)]
    fn to_i32(self) -> i32 {
        self as i32
    }
    #[inline(always)]
    fn is_oob(self, k: i32, plen: i32, tlen: i32) -> bool {
        // Exact as long as the true h and v are in (-65536 + max(plen,tlen), 65536), see Offset
        (self as u16) > tlen as u16 || (self.wrapping_sub(k as i16) as u16) > plen as u16
    }
}

/// Check if sequence lengths allow i16 offsets
///
/// Offsets are text positions (h), so valid offsets are at most text_length. Out-of-bounds
/// insertion offsets, and NULL-like offsets, grow by at most one per score, so with
/// i16_max_score() they stay below i16::MAX, and below zero, respectively. Diagonals k (from
/// -pattern_length-3 to text_length+3) must fit in i16, and for the unsigned bounds check, the
/// true value of v = h - k must stay in (-65536 + plen, 65536).
fn is_i16_compatible(plen: usize, tlen: usize) -> bool {
    plen <= 32_000 && tlen <= 32_000 && plen + tlen <= 48_000
}

/// Maximum score for i16 offsets, see is_i16_compatible
fn i16_max_score(tlen: i32) -> i32 {
    16_000.min(32_600 - tlen)
}

/// Penalties after WFA2's adjustment for a match score (Eizenga's formula)
#[derive(Clone, Copy)]
struct Penalties {
    match_: i32,
    mismatch: i32,
    gap_open1: i32,
    gap_ext1: i32,
    gap_open2: i32,
    gap_ext2: i32,
}

/// One wavefront component at one score
struct Wavefront<T> {
    /// False if the wavefront isn't allocated (WFA2 NULL wavefront pointer)
    allocated: bool,
    lo: i32,
    hi: i32,
    offsets: Vec<T>,
}

impl<T: Offset> Wavefront<T> {
    fn new(size: usize) -> Self {
        Self {
            allocated: false,
            lo: 1,
            hi: -1,
            offsets: vec![T::NULL; size],
        }
    }

    /// True for WFA2's null input wavefront: not allocated or trimmed to nothing
    fn is_null(&self) -> bool {
        !self.allocated || self.lo > self.hi
    }
}

/// Wavefront rings for one offset type
struct Rings<T> {
    m: Vec<Wavefront<T>>,
    i1: Vec<Wavefront<T>>,
    d1: Vec<Wavefront<T>>,
    i2: Vec<Wavefront<T>>,
    d2: Vec<Wavefront<T>>,
    null_wf: Wavefront<T>,
    /// Candidate (pre-extension, pre-clamp) M offsets of the most recent score, over
    /// cur_lo..=cur_hi
    pre: Vec<T>,
    /// Index of k=0 in each offsets array
    k_base: i32,
    /// Checkpoints of the ring state, in score order, see Affine2pAligner::checkpoint
    snaps: Vec<Snapshot<T>>,
    /// Budget bytes held for snaps
    snap_bytes: u64,
}

/// Saved state of one wavefront, with offsets only over lo..=hi, since no other offsets are read
/// before they are written (see pad in compute_score)
struct SnapWavefront<T> {
    allocated: bool,
    lo: i32,
    hi: i32,
    offsets: Vec<T>,
}

/// Saved state of all rings after the extension of one score
struct Snapshot<T> {
    score: i32,
    /// Wavefronts of rings m, i1, d1, i2, d2, in order
    wavefronts: Vec<SnapWavefront<T>>,
}

impl<T: Offset> Rings<T> {
    fn new() -> Self {
        Self {
            m: Vec::new(),
            i1: Vec::new(),
            d1: Vec::new(),
            i2: Vec::new(),
            d2: Vec::new(),
            null_wf: Wavefront::new(0),
            pre: Vec::new(),
            k_base: 0,
            snaps: Vec::new(),
            snap_bytes: 0,
        }
    }

    fn clear_snaps(&mut self) {
        self.snaps.clear();
        if self.snap_bytes > 0 {
            unsafe { mm_allocator_budget_release(self.snap_bytes) };
            self.snap_bytes = 0;
        }
    }

    /// Save the state of all rings, as of the given score
    fn save_snapshot(&mut self, score: i32) {
        let kb = self.k_base;
        let wavefronts = [&self.m, &self.i1, &self.d1, &self.i2, &self.d2]
            .into_iter()
            .flatten()
            .map(|w| SnapWavefront {
                allocated: w.allocated,
                lo: w.lo,
                hi: w.hi,
                offsets: if w.allocated && w.lo <= w.hi {
                    w.offsets[(w.lo + kb) as usize..=(w.hi + kb) as usize].to_vec()
                } else {
                    Vec::new()
                },
            })
            .collect::<Vec<_>>();
        let bytes = wavefronts
            .iter()
            .map(|x| (x.offsets.len() * std::mem::size_of::<T>()) as u64)
            .sum::<u64>();
        unsafe { mm_allocator_budget_acquire(bytes) };
        self.snap_bytes += bytes;
        self.snaps.push(Snapshot { score, wavefronts });
    }

    /// Restore the state of all rings from a snapshot
    fn restore_snapshot(&mut self, snap: &Snapshot<T>) {
        let kb = self.k_base;
        let wavefronts = [
            &mut self.m,
            &mut self.i1,
            &mut self.d1,
            &mut self.i2,
            &mut self.d2,
        ]
        .into_iter()
        .flatten();
        for (w, x) in wavefronts.zip(snap.wavefronts.iter()) {
            w.allocated = x.allocated;
            w.lo = x.lo;
            w.hi = x.hi;
            if !x.offsets.is_empty() {
                w.offsets[(x.lo + kb) as usize..=(x.hi + kb) as usize].copy_from_slice(&x.offsets);
            }
        }
    }

    fn setup(&mut self, pen: &Penalties, plen: usize, tlen: usize) {
        self.clear_snaps();
        let size = plen + tlen + 32;
        self.k_base = plen as i32 + 12;
        let m_n = [
            pen.mismatch,
            pen.gap_open1 + pen.gap_ext1,
            pen.gap_open2 + pen.gap_ext2,
        ]
        .into_iter()
        .max()
        .unwrap() as usize
            + 1;
        let i1_n = pen.gap_ext1 as usize + 1;
        let i2_n = pen.gap_ext2 as usize + 1;
        for (ring, n) in [
            (&mut self.m, m_n),
            (&mut self.i1, i1_n),
            (&mut self.d1, i1_n),
            (&mut self.i2, i2_n),
            (&mut self.d2, i2_n),
        ] {
            ring.truncate(n);
            ring.resize_with(n, || Wavefront::new(0));
            for w in ring.iter_mut() {
                w.allocated = false;
                w.lo = 1;
                w.hi = -1;
                if w.offsets.len() < size {
                    w.offsets.resize(size, T::NULL);
                }
            }
        }
        if self.null_wf.offsets.len() < size {
            self.null_wf = Wavefront::new(size);
        }
    }

    /// Release memory if the last alignment was much larger than typical
    fn shrink(&mut self, max_size: usize) {
        for ring in [
            &mut self.m,
            &mut self.i1,
            &mut self.d1,
            &mut self.i2,
            &mut self.d2,
        ] {
            if ring.first().is_some_and(|w| w.offsets.len() > max_size) {
                ring.clear();
            }
        }
        if self.null_wf.offsets.len() > max_size {
            self.null_wf = Wavefront::new(0);
        }
        if self.pre.len() > max_size {
            self.pre = Vec::new();
        }
    }
}

/// Backtrace data for one score
struct ScoreTrace {
    /// Lowest diagonal with trace data (all cells in lo..=hi have trace data)
    lo: i32,
    hi: i32,
    /// Location of this score's cells in the trace arena
    chunk: u32,
    start: u32,
    /// Diagonals of M cells set by ends-free initialization at this score
    injected: [i32; 2],
}

unsafe extern "C" {
    fn mm_allocator_budget_acquire(num_bytes: u64);
    fn mm_allocator_budget_release(num_bytes: u64);
}

/// Size of each trace arena chunk, in cells
const CHUNK_CELLS: usize = 1 << 20;

/// Default trace memory at which a checkpoint is taken
const DEFAULT_CHECKPOINT_BYTES: u64 = 128 << 20;

/// Inputs needed to recompute trace rows of a checkpointed alignment
#[derive(Clone, Copy)]
struct CheckpointState {
    is_i32: bool,
    plen: i32,
    tlen: i32,
    form: EndsFree,
    /// Trace rows are available for scores above this
    loaded_above: i32,
}

/// Backtrace data for all cells, in chunks so that memory grows without copying
///
/// Each cell has a choice byte and, if with_matches is set, a matches byte. Chunk memory is taken
/// from the WFA2 memory budget shared by all aligners.
#[derive(Default)]
struct TraceArena {
    /// Store matches (only needed for the MemoryHigh backtrace)
    with_matches: bool,
    choice: Vec<Box<[u8]>>,
    matches: Vec<Box<[u8]>>,
    /// Cells used in the last chunk
    used: usize,
    /// Trace bytes used in all chunks
    used_bytes: u64,
    /// Budget bytes held
    held: u64,
    /// Standard size chunks kept from earlier alignments, to avoid page faults on fresh memory
    pool_choice: Vec<Box<[u8]>>,
    pool_matches: Vec<Box<[u8]>>,
}

/// Bytes of trace chunks kept for reuse by each aligner
const POOL_BYTES: usize = 16 << 20;

impl TraceArena {
    /// Allocate n contiguous cells, returning (chunk, start)
    fn alloc(&mut self, n: usize) -> (u32, u32) {
        let last_len = self.choice.last().map(|x| x.len()).unwrap_or(0);
        if self.choice.is_empty() || self.used + n > last_len {
            let size = CHUNK_CELLS.max(n);
            let bytes = (1 + self.with_matches as u64) * size as u64;
            unsafe { mm_allocator_budget_acquire(bytes) };
            self.held += bytes;
            let reuse = size == CHUNK_CELLS;
            match reuse.then(|| self.pool_choice.pop()).flatten() {
                Some(c) => self.choice.push(c),
                None => self.choice.push(vec![0u8; size].into_boxed_slice()),
            }
            if self.with_matches {
                match reuse.then(|| self.pool_matches.pop()).flatten() {
                    Some(mut m) => {
                        m.fill(0);
                        self.matches.push(m);
                    }
                    None => self.matches.push(vec![0u8; size].into_boxed_slice()),
                }
            }
            self.used = 0;
        }
        let start = self.used;
        self.used += n;
        self.used_bytes += (1 + self.with_matches as u64) * n as u64;
        // Clear choices, since cells added by ends-free initialization aren't computed
        let c = self.choice.last_mut().unwrap();
        c[start..start + n].fill(BT_NONE);
        ((self.choice.len() - 1) as u32, start as u32)
    }

    fn clear(&mut self) {
        for c in self.choice.drain(..) {
            if c.len() == CHUNK_CELLS
                && (self.pool_choice.len() + 1) * CHUNK_CELLS <= POOL_BYTES / 2
            {
                self.pool_choice.push(c);
            }
        }
        for m in self.matches.drain(..) {
            if m.len() == CHUNK_CELLS
                && (self.pool_matches.len() + 1) * CHUNK_CELLS <= POOL_BYTES / 2
            {
                self.pool_matches.push(m);
            }
        }
        self.used = 0;
        self.used_bytes = 0;
        if self.held > 0 {
            unsafe { mm_allocator_budget_release(self.held) };
            self.held = 0;
        }
    }
}

impl Drop for TraceArena {
    fn drop(&mut self) {
        self.clear();
    }
}

/// Result of an alignment
pub struct Affine2pAlignment {
    pub score: i32,
    pub cigar: Vec<u8>,
}

/// Low-memory WFA2-equivalent aligner
pub struct Affine2pAligner {
    pen: Penalties,
    piggyback: bool,
    rings16: Rings<i16>,
    rings32: Rings<i32>,
    // Backtrace data
    traces: Vec<ScoreTrace>,
    arena: TraceArena,
    big_matches: HashMap<(i32, i32), i32>,
    /// Trace memory at which a checkpoint is taken, see checkpoint
    checkpoint_bytes: u64,
    /// Checkpoint state of the current alignment
    checkpoint: Option<CheckpointState>,
    /// Range of the most recent computed score
    cur_lo: i32,
    cur_hi: i32,
    // Sequences with distinct end sentinels and padding, for 8-byte block comparison
    pattern_buf: Vec<u8>,
    text_buf: Vec<u8>,
    /// Not Send: budget memory is accounted to the thread which allocates it
    _not_send: std::marker::PhantomData<*const ()>,
}

/// Index of score s in a ring of size n
#[inline(always)]
fn ring(s: i32, n: usize) -> usize {
    (s as usize) % n
}

/// Ends-free alignment form
#[derive(Clone, Copy)]
struct EndsFree {
    pattern_begin: i32,
    pattern_end: i32,
    text_begin: i32,
    text_end: i32,
}

impl Affine2pAligner {
    /// Create aligner with the same parameters as WFAlignerGapAffine2Pieces::new_with_match
    ///
    /// piggyback: true to match MemoryMed (or lower) backtrace, false to match MemoryHigh
    pub fn new(
        match_: i32,
        mismatch: i32,
        gap_opening1: i32,
        gap_extension1: i32,
        gap_opening2: i32,
        gap_extension2: i32,
        piggyback: bool,
    ) -> Self {
        assert!(match_ < 0, "Only negative match scores are supported");
        let pen = Penalties {
            match_,
            mismatch: 2 * mismatch - 2 * match_,
            gap_open1: 2 * gap_opening1,
            gap_ext1: 2 * gap_extension1 - match_,
            gap_open2: 2 * gap_opening2,
            gap_ext2: 2 * gap_extension2 - match_,
        };
        Self {
            pen,
            piggyback,
            rings16: Rings::new(),
            rings32: Rings::new(),
            traces: Vec::new(),
            arena: TraceArena::default(),
            big_matches: HashMap::new(),
            checkpoint_bytes: DEFAULT_CHECKPOINT_BYTES,
            checkpoint: None,
            cur_lo: 1,
            cur_hi: 0,
            pattern_buf: Vec::new(),
            text_buf: Vec::new(),
            _not_send: std::marker::PhantomData,
        }
    }

    /// Set the trace memory used at which a checkpoint is taken
    ///
    /// This doesn't change alignment results, only the memory and time used for long alignments.
    #[doc(hidden)]
    pub fn set_checkpoint_bytes(&mut self, checkpoint_bytes: u64) {
        self.checkpoint_bytes = checkpoint_bytes;
    }

    fn clear_traces(&mut self) {
        self.traces.clear();
        self.arena.clear();
        self.big_matches.clear();
        self.cur_lo = 1;
        self.cur_hi = 0;
        self.checkpoint = None;
        self.rings16.clear_snaps();
        self.rings32.clear_snaps();
    }

    /// Ends-free alignment of pattern to text, matching WFAligner::align_ends_free
    pub fn align_ends_free(
        &mut self,
        pattern: &[u8],
        text: &[u8],
        pattern_begin_free: i32,
        pattern_end_free: i32,
        text_begin_free: i32,
        text_end_free: i32,
    ) -> Affine2pAlignment {
        assert!(
            (0..=pattern.len() as i32).contains(&pattern_begin_free)
                && (0..=pattern.len() as i32).contains(&pattern_end_free)
                && (0..=text.len() as i32).contains(&text_begin_free)
                && (0..=text.len() as i32).contains(&text_end_free),
            "ends-free lengths must be between 0 and the sequence length"
        );
        debug_assert!(
            !pattern.contains(&b'!') && !text.contains(&b'?'),
            "sequences must not contain the padding sentinels '!' (pattern) and '?' (text)"
        );
        let plen = pattern.len() as i32;
        let tlen = text.len() as i32;
        let form = EndsFree {
            pattern_begin: pattern_begin_free,
            pattern_end: pattern_end_free,
            text_begin: text_begin_free,
            text_end: text_end_free,
        };
        self.pattern_buf.clear();
        self.pattern_buf.extend_from_slice(pattern);
        self.pattern_buf.resize(pattern.len() + 16, b'!');
        self.text_buf.clear();
        self.text_buf.extend_from_slice(text);
        self.text_buf.resize(text.len() + 16, b'?');
        self.arena.with_matches = !self.piggyback;

        // Compute wavefronts with i16 offsets if possible, else i32
        let mut result = None;
        if is_i16_compatible(pattern.len(), text.len()) {
            let mut rings = std::mem::replace(&mut self.rings16, Rings::new());
            rings.setup(&self.pen, pattern.len(), text.len());
            result = self.compute_wavefronts(&mut rings, plen, tlen, form, i16_max_score(tlen));
            self.rings16 = rings;
        }
        let (score, (end_k, end_offset)) = match result {
            Some(x) => x,
            None => {
                self.clear_traces();
                let mut rings = std::mem::replace(&mut self.rings32, Rings::new());
                rings.setup(&self.pen, pattern.len(), text.len());
                let x = self.compute_wavefronts(&mut rings, plen, tlen, form, i32::MAX);
                self.rings32 = rings;
                x.unwrap()
            }
        };

        let cigar = if self.piggyback {
            self.backtrace_piggyback(pattern, text, score, end_k, end_offset)
        } else {
            self.backtrace_high(plen, tlen, score, end_k, end_offset)
        };
        self.clear_traces();
        // Don't keep memory from unusually long alignments
        let max_kept_size = 100_000;
        self.rings16.shrink(max_kept_size);
        self.rings32.shrink(max_kept_size);
        if self.traces.capacity() > max_kept_size {
            self.traces = Vec::new();
        }

        let end_v = end_offset - end_k;
        let end_h = end_offset;
        let swg_match = -self.pen.match_;
        Affine2pAlignment {
            score: (swg_match * (end_v + end_h) - score) / 2,
            cigar,
        }
    }

    /// Set the backtrace data for score s, which is either the next score or a recomputed score
    fn set_trace(&mut self, s: i32, trace: ScoreTrace) {
        let s = s as usize;
        if s == self.traces.len() {
            self.traces.push(trace);
        } else {
            self.traces[s] = trace;
        }
    }

    /// Compute wavefronts until the end of the alignment is reached, storing backtrace data
    ///
    /// Returns the alignment score and end position (k, offset), or None if max_score is exceeded
    ///
    /// Backtrace data takes memory proportional to the number of cells over all scores, which can
    /// be several GB for long gaps. When the trace memory reaches checkpoint_bytes, the ring state
    /// is saved and the trace memory is released, so that the trace rows of each segment between
    /// checkpoints can be recomputed during the backtrace (see ensure_trace_row). Recomputation
    /// repeats exactly the same steps, so this has no effect on results.
    ///
    fn compute_wavefronts<T: ComputeKernel>(
        &mut self,
        rings: &mut Rings<T>,
        plen: i32,
        tlen: i32,
        form: EndsFree,
        max_score: i32,
    ) -> Option<(i32, (i32, i32))> {
        self.compute_score0(rings);
        let mut score = 0;
        loop {
            // Extend the M wavefront of this score, checking for the end of the alignment
            if let Some(end) = self.extend(rings, score, plen, tlen, form) {
                return Some((score, end));
            }

            if self.arena.used_bytes >= self.checkpoint_bytes {
                rings.save_snapshot(score);
                self.arena.clear();
                self.checkpoint = Some(CheckpointState {
                    is_i32: std::mem::size_of::<T>() == 4,
                    plen,
                    tlen,
                    form,
                    loaded_above: score,
                });
            }

            score += 1;
            if score > max_score {
                return None;
            }
            self.compute_score(rings, score, plen, tlen, form);
        }
    }

    /// Make the backtrace data for score s available, recomputing it if needed
    ///
    /// Scores must be requested in non-increasing order, as in the backtrace.
    fn ensure_trace_row(&mut self, s: i32) {
        let Some(checkpoint) = self.checkpoint else {
            return;
        };
        if s > checkpoint.loaded_above {
            return;
        }
        if checkpoint.is_i32 {
            let mut rings = std::mem::replace(&mut self.rings32, Rings::new());
            self.recompute_trace_rows(&mut rings, s, checkpoint);
            self.rings32 = rings;
        } else {
            let mut rings = std::mem::replace(&mut self.rings16, Rings::new());
            self.recompute_trace_rows(&mut rings, s, checkpoint);
            self.rings16 = rings;
        }
    }

    /// Recompute the trace rows of the segment between checkpoints which includes score s
    fn recompute_trace_rows<T: ComputeKernel>(
        &mut self,
        rings: &mut Rings<T>,
        s: i32,
        checkpoint: CheckpointState,
    ) {
        let CheckpointState {
            plen, tlen, form, ..
        } = checkpoint;
        let snaps = std::mem::take(&mut rings.snaps);
        // The segment is (snaps[j-1].score, snaps[j].score], or [0, snaps[0].score] for j = 0
        let j = snaps.partition_point(|x| x.score < s);
        let end = snaps[j].score;
        self.arena.clear();
        let start = if j == 0 {
            let snap_bytes = rings.snap_bytes;
            rings.snap_bytes = 0;
            rings.setup(&self.pen, plen as usize, tlen as usize);
            rings.snap_bytes = snap_bytes;
            self.compute_score0(rings);
            0
        } else {
            rings.restore_snapshot(&snaps[j - 1]);
            snaps[j - 1].score
        };
        if j == 0 {
            let x = self.extend(rings, 0, plen, tlen, form);
            debug_assert!(x.is_none());
        }
        for score in (start + 1)..=end {
            self.compute_score(rings, score, plen, tlen, form);
            let x = self.extend(rings, score, plen, tlen, form);
            debug_assert!(x.is_none());
        }
        rings.snaps = snaps;
        self.checkpoint.as_mut().unwrap().loaded_above = if j == 0 { -1 } else { start };
    }

    /// Initialize score 0
    fn compute_score0<T: Offset>(&mut self, rings: &mut Rings<T>) {
        let kb = rings.k_base;
        {
            let w = &mut rings.m[0];
            w.allocated = true;
            w.lo = 0;
            w.hi = 0;
            w.offsets[kb as usize] = T::from_i32(0);
        }
        let (chunk, start) = self.arena.alloc(1);
        self.set_trace(
            0,
            ScoreTrace {
                lo: 0,
                hi: 0,
                chunk,
                start,
                injected: [i32::MIN, i32::MIN],
            },
        );
        // The other components at score 0 are unallocated, as are all ring slots from setup
    }

    /// Compute the wavefronts of score s > 0 from the rings, before extension
    fn compute_score<T: ComputeKernel>(
        &mut self,
        rings: &mut Rings<T>,
        s: i32,
        plen: i32,
        tlen: i32,
        form: EndsFree,
    ) {
        let kb = rings.k_base;
        let pen = self.pen;
        let endsfree_required = form.text_begin != 0 || form.pattern_begin != 0;
        let m_n = rings.m.len();
        let i1_n = rings.i1.len();
        let i2_n = rings.i2.len();
        let null = T::NULL;
        {
            let fetch = |wfs: &Vec<Wavefront<T>>, n: usize, s2: i32| -> Option<usize> {
                if s2 < 0 {
                    None
                } else {
                    let i = ring(s2, n);
                    if wfs[i].is_null() {
                        None
                    } else {
                        Some(i)
                    }
                }
            };
            // Ring slots of the inputs, or None for a null wavefront
            let in_misms = fetch(&rings.m, m_n, s - pen.mismatch);
            let in_open1 = fetch(&rings.m, m_n, s - pen.gap_open1 - pen.gap_ext1);
            let in_i1 = fetch(&rings.i1, i1_n, s - pen.gap_ext1);
            let in_d1 = fetch(&rings.d1, i1_n, s - pen.gap_ext1);
            let in_open2 = fetch(&rings.m, m_n, s - pen.gap_open2 - pen.gap_ext2);
            let in_i2 = fetch(&rings.i2, i2_n, s - pen.gap_ext2);
            let in_d2 = fetch(&rings.d2, i2_n, s - pen.gap_ext2);

            let out_m = ring(s, m_n);
            let out_i1 = ring(s, i1_n);
            let out_i2 = ring(s, i2_n);

            let all_null = in_misms.is_none()
                && in_open1.is_none()
                && in_i1.is_none()
                && in_d1.is_none()
                && in_open2.is_none()
                && in_i2.is_none()
                && in_d2.is_none();
            if all_null {
                // wavefront_compute_allocate_output_null
                for w in [
                    &mut rings.i1[out_i1],
                    &mut rings.d1[out_i1],
                    &mut rings.i2[out_i2],
                    &mut rings.d2[out_i2],
                ] {
                    w.allocated = false;
                }
                let w = &mut rings.m[out_m];
                self.cur_lo = 1;
                self.cur_hi = 0;
                if endsfree_required {
                    // wavefront_compute_endsfree_allocate_null
                    let tbf = form.text_begin >= s;
                    let pbf = form.pattern_begin >= s;
                    let (lo, hi) = if tbf && pbf {
                        (-s, s)
                    } else if tbf {
                        (s, s)
                    } else if pbf {
                        (-s, -s)
                    } else {
                        (0, 0)
                    };
                    for k in lo..=hi {
                        w.offsets[(k + kb) as usize] = null;
                    }
                    let mut injected = [i32::MIN, i32::MIN];
                    if tbf {
                        w.offsets[(s + kb) as usize] = T::from_i32(s);
                        injected[0] = s;
                    }
                    if pbf {
                        w.offsets[(-s + kb) as usize] = T::from_i32(0);
                        injected[1] = -s;
                    }
                    w.allocated = true;
                    w.lo = lo;
                    w.hi = hi;
                    // No trim in this case
                    let (chunk, start) = self.arena.alloc((hi - lo + 1) as usize);
                    self.set_trace(
                        s,
                        ScoreTrace {
                            lo,
                            hi,
                            chunk,
                            start,
                            injected,
                        },
                    );
                } else {
                    w.allocated = false;
                    self.set_trace(
                        s,
                        ScoreTrace {
                            lo: 1,
                            hi: 0,
                            chunk: 0,
                            start: 0,
                            injected: [i32::MIN, i32::MIN],
                        },
                    );
                }
                return;
            }

            // wavefront_compute_limits_input (null inputs have lo=1, hi=-1)
            let lims = |slot: Option<usize>, wfs: &Vec<Wavefront<T>>| -> (i32, i32) {
                match slot {
                    Some(i) => (wfs[i].lo, wfs[i].hi),
                    None => (1, -1),
                }
            };
            let (l, h) = lims(in_misms, &rings.m);
            let mut lo = l;
            let mut hi = h;
            let mut upd = |l: i32, h: i32| {
                lo = lo.min(l);
                hi = hi.max(h);
            };
            let (l, h) = lims(in_open1, &rings.m);
            upd(l - 1, h + 1);
            let (l, h) = lims(in_i1, &rings.i1);
            upd(l + 1, h + 1);
            let (l, h) = lims(in_d1, &rings.d1);
            upd(l - 1, h - 1);
            let (l, h) = lims(in_open2, &rings.m);
            upd(l - 1, h + 1);
            let (l, h) = lims(in_i2, &rings.i2);
            upd(l + 1, h + 1);
            let (l, h) = lims(in_d2, &rings.d2);
            upd(l - 1, h - 1);

            // Fill input cells outside each input's limits with NULL over the range read
            // (wavefront_compute_init_ends)
            fn pad<T: Offset>(w: &mut Wavefront<T>, kb: i32, min_k: i32, max_k: i32) {
                for k in min_k..w.lo.min(max_k + 1) {
                    w.offsets[(k + kb) as usize] = T::NULL;
                }
                for k in (w.hi + 1).max(min_k)..=max_k {
                    w.offsets[(k + kb) as usize] = T::NULL;
                }
            }
            if let Some(i) = in_misms {
                pad(&mut rings.m[i], kb, lo, hi);
            }
            if let Some(i) = in_open1 {
                pad(&mut rings.m[i], kb, lo - 1, hi + 1);
            }
            if let Some(i) = in_open2 {
                pad(&mut rings.m[i], kb, lo - 1, hi + 1);
            }
            if let Some(i) = in_i1 {
                pad(&mut rings.i1[i], kb, lo - 1, hi);
            }
            if let Some(i) = in_d1 {
                pad(&mut rings.d1[i], kb, lo, hi + 1);
            }
            if let Some(i) = in_i2 {
                pad(&mut rings.i2[i], kb, lo - 1, hi);
            }
            if let Some(i) = in_d2 {
                pad(&mut rings.d2[i], kb, lo, hi + 1);
            }

            let alloc_i1 = in_open1.is_some() || in_i1.is_some();
            let alloc_d1 = in_open1.is_some() || in_d1.is_some();
            let alloc_i2 = in_open2.is_some() || in_i2.is_some();
            let alloc_d2 = in_open2.is_some() || in_d2.is_some();

            // Take output wavefronts out of their rings while computing
            let mut wm = std::mem::replace(&mut rings.m[out_m], Wavefront::new(0));
            let mut wi1 = std::mem::replace(&mut rings.i1[out_i1], Wavefront::new(0));
            let mut wd1 = std::mem::replace(&mut rings.d1[out_i1], Wavefront::new(0));
            let mut wi2 = std::mem::replace(&mut rings.i2[out_i2], Wavefront::new(0));
            let mut wd2 = std::mem::replace(&mut rings.d2[out_i2], Wavefront::new(0));

            let width = (hi - lo + 1) as usize;
            // Trace range, including cells which may be added by ends-free initialization
            let tlo = if endsfree_required && form.pattern_begin >= s {
                lo.min(-s)
            } else {
                lo
            };
            let thi = if endsfree_required && form.text_begin >= s {
                hi.max(s)
            } else {
                hi
            };
            let (chunk, tstart) = self.arena.alloc((thi - tlo + 1) as usize);
            if rings.pre.len() < width {
                rings.pre.resize(width, null);
            }
            self.cur_lo = lo;
            self.cur_hi = hi;
            {
                let nul = &rings.null_wf.offsets;
                let get = |slot: Option<usize>, wfs: &Vec<Wavefront<T>>| -> *const T {
                    match slot {
                        Some(i) => wfs[i].offsets.as_ptr(),
                        None => nul.as_ptr(),
                    }
                };
                let p_misms = get(in_misms, &rings.m);
                let p_open1 = get(in_open1, &rings.m);
                let p_open2 = get(in_open2, &rings.m);
                let p_i1 = get(in_i1, &rings.i1);
                let p_d1 = get(in_d1, &rings.d1);
                let p_i2 = get(in_i2, &rings.i2);
                let p_d2 = get(in_d2, &rings.d2);
                let c0 = tstart as usize + (lo - tlo) as usize;
                let choice = &mut self.arena.choice[chunk as usize][c0..c0 + width];
                let pre = &mut rings.pre[..width];
                unsafe {
                    T::compute_kernel(
                        lo,
                        hi,
                        kb,
                        plen,
                        tlen,
                        [p_misms, p_open1, p_open2, p_i1, p_d1, p_i2, p_d2],
                        [
                            wm.offsets.as_mut_ptr(),
                            wi1.offsets.as_mut_ptr(),
                            wd1.offsets.as_mut_ptr(),
                            wi2.offsets.as_mut_ptr(),
                            wd2.offsets.as_mut_ptr(),
                        ],
                        choice,
                        pre,
                    );
                }
            }
            wm.allocated = true;
            wm.lo = lo;
            wm.hi = hi;
            for (w, a) in [
                (&mut wi1, alloc_i1),
                (&mut wd1, alloc_d1),
                (&mut wi2, alloc_i2),
                (&mut wd2, alloc_d2),
            ] {
                w.allocated = a;
                w.lo = lo;
                w.hi = hi;
            }

            // wavefront_compute_process_ends: ends-free initialization
            let mut injected = [i32::MIN, i32::MIN];
            if endsfree_required {
                let ek = s;
                if form.text_begin >= ek {
                    let init = T::from_i32(ek);
                    if wm.hi >= ek {
                        if wm.offsets[(ek + kb) as usize] <= init {
                            wm.offsets[(ek + kb) as usize] = init;
                            injected[0] = ek;
                        }
                    } else {
                        for k in (wm.hi + 1)..ek {
                            wm.offsets[(k + kb) as usize] = null;
                        }
                        wm.offsets[(ek + kb) as usize] = init;
                        injected[0] = ek;
                        wm.hi = ek;
                    }
                }
                if form.pattern_begin >= ek {
                    let ek = -ek;
                    let init = T::from_i32(0);
                    if wm.lo <= ek {
                        if wm.offsets[(ek + kb) as usize] <= init {
                            wm.offsets[(ek + kb) as usize] = init;
                            injected[1] = ek;
                        }
                    } else {
                        wm.offsets[(ek + kb) as usize] = init;
                        injected[1] = ek;
                        for k in (ek + 1)..wm.lo {
                            wm.offsets[(k + kb) as usize] = null;
                        }
                        wm.lo = ek;
                    }
                }
            }
            self.set_trace(
                s,
                ScoreTrace {
                    lo: tlo,
                    hi: thi,
                    chunk,
                    start: tstart,
                    injected,
                },
            );

            // Trim ends
            trim(&mut wm, kb, plen, tlen);
            for (w, a) in [
                (&mut wi1, alloc_i1),
                (&mut wd1, alloc_d1),
                (&mut wi2, alloc_i2),
                (&mut wd2, alloc_d2),
            ] {
                if a {
                    trim(w, kb, plen, tlen);
                }
            }

            rings.m[out_m] = wm;
            rings.i1[out_i1] = wi1;
            rings.d1[out_i1] = wd1;
            rings.i2[out_i2] = wi2;
            rings.d2[out_i2] = wd2;
        }
    }

    /// Extend the M wavefront at score, returning the end position (k, offset) if the alignment
    /// end is reached
    fn extend<T: Offset>(
        &mut self,
        rings: &mut Rings<T>,
        score: i32,
        plen: i32,
        tlen: i32,
        form: EndsFree,
    ) -> Option<(i32, i32)> {
        let kb = rings.k_base;
        let mn = rings.m.len();
        let w = &mut rings.m[ring(score, mn)];
        if !w.allocated || w.lo > w.hi {
            return None;
        }
        let (lo, hi) = (w.lo, w.hi);
        let trace = &self.traces[score as usize];
        let chunk = trace.chunk as usize;
        let (cur_lo, cur_hi) = (self.cur_lo, self.cur_hi);
        let offsets = &mut w.offsets[(lo + kb) as usize..=(hi + kb) as usize];
        let choice = &self.arena.choice[chunk][trace.start as usize..];
        // Matches are only stored for the MemoryHigh backtrace
        let mut matches = if self.arena.with_matches {
            Some(&mut self.arena.matches[chunk][trace.start as usize..])
        } else {
            None
        };
        let pre = &rings.pre;
        let pattern_buf = &self.pattern_buf;
        let text_buf = &self.text_buf;

        for (j, off) in offsets.iter_mut().enumerate() {
            if *off == T::NULL {
                continue;
            }
            let offset = off.to_i32();
            let k = lo + j as i32;
            let new_offset = offset
                + extend_matches(
                    pattern_buf,
                    text_buf,
                    (offset - k) as usize,
                    offset as usize,
                );
            *off = T::from_i32(new_offset);

            // Store matches for the backtrace, relative to the candidate predecessor offset
            if let Some(matches) = matches.as_mut() {
                if cur_lo <= k && k <= cur_hi {
                    let ci = (k - trace.lo) as usize;
                    if choice[ci] & 0xf != BT_NONE {
                        let d = new_offset - pre[(k - cur_lo) as usize].to_i32();
                        if (0..255).contains(&d) {
                            matches[ci] = d as u8;
                        } else {
                            matches[ci] = 255;
                            self.big_matches.insert((score, k), d);
                        }
                    }
                }
            }

            // wavefront_termination_endsfree
            let h_pos = new_offset;
            let v_pos = new_offset - k;
            if h_pos >= tlen || v_pos >= plen {
                if h_pos >= tlen && plen - v_pos <= form.pattern_end {
                    return Some((k, new_offset));
                }
                if v_pos >= plen && tlen - h_pos <= form.text_end {
                    return Some((k, new_offset));
                }
            }
        }
        None
    }
}

/// Trim out-of-bounds cells from the wavefront ends (wavefront_compute_trim_ends)
fn trim<T: Offset>(w: &mut Wavefront<T>, kb: i32, plen: i32, tlen: i32) {
    let mut k = w.hi;
    while k >= w.lo {
        if !w.offsets[(k + kb) as usize].is_oob(k, plen, tlen) {
            break;
        }
        k -= 1;
    }
    w.hi = k;
    let mut k = w.lo;
    while k <= w.hi {
        if !w.offsets[(k + kb) as usize].is_oob(k, plen, tlen) {
            break;
        }
        k += 1;
    }
    w.lo = k;
}

/// Compute all components of one score over diagonals lo..=hi (wavefront_compute_affine2p_idm),
/// with the backtrace choice of each cell and the candidate M offset before clamping
///
/// inputs: misms, open1, open2, i1_ext, d1_ext, i2_ext, d2_ext
/// outputs: m, i1, d1, i2, d2
macro_rules! compute_kernel_body {
    ($T:ty, $lo:ident, $hi:ident, $kb:ident, $plen:ident, $tlen:ident, $inputs:ident, $outputs:ident,
     $choice:ident, $pre:ident) => {{
        let n = ($hi - $lo + 1) as usize;
        let base = ($lo + $kb) as isize;
        let [m_misms, m_open1, m_open2, i1_ext, d1_ext, i2_ext, d2_ext] = $inputs;
        // Slices over the cells read: k-1..=k+1 for each k in lo..=hi
        let sl = |p: *const $T| std::slice::from_raw_parts(p.offset(base - 1), n + 2);
        let misms_in = sl(m_misms);
        let open1 = sl(m_open1);
        let open2 = sl(m_open2);
        let i1e = sl(i1_ext);
        let d1e = sl(d1_ext);
        let i2e = sl(i2_ext);
        let d2e = sl(d2_ext);
        let [out_m, out_i1, out_d1, out_i2, out_d2] = $outputs;
        let slm = |p: *mut $T| std::slice::from_raw_parts_mut(p.offset(base), n);
        let om = slm(out_m);
        let oi1 = slm(out_i1);
        let od1 = slm(out_d1);
        let oi2 = slm(out_i2);
        let od2 = slm(out_d2);
        let pre = &mut $pre[..n];
        let null = <$T as Offset>::NULL;
        // Choices are computed with the offset type first, so that the main loop has a single
        // element width and can be vectorized
        let mut tmp: [$T; 256] = [0; 256];
        let mut j0 = 0;
        while j0 < n {
            let m = (n - j0).min(tmp.len());
            for jj in 0..m {
                let j = j0 + jj;
                let k = $lo + j as i32;
                // I1
                let ins1_o = *open1.get_unchecked(j);
                let ins1_e = *i1e.get_unchecked(j);
                let i1_is_ext = (ins1_e >= ins1_o) as $T;
                let ins1 = ins1_o.max(ins1_e) + 1;
                *oi1.get_unchecked_mut(j) = ins1;
                // I2
                let ins2_o = *open2.get_unchecked(j);
                let ins2_e = *i2e.get_unchecked(j);
                let i2_is_ext = (ins2_e >= ins2_o) as $T;
                let ins2 = ins2_o.max(ins2_e) + 1;
                *oi2.get_unchecked_mut(j) = ins2;
                // D1
                let del1_o = *open1.get_unchecked(j + 2);
                let del1_e = *d1e.get_unchecked(j + 2);
                let d1_is_ext = (del1_e >= del1_o) as $T;
                let del1 = del1_o.max(del1_e);
                *od1.get_unchecked_mut(j) = del1;
                // D2
                let del2_o = *open2.get_unchecked(j + 2);
                let del2_e = *d2e.get_unchecked(j + 2);
                let d2_is_ext = (del2_e >= del2_o) as $T;
                let del2 = del2_o.max(del2_e);
                *od2.get_unchecked_mut(j) = del2;
                // M
                let misms = *misms_in.get_unchecked(j + 1) + 1;
                let max = del1.max(del2).max(misms.max(ins1.max(ins2)));
                // Backtrace choice: highest type among candidates with the max offset, applied
                // from the lowest priority up (BT_xx_EXT = BT_xx_OPEN + 1)
                let mut c = BT_I1_OPEN as $T + i1_is_ext;
                c = if ins2 == max {
                    BT_I2_OPEN as $T + i2_is_ext
                } else {
                    c
                };
                c = if del1 == max {
                    BT_D1_OPEN as $T + d1_is_ext
                } else {
                    c
                };
                c = if del2 == max {
                    BT_D2_OPEN as $T + d2_is_ext
                } else {
                    c
                };
                c = if misms == max { BT_M as $T } else { c };
                c = if max < 0 { BT_NONE as $T } else { c };
                *tmp.get_unchecked_mut(jj) =
                    c | (i1_is_ext << 4) | (i2_is_ext << 5) | (d1_is_ext << 6) | (d2_is_ext << 7);
                *pre.get_unchecked_mut(j) = max;
                *om.get_unchecked_mut(j) = if max.is_oob(k, $plen, $tlen) {
                    null
                } else {
                    max
                };
            }
            let out = &mut $choice[j0..j0 + m];
            for jj in 0..m {
                *out.get_unchecked_mut(jj) = *tmp.get_unchecked(jj) as u8;
            }
            j0 += m;
        }
    }};
}

/// Offset types with a compute kernel
trait ComputeKernel: Offset {
    #[allow(clippy::too_many_arguments)]
    unsafe fn compute_kernel(
        lo: i32,
        hi: i32,
        kb: i32,
        plen: i32,
        tlen: i32,
        inputs: [*const Self; 7],
        outputs: [*mut Self; 5],
        choice: &mut [u8],
        pre: &mut [Self],
    );
}

macro_rules! impl_compute_kernel {
    ($T:ty, $avx2_fn:ident, $base_fn:ident) => {
        #[cfg(target_arch = "x86_64")]
        #[target_feature(enable = "avx2")]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $avx2_fn(
            lo: i32,
            hi: i32,
            kb: i32,
            plen: i32,
            tlen: i32,
            inputs: [*const $T; 7],
            outputs: [*mut $T; 5],
            choice: &mut [u8],
            pre: &mut [$T],
        ) {
            compute_kernel_body!($T, lo, hi, kb, plen, tlen, inputs, outputs, choice, pre)
        }

        #[allow(clippy::too_many_arguments)]
        unsafe fn $base_fn(
            lo: i32,
            hi: i32,
            kb: i32,
            plen: i32,
            tlen: i32,
            inputs: [*const $T; 7],
            outputs: [*mut $T; 5],
            choice: &mut [u8],
            pre: &mut [$T],
        ) {
            compute_kernel_body!($T, lo, hi, kb, plen, tlen, inputs, outputs, choice, pre)
        }

        impl ComputeKernel for $T {
            unsafe fn compute_kernel(
                lo: i32,
                hi: i32,
                kb: i32,
                plen: i32,
                tlen: i32,
                inputs: [*const Self; 7],
                outputs: [*mut Self; 5],
                choice: &mut [u8],
                pre: &mut [Self],
            ) {
                #[cfg(target_arch = "x86_64")]
                {
                    if std::is_x86_feature_detected!("avx2") {
                        return $avx2_fn(lo, hi, kb, plen, tlen, inputs, outputs, choice, pre);
                    }
                }
                $base_fn(lo, hi, kb, plen, tlen, inputs, outputs, choice, pre)
            }
        }
    };
}

impl_compute_kernel!(i16, compute_kernel_i16_avx2, compute_kernel_i16);
impl_compute_kernel!(i32, compute_kernel_i32_avx2, compute_kernel_i32);

/// Number of matching characters from (v,h), with sentinel-padded sequences
#[inline(always)]
fn extend_matches(pattern: &[u8], text: &[u8], v: usize, h: usize) -> i32 {
    let mut n = 0usize;
    loop {
        // The sentinels differ, so a mismatch is always found within the padding
        debug_assert!(v + n + 8 <= pattern.len() && h + n + 8 <= text.len());
        let p = unsafe { (pattern.as_ptr().add(v + n) as *const u64).read_unaligned() };
        let t = unsafe { (text.as_ptr().add(h + n) as *const u64).read_unaligned() };
        let x = p ^ t;
        if x != 0 {
            return (n + (x.trailing_zeros() / 8) as usize) as i32;
        }
        n += 8;
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Comp {
    M,
    I1,
    I2,
    D1,
    D2,
}

impl Affine2pAligner {
    fn trace_byte(&self, s: i32, k: i32) -> u8 {
        let t = &self.traces[s as usize];
        debug_assert!(t.lo <= k && k <= t.hi);
        self.arena.choice[t.chunk as usize][t.start as usize + (k - t.lo) as usize]
    }

    fn trace_matches(&self, s: i32, k: i32) -> i32 {
        let t = &self.traces[s as usize];
        let i = t.start as usize + (k - t.lo) as usize;
        let d = self.arena.matches[t.chunk as usize][i];
        if d == 255 {
            self.big_matches[&(s, k)]
        } else {
            d as i32
        }
    }

    /// Backtrace matching wavefront_backtrace_affine (MemoryHigh)
    fn backtrace_high(
        &mut self,
        plen: i32,
        tlen: i32,
        alignment_score: i32,
        alignment_k: i32,
        alignment_offset: i32,
    ) -> Vec<u8> {
        let pen = self.pen;
        // Operations in reverse order
        let mut ops: Vec<u8> = Vec::with_capacity((plen + tlen) as usize);
        let mut matrix = Comp::M;
        let mut score = alignment_score;
        let mut k = alignment_k;
        let mut offset = alignment_offset;
        let mut v = offset - k;
        let mut h = offset;
        for _ in v..plen {
            ops.push(b'D');
        }
        for _ in h..tlen {
            ops.push(b'I');
        }
        while v > 0 && h > 0 && score > 0 {
            let mismatch = score - pen.mismatch;
            let gap_open1 = score - pen.gap_open1 - pen.gap_ext1;
            let gap_open2 = score - pen.gap_open2 - pen.gap_ext2;
            let gap_extend1 = score - pen.gap_ext1;
            let gap_extend2 = score - pen.gap_ext2;
            self.ensure_trace_row(score);
            let byte = self.trace_byte(score, k);
            let bt = match matrix {
                Comp::M => {
                    let c = byte & 0xf;
                    if c == BT_NONE {
                        break;
                    }
                    let num_matches = self.trace_matches(score, k);
                    for _ in 0..num_matches {
                        ops.push(b'M');
                    }
                    offset -= num_matches;
                    v = offset - k;
                    h = offset;
                    if v <= 0 || h <= 0 {
                        break;
                    }
                    c
                }
                Comp::I1 => {
                    if byte & I1_EXT_FLAG != 0 {
                        BT_I1_EXT
                    } else {
                        BT_I1_OPEN
                    }
                }
                Comp::I2 => {
                    if byte & I2_EXT_FLAG != 0 {
                        BT_I2_EXT
                    } else {
                        BT_I2_OPEN
                    }
                }
                Comp::D1 => {
                    if byte & D1_EXT_FLAG != 0 {
                        BT_D1_EXT
                    } else {
                        BT_D1_OPEN
                    }
                }
                Comp::D2 => {
                    if byte & D2_EXT_FLAG != 0 {
                        BT_D2_EXT
                    } else {
                        BT_D2_OPEN
                    }
                }
            };
            match bt {
                BT_M => {
                    score = mismatch;
                    matrix = Comp::M;
                }
                BT_I1_OPEN => {
                    score = gap_open1;
                    matrix = Comp::M;
                }
                BT_I1_EXT => {
                    score = gap_extend1;
                    matrix = Comp::I1;
                }
                BT_I2_OPEN => {
                    score = gap_open2;
                    matrix = Comp::M;
                }
                BT_I2_EXT => {
                    score = gap_extend2;
                    matrix = Comp::I2;
                }
                BT_D1_OPEN => {
                    score = gap_open1;
                    matrix = Comp::M;
                }
                BT_D1_EXT => {
                    score = gap_extend1;
                    matrix = Comp::D1;
                }
                BT_D2_OPEN => {
                    score = gap_open2;
                    matrix = Comp::M;
                }
                BT_D2_EXT => {
                    score = gap_extend2;
                    matrix = Comp::D2;
                }
                _ => unreachable!(),
            }
            match bt {
                BT_M => {
                    ops.push(b'X');
                    offset -= 1;
                }
                BT_I1_OPEN | BT_I1_EXT | BT_I2_OPEN | BT_I2_EXT => {
                    ops.push(b'I');
                    k -= 1;
                    offset -= 1;
                }
                _ => {
                    ops.push(b'D');
                    k += 1;
                }
            }
            v = offset - k;
            h = offset;
        }
        if matrix == Comp::M {
            if v > 0 && h > 0 {
                let num_matches = v.min(h);
                for _ in 0..num_matches {
                    ops.push(b'M');
                }
                v -= num_matches;
                h -= num_matches;
            }
            for _ in 0..v.max(0) {
                ops.push(b'D');
            }
            for _ in 0..h.max(0) {
                ops.push(b'I');
            }
        } else {
            assert!(
                v == 0 && h == 0,
                "[Affine2pAligner] I?/D?-Beginning backtrace error"
            );
        }
        ops.reverse();
        ops
    }

    /// Backtrace matching the piggyback backtrace of MemoryMed (wavefront_backtrace_pcigar)
    fn backtrace_piggyback(
        &mut self,
        pattern: &[u8],
        text: &[u8],
        alignment_score: i32,
        alignment_k: i32,
        alignment_offset: i32,
    ) -> Vec<u8> {
        let pen = self.pen;
        let plen = pattern.len() as i32;
        let tlen = text.len() as i32;
        // Operations in reverse order, without matches
        let mut rev_ops: Vec<u8> = Vec::new();
        let mut comp = Comp::M;
        let mut s = alignment_score;
        let mut k = alignment_k;
        let (begin_v, begin_h) = loop {
            self.ensure_trace_row(s);
            match comp {
                Comp::M => {
                    let t = &self.traces[s as usize];
                    if t.injected[0] == k {
                        break (0, s);
                    }
                    if t.injected[1] == k {
                        break (s, 0);
                    }
                    if s == 0 {
                        assert_eq!(k, 0);
                        break (0, 0);
                    }
                    let c = self.trace_byte(s, k) & 0xf;
                    match c {
                        BT_M => {
                            rev_ops.push(b'X');
                            s -= pen.mismatch;
                        }
                        BT_I1_OPEN | BT_I1_EXT | BT_I2_OPEN | BT_I2_EXT | BT_D1_OPEN
                        | BT_D1_EXT | BT_D2_OPEN | BT_D2_EXT => {
                            // Gap close
                            rev_ops.push(b'X');
                            comp = match c {
                                BT_I1_OPEN | BT_I1_EXT => Comp::I1,
                                BT_I2_OPEN | BT_I2_EXT => Comp::I2,
                                BT_D1_OPEN | BT_D1_EXT => Comp::D1,
                                _ => Comp::D2,
                            };
                        }
                        _ => panic!("[Affine2pAligner] Missing piggyback backtrace source"),
                    }
                }
                Comp::I1 | Comp::I2 => {
                    let byte = self.trace_byte(s, k);
                    let (is_ext, ext, open) = if comp == Comp::I1 {
                        (
                            byte & I1_EXT_FLAG != 0,
                            pen.gap_ext1,
                            pen.gap_open1 + pen.gap_ext1,
                        )
                    } else {
                        (
                            byte & I2_EXT_FLAG != 0,
                            pen.gap_ext2,
                            pen.gap_open2 + pen.gap_ext2,
                        )
                    };
                    rev_ops.push(b'I');
                    k -= 1;
                    if is_ext {
                        s -= ext;
                    } else {
                        s -= open;
                        comp = Comp::M;
                    }
                }
                Comp::D1 | Comp::D2 => {
                    let byte = self.trace_byte(s, k);
                    let (is_ext, ext, open) = if comp == Comp::D1 {
                        (
                            byte & D1_EXT_FLAG != 0,
                            pen.gap_ext1,
                            pen.gap_open1 + pen.gap_ext1,
                        )
                    } else {
                        (
                            byte & D2_EXT_FLAG != 0,
                            pen.gap_ext2,
                            pen.gap_open2 + pen.gap_ext2,
                        )
                    };
                    rev_ops.push(b'D');
                    k += 1;
                    if is_ext {
                        s -= ext;
                    } else {
                        s -= open;
                        comp = Comp::M;
                    }
                }
            }
        };

        // Unpack forward, extending matches in the M state (wf_backtrace_buffer_unpack_cigar_affine)
        let mut cigar = Vec::with_capacity((plen + tlen) as usize);
        for _ in 0..begin_h {
            cigar.push(b'I');
        }
        for _ in 0..begin_v {
            cigar.push(b'D');
        }
        let mut v = begin_v;
        let mut h = begin_h;
        let mut in_m = true;
        for &op in rev_ops.iter().rev() {
            if in_m {
                while v < plen && h < tlen && pattern[v as usize] == text[h as usize] {
                    cigar.push(b'M');
                    v += 1;
                    h += 1;
                }
            }
            if !in_m && op == b'X' {
                in_m = true;
                continue;
            }
            cigar.push(op);
            match op {
                b'X' => {
                    v += 1;
                    h += 1;
                    in_m = true;
                }
                b'I' => {
                    h += 1;
                    in_m = false;
                }
                _ => {
                    v += 1;
                    in_m = false;
                }
            }
        }
        let end_v = alignment_offset - alignment_k;
        let end_h = alignment_offset;
        let num_matches = (end_v - v).min(end_h - h);
        for _ in 0..num_matches {
            cigar.push(b'M');
        }
        h += num_matches.max(0);
        v += num_matches.max(0);
        while h < tlen {
            cigar.push(b'I');
            h += 1;
        }
        while v < plen {
            cigar.push(b'D');
            v += 1;
        }
        cigar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aligner() -> Affine2pAligner {
        Affine2pAligner::new(-1, 4, 6, 2, 24, 1, false)
    }

    #[test]
    fn test_ends_free_at_sequence_length() {
        let mut aligner = aligner();
        let strict = aligner.align_ends_free(b"ACGTACGT", b"ACGTACGT", 0, 0, 0, 0).score;
        let free = aligner.align_ends_free(b"ACGTACGT", b"ACGTACGT", 8, 8, 8, 8).score;
        assert!(free <= strict);
    }

    #[test]
    #[should_panic(expected = "ends-free lengths")]
    fn test_ends_free_longer_than_sequence() {
        aligner().align_ends_free(b"ACGT", b"ACGTACGT", 0, 5, 0, 0);
    }

    #[test]
    #[should_panic(expected = "ends-free lengths")]
    fn test_ends_free_negative() {
        aligner().align_ends_free(b"ACGT", b"ACGTACGT", 0, 0, -1, 0);
    }
}
