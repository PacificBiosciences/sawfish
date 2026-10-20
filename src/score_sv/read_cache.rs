//! Cache of the reads used to score SVs, so that each read is decoded and filtered once
//!
//! Scoring reads the alignments overlapping each SV group region in each sample. The same region is
//! read again for every SV in the group, and the long reads overlapping one group often overlap the
//! next. Reading through the index decompresses the same BAM blocks many times, then each read's
//! sequence is decoded and its identity to the reference computed again.
//!
//! Instead, each worker keeps, for each sample, the scoring reads overlapping a few recently used
//! windows of the genome. Each window is read through the index once. The reads overlapping a
//! region within a window are found with the same overlap test and in the same order as an indexed
//! read of the region itself, so scoring sees exactly the same reads.
//!

use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use rust_vc_utils::bam_utils::filter_out_alignment_record;
use rust_vc_utils::cigar::get_cigar_ref_offset;
use rust_vc_utils::get_seq_order_read_split_segments;
use rust_vc_utils::{ChromList, GenomeSegment, IntRange};
use unwrap::unwrap;

use super::{BreakendMatchSplitReadInfo, ScoreSVSettings};
use crate::bam_utils::{bam_fetch_segment, is_split_read};
use crate::genome_ref_utils::{self, GenomeRef};

/// Window size used to read segments too large to cache
const READ_WINDOW_SIZE: i64 = 64 * 1024;

/// Planned segments are read in the same window if the gap between them is at most this
///
/// A shared window reads the long reads overlapping both segments once instead of twice, but
/// also reads those lying within the gap. Of 4, 8, 12, 16 and 32 kb, 8 kb was fastest on HG002.
const MAX_WINDOW_GAP: i64 = 8 * 1024;

/// Planned segments are only added to a window up to this size, to limit memory use
const MAX_WINDOW_SIZE: i64 = 256 * 1024;

/// Segments larger than this are read without caching, to limit memory use
const MAX_CACHED_SEGMENT_SIZE: i64 = 1024 * 1024;

/// Windows using more memory than this are only used for the segment they were read for
const MAX_KEPT_WINDOW_BYTES: usize = 4 * 1024 * 1024;

/// Number of windows kept for each sample
///
/// More than one allows the two breakend regions of a group to be read alternately.
const MAX_READ_WINDOWS: usize = 2;

/// Memory used by all of one worker's cached windows, over all samples, after each group
///
/// Without this, memory would grow with threads x samples (two windows of up to 4 MiB for each).
/// The least recently used windows are dropped first, which only costs a repeated read. At this
/// size the windows of one to three samples are never dropped.
pub const MAX_WORKER_CACHE_BYTES: usize = 64 * 1024 * 1024;

/// Orders window uses over all caches in the process, so that the least recently used window of
/// several caches can be found
static USE_STAMP: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

fn next_use_stamp() -> u64 {
    USE_STAMP.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

/// An alignment record that passes the scoring read filters, with what scoring needs from it
pub(super) struct ScoringRead {
    pub qname: String,
    /// Alignment start and end as defined by htslib for region queries
    begin: i64,
    end: i64,
    pub pos: i64,
    pub is_reverse: bool,
    pub cigar: Vec<Cigar>,
    /// Read sequence in BAM's 4-bit encoding, two bases per byte
    seq_4bit: Vec<u8>,
    seq_len: usize,
    /// Split read segments, if the record is part of a split read
    pub split_read_info: Option<BreakendMatchSplitReadInfo>,
}

impl ScoringRead {
    fn new(chrom_list: &ChromList, record: &bam::Record, begin: i64, end: i64) -> Self {
        let split_read_info = if is_split_read(record) {
            let mut seq_order_read_split_segments = Vec::new();
            let mut primary_split_segment_index = 0;
            for (split_segment_index, split_segment) in
                get_seq_order_read_split_segments(chrom_list, record)
                    .into_iter()
                    .enumerate()
            {
                if split_segment.from_primary_bam_record {
                    primary_split_segment_index = split_segment_index;
                }
                let end = split_segment.pos + get_cigar_ref_offset(&split_segment.cigar);
                let segment = GenomeSegment {
                    chrom_index: split_segment.chrom_index,
                    range: IntRange::from_pair(split_segment.pos, end),
                };
                seq_order_read_split_segments.push(segment);
            }
            Some(BreakendMatchSplitReadInfo {
                seq_order_read_split_segments,
                primary_split_segment_index,
            })
        } else {
            None
        };

        Self {
            qname: std::str::from_utf8(record.qname()).unwrap().to_string(),
            begin,
            end,
            pos: record.pos(),
            is_reverse: record.is_reverse(),
            cigar: record.cigar().take().0,
            seq_4bit: record.seq().encoded.to_vec(),
            seq_len: record.seq_len(),
            split_read_info,
        }
    }

    pub fn seq_len(&self) -> usize {
        self.seq_len
    }

    /// Get a range of the read sequence, with any base other than ACGTN converted to N, as from
    /// get_simplified_dna_seq
    pub fn get_seq(&self, range: std::ops::Range<usize>, seq: &mut Vec<u8>) {
        // BAM 4-bit base codes for "=ACMGRSVTWYHKDBN", simplified to ACGTN
        const DECODE: &[u8; 16] = b"NACNGNNNTNNNNNNN";
        assert!(range.end <= self.seq_len);
        seq.clear();
        seq.extend(range.map(|i| {
            let byte = self.seq_4bit[i / 2];
            let code = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
            DECODE[code as usize]
        }));
    }
}

/// Scoring reads overlapping one window of a chromosome, in BAM file order
struct ReadWindow {
    segment: GenomeSegment,
    reads: Vec<ScoringRead>,
    /// Approximate memory used by reads
    bytes: usize,
    last_use: u64,
}

impl ReadWindow {
    /// Call `f` for each read overlapping segment, which must be within the window
    fn for_each_read(&self, segment: &GenomeSegment, f: impl FnMut(&ScoringRead)) {
        // This is the overlap test of htslib's region iterator (hts_itr_next)
        let (begin, end) = (segment.range.start, segment.range.end);
        self.reads
            .iter()
            .take_while(|r| r.begin < end)
            .filter(|r| r.end > begin)
            .for_each(f);
    }
}

/// Scoring reads of one sample from recently read windows
pub struct SampleReadCache {
    windows: Vec<ReadWindow>,
    /// Segments expected to be read soon, sorted by chromosome and start
    plan: Vec<GenomeSegment>,
}

impl SampleReadCache {
    pub fn new() -> Self {
        Self {
            windows: Vec::new(),
            plan: Vec::new(),
        }
    }

    /// Memory used by the cached windows
    pub fn cached_bytes(&self) -> usize {
        self.windows.iter().map(|w| w.bytes).sum()
    }

    /// Drop the least recently used windows of all of a worker's caches until they use at most
    /// max_bytes
    pub fn trim_to_budget(caches: &mut [SampleReadCache], max_bytes: usize) {
        let mut total = caches.iter().map(|x| x.cached_bytes()).sum::<usize>();
        while total > max_bytes {
            let oldest = caches
                .iter()
                .enumerate()
                .flat_map(|(ci, c)| {
                    c.windows
                        .iter()
                        .enumerate()
                        .map(move |(wi, w)| (w.last_use, ci, wi))
                })
                .min();
            let Some((_, ci, wi)) = oldest else { break };
            total -= caches[ci].windows.swap_remove(wi).bytes;
        }
    }

    /// Set the segments expected to be read soon
    ///
    /// When a segment isn't in the cache, the window read for it is extended to cover any planned
    /// segments that follow it closely.
    pub fn set_plan(&mut self, mut plan: Vec<GenomeSegment>) {
        plan.sort_by_key(|x| (x.chrom_index, x.range.start, x.range.end));
        self.plan = plan;
    }

    /// End of the window to read for segment, extended over the planned segments that follow it
    fn get_window_end(&self, segment: &GenomeSegment) -> i64 {
        let key = (segment.chrom_index, segment.range.start);
        let first = self
            .plan
            .partition_point(|x| (x.chrom_index, x.range.start) < key);
        let mut window_end = segment.range.end;
        for p in self.plan[first..].iter() {
            if p.chrom_index != segment.chrom_index
                || p.range.start > window_end + MAX_WINDOW_GAP
                || p.range.end > segment.range.start + MAX_WINDOW_SIZE
            {
                break;
            }
            window_end = std::cmp::max(window_end, p.range.end);
        }
        window_end
    }

    /// Call `f` for each read that an indexed read of `segment` would return, in the same order,
    /// excluding those that fail the scoring read filters
    ///
    /// The filters are the alignment record filter, minimum MAPQ and minimum gap-compressed
    /// identity, none of which depend on the SV being scored.
    ///
    pub(super) fn for_each_read(
        &mut self,
        score_settings: &ScoreSVSettings,
        reference: &GenomeRef,
        chrom_list: &ChromList,
        bam_reader: &mut bam::IndexedReader,
        segment: &GenomeSegment,
        f: impl FnMut(&ScoringRead),
    ) {
        if verify_read_cache() {
            let mut expected = Vec::new();
            stream_reads(
                score_settings,
                reference,
                chrom_list,
                bam_reader,
                segment,
                usize::MAX,
                |r| expected.push((r.qname.clone(), r.begin, r.end, r.is_reverse)),
            );
            let mut found = Vec::new();
            self.for_each_read_cached(
                score_settings,
                reference,
                chrom_list,
                bam_reader,
                segment,
                |r| found.push((r.qname.clone(), r.begin, r.end, r.is_reverse)),
            );
            assert_eq!(
                found, expected,
                "Read cache differs from an indexed read of {segment:?}"
            );
        }
        self.for_each_read_cached(
            score_settings,
            reference,
            chrom_list,
            bam_reader,
            segment,
            f,
        );
    }

    fn for_each_read_cached(
        &mut self,
        score_settings: &ScoreSVSettings,
        reference: &GenomeRef,
        chrom_list: &ChromList,
        bam_reader: &mut bam::IndexedReader,
        segment: &GenomeSegment,
        f: impl FnMut(&ScoringRead),
    ) {
        if segment.range.size() > MAX_CACHED_SEGMENT_SIZE {
            stream_reads(
                score_settings,
                reference,
                chrom_list,
                bam_reader,
                segment,
                READ_WINDOW_SIZE as usize,
                f,
            );
            return;
        }

        let covers = |w: &ReadWindow| {
            w.segment.chrom_index == segment.chrom_index
                && w.segment.range.start <= segment.range.start
                && segment.range.end <= w.segment.range.end
        };
        let window_index = match self.windows.iter().position(covers) {
            Some(i) => i,
            None => {
                let window_end = self.get_window_end(segment);
                let window = read_window(
                    score_settings,
                    reference,
                    chrom_list,
                    bam_reader,
                    segment,
                    window_end,
                );
                if window.bytes > MAX_KEPT_WINDOW_BYTES {
                    window.for_each_read(segment, f);
                    return;
                }
                if self.windows.len() < MAX_READ_WINDOWS {
                    self.windows.push(window);
                    self.windows.len() - 1
                } else {
                    let (i, _) = self
                        .windows
                        .iter()
                        .enumerate()
                        .min_by_key(|(_, w)| w.last_use)
                        .unwrap();
                    self.windows[i] = window;
                    i
                }
            }
        };
        let window = &mut self.windows[window_index];
        window.last_use = next_use_stamp();
        window.for_each_read(segment, f);
    }
}

/// If true, check the reads of each segment against an indexed read of the segment, and abort on
/// any difference
fn verify_read_cache() -> bool {
    static VERIFY: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *VERIFY.get_or_init(|| std::env::var("SAWFISH_VERIFY_READ_CACHE").is_ok())
}

/// Call `f` for each read of a segment without caching them, reading the segment in windows of
/// window_size to limit memory use
///
/// Each window's reads overlap that window, and those that start before it also overlap the
/// previous window. So the reads overlapping the segment, in file order, are those of the first
/// window followed by those starting within each later window.
///
fn stream_reads(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    segment: &GenomeSegment,
    window_size: usize,
    mut f: impl FnMut(&ScoringRead),
) {
    let mut start = segment.range.start;
    while start < segment.range.end {
        let window_size = i64::try_from(window_size).unwrap_or(i64::MAX);
        let end = std::cmp::min(start.saturating_add(window_size), segment.range.end);
        let window_start = GenomeSegment {
            chrom_index: segment.chrom_index,
            range: IntRange::from_pair(start, end),
        };
        let window = read_window(
            score_settings,
            reference,
            chrom_list,
            bam_reader,
            &window_start,
            end,
        );
        let is_first_window = start == segment.range.start;
        window
            .reads
            .iter()
            .filter(|r| is_first_window || r.begin >= start)
            .for_each(&mut f);
        start = end;
    }
}

/// Read the scoring reads overlapping the window from the start of segment to window_end, limited
/// to the chromosome end
fn read_window(
    score_settings: &ScoreSVSettings,
    reference: &GenomeRef,
    chrom_list: &ChromList,
    bam_reader: &mut bam::IndexedReader,
    segment: &GenomeSegment,
    window_end: i64,
) -> ReadWindow {
    let chrom_size = chrom_list.data[segment.chrom_index].length as i64;
    let window_end = std::cmp::max(segment.range.end, std::cmp::min(window_end, chrom_size));
    let window_segment = GenomeSegment {
        chrom_index: segment.chrom_index,
        range: IntRange::from_pair(segment.range.start, window_end),
    };

    let mut reads = Vec::new();
    bam_fetch_segment(bam_reader, &window_segment);
    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        if filter_out_alignment_record(&record) {
            continue;
        }

        if record.mapq() < score_settings.min_sv_mapq as u8 {
            continue;
        }

        if genome_ref_utils::get_gap_compressed_identity(&record, reference, segment.chrom_index)
            < score_settings.min_gap_compressed_identity
        {
            continue;
        }

        let begin = record.pos();
        let end = unsafe { htslib::bam_endpos(record.inner()) };
        reads.push(ScoringRead::new(chrom_list, &record, begin, end));
    }

    let bytes = reads.capacity() * std::mem::size_of::<ScoringRead>()
        + reads
            .iter()
            .map(|r| {
                r.qname.capacity()
                    + r.seq_4bit.capacity()
                    + r.cigar.capacity() * std::mem::size_of::<Cigar>()
                    + r.split_read_info.as_ref().map_or(0, |x| {
                        x.seq_order_read_split_segments.capacity()
                            * std::mem::size_of::<GenomeSegment>()
                    })
            })
            .sum::<usize>();
    ReadWindow {
        segment: window_segment,
        reads,
        bytes,
        last_use: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// Write a sorted, indexed BAM of random long reads on one chromosome, some of which fail the
    /// scoring read filters
    fn write_test_bam(path: &str, ref_seq: &[u8], rand: &mut impl FnMut(usize) -> usize) {
        let mut header = bam::Header::new();
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", ref_seq.len());
        header.push_record(&sq);
        let mut records = Vec::new();
        for i in 0..800 {
            let len = 1 + rand(if i % 5 == 0 { 30_000 } else { 3000 });
            let pos = rand(ref_seq.len() - len);
            // Read the reference with a deletion, or as an exact copy
            let (cigar, seq) = if rand(3) == 0 && len > 20 {
                let a = rand(len - 10) + 1;
                let d = 1 + rand(9);
                let mut seq = ref_seq[pos..pos + a].to_vec();
                seq.extend_from_slice(&ref_seq[pos + a + d..pos + len]);
                (
                    vec![
                        Cigar::Match(a as u32),
                        Cigar::Del(d as u32),
                        Cigar::Match((len - a - d) as u32),
                    ],
                    seq,
                )
            } else {
                (
                    vec![Cigar::Match(len as u32)],
                    ref_seq[pos..pos + len].to_vec(),
                )
            };
            let mut seq = seq;
            // Some reads have too low an identity to the reference
            if rand(10) == 0 {
                for b in seq.iter_mut() {
                    *b = b"ACGT"[rand(4)];
                }
            }
            let mut record = bam::Record::new();
            let qual = vec![30u8; seq.len()];
            record.set(
                format!("read{i}").as_bytes(),
                Some(&CigarString(cigar)),
                &seq,
                &qual,
            );
            record.set_tid(0);
            record.set_pos(pos as i64);
            record.set_mapq(if rand(10) == 0 { 0 } else { 60 });
            record.set_flags(match rand(20) {
                0 => htslib::BAM_FDUP as u16,
                1 => htslib::BAM_FSECONDARY as u16,
                2 => htslib::BAM_FREVERSE as u16,
                _ => 0,
            });
            records.push(record);
        }
        records.sort_by_key(|x| x.pos());
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
            for r in records.iter() {
                writer.write(r).unwrap();
            }
        }
        bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
    }

    /// Check cached and streamed reads against indexed reads of each segment
    #[test]
    fn test_read_cache_matches_indexed_reads() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut rand = move |n: usize| -> usize {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let chrom_len = 300_000;
        let ref_seq = (0..chrom_len).map(|_| b"ACGT"[rand(4)]).collect::<Vec<_>>();
        let path = std::env::temp_dir().join(format!(
            "sawfish_read_cache_test_{}.bam",
            std::process::id()
        ));
        let path = path.to_str().unwrap();
        write_test_bam(path, &ref_seq, &mut rand);

        let mut chrom_list = ChromList::default();
        chrom_list.add_chrom("chr1", chrom_len as u64);
        let reference = GenomeRef::from_chrom_seqs(vec![ref_seq]);
        let score_settings = ScoreSVSettings::new(10, 0.9, 60, false, true, false);
        let mut bam_reader = bam::IndexedReader::from_path(path).unwrap();

        let reads = |cache: Option<&mut SampleReadCache>,
                     bam_reader: &mut bam::IndexedReader,
                     segment: &GenomeSegment,
                     window_size: usize| {
            let mut x = Vec::new();
            let f = |r: &ScoringRead| x.push((r.qname.clone(), r.begin, r.end, r.seq_4bit.clone()));
            match cache {
                Some(cache) => cache.for_each_read_cached(
                    &score_settings,
                    &reference,
                    &chrom_list,
                    bam_reader,
                    segment,
                    f,
                ),
                None => stream_reads(
                    &score_settings,
                    &reference,
                    &chrom_list,
                    bam_reader,
                    segment,
                    window_size,
                    f,
                ),
            }
            x
        };

        for round in 0..20 {
            let segments = (0..200)
                .map(|_| {
                    let max_size = if rand(20) == 0 { 2_000_000 } else { 5000 };
                    let size = 1 + rand(max_size);
                    let start = rand(chrom_len) as i64;
                    let end = std::cmp::min(start + size as i64, chrom_len as i64);
                    GenomeSegment {
                        chrom_index: 0,
                        range: IntRange::from_pair(start, std::cmp::max(end, start + 1)),
                    }
                })
                .collect::<Vec<_>>();
            let mut cache = SampleReadCache::new();
            if round % 2 == 0 {
                cache.set_plan(segments.clone());
            }
            let mut order = segments.clone();
            if round % 4 < 2 {
                order.sort_by_key(|x| x.range.start);
            }
            let mut total = 0;
            for segment in order.iter() {
                let expected = reads(None, &mut bam_reader, segment, usize::MAX);
                assert_eq!(
                    reads(Some(&mut cache), &mut bam_reader, segment, 0),
                    expected
                );
                assert_eq!(
                    reads(None, &mut bam_reader, segment, 1 + rand(20_000)),
                    expected
                );
                total += expected.len();
            }
            assert!(total > 0);
        }
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(format!("{path}.bai"));
    }

    fn test_window(bytes: usize, last_use: u64) -> ReadWindow {
        ReadWindow {
            segment: GenomeSegment {
                chrom_index: 0,
                range: IntRange::from_pair(0, 1),
            },
            reads: Vec::new(),
            bytes,
            last_use,
        }
    }

    /// The least recently used windows over all samples are dropped first, and nothing is dropped
    /// within the budget
    #[test]
    fn test_trim_to_budget_drops_least_recently_used() {
        let mut caches = (0..3).map(|_| SampleReadCache::new()).collect::<Vec<_>>();
        caches[0].windows = vec![test_window(10, 5), test_window(10, 1)];
        caches[1].windows = vec![test_window(10, 3)];
        caches[2].windows = vec![test_window(10, 4), test_window(10, 2)];
        let uses = |caches: &[SampleReadCache]| {
            let mut x = caches
                .iter()
                .flat_map(|c| c.windows.iter().map(|w| w.last_use))
                .collect::<Vec<_>>();
            x.sort();
            x
        };

        SampleReadCache::trim_to_budget(&mut caches, 50);
        assert_eq!(uses(&caches), vec![1, 2, 3, 4, 5]);
        SampleReadCache::trim_to_budget(&mut caches, 35);
        assert_eq!(uses(&caches), vec![3, 4, 5]);
        SampleReadCache::trim_to_budget(&mut caches, 10);
        assert_eq!(uses(&caches), vec![5]);
        SampleReadCache::trim_to_budget(&mut caches, 0);
        assert!(uses(&caches).is_empty());
        assert_eq!(caches.iter().map(|x| x.cached_bytes()).sum::<usize>(), 0);
    }
}
