//! Read ahead of the genotyping workers in each alignment file
//!
//! Genotyping reads the regions of each chunk of SV groups through the alignment file index. Read
//! this way, the OS reads ahead around each
//! scattered request, so that a cold run of HG002 read 51 GB of its 73 GB BAM, about twice the data
//! needed, mostly by seeking on a hard disk. Instead, this asks the OS to read the file ranges of
//! upcoming chunks, from the alignment file indexes, in genome order. The OS then reads only those
//! ranges, with many requests queued in file order, and the workers find their data in the page
//! cache. It only changes when data is read from disk, never what is read.
//!
//! Reading the reference ranges ahead in the same way made cold runs slower, from the seeks between
//! the files.
//!

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use rust_vc_utils::GenomeSegment;

/// Prefetch at most this many chunks ahead of the lowest chunk not yet scored
///
/// A lead of 2 chunks per thread left workers waiting on cold reads when one chunk ran long (cold
/// HPC runs were up to 35% slower); 8 per thread matched the unbounded read-ahead.
pub fn get_prefetch_chunk_lead(thread_count: usize) -> usize {
    8 * thread_count
}

/// Tracks the lowest chunk that is not yet scored, while workers finish chunks in any order
pub struct ChunkProgress {
    finished: Vec<AtomicBool>,
    lowest_unfinished: AtomicUsize,
}

impl ChunkProgress {
    pub fn new(chunk_count: usize) -> Self {
        Self {
            finished: (0..chunk_count).map(|_| AtomicBool::new(false)).collect(),
            lowest_unfinished: AtomicUsize::new(0),
        }
    }

    /// Index of the lowest chunk not yet finished, or the chunk count if all are
    #[cfg(any(target_os = "linux", test))]
    pub fn lowest_unfinished(&self) -> usize {
        self.lowest_unfinished.load(Ordering::Relaxed)
    }

    /// Record that a chunk is scored
    ///
    /// The thread that finishes the last chunk below a run of finished chunks sees the whole run,
    /// because every thread marks its chunk before it looks.
    pub fn mark_finished(&self, chunk_index: usize) {
        self.finished[chunk_index].store(true, Ordering::Release);
        loop {
            let lowest = self.lowest_unfinished.load(Ordering::Relaxed);
            if lowest >= self.finished.len() || !self.finished[lowest].load(Ordering::Acquire) {
                break;
            }
            let _ = self.lowest_unfinished.compare_exchange(
                lowest,
                lowest + 1,
                Ordering::Relaxed,
                Ordering::Relaxed,
            );
        }
    }
}

/// Request the alignment file ranges of each chunk's segments, staying within `chunk_lead` chunks of
/// the lowest unfinished chunk in `progress`, until all chunks are requested or `stop` is set
///
/// The lead is measured from the lowest chunk that is not finished, not from the last chunk started,
/// so that a slow chunk keeps the read-ahead from running far ahead and evicting its data from the
/// page cache.
///
/// # Arguments
/// * `chunk_segments` - Segments read by each chunk, in the order chunks are scored
///
pub fn prefetch_chunks(
    bam_filenames: &[&str],
    chunk_segments: &[Vec<GenomeSegment>],
    progress: &ChunkProgress,
    chunk_lead: usize,
    stop: &AtomicBool,
) {
    #[cfg(target_os = "linux")]
    {
        use rust_htslib::bam;
        use rust_htslib::htslib;

        // Size limit of a BGZF block, so that a range ending at a block start covers that block
        const MAX_BGZF_BLOCK_SIZE: u64 = 64 * 1024;

        let readers = bam_filenames
            .iter()
            .map(|x| bam::IndexedReader::from_path(x).unwrap())
            .collect::<Vec<_>>();
        let bam_files = bam_filenames
            .iter()
            .map(|x| std::fs::File::open(x).unwrap())
            .collect::<Vec<_>>();

        for (chunk_index, segments) in chunk_segments.iter().enumerate() {
            while chunk_index > progress.lowest_unfinished() + chunk_lead {
                if stop.load(Ordering::Relaxed) {
                    return;
                }
                std::thread::sleep(std::time::Duration::from_millis(20));
            }
            if stop.load(Ordering::Relaxed) {
                return;
            }

            for (reader, file) in readers.iter().zip(bam_files.iter()) {
                let mut ranges = Vec::new();
                for segment in segments.iter() {
                    let itr = unsafe {
                        htslib::sam_itr_queryi(
                            reader.index().inner_ptr(),
                            segment.chrom_index as i32,
                            segment.range.start,
                            segment.range.end,
                        )
                    };
                    if itr.is_null() {
                        continue;
                    }
                    let n_off = unsafe { (*itr).n_off } as usize;
                    let off = unsafe { (*itr).off };
                    for i in 0..n_off {
                        let chunk = unsafe { &*off.add(i) };
                        ranges.push((chunk.u >> 16, (chunk.v >> 16) + MAX_BGZF_BLOCK_SIZE));
                    }
                    unsafe { htslib::hts_itr_destroy(itr) };
                }
                advise_will_need(file, ranges);
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    let _ = (bam_filenames, chunk_segments, progress, chunk_lead, stop);
}

/// Ask the OS to read file ranges into the page cache, in file order
#[cfg(target_os = "linux")]
fn advise_will_need(file: &std::fs::File, mut ranges: Vec<(u64, u64)>) {
    use std::os::unix::io::AsRawFd;
    ranges.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for (start, end) in ranges {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    for (start, end) in merged {
        unsafe {
            libc::posix_fadvise(
                file.as_raw_fd(),
                start as libc::off_t,
                (end - start) as libc::off_t,
                libc::POSIX_FADV_WILLNEED,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunk_progress_in_order() {
        let progress = ChunkProgress::new(3);
        assert_eq!(progress.lowest_unfinished(), 0);
        progress.mark_finished(0);
        assert_eq!(progress.lowest_unfinished(), 1);
        progress.mark_finished(1);
        progress.mark_finished(2);
        assert_eq!(progress.lowest_unfinished(), 3);
    }

    #[test]
    fn test_chunk_progress_out_of_order() {
        let progress = ChunkProgress::new(5);
        progress.mark_finished(2);
        progress.mark_finished(1);
        // Chunk 0 is still running
        assert_eq!(progress.lowest_unfinished(), 0);
        progress.mark_finished(4);
        progress.mark_finished(0);
        // The run 0..=2 is finished, chunk 3 isn't
        assert_eq!(progress.lowest_unfinished(), 3);
        progress.mark_finished(3);
        assert_eq!(progress.lowest_unfinished(), 5);
    }

    #[test]
    fn test_chunk_progress_empty() {
        assert_eq!(ChunkProgress::new(0).lowest_unfinished(), 0);
    }

    #[test]
    fn test_chunk_progress_concurrent() {
        // Finish every chunk from several threads in a scrambled order, many times
        for _ in 0..200 {
            let chunk_count = 64;
            let progress = ChunkProgress::new(chunk_count);
            let next = AtomicUsize::new(0);
            std::thread::scope(|scope| {
                for _ in 0..8 {
                    let (progress, next) = (&progress, &next);
                    scope.spawn(move || {
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= chunk_count {
                                break;
                            }
                            // 37 is coprime with 64, so this visits every chunk once, out of order
                            progress.mark_finished(i * 37 % chunk_count);
                        }
                    });
                }
            });
            assert_eq!(progress.lowest_unfinished(), chunk_count);
        }
    }
}
