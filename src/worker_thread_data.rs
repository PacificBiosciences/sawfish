use std::sync::{Arc, Mutex};

use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;

use crate::cli;

/// For worker threads making indexed bam reads, this provides a persistent worker specific reader
/// for each bam file
pub struct BamReaderWorkerThreadData {
    pub bam_readers: Vec<bam::IndexedReader>,
}

impl BamReaderWorkerThreadData {
    pub fn new(bam_filenames: &[&str]) -> Self {
        let bam_readers = bam_filenames
            .iter()
            .map(|x| bam::IndexedReader::from_path(x).unwrap())
            .collect();
        Self { bam_readers }
    }

    /// Setup shared reference data structure for all file reads, in case any of the input alignment files are in CRAM format.
    ///
    /// CRAM files will usually need a pointer to the reference filename (or at least they will be faster than they would pulling
    /// the reference from a URL).
    ///
    /// Additionally, each CRAM file pointer would independently cache reference information to decode all sequence matches to
    /// reference, which can start to make a huge number of copied reference caches when many file pointers are opened. This
    /// can be especially tricky at high-sample counts.
    ///
    /// Here we prepare any CRAM-based file pointers to (1) point to the local reference sequence, and enable a shared reference
    /// pointer between all file pointers.
    ///
    /// Note this logic should all be ignored by htslib for any BAM-formatted input streams
    ///
    /// # Arguments
    /// * ref_filename - the local reference path to be used by CRAM files for sequence decoding
    /// * refs - pointer the the htslib refs datastructure, which can be shared between CRAM file pointers, pass in a null ptr
    ///   if the first refs data structure has not been identified yet.
    ///
    /// Returns the current refs data structure pointer, or null if none were identified
    ///
    pub fn enable_cram_shared_reference(
        &mut self,
        ref_filename: &str,
        mut refs: *mut htslib::refs_t,
    ) -> *mut htslib::refs_t {
        for x in self.bam_readers.iter_mut() {
            if refs.is_null() {
                x.set_reference(ref_filename).unwrap();
                refs = unsafe { htslib::cram_get_refs(x.htsfile()) };
            } else {
                unsafe {
                    htslib::hts_set_opt(
                        x.htsfile(),
                        htslib::hts_fmt_option_CRAM_OPT_SHARED_REF,
                        refs,
                    );
                }
            }
        }

        refs
    }
}

pub type BamReaderWorkerThreadDataSet = Arc<Vec<Mutex<BamReaderWorkerThreadData>>>;

pub fn get_bam_reader_worker_thread_data(
    shared_settings: &cli::SharedSettings,
    ref_filename: &str,
    bam_filenames: &[&str],
) -> BamReaderWorkerThreadDataSet {
    assert!(!bam_filenames.is_empty());

    let mut worker_thread_data = Vec::new();

    // refs enables all cram file pointers across all threads to share the same reference
    // data structure to save memory
    //
    let mut refs: *mut htslib::refs_t = std::ptr::null_mut();

    let orig_hts_log_level = unsafe { htslib::hts_get_log_level() };

    for thread_index in 0..shared_settings.thread_count {
        // Turn htslib warnings off after first bam opening iteration, so that warnings are reported once but aren't
        // repeated thread_count times. Most often we see this repetition for the index timestamp warning:
        //
        // "[W::hts_idx_load3] The index file is older than the data file: ..."
        //
        if thread_index == 1 {
            unsafe {
                htslib::hts_set_log_level(htslib::htsLogLevel_HTS_LOG_ERROR);
            };
        }

        let mut x = BamReaderWorkerThreadData::new(bam_filenames);
        refs = x.enable_cram_shared_reference(ref_filename, refs);
        worker_thread_data.push(Mutex::new(x));
    }

    if shared_settings.thread_count > 1 {
        unsafe {
            htslib::hts_set_log_level(orig_hts_log_level);
        };
    }

    Arc::new(worker_thread_data)
}
