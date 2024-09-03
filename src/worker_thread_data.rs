use std::sync::{Arc, Mutex};

use rust_htslib::bam;

use crate::cli;

/// For worker threads making indexed bam reads, this provides a persistent worker specific reader
/// for each bam file
pub struct BamReaderWorkerThreadData {
    pub bam_readers: Vec<bam::IndexedReader>,
}

impl BamReaderWorkerThreadData {
    pub fn new(bam_filenames: &[&str], ref_filename: &str) -> Self {
        let mut bam_readers = Vec::new();
        for bam_filename in bam_filenames {
            let mut x = bam::IndexedReader::from_path(bam_filename).unwrap();
            x.set_reference(ref_filename).unwrap();
            bam_readers.push(x);
        }
        Self { bam_readers }
    }
}

pub type BamReaderWorkerThreadDataSet = Arc<Vec<Mutex<BamReaderWorkerThreadData>>>;

pub fn get_bam_reader_worker_thread_data(
    shared_settings: &cli::SharedSettings,
    discovery_settings: &[&cli::DiscoverSettings],
) -> BamReaderWorkerThreadDataSet {
    assert!(!discovery_settings.is_empty());
    let bam_filenames = discovery_settings
        .iter()
        .map(|x| x.bam_filename.as_str())
        .collect::<Vec<_>>();
    let ref_filename = &discovery_settings[0].ref_filename;
    let mut worker_thread_data = Vec::new();
    for _ in 0..shared_settings.thread_count {
        worker_thread_data.push(Mutex::new(BamReaderWorkerThreadData::new(
            &bam_filenames,
            ref_filename,
        )));
    }
    Arc::new(worker_thread_data)
}
