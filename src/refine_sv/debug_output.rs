use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{
    Arc, Mutex,
    atomic::{self, AtomicBool},
};
use std::time::{Duration, Instant};

use camino::Utf8Path;
use unwrap::unwrap;

use crate::breakpoint::{Breakpoint, ConsolidatedAssemblySegmentInfo};

#[derive(Default)]
pub struct ClusterDurationInfo {
    pub assembly: Duration,
    pub alignment: Duration,
    pub scoring: Duration,
    pub total: Duration,
}

pub struct ClusterDebugInfo {
    pub duration: ClusterDurationInfo,
    pub bp: Breakpoint,
}

/// Sort by processing duration and report out the N slowest clusters:
pub fn write_debug_cluster_info(output_dir: &Utf8Path, mut debug_info: Vec<ClusterDebugInfo>) {
    let cluster_debug_report_cases = 100;
    debug_info.sort_by_key(|x| std::cmp::Reverse(x.duration.total));
    debug_info.truncate(cluster_debug_report_cases);

    let filename = output_dir.join("debug.cluster.refinement.txt");
    let f = unwrap!(
        File::create(&filename),
        "Unable to create cluster refinement debug file: '{filename}'"
    );
    let mut f = BufWriter::new(f);

    for x in debug_info {
        writeln!(
            f,
            "CLUSTER total: {:?} asm: {:?} aln: {:?} score: {:?} bp: {:?}",
            x.duration.total, x.duration.assembly, x.duration.alignment, x.duration.scoring, x.bp
        )
        .unwrap();
    }
}

#[derive(Clone)]
pub struct WorkerDebugInfo {
    pub cluster_start_time: Instant,
    pub bp: Breakpoint,
    pub car: Option<ConsolidatedAssemblySegmentInfo>,
}

pub type AllWorkerDebugInfo = Arc<Vec<Mutex<Option<WorkerDebugInfo>>>>;

fn worker_status_report(worker_debug_info: &AllWorkerDebugInfo) {
    let now = chrono::Local::now();
    eprintln!("[{}] Refinement worker thread status", now);
    for (worker_id, x) in worker_debug_info.iter().enumerate() {
        let worker_id_debug_info = &*x.lock().unwrap();
        eprint!("Worker {} Status: ", worker_id);
        match worker_id_debug_info {
            Some(x) => {
                let duration = x.cluster_start_time.elapsed();
                eprintln!("Job Duration: {duration:?}");
                eprintln!("\tJob Details: {:?}", x.bp);
                if let Some(car) = &x.car {
                    eprintln!("\tconsolidated_region: {:?}", car.segment);
                };
            }
            None => {
                eprintln!("No Job");
            }
        }
    }
}

pub fn worker_status_reporter(
    run_worker_status_reporter: Arc<AtomicBool>,
    worker_debug_info: AllWorkerDebugInfo,
) {
    // Report interval in seconds:
    let worker_status_reporter_interval = 5 * 60;
    loop {
        std::thread::sleep(Duration::from_secs(worker_status_reporter_interval));
        if !run_worker_status_reporter.load(atomic::Ordering::Relaxed) {
            return;
        }
        worker_status_report(&worker_debug_info);
    }
}
