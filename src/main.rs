mod bam_sa_parser;
mod bam_scanner;
mod bam_utils;
mod bio_align_utils;
mod breakpoint;
mod cli;
mod cluster_breakpoints;
mod cnv_output;
mod contig_output;
mod copy_number_segmentation;
mod depth_bins;
mod discover;
mod expected_ploidy;
mod gc_correction;
mod genome_ref_utils;
mod genome_regions;
mod genome_segment;
mod int_range;
mod joint_call;
mod large_variant_output;
mod log_utils;
mod maf_utils;
mod prob_utils;
mod refine_sv;
mod run_stats;
mod score_sv;
mod simple_alignment;
mod spoa_utils;
mod sv_group;
mod sv_id;
mod utils;
mod vcf_utils;
mod version;
mod wfa2_utils;
mod worker_thread_data;

use std::path::Path;
use std::{error, process};

use hhmmss::Hhmmss;
use log::info;

use crate::cli::Commands;
use crate::copy_number_segmentation::*;
use crate::discover::run_discover;
use crate::joint_call::run_joint_call;
use crate::version::SAWFISH_VERSION;

static PROG_NAME: &str = env!("CARGO_PKG_NAME");

fn setup_logger(output_dir: Option<&Path>, debug: bool) -> Result<(), fern::InitError> {
    let level = if debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };
    let logger = fern::Dispatch::new()
        .format(|out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                PROG_NAME,
                record.level(),
                message
            ))
        })
        .level(level)
        .chain(std::io::stderr());

    let logger = if let Some(output_dir) = output_dir {
        let log_filename = output_dir.join(PROG_NAME.to_string() + ".log");
        logger.chain(fern::log_file(log_filename)?)
    } else {
        logger
    };

    logger.apply()?;
    Ok(())
}

/// Check and create output directory, then setup logger to write there
///
/// All error messaging in this method needs to account for no logger being setup yet.
///
fn setup_output_dir_and_logger(output_dir: &Path, clobber: bool, debug: bool) {
    let mut output_dir_exists = false;
    if let Err(msg) = cli::check_novel_dirname(output_dir, "Output directory") {
        if clobber && output_dir.is_dir() {
            output_dir_exists = true;
        } else {
            eprintln!("Invalid command-line setting: {}", msg);
            std::process::exit(exitcode::USAGE);
        }
    };
    if !output_dir_exists {
        match std::fs::create_dir_all(output_dir) {
            Ok(_) => {}
            Err(e) => {
                panic!(
                    "Can't create new output directory at '{}': {}",
                    output_dir.display(),
                    e
                );
            }
        }
    }
    setup_logger(Some(output_dir), debug).unwrap();
}

fn run(settings: &cli::Settings) -> Result<(), Box<dyn error::Error>> {
    info!("Starting {PROG_NAME} {SAWFISH_VERSION}");
    info!(
        "cmdline: {}",
        std::env::args().collect::<Vec<_>>().join(" ")
    );
    info!("Running on {} threads", settings.shared.thread_count);

    let start = std::time::Instant::now();

    match &settings.command {
        Commands::Discover(x) => {
            run_discover(&settings.shared, x);
        }
        Commands::JointCall(x) => {
            run_joint_call(&settings.shared, x);
        }
    }

    info!(
        "{PROG_NAME} completed. Total Runtime: {}",
        start.elapsed().hhmmssxxx()
    );
    Ok(())
}

fn main() {
    let settings = cli::parse_settings();

    // Validation of output_dir needs to be handled separately so that we don't log error messages
    // before logging is setup.
    setup_output_dir_and_logger(
        settings.get_output_dir(),
        settings.shared.clobber,
        settings.shared.debug,
    );

    let settings = cli::validate_and_fix_settings(settings);

    if let Err(err) = run(&settings) {
        eprintln!("{}", err);
        process::exit(2);
    }
}
