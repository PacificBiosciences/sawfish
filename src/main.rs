mod bam_sa_parser;
mod bam_scanner;
mod bam_utils;
mod bio_align_utils;
mod breakpoint;
mod cli;
mod cluster_breakpoints;
mod contig_output;
mod copy_number_segmentation;
mod depth_bins;
mod discover;
mod expected_ploidy;
mod gc_correction;
mod genome_ref_utils;
mod genome_regions;
mod genome_segment;
mod globals;
mod int_range;
mod joint_call;
mod large_variant_output;
mod log_utils;
mod logger;
mod maf_utils;
mod os_utils;
mod prob_utils;
mod refine_sv;
mod refined_cnv;
mod run_stats;
mod score_sv;
mod simple_alignment;
mod spoa_utils;
mod sv_group;
mod sv_id;
mod utils;
mod vcf_utils;
mod wfa2_utils;
mod worker_thread_data;

use std::{error, process};

use hhmmss::Hhmmss;
use log::info;

use crate::cli::Commands;
use crate::discover::run_discover;
use crate::globals::{PROGRAM_NAME, PROGRAM_VERSION};
use crate::joint_call::run_joint_call;
use crate::logger::setup_output_dir_and_logger;

/// Run system configuration steps prior to starting any other program logic
///
fn system_configuration_prelude() {
    os_utils::attempt_max_open_file_limit();
}

fn run(settings: &cli::Settings) -> Result<(), Box<dyn error::Error>> {
    info!("Starting {PROGRAM_NAME} {PROGRAM_VERSION}");
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
        "{PROGRAM_NAME} completed. Total Runtime: {}",
        start.elapsed().hhmmssxxx()
    );
    Ok(())
}

fn main() {
    system_configuration_prelude();

    let settings = cli::validate_and_fix_settings(cli::parse_settings());

    // Setup logger, including creation of the output directory for the log file:
    setup_output_dir_and_logger(
        settings.get_output_dir(),
        settings.shared.clobber,
        settings.shared.debug,
    );

    if let Err(err) = run(&settings) {
        eprintln!("{err}");
        process::exit(2);
    }
}
