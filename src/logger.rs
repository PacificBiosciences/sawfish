//! Methods specific to the sawfish logger
//!

use std::path::Path;

use crate::cli;
use crate::globals::PROGRAM_NAME;

/// If debug is true set the default logger to the more verbose debug level
///
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
                PROGRAM_NAME,
                record.level(),
                message
            ))
        })
        .level(level)
        .chain(std::io::stderr());

    let logger = if let Some(output_dir) = output_dir {
        let log_filename = output_dir.join(PROGRAM_NAME.to_string() + ".log");
        logger.chain(fern::log_file(log_filename)?)
    } else {
        logger
    };

    logger.apply()?;
    Ok(())
}

/// Check and create output directory, then setup logger to write there
///
/// #Arguments
/// * `debug` - If true use debug log level, and info level otherwise
///
pub fn setup_output_dir_and_logger(output_dir: &Path, clobber: bool, debug: bool) {
    // All error messaging in this method needs to account for no logger being setup yet.
    //
    // We try to match the pre-logging error pattern used in the command-line settings verification methods
    //

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
