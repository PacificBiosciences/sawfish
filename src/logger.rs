//! Methods specific to the sawfish logger
//!

use camino::Utf8Path;

use crate::cli;
use crate::globals::PROGRAM_NAME;
use crate::os_utils::create_dir_all;

/// If debug is true set the default logger to the more verbose debug level
///
fn setup_logger(output_dir: Option<&Utf8Path>, debug: bool) -> Result<(), fern::InitError> {
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
pub fn setup_output_dir_and_logger(output_dir: &Utf8Path, clobber: bool, debug: bool) {
    // All error messaging in this method needs to account for no logger being setup yet.
    //
    // We try to match the pre-logging error pattern used in the command-line settings verification methods
    //

    if let Err(msg) = cli::check_novel_dirname(output_dir, "Output directory") {
        if !(clobber || output_dir.is_dir()) {
            eprintln!("Invalid command-line setting: {}", msg);
            std::process::exit(exitcode::USAGE);
        }
    };
    create_dir_all(output_dir, "output");
    setup_logger(Some(output_dir), debug).unwrap();
}
