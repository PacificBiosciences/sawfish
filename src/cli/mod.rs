mod discover;
mod joint_call;
mod shared;

use camino::Utf8Path;
use chrono::Datelike;
use clap::{Parser, Subcommand};
use simple_error::{SimpleResult, bail};

use self::discover::validate_and_fix_discovery_settings;
pub use self::discover::{
    DiscoverSettings, read_discover_settings, validate_discover_settings_data,
    write_discover_settings,
};
pub use self::joint_call::JointCallSettings;
use self::joint_call::validate_and_fix_joint_call_settings;
pub use self::shared::SharedSettings;
use self::shared::validate_and_fix_shared_settings;
use crate::globals::PROGRAM_VERSION;

#[derive(Subcommand)]
pub enum Commands {
    /// Discover SV candidate alleles in one sample
    Discover(DiscoverSettings),

    /// Merge and genotype SVs from one or more samples, given the discover command results from each
    JointCall(JointCallSettings),
}

#[derive(Parser)]
#[command(
    author,
    version = PROGRAM_VERSION,
    about,
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
    help_template = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}"
)]
#[clap(propagate_version = true, rename_all = "kebab_case")]
pub struct Settings {
    #[command(flatten)]
    pub shared: SharedSettings,

    #[command(subcommand)]
    pub command: Commands,
}

impl Settings {
    pub fn get_output_dir(&self) -> &Utf8Path {
        match &self.command {
            Commands::Discover(x) => &x.output_dir,
            Commands::JointCall(x) => &x.output_dir,
        }
    }
}

/// Checks if a directory does not exist
///
pub fn check_novel_dirname(dirname: &Utf8Path, label: &str) -> SimpleResult<()> {
    if dirname.exists() {
        bail!("{label} already exists: \"{dirname}\"");
    }
    Ok(())
}

/// Validate settings and update settings that can't be processed by clap
///
/// Assumes that the logger is not setup
///
pub fn validate_and_fix_settings(settings: Settings) -> Settings {
    fn validate_and_fix_settings_impl(mut settings: Settings) -> SimpleResult<Settings> {
        settings.shared = validate_and_fix_shared_settings(settings.shared)?;

        settings.command = match settings.command {
            Commands::Discover(x) => {
                let x = validate_and_fix_discovery_settings(x)?;
                Commands::Discover(x)
            }
            Commands::JointCall(x) => {
                let x = validate_and_fix_joint_call_settings(x)?;
                Commands::JointCall(x)
            }
        };

        Ok(settings)
    }

    match validate_and_fix_settings_impl(settings) {
        Ok(x) => x,
        Err(msg) => {
            eprintln!("Invalid command-line setting: {msg}");
            std::process::exit(exitcode::USAGE);
        }
    }
}

pub fn parse_settings() -> Settings {
    Settings::parse()
}

// Shared default parameter values
pub mod defaults {
    pub const MIN_GAP_COMPRESSED_IDENTITY: f64 = 0.97;
    pub const MIN_SV_MAPQ: u32 = 10;
}
