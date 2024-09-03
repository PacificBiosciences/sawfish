mod discover;
mod joint_call;
mod shared;

use std::path::Path;

use chrono::Datelike;
use clap::{Parser, Subcommand};
use simple_error::{bail, SimpleResult};

use self::discover::validate_and_fix_discovery_settings;
pub use self::discover::{
    read_discover_settings, validate_discover_settings_data, write_discover_settings,
    DiscoverSettings,
};
use self::joint_call::validate_and_fix_joint_call_settings;
pub use self::joint_call::JointCallSettings;
use self::shared::validate_and_fix_shared_settings;
pub use self::shared::SharedSettings;

#[derive(Subcommand)]
pub enum Commands {
    /// Discover SV candidate alleles in one sample
    Discover(DiscoverSettings),

    /// Call and genotype SVs in one to many samples, given the discover command results from each
    JointCall(JointCallSettings),
}

#[derive(Parser)]
#[command(
    author,
    version,
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
    pub fn get_output_dir(&self) -> &Path {
        match &self.command {
            Commands::Discover(x) => &x.output_dir,
            Commands::JointCall(x) => &x.output_dir,
        }
    }
}

/// Checks if a directory does not exist
///
pub fn check_novel_dirname(dirname: &Path, label: &str) -> SimpleResult<()> {
    if dirname.exists() {
        bail!("{} already exists: \"{}\"", label, dirname.display());
    }
    Ok(())
}

/// Validate settings and update parameters that can't be processed by clap
///
/// Parts of this process assume logging is already setup
///
pub fn validate_and_fix_settings_impl(mut settings: Settings) -> SimpleResult<Settings> {
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

/// Validate settings and update to parameters that can't be processed automatically by clap.
///
pub fn validate_and_fix_settings(settings: Settings) -> Settings {
    match validate_and_fix_settings_impl(settings) {
        Ok(x) => x,
        Err(msg) => {
            eprintln!("Invalid command-line setting: {}", msg);
            std::process::exit(exitcode::USAGE);
        }
    }
}

pub fn parse_settings() -> Settings {
    Settings::parse()
}
