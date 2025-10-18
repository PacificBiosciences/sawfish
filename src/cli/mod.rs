mod discover;
mod joint_call;
mod shared;
mod utils;

use camino::Utf8Path;
use chrono::Datelike;
use clap::{Parser, Subcommand};
use simple_error::{SimpleResult, bail};

use self::discover::validate_and_fix_discovery_settings;
pub use self::discover::{
    DiscoverSettings, read_discover_settings, validate_discover_settings_data,
    write_discover_settings,
};
use self::joint_call::validate_and_fix_joint_call_settings;
pub use self::joint_call::{InputSampleData, JointCallDerivedSettings, JointCallSettings};
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

pub enum CommandDerivedSettings {
    Discover,
    JointCall(JointCallDerivedSettings),
}

pub struct DerivedSettings {
    pub command: CommandDerivedSettings,
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
pub fn validate_and_fix_settings(settings: Settings) -> (Settings, DerivedSettings) {
    fn validate_and_fix_settings_impl(
        mut settings: Settings,
    ) -> SimpleResult<(Settings, DerivedSettings)> {
        settings.shared = validate_and_fix_shared_settings(settings.shared)?;

        let (command, derived_command) = match settings.command {
            Commands::Discover(x) => {
                let x = validate_and_fix_discovery_settings(x)?;
                (Commands::Discover(x), CommandDerivedSettings::Discover)
            }
            Commands::JointCall(x) => {
                let derived_settings = validate_and_fix_joint_call_settings(&x)?;
                (
                    Commands::JointCall(x),
                    CommandDerivedSettings::JointCall(derived_settings),
                )
            }
        };

        settings.command = command;

        let derived_settings = DerivedSettings {
            command: derived_command,
        };

        Ok((settings, derived_settings))
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
    pub const MIN_SV_MAPQ: u32 = 5;
}
