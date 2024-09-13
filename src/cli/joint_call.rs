use std::collections::HashSet;
use std::path::PathBuf;

use clap::Args;
use const_format::concatcp;
use simple_error::{bail, SimpleResult};

use super::defaults::{MIN_GAP_COMPRESSED_IDENTITY, MIN_SV_MAPQ};

#[derive(Args)]
pub struct JointCallSettings {
    /// Directory for all joint-call command output (must not already exist)
    #[arg(long, value_name = "DIR", default_value = concatcp!(env!("CARGO_PKG_NAME"), "_joint-call_output"))]
    pub output_dir: PathBuf,

    /// Sample discover-mode results directory. Can be specified multiple times
    /// to joint call over multiple samples.
    ///
    #[arg(long, value_name = "DIR")]
    pub sample: Vec<PathBuf>,

    /// Minimum QUAL score below which the VCF record is marked as filtered
    #[arg(hide = true, long, default_value_t = 10)]
    pub min_qual: i32,

    /// Min mapq value for reads to be used in SV breakend finding. This does not change depth
    /// analysis.
    #[arg(long, default_value_t = MIN_SV_MAPQ)]
    pub min_sv_mapq: u32,

    /// Threshold for the gap-compressed identity filter, to filter out reads with identity to the
    /// reference so low that they are likely to reflect a reference compression or other form
    /// of mismapping.
    #[arg(hide = true, long, default_value_t = MIN_GAP_COMPRESSED_IDENTITY)]
    pub min_gap_compressed_identity: f64,

    /// Disable end-stage VCF duplicate record filter
    ///
    /// This is only intended for internal debugging cases
    ///
    #[arg(hide = true, long)]
    pub no_vcf_dedup: bool,

    /// Create a JSON output listing the reads supporting each variant
    ///
    /// This is only intended for internal debugging use cases
    ///
    #[arg(hide = true, long)]
    pub report_supporting_reads: bool,
}

/// Validate settings and update to parameters that can't be processed automatically by clap.
///
pub fn validate_and_fix_joint_call_settings(
    mut settings: JointCallSettings,
) -> SimpleResult<JointCallSettings> {
    // Check that all sample paths exist and correspond to directories:
    for sample in settings.sample.iter() {
        if !sample.exists() {
            bail!("--sample argument does not exist: '{}'", sample.display());
        }
        if !sample.is_dir() {
            bail!(
                "--sample argument is not a directory: '{}'",
                sample.display()
            );
        }
    }

    // Canonicalize sample paths
    settings.sample = settings
        .sample
        .iter_mut()
        .map(|x| x.canonicalize().unwrap())
        .collect();

    // Check for repeated entries
    let mut check_paths = HashSet::new();
    for sample in settings.sample.iter() {
        if !check_paths.insert(sample.to_str().unwrap()) {
            bail!("Duplicated --sample path: '{}'", sample.display());
        }
    }

    Ok(settings)
}
