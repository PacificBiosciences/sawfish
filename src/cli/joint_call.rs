use std::collections::HashSet;

use camino::Utf8PathBuf;
use clap::Args;
use const_format::concatcp;
use simple_error::{bail, SimpleResult};

use crate::run_stats::read_discover_run_stats;

use super::defaults::{MIN_GAP_COMPRESSED_IDENTITY, MIN_SV_MAPQ};

#[derive(Args)]
pub struct JointCallSettings {
    /// Directory for all joint-call command output (must not already exist)
    #[arg(long, value_name = "DIR", default_value = concatcp!(env!("CARGO_PKG_NAME"), "_joint-call_output"))]
    pub output_dir: Utf8PathBuf,

    /// Sample discover-mode results directory (required). Can be specified multiple times
    /// to joint call over multiple samples.
    ///
    #[arg(long, value_name = "DIR", required = true)]
    pub sample: Vec<Utf8PathBuf>,

    /// Minimum QUAL score below which the VCF record is marked as filtered
    #[arg(hide = true, long, default_value_t = 10)]
    pub min_qual: i32,

    /// Minimum MAPQ value for reads to be used in joint-genotyping.
    #[arg(long, default_value_t = MIN_SV_MAPQ)]
    pub min_sv_mapq: u32,

    /// Threshold for the gap-compressed identity filter, to filter out reads with identity to the
    /// reference so low that they are likely to reflect a reference compression or other form
    /// of mismapping.
    #[arg(hide = true, long, default_value_t = MIN_GAP_COMPRESSED_IDENTITY)]
    pub min_gap_compressed_identity: f64,

    /// Max quality score for all quality outputs
    #[arg(hide = true, long, default_value_t = 999)]
    pub max_qscore: u32,

    /// Disable end-stage VCF duplicate record filter
    ///
    /// This is only intended for internal debugging cases
    ///
    #[arg(hide = true, long)]
    pub no_vcf_dedup: bool,

    /// Create a JSON output file listing the reads supporting each variant
    ///
    #[arg(long)]
    pub report_supporting_reads: bool,

    /// Disable filtering out unused contigs and trimming alignment artifacts in
    /// the aligned contig output.
    ///
    /// This is only intended for internal debugging cases
    ///
    #[arg(hide = true, long)]
    pub no_contig_cleanup: bool,

    /// By default the SV caller ignores the expected copy number values ("--expected-cn" argument in discover).
    ///
    /// With this option, all SVs in single copy regions are genotyped as haploid.
    ///
    #[arg(long)]
    pub treat_single_copy_as_haploid: bool,

    /// Skip checks on input sample discover directory contents
    #[arg(long)]
    pub skip_sample_input_check: bool,
}

/// Check that every sample path represents a completed sawfish discover run
///
/// To do this we rely on the last file written in the discover step 'run.stats.json'
///
fn check_valid_sample_discover_paths(settings: &JointCallSettings) -> SimpleResult<()> {
    use crate::discover::RUN_STATS_FILENAME;

    for sample in settings.sample.iter() {
        let sample_discover_run_stats_filename = sample.join(RUN_STATS_FILENAME);
        if !sample_discover_run_stats_filename.exists()
            || !sample_discover_run_stats_filename.is_file()
        {
            bail!(
                "Sample path: '{sample}' does not contain completed sawfish discover step output"
            );
        }

        let discover_stats = match read_discover_run_stats(sample) {
            Ok(x) => x,
            Err(_) => {
                bail!(
                    "Sample path: '{sample}' does not contain completed sawfish discover step output"
                );
            }
        };

        if discover_stats.run_step.name != "discover" {
            bail!(
                "Sample path: '{sample}' does not contain completed sawfish discover step output"
            );
        }

        // add check on discover step version here later:
        //discover_stats.run_step.version
    }

    Ok(())
}

/// Validate settings and update to parameters that can't be processed automatically by clap.
///
/// Assumes that the logger is not setup
///
pub fn validate_and_fix_joint_call_settings(
    mut settings: JointCallSettings,
) -> SimpleResult<JointCallSettings> {
    // Check that all sample paths exist and correspond to directories:
    for sample in settings.sample.iter() {
        if !sample.exists() {
            bail!("sample argument does not exist: '{sample}'");
        }
        if !sample.is_dir() {
            bail!("sample argument is not a directory: '{sample}'");
        }
    }

    // Canonicalize sample paths
    settings.sample = settings
        .sample
        .into_iter()
        .map(|x| x.canonicalize_utf8().unwrap())
        .collect();

    // Check for repeated entries
    let mut check_paths = HashSet::new();
    for sample in settings.sample.iter() {
        if !check_paths.insert(sample.as_str()) {
            bail!("Duplicated sample path: '{sample}'");
        }
    }

    if !settings.skip_sample_input_check {
        check_valid_sample_discover_paths(&settings)?;
    };

    Ok(settings)
}
