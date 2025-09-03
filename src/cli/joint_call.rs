use std::collections::HashSet;

use camino::Utf8PathBuf;
use clap::Args;
use const_format::concatcp;
use csv::{ReaderBuilder, Trim};
use regex::Regex;
use simple_error::{SimpleResult, bail, map_err_with};
use unwrap::unwrap;

use super::utils::check_optional_filename;
use crate::cli::utils::check_required_filename;
use crate::run_stats::read_discover_run_stats;
use crate::utils::canonicalize_string_path;

use super::defaults::{MIN_GAP_COMPRESSED_IDENTITY, MIN_SV_MAPQ};
use super::discover::read_discover_settings;

#[derive(Args)]
#[group(required = true, multiple = false)]
pub struct SampleInputGroup {
    /// Sample discover-mode results directory (required). Can be specified multiple times
    /// to joint call over multiple samples.
    ///
    #[arg(long, value_name = "DIR")]
    pub sample: Vec<Utf8PathBuf>,

    /// Sample csv file
    #[arg(long, value_name = "FILE")]
    pub sample_csv: Option<Utf8PathBuf>,
}

#[derive(Args)]
pub struct JointCallSettings {
    /// Directory for all joint-call command output (must not already exist)
    #[arg(long, value_name = "DIR", default_value = concatcp!(env!("CARGO_PKG_NAME"), "_joint-call_output"))]
    pub output_dir: Utf8PathBuf,

    #[command(flatten)]
    pub sample_input_group: SampleInputGroup,

    /// Genome reference in FASTA format. If not specified, the reference path from the first sample's discover command
    /// will be automatically selected instead.
    ///
    #[arg(long = "ref", value_name = "FILE")]
    pub ref_filename: Option<String>,

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

    /// Regex used to select chromosomes where max depth filtration is disabled
    ///
    /// Default value is intended to disable high depth filtration on mitochondria
    ///
    #[arg(
        long,
        value_name = "REGEX",
        default_value = r"(?i)^(chr)?(m|mt|mito|mitochondria)$"
    )]
    pub disable_max_dapth_chrom_regex: String,

    /// Disable SV/CNV merge, for internal debug only
    #[arg(hide = true, long)]
    pub disable_sv_cnv_merge: bool,

    /// Disable any additional CNV segmentation and output. Note compared to the discover step
    /// option of the same name, this disables CNV in all samples
    #[arg(hide = true, long)]
    pub disable_cnv: bool,
}

/// All per-sample input values to joint-call
pub struct InputSampleData {
    /// Sawfish discover output directory of the sample:
    pub discover_dir: Utf8PathBuf,

    /// Optional bam filename specification. If given, this value will be used instead of the one
    /// found in the discover_dir configuration files
    pub bam_filename: Option<String>,
}

pub struct JointCallDerivedSettings {
    pub all_input_sample_data: Vec<InputSampleData>,
}

/// Build input sample data structure from either the --sample or --sample-csv options
///
/// This routine does not QC the discover dir, etc. but only the format of the CSV file itself. Other QC checks should
/// be run downstream from this routine.
///
fn get_all_input_sample_data(
    sample_input_group: &SampleInputGroup,
) -> SimpleResult<Vec<InputSampleData>> {
    // clap should already be restricting input to one of these two options, but double-check here:
    assert!(sample_input_group.sample.is_empty() ^ sample_input_group.sample_csv.is_none());

    let all_input_sample_data = if let Some(sample_csv_filename) = &sample_input_group.sample_csv {
        if !sample_csv_filename.exists() || !sample_csv_filename.is_file() {
            bail!("Input sample csv path does not exist or is not a file: '{sample_csv_filename}'");
        }
        let mut rdr = ReaderBuilder::new()
            .has_headers(false)
            .flexible(true)
            .trim(Trim::All)
            .comment(Some(b'#'))
            .delimiter(b',')
            .from_path(sample_csv_filename)
            .unwrap();

        let mut sample_data = Vec::new();
        for (line_index, result) in rdr.records().enumerate() {
            let record = unwrap!(
                result,
                "Failed to parse csv record line {} from sample csv file: '{sample_csv_filename}'",
                line_index + 1,
            );

            let discover_dir = unwrap!(
                record.get(0),
                "Missing required first column in csv record line {} from sample csv file: '{sample_csv_filename}'",
                line_index + 1,
            );

            let bam_filename = match record.get(1) {
                Some(x) => {
                    if x.is_empty() {
                        None
                    } else {
                        Some(x.to_owned())
                    }
                }
                None => None,
            };

            sample_data.push(InputSampleData {
                discover_dir: discover_dir.into(),
                bam_filename,
            });
        }
        sample_data
    } else {
        sample_input_group
            .sample
            .iter()
            .map(|x| InputSampleData {
                discover_dir: x.clone(),
                bam_filename: None,
            })
            .collect()
    };

    Ok(all_input_sample_data)
}

/// Check that every sample discover dir represents a completed sawfish discover run
///
/// To do this we rely on the last file written in the discover step 'run.stats.json'
///
/// Also check that the ref_filename and bam_filename in the discover settings json file are valid
/// if no alternate ref and bam filenames have been provided.
///
fn check_valid_sample_discover_dirs(
    settings: &JointCallSettings,
    all_input_sample_data: &[InputSampleData],
) -> SimpleResult<()> {
    use crate::filenames::RUN_STATS_FILENAME;

    for input_sample_data in all_input_sample_data.iter() {
        let discover_dir = &input_sample_data.discover_dir;
        let discover_run_stats_filename = discover_dir.join(RUN_STATS_FILENAME);
        if !discover_run_stats_filename.exists() || !discover_run_stats_filename.is_file() {
            bail!(
                "Input sample discover dir: '{discover_dir}' does not contain completed sawfish discover step output"
            );
        }

        let discover_stats = match read_discover_run_stats(discover_dir) {
            Ok(x) => x,
            Err(_) => {
                bail!(
                    "Input sample discover dir: '{discover_dir}' does not contain completed sawfish discover step output"
                );
            }
        };

        if discover_stats.run_step.name != "discover" {
            bail!(
                "Input sample discover dir: '{discover_dir}' does not contain completed sawfish discover step output"
            );
        }

        // add check on discover step version here later:
        //discover_stats.run_step.version

        // Check ref and bam filenames in the discover settings json file
        let discover_settings = read_discover_settings(discover_dir);

        if settings.ref_filename.is_none() {
            match check_required_filename(&discover_settings.ref_filename, "reference") {
                Ok(x) => x,
                Err(msg) => {
                    bail!(
                        "For reference file specified in input sample discover directory '{discover_dir}': {msg}"
                    );
                }
            }
        }

        if input_sample_data.bam_filename.is_none() {
            match check_required_filename(&discover_settings.bam_filename, "alignment") {
                Ok(x) => x,
                Err(msg) => {
                    bail!(
                        "For alignment file specified in input sample discover directory '{discover_dir}': {msg}"
                    );
                }
            }
        }
    }

    Ok(())
}

/// Validate settings and update to parameters that can't be processed automatically by clap.
///
/// Assumes that the logger is not setup
///
pub fn validate_and_fix_joint_call_settings(
    settings: &JointCallSettings,
) -> SimpleResult<JointCallDerivedSettings> {
    let mut all_input_sample_data = get_all_input_sample_data(&settings.sample_input_group)?;

    // Check that all sample discover paths exist and correspond to directories:
    for (sample_index, input_sample_data) in all_input_sample_data.iter().enumerate() {
        let discover_dir = &input_sample_data.discover_dir;
        if !discover_dir.exists() {
            bail!(
                "For input sample {}, discover path does not exist: '{discover_dir}'",
                sample_index + 1
            );
        }
        if !discover_dir.is_dir() {
            bail!(
                "For input sample {}, discover path is not a directory: '{discover_dir}'",
                sample_index + 1
            );
        }
    }

    // Check that any bam filenames from csv input exist and correspond to files:
    for (sample_index, input_sample_data) in all_input_sample_data.iter().enumerate() {
        let bam_filename = &input_sample_data.bam_filename;
        if let Some(bam_filename) = bam_filename {
            let bam_path = std::path::Path::new(bam_filename);
            if !bam_path.exists() {
                bail!(
                    "For input sample {}, the bam path input from the sample csv does not exist: '{bam_filename}'",
                    sample_index + 1
                );
            }
            if !bam_path.is_file() {
                bail!(
                    "For input sample {}, the bam path input from the sample csv is not a file: '{bam_filename}'",
                    sample_index + 1
                );
            }
        }
    }

    // Canonicalize sample paths
    for input_sample_data in all_input_sample_data.iter_mut() {
        input_sample_data.discover_dir =
            input_sample_data.discover_dir.canonicalize_utf8().unwrap();
        input_sample_data.bam_filename = input_sample_data
            .bam_filename
            .as_ref()
            .map(|x| canonicalize_string_path(x));
    }

    // Check for repeated entries
    let mut check_dirs = HashSet::new();
    for input_sample_data in all_input_sample_data.iter() {
        if !check_dirs.insert(input_sample_data.discover_dir.as_str()) {
            bail!(
                "Duplicated input sample discover dir: '{}'",
                input_sample_data.discover_dir
            );
        }
    }

    if !settings.skip_sample_input_check {
        check_valid_sample_discover_dirs(settings, &all_input_sample_data)?;
    };

    check_optional_filename(settings.ref_filename.as_ref(), "reference")?;

    // Check that regex is valid
    let _ = map_err_with!(
        Regex::new(&settings.disable_max_dapth_chrom_regex),
        "Invalid regex for --disable-max-dapth-chrom-regex"
    )?;

    let derived_settings = JointCallDerivedSettings {
        all_input_sample_data,
    };

    Ok(derived_settings)
}
