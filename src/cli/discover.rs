use camino::{Utf8Path, Utf8PathBuf};
use clap::Args;
use const_format::concatcp;
use rust_vc_utils::ChromList;
use serde::{Deserialize, Serialize};
use simple_error::{bail, SimpleResult};
use unwrap::unwrap;

use super::defaults::{MIN_GAP_COMPRESSED_IDENTITY, MIN_SV_MAPQ};
use crate::discover::SETTINGS_FILENAME;

#[derive(Args, Default, Deserialize, Serialize)]
pub struct DiscoverSettings {
    /// Directory for all discover command output (must not already exist)
    #[arg(long, value_name = "DIR", default_value = concatcp!(env!("CARGO_PKG_NAME"), "_discover_output"))]
    pub output_dir: Utf8PathBuf,

    /// Alignment file for query sample in BAM or CRAM format.
    #[arg(long = "bam", value_name = "FILE")]
    pub bam_filename: String,

    /// Genome reference in FASTA format
    #[arg(long = "ref", value_name = "FILE")]
    pub ref_filename: String,

    /// Expected copy number values by genomic interval, in BED format.
    ///
    /// Copy number will be read from column 5 of the input BED file. Column 4 is ignored and can
    /// be used as a region label. The default copy number is 2 for unspecified regions. These copy
    /// number values will be stored in discover output for each samplme and automatically selected
    /// for the given sample during joint-calling.
    ///
    /// Note this option is especially useful to clarify expected sex chromosome copy number in the sample.
    ///
    /// By default, the expected copy number only affects CNV calling behavior. If the option
    /// "treat-single-copy-as-haploid" is given during joint calling, then SV genotyping is updated for
    /// single copy regions as well.
    ///
    #[arg(long = "expected-cn", value_name = "FILE")]
    pub expected_copy_number_filename: Option<String>,

    /// Regions of the genome to exclude from CNV analysis, in BED format.
    ///
    /// Note this does not impact SV breakpoint analysis
    ///
    #[arg(long = "cnv-excluded-regions", value_name = "FILE")]
    pub cnv_excluded_regions_filename: Option<String>,

    /// Variant file used to generate minor allele frequency track for this sample, in VCF or BCF format.
    ///
    /// The bigwig track file will be output for each sample in the joint-call step as 'maf.bw' in each
    /// sample directory. The frequences may also be used to improve copy-number segmentation in a future
    /// update.
    ///
    #[arg(long = "maf", value_name = "FILE")]
    pub maf_filename: Option<String>,

    /// Size of bins used for CNV depth track. The segmentation parameters are highly dependent on
    /// this value so it cannot be changed independently of the transition probability and
    /// dependency correction factors
    #[arg(hide = true, long, default_value_t = 1000)]
    pub depth_bin_size: u32,

    /// The transition probability for going away from the current copy number.
    /// Stay probability is derived from this value.
    #[arg(hide = true, long = "transition-probability", default_value_t = 0.02)]
    pub transition_prob: f64,

    /// Depth bin Poisson emision probabilities are raised to this power to help correct
    /// for their having dependencies on neighboring bins, even though they ar treated as
    /// independent for segmentation and quality scoring.
    #[arg(hide = true, long = "dependency-correction", default_value_t = 0.02)]
    pub bin_dependency_correction_factor: f64,

    /// Threshold for the gap-compressed identity filter, to filter out reads with identity to the
    /// reference so low that they are likely to reflect a reference compression or other form
    /// of mismapping.
    #[arg(hide = true, long, default_value_t = MIN_GAP_COMPRESSED_IDENTITY)]
    pub min_gap_compressed_identity: f64,

    /// Regex used to select chromosomes for mean haploid coverage estimation. All selected
    /// chromosomes are assumed diploid.
    #[arg(
        long = "cov-regex",
        value_name = "REGEX",
        default_value = r"^(chr)?\d{1,2}$"
    )]
    pub coverage_est_regex: String,

    /// Number of equal divisions of the GC frequency range to use in the GC correction process
    ///
    /// Depth will be estimated for each GC division, or 'level', which has enough support in the
    /// genome, and the relative depth of each gc level will be used to implement the coverage
    /// correction.
    #[arg(hide = true, long, default_value_t = 40)]
    pub gc_level_count: usize,

    /// This flag enables additional GC correction tracks and files useful for debugging
    #[arg(hide = true, long)]
    pub debug_gc_correction: bool,

    /// Each depth bin is GC corrected based on the GC content of a genomic segment of this size,
    /// centered on the bin.
    #[arg(hide = true, long, default_value_t = 20000)]
    pub gc_genome_window_size: u32,

    /// Co-linear SVs must have either an insertion or deletion of this size or greater to be
    /// included in the output. All other SV evidence patterns such as those consistent with
    /// duplications, inversions and translocations will always be included in the output.
    ///
    #[arg(long, default_value_t = 35)]
    pub min_indel_size: u32,

    /// Specify how much noise is expected in the obs indel size near the minimum size threshold
    ///
    /// This is used to find the minimum evidence indel size, the size used to determine if an
    /// assembly is triggered and stored for a locus. Allowing this trigger for indels smaller than
    /// the reporting size serves several purposes:
    /// 1. Adjacent small indels may merge into a larger indel.
    /// 2. Smaller indels may form a complex haplotype that should be considered in the genotyping
    ///    process, this will be a more accurate non-SV haplotype model then simply using the reference
    ///    sequence
    ///
    #[arg(hide = true, long, default_value_t = 10)]
    pub min_indel_size_noise_margin: u32,

    /// Minimum MAPQ value for reads to be used in SV breakend finding. This does not change depth
    /// analysis.
    ///
    #[arg(long, default_value_t = MIN_SV_MAPQ)]
    pub min_sv_mapq: u32,

    /// Reduce overlapping SV alleles to a single copy
    ///
    /// Certain benchmarks like GIAB SV v0.6 tend to compress both alleles at a VNTR to a single
    /// allele, often genotyped as homozygous. Using this flag will make the SV caller's output
    /// more directly comparable to such representations, but should not be used otherwise.
    ///
    #[arg(long)]
    pub reduce_overlapping_sv_alleles: bool,

    /// Don't canonicalize input file paths
    ///
    /// By default, all file paths input to the discover step are canonicalized and stored in the
    /// discover output directory for use in follow-on joint-call steps. This flag disables all
    /// canonicalization, which allows all paths to be stored as-is, including as relative paths. This
    /// may be useful for situations where the sample discover and joint-call steps are run in
    /// different directory structures.
    ///
    #[arg(long)]
    pub disable_path_canonicalization: bool,

    /// Disable large insertion assembly to save memory/runtime.
    ///
    /// This is currently for debug only.
    ///
    #[arg(hide = true, long)]
    pub disable_large_insertions: bool,

    /// Minimum QUAL score below which the output SV VCF record is marked as filtered
    ///
    /// This shouldn't normally be used in the new discover/joint-call two step mode.
    ///
    #[arg(hide = true, long, default_value_t = 10)]
    pub min_qual: i32,

    /// Maximum distance between breakends for them to be treated as 'close' for the porpose of modifying
    /// the alt haplotype model
    ///
    #[arg(hide = true, long, default_value_t = 2000)]
    pub max_close_breakend_distance: usize,

    /// Remove all clusters except for the given index and turn on high level refinement debug output.
    ///
    #[arg(hide = true, long)]
    pub target_cluster_index: Option<usize>,
}

impl DiscoverSettings {
    pub fn get_min_evidence_indel_size(&self) -> u32 {
        std::cmp::max(self.min_indel_size, self.min_indel_size_noise_margin)
            - self.min_indel_size_noise_margin
    }
}

/// Validate settings and update to parameters that can't be processed automatically by clap.
///
/// Assumes that the logger is not setup
///
pub fn validate_and_fix_discovery_settings(
    settings: DiscoverSettings,
) -> SimpleResult<DiscoverSettings> {
    fn check_required_filename(filename: &str, label: &str) -> SimpleResult<()> {
        if filename.is_empty() {
            bail!("Must specify {label} file");
        }
        if !std::path::Path::new(&filename).exists() {
            bail!("Can't find specified {label} file: '{filename}'");
        }
        Ok(())
    }

    fn check_optional_filename(filename_opt: Option<&String>, label: &str) -> SimpleResult<()> {
        if let Some(filename) = filename_opt {
            if !std::path::Path::new(&filename).exists() {
                bail!("Can't find specified {label} file: '{filename}'");
            }
        }
        Ok(())
    }

    check_required_filename(&settings.ref_filename, "reference")?;

    check_required_filename(&settings.bam_filename, "alignment")?;

    check_optional_filename(settings.maf_filename.as_ref(), "minor allele frequency")?;

    check_optional_filename(
        settings.cnv_excluded_regions_filename.as_ref(),
        "CNV excluded regions",
    )?;

    check_optional_filename(
        settings.expected_copy_number_filename.as_ref(),
        "expected copy number",
    )?;

    if settings.gc_level_count == 0 {
        bail!("--gc-level-count argument must be greater than 0");
    }

    if settings.gc_genome_window_size < settings.depth_bin_size {
        bail!(
            "--gc-genome-window-size is set below the depth bin size of {}",
            settings.depth_bin_size
        );
    }

    // Canonicalize file paths:
    fn canonicalize_string_path(s: &str) -> String {
        Utf8PathBuf::from(s)
            .canonicalize_utf8()
            .unwrap()
            .to_string()
    }

    let mut settings = settings;
    if !settings.disable_path_canonicalization {
        settings.ref_filename = canonicalize_string_path(&settings.ref_filename);
        settings.bam_filename = canonicalize_string_path(&settings.bam_filename);

        settings.cnv_excluded_regions_filename = settings
            .cnv_excluded_regions_filename
            .map(|x| canonicalize_string_path(&x));
        settings.expected_copy_number_filename = settings
            .expected_copy_number_filename
            .map(|x| canonicalize_string_path(&x));
        settings.maf_filename = settings.maf_filename.map(|x| canonicalize_string_path(&x));
    }

    Ok(settings)
}

#[derive(Debug, PartialEq)]
enum SettingValidationError {
    NotFound,
    UnMapped,
    NoChromMatch,
}

fn validate_discover_settings_data_impl(
    settings: &DiscoverSettings,
) -> Result<(), SettingValidationError> {
    use log::error;
    use regex::Regex;
    use rust_htslib::bam::{self, Read};

    // Pull chromosome list from alignment file header, and also check that htslib recognizes the index
    let chrom_list = {
        // Note that IndexedReader is not required here just for the header lookup, but used
        // to additionally test for an index recognized by htslib.
        //
        let bam_reader = match bam::IndexedReader::from_path(&settings.bam_filename) {
            Ok(x) => x,
            Err(error) => {
                error!("Failed to open input alignment file: {}", error);
                return Err(SettingValidationError::NotFound);
            }
        };
        ChromList::from_bam_header(bam_reader.header())
    };

    // Check for unmapped input
    if chrom_list.data.is_empty() {
        error!(
            "Input alignment file is not mapped: '{}'",
            &settings.bam_filename
        );
        return Err(SettingValidationError::UnMapped);
    }

    // Check chromosome regex against alignment file chromosomes to make sure at least one matches:
    let chrom_include_regex = Regex::new(&settings.coverage_est_regex).unwrap();
    let is_any_match = chrom_list
        .data
        .iter()
        .any(|x| chrom_include_regex.is_match(x.label.as_str()));

    if !is_any_match {
        error!(
            "Diploid chromosome regex '{}' does not match any chromosome names in the input alignment file, use '--cov-regex \".\"' to match all available chromosomes.",
            settings.coverage_est_regex
            );
        return Err(SettingValidationError::NoChromMatch);
    }

    Ok(())
}

/// Extended input data/settings validation that's too complex/slow to put in the cmdline parser
///
/// Assumes that the logger is setup
///
pub fn validate_discover_settings_data(settings: &DiscoverSettings) {
    if let Err(err) = validate_discover_settings_data_impl(settings) {
        match err {
            SettingValidationError::NotFound => std::process::exit(exitcode::USAGE),
            _ => std::process::exit(exitcode::DATAERR),
        }
    }
}

/// Write discover settings out in json format
pub fn write_discover_settings(output_dir: &Utf8Path, settings: &DiscoverSettings) {
    use log::info;

    let filename = output_dir.join(SETTINGS_FILENAME);

    info!("Writing discover settings to file: '{filename}'");

    let f = unwrap!(
        std::fs::File::create(&filename),
        "Unable to create discover settings json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &settings).unwrap();
}

pub fn read_discover_settings(discover_dir: &Utf8Path) -> DiscoverSettings {
    use std::fs::File;
    use std::io::BufReader;

    let filename = discover_dir.join(SETTINGS_FILENAME);
    let file = unwrap!(
        File::open(&filename),
        "Unable to read discover-mode settings json file: `{filename}`"
    );
    let reader = BufReader::new(file);
    unwrap!(
        serde_json::from_reader(reader),
        "Unable to parse discover-mode settings from json file: `{filename}`"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_default_settings() -> DiscoverSettings {
        DiscoverSettings::default()
    }

    #[test]
    fn test_not_there() {
        let mut test_settings = get_default_settings();
        test_settings.bam_filename = "./test_data/not_there.bam".to_string();
        let rc = validate_discover_settings_data_impl(&test_settings);
        assert_eq!(rc, Err(SettingValidationError::NotFound))
    }

    #[test]
    fn test_no_index() {
        let mut test_settings = get_default_settings();
        test_settings.bam_filename = "./test_data/no_index.bam".to_string();
        let rc = validate_discover_settings_data_impl(&test_settings);
        assert_eq!(rc, Err(SettingValidationError::NotFound))
    }

    #[test]
    fn test_unmapped() {
        let mut test_settings = get_default_settings();
        test_settings.bam_filename = "./test_data/minimal.bam".to_string();
        let rc = validate_discover_settings_data_impl(&test_settings);
        assert_eq!(rc, Err(SettingValidationError::UnMapped))
    }

    #[test]
    fn test_reference_mismatch() {
        let mut test_settings = get_default_settings();
        test_settings.bam_filename = "./test_data/header_only.bam".to_string();
        test_settings.coverage_est_regex = "test".to_string();
        let rc = validate_discover_settings_data_impl(&test_settings);
        assert_eq!(rc, Err(SettingValidationError::NoChromMatch))
    }
}
