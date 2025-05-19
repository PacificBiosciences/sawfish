use clap::Args;
use simple_error::{bail, SimpleResult};

#[derive(Args)]
pub struct SharedSettings {
    /// Number of threads to use. Defaults to all logical cpus detected.
    #[arg(long = "threads", global = true, value_name = "THREAD_COUNT")]
    thread_count_option: Option<usize>,

    /// This value will be filled in by thread_count_option
    #[arg(hide = true, default_value_t = 0)]
    pub thread_count: usize,

    /// Overwrite an existing output directory
    #[arg(long, global = true)]
    pub clobber: bool,

    /// Turn on extra debug logging
    ///
    /// This option enables extra logging intended for debugging only. It is highly
    /// recommended (but not required) to set --threads to 1 when this is enabled.
    ///
    #[arg(long, global = true)]
    pub debug: bool,

    /// Specify one or more target regions for SV calling
    ///
    /// This option is provided strictly for debugging at this point. It is the caller's
    /// responsibility to ensure that regions do not overlap.
    ///
    #[arg(hide = true, long = "target-region", global = true)]
    pub target_region_list: Vec<String>,

    /// Disable small indels to save memory/runtime.
    ///
    /// This is currently for debug only. Useful when quickly analyzing large and complex SVs at the genome level.
    ///
    /// Note that "small" is defined by sawfish mechanics rather than size, so this will disable all single region assemblies
    ///
    #[arg(hide = true, long, global = true)]
    pub disable_small_indels: bool,
}

pub fn validate_and_fix_shared_settings(
    mut settings: SharedSettings,
) -> SimpleResult<SharedSettings> {
    settings.thread_count = match settings.thread_count_option {
        Some(count) => {
            if count == 0 {
                bail!("--threads argument must be greater than 0");
            }
            count
        }
        None => num_cpus::get(),
    };

    Ok(settings)
}
