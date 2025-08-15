use simple_error::{SimpleResult, bail};

/// Check a required input filename
///
/// Assumes no logger has been configured yet
///
pub fn check_required_filename(filename: &str, label: &str) -> SimpleResult<()> {
    if filename.is_empty() {
        bail!("Must specify {label} file");
    }
    let path = std::path::Path::new(&filename);
    if !path.exists() {
        bail!("Can't find specified {label} file: '{filename}'");
    }
    if !path.is_file() {
        bail!("Specified {label} file path does not appear to be a file: '{filename}'");
    }
    Ok(())
}

/// Check an optional input filename
///
/// Assumes no logger has been configured yet
///
pub fn check_optional_filename(filename_opt: Option<&String>, label: &str) -> SimpleResult<()> {
    if let Some(filename) = filename_opt {
        let path = std::path::Path::new(&filename);
        if !path.exists() {
            bail!("Can't find specified {label} file: '{filename}'");
        }
        if !path.is_file() {
            bail!("Specified {label} file path does not appear to be a file: '{filename}'");
        }
    }
    Ok(())
}
