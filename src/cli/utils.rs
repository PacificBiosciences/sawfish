use simple_error::{SimpleResult, bail};

/// Remote paths (URLs) are not supported, since sawfish is built without remote file access
///
fn bail_if_remote_path(filename: &str, label: &str) -> SimpleResult<()> {
    const REMOTE_PREFIXES: [&str; 7] = [
        "http://", "https://", "ftp://", "s3://", "gs://", "gcs://", "az://",
    ];
    if REMOTE_PREFIXES.iter().any(|x| filename.starts_with(x)) {
        bail!(
            "Remote {label} paths are not supported, localize the file before running sawfish: '{filename}'"
        );
    }
    Ok(())
}

/// Check a required input filename
///
/// Assumes no logger has been configured yet
///
pub fn check_required_filename(filename: &str, label: &str) -> SimpleResult<()> {
    if filename.is_empty() {
        bail!("Must specify {label} file");
    }
    bail_if_remote_path(filename, label)?;
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
        bail_if_remote_path(filename, label)?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remote_paths_are_rejected() {
        for x in ["https://host/a.bam", "s3://bucket/a.bam", "gs://b/ref.fa"] {
            let err = check_required_filename(x, "input").unwrap_err().to_string();
            assert!(
                err.contains("Remote input paths are not supported"),
                "{err}"
            );
            let err = check_optional_filename(Some(&x.to_string()), "input")
                .unwrap_err()
                .to_string();
            assert!(
                err.contains("Remote input paths are not supported"),
                "{err}"
            );
        }
        let err = check_required_filename("/no/such/file.bam", "input")
            .unwrap_err()
            .to_string();
        assert!(err.contains("Can't find"), "{err}");
    }
}
