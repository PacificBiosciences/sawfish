//! Utilities pertaining to filesystem and other os-level settings
//!

#![allow(dead_code)]

use camino::Utf8Path;

/// Create a novel directory path if it does not exist already
///
/// If the directory already exists no operations are performed
///
/// * `label` - used to describe the error directory in an error message
///
pub fn create_dir_all(dir: &Utf8Path, label: &str) {
    if !dir.is_dir() {
        match std::fs::create_dir_all(dir) {
            Ok(_) => {}
            Err(e) => {
                panic!("Can't create new {} directory at '{}': {}", label, dir, e);
            }
        }
    }
}

/// Attempt to increase open file limit to the system's hard limit on *nix-like systems
///
/// This is an optional increase so continue through all failure cases without error.
///
pub fn attempt_max_open_file_limit() {
    use rlimit::Resource;

    let (soft, hard) = match Resource::NOFILE.get() {
        Ok(x) => x,
        Err(_) => return,
    };

    if soft < hard {
        rlimit::setrlimit(Resource::NOFILE, hard, hard).unwrap_or_default();
    }
}
