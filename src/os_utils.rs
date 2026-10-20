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
                panic!("Can't create new {label} directory at '{dir}': {e}");
            }
        }
    }
}

/// Attempt to increase open file limit to the system's hard limit on *nix-like systems
///
/// This is an optional increase so continue through all failure cases without error. Returns the
/// soft limit in effect afterwards, if it could be determined.
///
/// On macOS the hard limit is reported as unlimited but the kernel rejects a soft limit above
/// OPEN_MAX, so `rlimit::increase_nofile_limit` is used to find a value the system accepts.
///
pub fn attempt_max_open_file_limit() -> Option<u64> {
    rlimit::increase_nofile_limit(u64::MAX).ok()
}

/// Set how long the allocator keeps freed memory before returning it to the OS
///
/// jemalloc's default of 10 s is fine for a long run, but a step that finishes in a few seconds
/// then reports most of its freed memory in its peak RSS.
///
pub fn set_allocator_dirty_decay_ms(decay_ms: isize) {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    {
        use std::ffi::{CString, c_void};
        let set = |name: &str| {
            let name = CString::new(name).unwrap();
            let mut value = decay_ms;
            unsafe {
                tikv_jemalloc_sys::mallctl(
                    name.as_ptr(),
                    std::ptr::null_mut(),
                    std::ptr::null_mut(),
                    &mut value as *mut isize as *mut c_void,
                    std::mem::size_of::<isize>(),
                )
            }
        };
        // Arenas initialized later
        if set("arenas.dirty_decay_ms") != 0 {
            log::warn!("Unable to set allocator dirty page decay time");
            return;
        }
        // Arenas already initialized. This fails for those that aren't, which is fine.
        let mut narenas: u32 = 0;
        let mut len = std::mem::size_of::<u32>();
        let ret = unsafe {
            tikv_jemalloc_sys::mallctl(
                c"arenas.narenas".as_ptr(),
                &mut narenas as *mut u32 as *mut c_void,
                &mut len,
                std::ptr::null_mut(),
                0,
            )
        };
        if ret == 0 {
            for i in 0..narenas {
                set(&format!("arena.{i}.dirty_decay_ms"));
            }
        }
    }
    #[cfg(not(all(target_os = "linux", target_env = "gnu")))]
    let _ = decay_ms;
}
