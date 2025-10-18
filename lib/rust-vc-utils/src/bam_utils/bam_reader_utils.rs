use rust_htslib::bam;

pub fn get_bam_reader_filename<T: bam::Read>(reader: &T) -> String {
    use std::ffi::CStr;

    unsafe {
        let c_fn = (*reader.htsfile()).fn_;
        if c_fn.is_null() {
            String::new()
        } else {
            let u8_fn = CStr::from_ptr(c_fn).to_bytes();
            std::str::from_utf8(u8_fn).unwrap().to_string()
        }
    }
}

/// Assert that a BAM/CRAM file has expected EOF marker
///
pub fn assert_bam_eof<T: bam::Read>(reader: &T) {
    use rust_htslib::htslib;

    // Return value info from htslib:
    //    3 for a non-EOF checkable filetype;
    //    2 for an unseekable file type where EOF cannot be checked;
    //    1 for a valid EOF block;
    //    0 for if the EOF marker is absent when it should be present;
    //   -1 (with errno set) on failure
    //
    let eof_check = unsafe { htslib::hts_check_EOF(reader.htsfile()) };
    if eof_check != 1 {
        if eof_check == 0 {
            panic!(
                "Missing EOF marker in BAM/CRAM file: '{}'",
                get_bam_reader_filename(reader)
            );
        } else {
            panic!(
                "Unexpected error ({eof_check}) while checking for EOF marker in BAM/CRAM file: '{}'",
                get_bam_reader_filename(reader)
            );
        }
    }
}
