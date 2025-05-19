use rust_htslib::bcf::header::Header;
use rust_htslib::{htslib, utils};
use rust_vc_utils::ChromList;

use crate::globals::PROGRAM_VERSION;

// Imported non-public constants from rust-htslib
pub const VECTOR_END_INTEGER: i32 = i32::MIN + 1;

/*
lazy_static! {
    pub static ref MISSING_FLOAT: f32 = Ieee754::from_bits(0x7F80_0001);
    pub static ref VECTOR_END_FLOAT: f32 = Ieee754::from_bits(0x7F80_0002);
}
 */

/// Get a new bcf header which is actually empty so that we can set our own version number
///
pub fn get_empty_bcf_header() -> Header {
    // Give the incorrect mode to htslib to prevent it from writing the wrong VCF version number.
    // As of 202409 htslib doesn't store the mode argument so this shouldn't break the output.
    //
    let mode = std::ffi::CString::new("r").unwrap();
    Header {
        inner: unsafe { htslib::bcf_hdr_init(mode.as_ptr()) },
        subset: None,
    }
}

/// Builds common fields into a VCF header, upon which more app specific details can be added
///
pub fn get_basic_vcf_header(
    ref_filename: &str,
    chrom_list: &ChromList,
    sample_names: &[&str],
) -> Header {
    let mut header = get_empty_bcf_header();
    header.push_record(b"##fileformat=VCFv4.4");

    let date_string = chrono::Local::now().format("%Y%m%d").to_string();
    header.push_record(format!("##fileDate={date_string}").as_bytes());
    header.push_record(format!("##reference=file://{ref_filename}").as_bytes());
    let prog_name = env!("CARGO_PKG_NAME");
    header.push_record(format!("##source=\"{prog_name} {}\"", PROGRAM_VERSION).as_bytes());
    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    header.push_record(format!("##{prog_name}_cmdline=\"{cmdline}\"").as_bytes());

    // Add contig records
    for chrom_info in chrom_list.data.iter() {
        let header_contig_line = format!(
            "##contig=<ID={},length={}>",
            chrom_info.label, chrom_info.length
        );
        header.push_record(header_contig_line.as_bytes());
    }

    // Add sample names
    for sample_name in sample_names {
        header.push_sample(sample_name.as_bytes());
    }

    header
}

#[derive(Debug)]
pub struct BcfBuildError {
    pub msg: String,
}

impl BcfBuildError {
    pub fn error_message(error: i32) -> &'static str {
        match error {
            -1 => "indexing failed",
            -2 => "opening @fn failed",
            -3 => "format not indexable",
            -4 => "failed to create and/or save the index",
            _ => "unknown error",
        }
    }
}

impl std::error::Error for BcfBuildError {}

impl std::fmt::Display for BcfBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BcfBuildError{{msg: {}}}", self.msg)
    }
}

/// Build index for bcf or vcf.gz file
///
/// # Arguments
/// * `bcf_path` - Path to bcf/vcf file for indexing
/// * `build_tbi` - If true build older tbi style index, otherwise build csi index
///
pub fn build_bcf_index<P: AsRef<std::path::Path>>(
    bcf_path: P,
    n_threads: usize,
    build_tbi: bool,
) -> Result<(), BcfBuildError> {
    let min_shift = if build_tbi { 0 } else { 14 };
    let idx_path_ptr = std::ptr::null();
    let ret = unsafe {
        /*
         *  bcf_index_build3() - Generate and save an index to a specific file
         *  @fn:         Input VCF/BCF filename
         *  @fnidx:      Output filename, or NULL to add .csi/.tbi to @fn
         *  @min_shift:  Positive to generate CSI, or 0 to generate TBI
         *  @n_threads:  Number of VCF/BCF decoder threads
         *
         *  Returns 0 if successful, or negative if an error occurred.
         *
         *  List of error codes:
         *      -1 .. indexing failed
         *      -2 .. opening @fn failed
         *      -3 .. format not indexable
         *      -4 .. failed to create and/or save the index
         */
        htslib::bcf_index_build3(
            utils::path_to_cstring(&bcf_path).unwrap().as_ptr(),
            idx_path_ptr,
            min_shift,
            n_threads as i32,
        )
    };
    match ret {
        0 => Ok(()),
        e => Err(BcfBuildError {
            msg: format!(
                "Failed to build  bcf index. Error: {e:?}/{}",
                BcfBuildError::error_message(e)
            ),
        }),
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf;
    use rust_htslib::bcf::Read;

    /// Test custom header generation:
    #[test]
    fn vcf_from_string() {
        use tempfile::NamedTempFile;

        let mut file = NamedTempFile::new().unwrap();

        writeln!(file, "##fileformat=VCFv4.3").unwrap();
        writeln!(
            file,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        );
        let vcf = bcf::Reader::from_path(file.path()).unwrap();

        let mut header = Header::from_template(vcf.header());
        vcf.close();
    }
}
*/
