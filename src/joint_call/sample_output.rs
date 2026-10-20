use camino::{Utf8Path, Utf8PathBuf};
use rust_vc_utils::ChromList;

use super::read_sample_data::SampleJointCallData;

use crate::depth_bins::write_depth_bigwig_file;
use crate::filenames::DEPTH_BINS_BIGWIG_FILENAME;
use crate::gc_correction::write_gc_scaled_depth_track_files;
use crate::joint_call::SharedJointCallData;
use crate::maf_utils::write_maf_track_files;
use crate::os_utils::create_dir_all;

pub fn get_sample_output_dir(
    output_dir: &Utf8Path,
    sample_index: usize,
    sample_name: &str,
) -> Utf8PathBuf {
    let sample_label = format!("sample{:04}_{}", sample_index + 1, sample_name);

    output_dir.to_path_buf().join("samples").join(sample_label)
}

fn setup_sample_output_directories(output_dir: &Utf8Path, all_sample_data: &[SampleJointCallData]) {
    for sample_data in all_sample_data.iter() {
        let sample_dir = get_sample_output_dir(
            output_dir,
            sample_data.sample_index,
            &sample_data.sample_name,
        );
        create_dir_all(&sample_dir, "sample output");
    }
}

/// Create directories for per-sample output, and write any per-sample output files that can be
/// completed at the start of the run
///
pub fn setup_sample_output(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    shared_data: &SharedJointCallData,
    all_sample_data: &[SampleJointCallData],
    thread_count: usize,
) {
    setup_sample_output_directories(output_dir, all_sample_data);

    // Each sample writes its depth and gc scaled depth tracks at the same time, on 2 threads, so
    // this many samples are written at once
    let batch_size = std::cmp::max(thread_count / 2, 1);
    for batch in all_sample_data.chunks(batch_size) {
        std::thread::scope(|scope| {
            for sample_data in batch {
                scope.spawn(move || {
                    write_sample_tracks(output_dir, chrom_list, shared_data, sample_data)
                });
            }
        });
    }
}

fn write_sample_tracks(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    shared_data: &SharedJointCallData,
    sample_data: &SampleJointCallData,
) {
    let sample_index = sample_data.sample_index;
    let sample_dir = get_sample_output_dir(output_dir, sample_index, &sample_data.sample_name);
    let depth_bigwig_filename = &sample_dir.join(DEPTH_BINS_BIGWIG_FILENAME);

    std::thread::scope(|scope| {
        if let Some(sample_gc_bias_data) = sample_data.sample_gc_bias_data.as_ref() {
            let sample_dir = &sample_dir;
            scope.spawn(move || {
                write_gc_scaled_depth_track_files(
                    sample_dir,
                    chrom_list,
                    sample_gc_bias_data,
                    &shared_data.genome_gc_levels,
                    &sample_data.genome_depth_bins,
                )
            });
        }
        // Write depth bins out to bigwig track (for reading in IGV)
        write_depth_bigwig_file(
            depth_bigwig_filename,
            sample_data.genome_depth_bins.bin_size,
            &sample_data.genome_depth_bins.depth_bins,
            chrom_list,
            "depth",
        );
    });

    if let Some(maf_data) = &sample_data.maf_data {
        write_maf_track_files(&sample_dir, chrom_list, maf_data, &sample_data.sample_name);
    }
}
