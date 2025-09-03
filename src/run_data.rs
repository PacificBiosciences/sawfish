//! Track misc data for sawfish run
//!
//! Run data is used to capture various values from the sawfish run which don't have an obvious location in
//! any other output file.
//!

use std::fs::File;

use camino::Utf8Path;
use log::info;
use serde::{Deserialize, Serialize};
use simple_error::{SimpleResult, try_with};
use unwrap::unwrap;

use crate::filenames::RUN_DATA_FILENAME;
use crate::run_stats::RunStep;

#[derive(Deserialize, Serialize)]
pub struct DiscoverRunData {
    pub run_step: RunStep,
    pub sample_name: String,
    pub max_sv_scoring_depth: u32,
}

#[derive(Deserialize, Serialize)]
pub struct JointCalSampleRunData {
    pub run_step: RunStep,
    pub sample_name: String,
}

/// Write run data structure out in json format
pub fn write_discover_run_data(discover_dir: &Utf8Path, run_data: &DiscoverRunData) {
    let filename = discover_dir.join(RUN_DATA_FILENAME);

    info!("Writing run data to file: '{filename}'");

    let f = unwrap!(
        File::create(&filename),
        "Unable to create run data json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &run_data).unwrap();
}

#[allow(unused)]
pub fn read_discover_run_data(discover_dir: &Utf8Path) -> SimpleResult<DiscoverRunData> {
    use std::io::BufReader;

    let filename = discover_dir.join(RUN_DATA_FILENAME);
    let file = try_with!(
        File::open(&filename),
        "Unable to read discover-step run data json file: '{filename}'"
    );

    let reader = BufReader::new(file);
    let run_stats = try_with!(
        serde_json::from_reader(reader),
        "Unable to parse discover-step run data from json file: '{filename}'"
    );

    Ok(run_stats)
}

/*
/// Write run_stats structure out in json format
pub fn write_joint_call_run_stats(output_dir: &Utf8Path, run_stats: &JointCallRunStats) {
    let filename = output_dir.join(RUN_STATS_FILENAME);

    info!("Writing run data to file: '{filename}'");

    let f = unwrap::unwrap!(
        std::fs::File::create(&filename),
        "Unable to create run statistics json file: '{filename}'"
    );

    serde_json::to_writer_pretty(&f, &run_stats).unwrap();
}

/// Delete run stats file
///
/// Delete a run stats file if one exists. This is typically done during a clobber run,
/// to prevent the old run stats file from being misinterpreted in the event of a crash.
///
pub fn delete_run_stats(output_dir: &Utf8Path) {
    let filename = output_dir.join(RUN_STATS_FILENAME);

    if filename.exists() {
        unwrap!(
            remove_file(&filename),
            "Can't remove original run statistics json file: '{filename}'"
        );
    }
}
*/
