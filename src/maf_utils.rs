use std::collections::{BTreeMap, HashMap};

use camino::Utf8Path;
use itertools::Itertools;
use log::info;
use rust_htslib::bcf::{self, Read};
use rust_vc_utils::{ChromList, bigwig_utils};
use serde::{Deserialize, Serialize};
use unwrap::unwrap;

use crate::discover::{MAF_BIGWIG_FILENAME, MAF_MESSAGEPACK_FILENAME};

#[derive(Clone, Deserialize, Serialize)]
pub struct SampleMafData {
    pub sample_name: String,

    /// A vector of biallelic counts for each chromosome
    ///
    /// Chromosome order in the vector follows that of the ChromList used during file scanning
    ///
    /// For each chrom, the map key is zero-indexed position, the value is a tuple of the two allele counts
    pub data: Vec<BTreeMap<i64, (i32, i32)>>,
}

/// Minor allele frequency results from all samples
///
#[derive(Default, Deserialize, Serialize)]
pub struct MafData {
    pub samples: Vec<SampleMafData>,
}

/// Scan VCF/BCF file to create minor allele frequency track for the genome
///
/// # Arguments
/// * `chrom_list` - All VCF chromosome orders will be arranged to fit the chromosome order in chrom_list
/// * `sample_names` - If empty, then read all samples in the VCF, otherwise only read the
///   specified samples, and panic if these are not found. Order of samples in the result
///   will match those in sample_names.
///
pub fn scan_maf_file(filename: &str, chrom_list: &ChromList, sample_names: &[&str]) -> MafData {
    assert!(!filename.is_empty());

    info!("Scanning minor allele frequency data from file '{filename}'");

    let mut maf_data = MafData::default();

    let mut reader = bcf::Reader::from_path(filename).unwrap();
    let header = reader.header();

    // Translate chrom order from VCF (rid) to match that in chrom_list
    let vcf_to_output_chrom_index_map = {
        let mut x = Vec::new();
        for vcf_chrom_index in 0..header.contig_count() {
            let vcf_chrom_bytes = header.rid2name(vcf_chrom_index).unwrap();
            let vcf_chrom_name = std::str::from_utf8(vcf_chrom_bytes).unwrap().to_string();
            match chrom_list.label_to_index.get(&vcf_chrom_name) {
                Some(&output_chrom_index) => {
                    x.push(output_chrom_index);
                }
                None => {
                    panic!(
                        "Input alignment file does not includes chromosome '{vcf_chrom_name}' from minor allele frequency file '{filename}'"
                    );
                }
            }
        }
        x
    };

    // Get sample names from header, and build a map of VCF sample index values to output sample index values
    let vcf_to_output_sample_index_map = {
        let mut x = Vec::new();
        let output_chrom_count = chrom_list.data.len();
        if sample_names.is_empty() {
            for (vcf_sample_index, &sample_bytes) in header.samples().iter().enumerate() {
                let sample_name = std::str::from_utf8(sample_bytes).unwrap().to_string();
                maf_data.samples.push(SampleMafData {
                    sample_name,
                    data: vec![BTreeMap::new(); output_chrom_count],
                });
                x.push(Some(vcf_sample_index));
            }
        } else {
            let sample_map: HashMap<String, usize> = sample_names
                .iter()
                .enumerate()
                .map(|(i, &x)| (x.to_string(), i))
                .collect();
            maf_data.samples = sample_names
                .iter()
                .map(|&x| SampleMafData {
                    sample_name: x.to_string(),
                    data: vec![BTreeMap::new(); output_chrom_count],
                })
                .collect();
            for sample_bytes in header.samples() {
                let sample_name = std::str::from_utf8(sample_bytes).unwrap().to_string();
                x.push(sample_map.get(&sample_name).copied());
            }

            // Check that all output samples have been found:
            let mut testx = x.iter().filter_map(|&x| x).sorted().collect::<Vec<_>>();
            testx.sort();
            for (sample_index, sample_name) in sample_names.iter().enumerate() {
                assert_eq!(
                    testx[sample_index], sample_index,
                    "Can't find sample '{sample_name}' in minor allele frequency input file: '{filename}'"
                );
            }
        }
        x
    };

    let pass_filter_id = header.name_to_id("PASS".as_bytes());

    let mut rec = reader.empty_record();
    while let Some(r) = reader.read(&mut rec) {
        unwrap!(r, "Failed to parse variant record");

        // Only use passing variant calls
        let is_pass = match pass_filter_id {
            Ok(id) => rec.has_filter(&id),
            Err(_) => rec.filters().next().is_none(),
        };
        if !is_pass {
            continue;
        }

        // Test if variant is biallelic
        if rec.allele_count() != 2 {
            continue;
        }

        // Test if AD field exists
        let ad_values = rec.format(b"AD").integer();
        if ad_values.is_err() {
            continue;
        }
        let ad_values = ad_values.unwrap();

        let rid = rec.rid().unwrap() as usize;
        let output_chrom_index = vcf_to_output_chrom_index_map[rid];
        let pos = rec.pos();

        for (vcf_sample_index, v) in ad_values.iter().enumerate() {
            if let Some(output_sample_index) = vcf_to_output_sample_index_map[vcf_sample_index] {
                // VCF records have already been filtered for exactly two alleles, so AD should
                // have 2 values unless it is given a single unknown entry ('.')
                let missing_int: i32 = bcf::record::Numeric::missing();
                if v.len() != 2 {
                    if !(v.len() == 1 && v[0] == missing_int) {
                        let rstr = rec.to_vcf_string().unwrap();
                        panic!(
                            "Failed to process allele depth from minor allele frequency file record:\n\n{rstr}"
                        );
                    }
                } else if v.iter().all(|&x| x != missing_int) {
                    maf_data.samples[output_sample_index].data[output_chrom_index]
                        .insert(pos, (v[0], v[1]));
                }
            }
        }
    }

    maf_data
}

pub fn write_maf_bigwig_file(
    filename: &Utf8Path,
    sample_maf_data: &SampleMafData,
    chrom_list: &ChromList,
) {
    info!("Writing bigwig minor allele frequency track to file: '{filename}'");

    let mut bigwig_writer = bigwig_utils::get_new_writer(filename.as_str(), chrom_list);

    for (chrom_index, chrom_name) in chrom_list.data.iter().map(|x| &x.label).enumerate() {
        let chrom_biallelic_counts = &sample_maf_data.data[chrom_index];

        // Get minor allele frequency from AD counts when sufficient counts are available
        let mut start_pos: Vec<u32> = Vec::new();
        let mut mafs: Vec<f32> = Vec::new();
        for (pos, pair_count) in chrom_biallelic_counts.iter() {
            let total = pair_count.0 + pair_count.1;
            if total < 8 {
                continue;
            }
            let maf = std::cmp::min(pair_count.0, pair_count.1) as f32 / total as f32;
            start_pos.push(*pos as u32);
            mafs.push(maf);
        }
        bigwig_writer
            .add_interval_spans(chrom_name, &mut start_pos, 1, &mut mafs)
            .unwrap();
    }
}

/// Write out a bigwig maf track file for the named sample
///
pub fn write_maf_track_files(
    output_dir: &Utf8Path,
    chrom_list: &ChromList,
    maf_data: &MafData,
    sample_name: &str,
) {
    for sample_maf_data in maf_data.samples.iter() {
        if sample_name != sample_maf_data.sample_name {
            continue;
        }
        let bigwig_filename = output_dir.join(MAF_BIGWIG_FILENAME);
        write_maf_bigwig_file(&bigwig_filename, sample_maf_data, chrom_list);
    }
}

pub fn serialize_maf_data(discover_dir: &Utf8Path, maf_data: &MafData) {
    let mut buf = Vec::new();
    maf_data
        .serialize(&mut rmp_serde::Serializer::new(&mut buf))
        .unwrap();

    let filename = discover_dir.join(MAF_MESSAGEPACK_FILENAME);

    info!("Writing minor allele frequency binary file: '{filename}'");

    unwrap!(
        std::fs::write(&filename, buf.as_slice()),
        "Unable to open and write minor allele frequency binary file: '{filename}'"
    );
}

pub fn deserialize_maf_data(discover_dir: &Utf8Path) -> MafData {
    let filename = discover_dir.join(MAF_MESSAGEPACK_FILENAME);
    let buf = unwrap!(
        std::fs::read(&filename),
        "Unable to open and read minor allele frequency binary file: '{filename}'"
    );
    rmp_serde::from_slice(&buf).unwrap()
}
