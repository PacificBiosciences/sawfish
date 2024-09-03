use std::collections::BTreeMap;
use std::path::Path;

use log::info;
use rust_htslib::bcf::{self, Read};
use rust_vc_utils::{bigwig_utils, ChromList};
use unwrap::unwrap;

use crate::genome_regions::GenomeRegions;

#[derive(Clone)]
pub struct GenomeBiallelicCounts {
    /// A vector of biallelic counts for each chromosome
    ///
    /// Each chrom is indexed in the vector by its vcf chrom index (called `rid` within htslib)
    ///
    /// For each chrom, map key is zero-indexed position, value is a tuple of the two allele counts
    pub data: Vec<BTreeMap<i64, (i32, i32)>>,
}

/// Minor allele frequency results from all samples
///
pub struct MultiSampleMafScanResult {
    pub rid_to_chrom_name: Vec<String>,
    pub sid_to_sample_name: Vec<String>,
    pub samples: Vec<GenomeBiallelicCounts>,
}

impl MultiSampleMafScanResult {
    pub fn new() -> Self {
        Self {
            rid_to_chrom_name: Vec::new(),
            sid_to_sample_name: Vec::new(),
            samples: Vec::new(),
        }
    }
}

/// Scan VCF or BCF file to create minor allele frequency track for the genome
pub fn scan_maf_file(
    filename: &str,
    excluded_regions: &GenomeRegions,
) -> Option<MultiSampleMafScanResult> {
    if filename.is_empty() {
        return None;
    }

    info!(
        "Scanning minor allele frequency data from file '{}'",
        filename
    );

    let mut maf_result = MultiSampleMafScanResult::new();

    let mut reader = bcf::Reader::from_path(filename).unwrap();
    let header = reader.header();

    // Get chrom names from header:
    let chrom_count = header.contig_count();
    for rid in 0..chrom_count {
        let chrom_bytes = header.rid2name(rid).unwrap();
        let chrom_name = std::str::from_utf8(chrom_bytes).unwrap().to_string();
        maf_result.rid_to_chrom_name.push(chrom_name);
    }

    // Get sample names from header:
    for sample_bytes in header.samples() {
        let sample_name = std::str::from_utf8(sample_bytes).unwrap().to_string();
        maf_result.sid_to_sample_name.push(sample_name);
    }
    maf_result.samples.resize(
        maf_result.sid_to_sample_name.len(),
        GenomeBiallelicCounts {
            data: vec![BTreeMap::new(); chrom_count as usize],
        },
    );

    // pre-map excluded regions to rid before scanning vcf
    let rid_to_excluded_regions = maf_result
        .rid_to_chrom_name
        .iter()
        .map(|x| excluded_regions.chroms.get(x))
        .collect::<Vec<_>>();

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
        let pos = rec.pos();

        // Test to see if record is excluded
        if let Some(regions) = rid_to_excluded_regions[rid] {
            if regions.intersect_pos(pos) {
                continue;
            }
        }

        for (sample_index, v) in ad_values.iter().enumerate() {
            assert_eq!(v.len(), 2);
            maf_result.samples[sample_index].data[rid].insert(pos, (v[0], v[1]));
        }
    }
    Some(maf_result)
}

pub fn write_maf_bigwig_file(
    filename: &Path,
    vcf_chrom_index_to_chrom_name: &[String],
    genome_biallelic_counts: &GenomeBiallelicCounts,
    chrom_list: &ChromList,
) {
    info!("Writing bigwig maf track to file: '{}'", filename.display());

    let mut bigwig_writer = bigwig_utils::get_new_writer(filename.to_str().unwrap(), chrom_list);

    // translate chrom order from VCF rid to match that in chrom_list
    let mut chrom_index_to_vcf_chrom_index = vec![None; chrom_list.data.len()];
    for (vcf_chrom_index, chrom_name) in vcf_chrom_index_to_chrom_name.iter().enumerate() {
        let val = chrom_list.label_to_index.get(chrom_name);
        if val.is_none() {
            panic!("Minor allele frequency file includes chromosome '{}', which is not found in the alignment file input", chrom_name);
        }
        let val = val.unwrap();
        assert!(chrom_index_to_vcf_chrom_index[*val].is_none());
        chrom_index_to_vcf_chrom_index[*val] = Some(vcf_chrom_index);
    }
    let chrom_index_to_vcf_chrom_index = chrom_index_to_vcf_chrom_index;

    for (chrom_index, chrom_lookup) in chrom_index_to_vcf_chrom_index.iter().enumerate() {
        let vcf_chrom_index: usize = match chrom_lookup {
            Some(x) => *x,
            None => {
                continue;
            }
        };
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        let chrom_biallelic_counts = &genome_biallelic_counts.data[vcf_chrom_index];

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

/// Write out a bigwig maf track file for every sample
pub fn write_maf_track_files(
    output_dir: &Path,
    chrom_list: &ChromList,
    maf_scan_result: &MultiSampleMafScanResult,
) {
    for (sample_index, sample_result) in maf_scan_result.samples.iter().enumerate() {
        let sample_name = maf_scan_result.sid_to_sample_name[sample_index].as_str();
        let bigwig_filename = output_dir.join(sample_name.to_owned() + ".maf.bw");
        write_maf_bigwig_file(
            &bigwig_filename,
            &maf_scan_result.rid_to_chrom_name,
            sample_result,
            chrom_list,
        );
    }
}
