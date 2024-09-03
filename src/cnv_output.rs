use std::collections::HashMap;
use std::path::Path;

use log::info;
use rust_vc_utils::ChromList;
use unwrap::unwrap;

use crate::cli;
use crate::copy_number_segmentation::{CNState, SampleCopyNumberSegments};
use crate::genome_regions::GenomeRegions;
use crate::vcf_utils;

/// Write out a bedgraph track for copy number segments
pub fn write_copy_number_segment_file(
    output_dir: &Path,
    sample_cn_segments: &SampleCopyNumberSegments,
) {
    use std::fs::File;
    use std::io::{BufWriter, Write};

    let filename = output_dir.join("copynum.bedgraph");

    info!(
        "Writing bedgraph copy number track to file: '{}'",
        filename.display()
    );

    let f = unwrap!(
        File::create(&filename),
        "Unable to create bedgraph copy number track file: '{}'",
        filename.display()
    );
    let mut f = BufWriter::new(f);

    let chrom_list = &sample_cn_segments.chrom_list;
    let chrom_count = sample_cn_segments.chrom_list.data.len();

    for chrom_index in 0..chrom_count {
        let chrom_label = &chrom_list.data[chrom_index].label;
        let chrom_cn_segments = &sample_cn_segments.data[chrom_index];
        for s in chrom_cn_segments.iter() {
            if s.state == CNState::Unknown {
                continue;
            }
            writeln!(
                f,
                "{}\t{}\t{}\t{}",
                chrom_label, s.begin, s.end, s.state as u32
            )
            .unwrap();
        }
    }
}

#[allow(dead_code)]
// this is just temporary until we get a system for denoting sex expected CN
fn estimate_sex(sample_cn_segments: &SampleCopyNumberSegments, chrom_list: &ChromList) -> u32 {
    let chrom_count = chrom_list.data.len();
    for chrom_index in 0..chrom_count {
        let chrom_label = &chrom_list.data[chrom_index].label;
        if chrom_label == "chrX" {
            let chrom_cn_segments = &sample_cn_segments.data[chrom_index];
            let state_count = CNState::Unknown as usize;
            let mut bin_counts = vec![0; state_count];
            for s in chrom_cn_segments.iter() {
                if s.state == CNState::Unknown {
                    continue;
                }
                bin_counts[s.state as usize] += s.end - s.begin;
            }
            if bin_counts[2] >= bin_counts[1] {
                return 2;
            } else {
                return 1;
            }
        }
    }
    // no chrX detected, so we will just default to 2
    2
}

fn get_expected_cn(expected_cn_regions: &GenomeRegions, chrom: &str, begin: i64, end: i64) -> u32 {
    let default_expected_cn: u32 = 2;
    match expected_cn_regions.find_overlaps(chrom, begin, end) {
        Some(all_overlaps) => {
            //start with assuming everything is 2 copy, then remove non-2s
            let mut cn_lookup: HashMap<u32, i64> = HashMap::<u32, i64>::new();
            cn_lookup.insert(2, end - begin);

            for overlaps in all_overlaps {
                let cn = overlaps.data();
                if *cn as u32 != default_expected_cn {
                    let interval = overlaps.interval();
                    let max_start = std::cmp::max(begin, interval.start);
                    let min_end = std::cmp::min(end, interval.end);
                    let overlap = min_end - max_start;
                    assert!(overlap >= 0);

                    *cn_lookup.entry(2).or_insert(0) -= overlap;
                    *cn_lookup.entry(*cn as u32).or_insert(0) += overlap;
                }
            }

            //iterate over the pairs, pick the maximum by the value, return the CN value "k"
            cn_lookup
                .iter()
                .max_by(|a, b| a.1.cmp(b.1))
                .map(|(&k, _v)| k)
                .unwrap()
        }
        None => default_expected_cn,
    }
}

fn write_cnv_vcf_file(
    settings: &cli::DiscoverSettings,
    filename: &Path,
    sample_cn_segments: &SampleCopyNumberSegments,
    expected_cn_regions: &GenomeRegions,
) {
    use rust_htslib::bcf::record::{GenotypeAllele, Record};
    use rust_htslib::bcf::{Format, Writer};

    let bin_size = settings.depth_bin_size;
    let chrom_list = &sample_cn_segments.chrom_list;

    //first we need to build out the VCF header
    let mut header = vcf_utils::get_basic_vcf_header(
        &settings.ref_filename,
        chrom_list,
        &[&sample_cn_segments.sample_name],
    );

    //add any relevant filters
    header.push_record(r#"##FILTER=<ID=PASS,Description="All filters passed">"#.as_bytes());
    let target_size = 100000;
    header.push_record(
        r#"##FILTER=<ID=TARGET_SIZE,Description="Call is smaller than the target size of 100kbp">"#
            .as_bytes(),
    );

    //add in all our header info
    //TODO: add reference file somewhere? format below
    //header.push_record(r#"##reference=file://reference/human_GRCh38_no_alt_analysis_set.fasta"#.as_bytes());
    header.push_record(
        r#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">"#
            .as_bytes(),
    );
    header.push_record(
        r#"##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">"#
            .as_bytes(),
    );
    header.push_record(r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">"#.as_bytes());
    header.push_record(
        r#"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">"#
            .as_bytes(),
    );
    header.push_record(
        r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">"#.as_bytes(),
    );
    header.push_record(
        r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">"#.as_bytes(),
    );
    header.push_record(r#"##ALT=<ID=DEL,Description="Deletion">"#.as_bytes());
    header.push_record(r#"##ALT=<ID=DUP,Description="Duplication">"#.as_bytes());

    header.push_record(r#"##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number genotype for imprecise events">"#.as_bytes());
    header
        .push_record(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#.as_bytes());

    // Write uncompressed VCF to stdout with above header and get an empty record
    let mut vcf: Writer = Writer::from_path(filename, &header, false, Format::Vcf).unwrap();

    /*
    // old original approach for deriving sex
    let default_expected_cn = 2;
    let expected_sex_cn = estimate_sex(sample_cn_segments, chrom_list);
    */

    //get filter ids ready for use
    let pass_id = vcf.header().name_to_id("PASS".as_bytes()).unwrap();
    let filter_target_size_id = vcf.header().name_to_id("TARGET_SIZE".as_bytes()).unwrap();

    let chrom_count = chrom_list.data.len();
    for chrom_index in 0..chrom_count {
        let chrom_label = &chrom_list.data[chrom_index].label;
        let chrom_cn_segments = &sample_cn_segments.data[chrom_index];
        let chrom_rid = vcf.header().name2rid(chrom_label.as_bytes()).unwrap();

        /*
        // this was an old approach based on sex estimates
        let expected_cn: u32 = if chrom_label == "chrX" || chrom_label == "chrY" {
            expected_sex_cn
        } else {
            default_expected_cn
        };
        */

        for s in chrom_cn_segments.iter() {
            if s.state == CNState::Unknown {
                continue;
            }
            //TODO: this occasionally happens, I think for chromosomes shorter than the bin size; should verify
            if s.begin == s.end {
                continue;
            }

            //calculate expected copy number for the region
            let expected_cn = get_expected_cn(expected_cn_regions, chrom_label, s.begin, s.end);
            if s.state as u32 == expected_cn {
                continue;
            }

            let mut record: Record = vcf.empty_record();
            record.set_rid(Some(chrom_rid));
            record.set_pos(s.begin);
            let sv_len: i64 = s.end - s.begin;

            //set REF and ALT
            let svtype = if (s.state as u32) < expected_cn {
                "DEL"
            } else {
                "DUP"
            };
            record
                .set_alleles(&["N".as_bytes(), format!("<{}>", svtype).as_bytes()])
                .unwrap();

            //TODO: placeholder fixed value, make QUAL meaningful long-term if we can
            let max_qual: f32 = 100.0;
            let quality: f32 = s.qual.floor().min(max_qual);
            record.set_qual(quality);

            /*
            TODO: do we have any filters we want to add? some ideas:
            - mostly excluded, e.g. if the event has significant overlap with an exclude region
            - common calls - would likely require additional input; probably better separate IMO
            */
            let mut filters = Vec::new();
            //check if it's small
            if sv_len < target_size {
                filters.push(&filter_target_size_id);
            }
            if filters.is_empty() {
                //if no other filters are used, add in PASS
                filters.push(&pass_id);
            }
            record.set_filters(&filters).unwrap();

            record.push_info_flag("IMPRECISE".as_bytes()).unwrap();
            record
                .push_info_string("SVTYPE".as_bytes(), &[svtype.as_bytes()])
                .unwrap();
            record
                .push_info_integer("END".as_bytes(), &[s.end as i32])
                .unwrap();
            record
                .push_info_integer("SVLEN".as_bytes(), &[sv_len as i32])
                .unwrap();

            //for now, set the confidence intervals as +- bin size
            record
                .push_info_string(
                    "CIPOS".as_bytes(),
                    &[format!("-{},{}", bin_size, bin_size).as_bytes()],
                )
                .unwrap();
            record
                .push_info_string(
                    "CIEND".as_bytes(),
                    &[format!("-{},{}", bin_size, bin_size).as_bytes()],
                )
                .unwrap();

            //TODO: do we care about genotypes other than 0/1?
            let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];
            record.push_genotypes(alleles).unwrap();
            record
                .push_format_integer("CN".as_bytes(), &[s.state as i32])
                .unwrap();

            //save everything
            vcf.write(&record).unwrap()
        }
    }
}

/// Generates a single-sample CNV VCF
///
/// # Arguments
/// * `output_prefix` - file is saved to "{prefix}.{sample_name}.vcf.gz"
/// * `copy_number_segment_result` - the copy number segments for the sample
///
pub fn write_indexed_cnv_vcf_file(
    shared_settings: &cli::SharedSettings,
    settings: &cli::DiscoverSettings,
    copy_number_segment_result: &SampleCopyNumberSegments,
    expected_cn: &GenomeRegions,
) {
    // Write VCF file
    let filename = settings.output_dir.join("cnv.vcf.gz");
    info!(
        "Writing copy number variants to file: '{}'",
        filename.display()
    );
    write_cnv_vcf_file(settings, &filename, copy_number_segment_result, expected_cn);

    // Index VCF file
    vcf_utils::build_bcf_index(&filename, shared_settings.thread_count, true).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_expected_cn() {
        //no overlaps allowed for these
        let mut genome_regions = GenomeRegions::new(false);
        genome_regions.add_region_value("chr1", 10, 20, 1);
        genome_regions.add_region_value("chr1", 30, 40, 3);
        genome_regions.add_region_value("chr1", 50, 60, 0);
        genome_regions.add_region_value("chr1", 60, 70, 1);

        //full overlaps
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 0, 10), 2);
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 10, 20), 1);
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 30, 40), 3);

        //straddling an implicit border
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 5, 12), 2);
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 8, 15), 1);

        //straddling an explicit border
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 55, 64), 0);
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 56, 65), 1);

        //out of bounds or a completely un-set chromosome
        assert_eq!(get_expected_cn(&genome_regions, "chr1", 100, 120), 2);
        assert_eq!(get_expected_cn(&genome_regions, "chr2", 100, 120), 2);
    }
}
