use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use log::info;
use rust_vc_utils::ChromList;
use unwrap::unwrap;

use crate::breakpoint::BreakpointCluster;
use crate::discover::ASSEMBLY_REGIONS_FILENAME;
use crate::genome_segment::{GenomeSegment, IntRange};

fn get_max_insert_from_breakpoint_cluster(bpc: &BreakpointCluster) -> usize {
    if let Some(consolidated_region) = &bpc.consolidated_assembly_segment {
        consolidated_region.insert_info.max_size()
    } else {
        bpc.breakpoint.insert_info.max_size()
    }
}

pub fn get_assembly_regions_from_breakpoint_cluster(bpc: &BreakpointCluster) -> Vec<GenomeSegment> {
    if let Some(consolidated_region) = &bpc.consolidated_assembly_segment {
        vec![consolidated_region.segment.clone()]
    } else if let Some(large_insertion_info) = &bpc.large_insertion_info {
        vec![large_insertion_info.segment.clone()]
    } else {
        bpc.assembly_segments.clone()
    }
}

struct AssemblyRegionRecord {
    segment: GenomeSegment,
    cluster_index: usize,
    assembly_region_id: usize,
    consolidated_ids: Vec<usize>,

    /// A record that has been consolidated and will therefore be skipped as an individual item
    is_consolidated: bool,

    /// Largest insertion on this or the consolidated cluster set
    max_insert: usize,

    /// True if this assembly region is associated with an unpaired breakend
    unpaired_breakend: bool,
}

pub fn write_assembly_regions_to_bed(
    output_dir: &Path,
    chrom_list: &ChromList,
    clusters: &[BreakpointCluster],
) {
    // First pull all the records together for sorting:
    let mut arecs = Vec::new();

    for (cluster_index, bpc) in clusters.iter().enumerate() {
        let assembly_regions = get_assembly_regions_from_breakpoint_cluster(bpc);
        let max_insert = get_max_insert_from_breakpoint_cluster(bpc);
        let unpaired_breakend = bpc.breakpoint.breakend2.is_none();
        for (segment_index, segment) in assembly_regions.into_iter().enumerate() {
            arecs.push(AssemblyRegionRecord {
                segment,
                cluster_index,
                assembly_region_id: segment_index,
                consolidated_ids: if let Some(x) = &bpc.consolidated_assembly_segment {
                    x.cluster_ids.clone()
                } else {
                    Vec::new()
                },
                is_consolidated: bpc.consolidated_cluster_index.is_some(),
                max_insert,
                unpaired_breakend,
            });
        }
    }

    arecs.sort_by(|a, b| a.segment.partial_cmp(&b.segment).unwrap());

    let filename = output_dir.join(ASSEMBLY_REGIONS_FILENAME);

    info!(
        "Writing debug assembly region bed file: '{}'",
        filename.display()
    );

    let f = unwrap!(
        File::create(&filename),
        "Unable to create debug assembly region bed file: '{}'",
        filename.display()
    );
    let mut f = BufWriter::new(f);
    writeln!(f, "#gffTags").unwrap();
    for arec in arecs {
        let chrom_label = &chrom_list.data[arec.segment.chrom_index].label;
        let range = &arec.segment.range;
        writeln!(
            f,
            "{chrom_label}\t{}\t{}\tName={};AssemblyRegionID={};ConsolidatedIds={:?};IsConsolidated={};MaxInsert={}{}",
            range.start,
            range.end,
            arec.cluster_index,
            arec.assembly_region_id,
            arec.consolidated_ids,
            arec.is_consolidated,
            arec.max_insert,
            if arec.unpaired_breakend { ";UnpairedBND" } else { "" },
        )
            .unwrap();
    }
}

/// Return a vector with index matched to cluster_index, and value matching the
/// get_assembly_regions_from_breakpoint_cluster() function for the given cluster index
///
pub fn read_assembly_regions_from_bed(
    discover_dir: &Path,
    chrom_list: &ChromList,
) -> Vec<Vec<GenomeSegment>> {
    use rust_htslib::bgzf;
    use std::io::Read;

    let filename = discover_dir.join(ASSEMBLY_REGIONS_FILENAME);

    /*
    info!(
        "Reading assembly regions from bed file: '{}'",
        filename.display()
    );
     */

    let mut reader = unwrap!(
        bgzf::Reader::from_path(&filename),
        "Unable to open assembly regions file: '{}'",
        filename.display()
    );

    let mut content = String::new();
    unwrap!(
        reader.read_to_string(&mut content),
        "Can't parse text from assembly regions file: '{}'",
        filename.display()
    );

    let mut regions = Vec::new();

    let mut last_chrom_index = 0;
    let mut last_chrom = "";
    for line in content.split('\n') {
        // The last line is expected to be empty
        if line.is_empty() {
            break;
        }

        if line.starts_with('#') {
            continue;
        }

        let words = line.split('\t').collect::<Vec<_>>();

        assert!(words.len() >= 4);
        let chrom = words[0];
        let start = words[1].parse::<i64>().unwrap();
        let end = words[2].parse::<i64>().unwrap();
        let keyval_string = words[3];

        if chrom != last_chrom {
            last_chrom = chrom;
            last_chrom_index = *unwrap!(
                chrom_list.label_to_index.get(chrom),
                "assembly regions include unknown chromosome name` `{chrom}` in file: '{}'",
                filename.display()
            );
        }

        let mut cluster_index = None;
        let keyvals = keyval_string.split(';').collect::<Vec<_>>();
        for keyval in keyvals {
            let x = keyval.split('=').collect::<Vec<_>>();
            assert_eq!(x.len(), 2);
            if x[0] == "Name" {
                cluster_index = Some(x[1].parse::<usize>().unwrap());
                break;
            }
        }
        let cluster_index = cluster_index.unwrap();

        let assembly_region = GenomeSegment {
            chrom_index: last_chrom_index,
            range: IntRange::from_pair(start, end),
        };

        if cluster_index >= regions.len() {
            regions.resize(cluster_index + 1, Vec::new());
        }

        regions[cluster_index].push(assembly_region);
    }

    regions
}
