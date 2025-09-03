use rust_htslib::bam;

use camino::Utf8Path;
use rust_htslib::bcf::{self, Read};
use rust_vc_utils::aux::{get_string_aux_tag, is_aux_tag_found};
use rust_vc_utils::{ChromList, rev_comp_in_place};
use scanf::sscanf;
use unwrap::unwrap;

use crate::bam_sa_parser::get_seq_order_read_split_segments;
use crate::bam_utils::get_simplified_dna_seq;
use crate::breakpoint::{Breakend, BreakendDirection, BreakendNeighbor, Breakpoint, InsertInfo};
use crate::contig_output::{CONTIG_AUX_TAG, SA_AUX_TAG};
use crate::expected_ploidy::get_max_haplotype_count_for_regions;
use crate::filenames::{CANDIDATE_SV_FILENAME, CONTIG_ALIGNMENT_FILENAME};
use crate::genome_regions::{GenomeRegions, GenomeRegionsByChromIndex};
use crate::genome_segment::GenomeSegment;
use crate::int_range::IntRange;
use crate::large_variant_output::{
    BND0_NEIGHBOR_INFO_KEY, BND1_NEIGHBOR_INFO_KEY, CONTIG_POS_INFO_KEY, OVERLAP_ASM_INFO_KEY,
};
use crate::refine_sv::RefinedSV;
use crate::refine_sv::assembly_regions::read_assembly_regions_from_bed;
use crate::simple_alignment::SimpleAlignment;
use crate::sv_group::{
    ClusterAssembly, ClusterAssemblyAlignment, GroupHaplotypeId, SVGroup, SVGroupHaplotype,
};
use crate::sv_id::get_sv_id_from_label;

/// All contigs associated with the same cluster_index
type ClusterAssemblies = Vec<ClusterAssembly>;

/// All contigs from a sample
type SampleAssemblies = Vec<ClusterAssemblies>;

/// Extract supporting read count from an assembly contig following standard sawfish contig BAM
/// formatting conventions.
///
fn get_supporting_read_count(record: &bam::Record) -> usize {
    let contig_tag = get_string_aux_tag(record, CONTIG_AUX_TAG);
    let words = contig_tag.split(';').collect::<Vec<_>>();
    let mut x = 0;
    for word in words {
        if word.is_empty() {
            continue;
        }
        let kv = word.split(':').collect::<Vec<_>>();
        assert_eq!(kv.len(), 2);
        if kv[0] == "n_reads" {
            x = kv[1].parse::<usize>().unwrap();
            break;
        }
    }
    x
}

fn get_high_quality_contig_range(record: &bam::Record) -> Option<IntRange> {
    let contig_tag = get_string_aux_tag(record, CONTIG_AUX_TAG);
    let words = contig_tag.split(';').collect::<Vec<_>>();
    let mut x = None;
    for word in words {
        if word.is_empty() {
            continue;
        }
        let kv = word.split(':').collect::<Vec<_>>();
        assert_eq!(kv.len(), 2);
        if kv[0] == "hq_range" {
            let range_vals = kv[1].split('-').collect::<Vec<_>>();
            assert_eq!(range_vals.len(), 2);
            let start = range_vals[0].parse::<i64>().unwrap();
            let end = range_vals[1].parse::<i64>().unwrap();
            x = Some(IntRange::from_pair(start, end));
            break;
        }
    }
    x
}

/// Parse BAM record from sawfish contig alignment file into the sample_assembles structure to use in haplotype merging
/// and joint GT steps
///
/// Results are accumulated in the sample_assemblies structure. The primary index for this struct is cluster_index, so
/// that it can easily be looked up from any SV cluster_index.
///
fn update_sample_assemblies_from_bam_record(
    chrom_list: &ChromList,
    record: &bam::Record,
    sample_assemblies: &mut SampleAssemblies,
) {
    // Only read the primary alignment from each split record
    //
    // The convention used in the contig bam is that the primary alignment corresponds to breakend1 and the
    // supplementary alignment is used for breakend2. We can setup the contig alignment information for breakend2
    // from the SA tags on the primary read, so we shouldn't need to directly read in the supplementary alignment.
    //
    if ((record.flags() as u32) & rust_htslib::htslib::BAM_FSUPPLEMENTARY) != 0 {
        return;
    }

    // Parse indexes from the qname:
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    let words = qname.split(':').collect::<Vec<_>>();
    assert!(words.len() >= 4);

    // Ignore sample index. The joint call routine will add its own sample indexes based on the input order
    // of the discover directories
    //
    //let sample_index = words[1].parse::<usize>().unwrap();

    let cluster_index = words[2].parse::<usize>().unwrap();
    let assembly_index = words[3].parse::<usize>().unwrap();

    if sample_assemblies.len() <= cluster_index {
        sample_assemblies.resize(cluster_index + 1, vec![ClusterAssembly::default()]);
    }

    let cluster_assemblies = &mut sample_assemblies[cluster_index];

    if cluster_assemblies.len() <= assembly_index {
        cluster_assemblies.resize(assembly_index + 1, ClusterAssembly::default());
    }

    let cluster_assembly = &mut cluster_assemblies[assembly_index];
    cluster_assembly.supporting_read_count = get_supporting_read_count(record);

    let high_quality_contig_range = get_high_quality_contig_range(record).unwrap();

    let primary_segment_alignment = {
        let contig_seq = get_simplified_dna_seq(record);
        let contig_alignment = SimpleAlignment {
            ref_offset: record.pos(),
            cigar: record.cigar().to_vec(),
        };
        ClusterAssemblyAlignment {
            contig_seq,
            chrom_index: record.tid() as usize,
            contig_alignment,
            high_quality_contig_range: high_quality_contig_range.clone(),
            is_fwd_strand: true, // The contig output format requires fwd-strand for the primary alignment
        }
    };

    cluster_assembly
        .contig_alignments
        .push(primary_segment_alignment);

    if is_aux_tag_found(record, SA_AUX_TAG) {
        for segment in get_seq_order_read_split_segments(chrom_list, record)
            .into_iter()
            .filter(|x| !x.from_primary_bam_record)
        {
            let mut contig_seq = get_simplified_dna_seq(record);
            let mut high_quality_contig_range = high_quality_contig_range.clone();
            if !segment.is_fwd_strand {
                rev_comp_in_place(&mut contig_seq);
                high_quality_contig_range.reverse(contig_seq.len() as i64);
            }

            let contig_alignment = SimpleAlignment {
                ref_offset: segment.pos,
                cigar: segment.cigar.to_vec(),
            };

            let supp_segment_alignment = ClusterAssemblyAlignment {
                contig_seq,
                chrom_index: segment.chrom_index,
                contig_alignment,
                high_quality_contig_range,
                is_fwd_strand: segment.is_fwd_strand,
            };
            cluster_assembly
                .contig_alignments
                .push(supp_segment_alignment);
        }
    }
}

fn get_candidate_sv_contig_info(
    discover_dir: &Utf8Path,
    chrom_list: &ChromList,
) -> SampleAssemblies {
    use rust_htslib::bam::Read;

    let contig_filename = discover_dir.join(CONTIG_ALIGNMENT_FILENAME);
    if !contig_filename.exists() {
        panic!("Can't find contig alignment file expected at {contig_filename}");
    }

    let mut sample_assemblies = Vec::new();
    let mut bam_reader = bam::Reader::from_path(contig_filename).unwrap();
    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");
        update_sample_assemblies_from_bam_record(chrom_list, &record, &mut sample_assemblies);
    }
    sample_assemblies
}

struct BndAltInfo {
    breakend1_direction: BreakendDirection,
    breakend2_chrom_index: usize,
    breakend2_pos: i64,
    breakend2_direction: BreakendDirection,
    insert_seq: Vec<u8>,
}

/// Parse a VCF breakend record alt-allele string into component parts describing the breakpoint
///
/// Example BND alt allele string is "TCT]chr1:500]", see VCF spec for full description.
///
/// # Arguments
/// * `homlen` - this is used to correct the remote breakend into a left-shifted position when the breakpoint is inverted
///
fn parse_bnd_alt_allele(chrom_list: &ChromList, alt_allele: &[u8], homlen: i64) -> BndAltInfo {
    let alt_allele = std::str::from_utf8(alt_allele).unwrap();

    let breakend2_direction = if alt_allele.contains(']') {
        BreakendDirection::LeftAnchor
    } else if alt_allele.contains('[') {
        BreakendDirection::RightAnchor
    } else {
        panic!("Invalid BND alt alt allele '{alt_allele}'");
    };

    let split_char = match breakend2_direction {
        BreakendDirection::LeftAnchor => ']',
        BreakendDirection::RightAnchor => '[',
    };

    let words = alt_allele.split(split_char).collect::<Vec<_>>();
    assert_eq!(
        words.len(),
        3,
        "Unexpected BND alt allele format `{alt_allele}`"
    );

    let (breakend2_chrom_index, breakend2_raw_pos) = {
        // Note that rsplitn orders words in reverse order compared to how they appear in the string:
        let mate_location = words[1].rsplitn(2, ':').collect::<Vec<_>>();
        assert_eq!(
            mate_location.len(),
            2,
            "Unexpected BND alt allele format `{alt_allele}`"
        );

        let breakend2_chrom_name = mate_location[1];
        let breakend2_chrom_index = *chrom_list.label_to_index.get(breakend2_chrom_name).unwrap();
        let breakend2_raw_pos = mate_location[0].parse::<i64>().unwrap() - 1;
        (breakend2_chrom_index, breakend2_raw_pos)
    };

    let breakend1_direction = if !words[0].is_empty() && words[2].is_empty() {
        BreakendDirection::LeftAnchor
    } else if words[0].is_empty() && !words[2].is_empty() {
        BreakendDirection::RightAnchor
    } else {
        panic!("Invalid BND alt alt allele '{alt_allele}'");
    };

    let insert_seq = match breakend1_direction {
        BreakendDirection::LeftAnchor => words[0].as_bytes()[1..].to_vec(),
        BreakendDirection::RightAnchor => {
            let len = words[2].len();
            words[2].as_bytes()[..len - 1].to_vec()
        }
    };

    // Correct breakend2 raw pos
    let breakend2_pos = {
        let mut x = breakend2_raw_pos;
        let is_same_breakend_orientation = breakend1_direction == breakend2_direction;
        if is_same_breakend_orientation {
            x -= homlen;
        }
        if breakend2_direction == BreakendDirection::RightAnchor {
            x -= 1;
        }
        x
    };

    BndAltInfo {
        breakend1_direction,
        breakend2_chrom_index,
        breakend2_pos,
        breakend2_direction,
        insert_seq,
    }
}

/// RefinedSV plus auxiliary data which will be added to the SVGroup
struct AnnotatedRefinedSV {
    refined_sv: RefinedSV,

    /// Assembly index or indexes used to associate the contig back to the SV call
    ///
    /// This provides a way to clarify which contigs are associated with the variant in overlapping contig regions.
    ///
    contig_map: Vec<usize>,
}

fn parse_breakend_neighbor_info(rec: &bcf::Record, info_key: &str) -> Option<BreakendNeighbor> {
    // A BcfUndefinedTag error on this INFO key lookup is expected from older versions of discover
    //
    let x = match rec.info(info_key.as_bytes()).string() {
        Ok(x) => x,
        Err(err) => match &err {
            rust_htslib::errors::Error::BcfUndefinedTag { .. } => None,
            _ => panic!(
                "Unexpected error reading INFO key '{}': {:?}",
                info_key, err
            ),
        },
    }?;

    let val_bytes = x.first().unwrap();
    let val_str = std::str::from_utf8(val_bytes).unwrap();

    let mut cluster_index = 0usize;
    let mut breakend_index = 0usize;
    unwrap!(
        sscanf!(val_str, "{cluster_index}:{breakend_index}"),
        "Unexpected value for candidate VCF INFO key '{info_key}': '{val_str}', from VCF record:\n{}",
        rec.to_vcf_string().unwrap_or(String::from("ERROR"))
    );
    Some(BreakendNeighbor {
        sample_index: 0,
        cluster_index,
        breakend_index,
    })
}

/// Process a single SV candidate record from VCF into a refined SV -- this is the first step of converting it
/// into an SVGroup
///
fn process_bcf_record_to_anno_refine_sv(
    chrom_list: &ChromList,
    all_cluster_assembly_regions: &[Vec<GenomeSegment>],
    contigs: &SampleAssemblies,
    rid_to_chrom_index: &[usize],
    rec: &bcf::Record,
) -> Option<AnnotatedRefinedSV> {
    let rid = rec.rid().unwrap() as usize;
    let chrom_index = rid_to_chrom_index[rid];
    let pos = rec.pos();
    let (id, is_breakend1) =
        get_sv_id_from_label(std::str::from_utf8(rec.id().as_slice()).unwrap());

    // Assert that locus is biallelic
    let allele_count = rec.allele_count();
    if allele_count != 2 {
        panic!("Unexpected allele count in sample candidate SV input {allele_count}",);
    }

    let alt_allele = rec.alleles()[1];

    let svtype = std::str::from_utf8(
        rec.info(b"SVTYPE")
            .string()
            .unwrap()
            .unwrap()
            .first()
            .unwrap(),
    )
    .unwrap()
    .to_string();

    let homlen = {
        let x = rec.info(b"HOMLEN").integer().unwrap();
        if let Some(x) = x {
            x.first().cloned().unwrap()
        } else {
            0
        }
    } as i64;

    let homology_seq = {
        let x = rec.info(b"HOMSEQ").string().unwrap();
        if let Some(x) = x {
            x.first().unwrap().to_vec()
        } else {
            Vec::new()
        }
    };

    assert_eq!(homlen, homology_seq.len() as i64);

    let contig_pos = *rec
        .info(CONTIG_POS_INFO_KEY.as_bytes())
        .integer()
        .unwrap()
        .unwrap()
        .first()
        .unwrap() as usize;

    let cluster_assembly = &contigs[id.cluster_index][id.assembly_index];

    let assembly_regions = all_cluster_assembly_regions[id.cluster_index].clone();

    let overlapping_assembly_indexes = {
        let x = rec.info(OVERLAP_ASM_INFO_KEY.as_bytes()).integer().unwrap();
        if let Some(x) = x {
            x.iter().map(|&x| x as usize).collect()
        } else {
            Vec::new()
        }
    };

    let contig_map = {
        if overlapping_assembly_indexes.is_empty() {
            vec![id.assembly_index]
        } else {
            overlapping_assembly_indexes
        }
    };

    if svtype == "DEL" || svtype == "INS" {
        assert!(is_breakend1.is_none());

        let end = *rec
            .info(b"END")
            .integer()
            .unwrap()
            .unwrap()
            .first()
            .unwrap() as i64
            - 1;

        let breakend1 = Breakend {
            segment: GenomeSegment {
                chrom_index,
                range: IntRange::from_pair(pos, pos + homlen + 1),
            },
            dir: BreakendDirection::LeftAnchor,
        };
        let breakend2 = Some(Breakend {
            segment: GenomeSegment {
                chrom_index,
                range: IntRange::from_pair(end, end + homlen + 1),
            },
            dir: BreakendDirection::RightAnchor,
        });

        let insert_seq = if alt_allele == b"<DEL>" {
            let x = rec.info(b"INSSEQ").string().unwrap();
            if let Some(x) = x {
                x.first().unwrap().to_vec()
            } else {
                Vec::new()
            }
        } else {
            alt_allele[1..].to_vec()
        };

        let insert_info = if insert_seq.is_empty() {
            InsertInfo::NoInsert
        } else {
            InsertInfo::Seq(insert_seq)
        };

        let bp = Breakpoint {
            breakend1,
            breakend2,
            insert_info,
            is_precise: true,
            ..Default::default()
        };

        let refined_sv = RefinedSV {
            id,
            bp,
            breakend1_homology_seq: homology_seq,
            single_region_refinement: cluster_assembly.single_region_refinement(),
            contig_pos_before_breakend1: contig_pos,
            assembly_regions,
            ext: Default::default(),

            // All items below are set by the scoring routine, so we just default them here
            //
            score: Default::default(),
        };

        Some(AnnotatedRefinedSV {
            refined_sv,
            contig_map,
        })
    } else if svtype == "BND" {
        // In theory we should be able to pull the full rsv out of breakend1
        assert!(is_breakend1.is_some());
        if !is_breakend1.unwrap() {
            return None;
        }

        // parse the alt allele
        let bnd_alt_info = parse_bnd_alt_allele(chrom_list, alt_allele, homlen);

        let breakend1_pos = match bnd_alt_info.breakend1_direction {
            BreakendDirection::LeftAnchor => pos,
            BreakendDirection::RightAnchor => pos - 1,
        };
        let breakend1 = Breakend {
            segment: GenomeSegment {
                chrom_index,
                range: IntRange::from_pair(breakend1_pos, breakend1_pos + homlen + 1),
            },
            dir: bnd_alt_info.breakend1_direction,
        };

        let breakend2_pos = bnd_alt_info.breakend2_pos;
        let breakend2 = Some(Breakend {
            segment: GenomeSegment {
                chrom_index: bnd_alt_info.breakend2_chrom_index,
                range: IntRange::from_pair(breakend2_pos, breakend2_pos + homlen + 1),
            },
            dir: bnd_alt_info.breakend2_direction,
        });

        let insert_info = if bnd_alt_info.insert_seq.is_empty() {
            InsertInfo::NoInsert
        } else {
            InsertInfo::Seq(bnd_alt_info.insert_seq)
        };

        let breakend1_neighbor = parse_breakend_neighbor_info(rec, BND0_NEIGHBOR_INFO_KEY);
        let breakend2_neighbor = parse_breakend_neighbor_info(rec, BND1_NEIGHBOR_INFO_KEY);

        let bp = Breakpoint {
            breakend1,
            breakend2,
            insert_info,
            is_precise: true,
            breakend1_neighbor,
            breakend2_neighbor,
        };

        let refined_sv = RefinedSV {
            id,
            bp,
            breakend1_homology_seq: homology_seq,
            single_region_refinement: cluster_assembly.single_region_refinement(),
            contig_pos_before_breakend1: contig_pos,
            assembly_regions,
            ext: Default::default(),

            // All items below are set by the scoring routine, so we just default them here
            //
            score: Default::default(),
        };
        Some(AnnotatedRefinedSV {
            refined_sv,
            contig_map,
        })
    } else {
        panic!("Unexpected SV type: {svtype}");
    }
}

/// Attempt to convert the full set of AnnotatedRefinedSVs from the same cluster_index into a single-sample SVGroup
///
/// The group_arsvs will always be cleared by this method.
///
fn process_rsv_cluster_to_sv_group(
    treat_single_copy_as_haploid: bool,
    contigs: &SampleAssemblies,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    group_arsvs: &mut Vec<AnnotatedRefinedSV>,
) -> Option<SVGroup> {
    assert!(!group_arsvs.is_empty());

    let first_arsv = &group_arsvs[0];
    let first_rsv = &first_arsv.refined_sv;
    let group_regions = first_rsv.assembly_regions.clone();
    let cluster_index = first_rsv.id.cluster_index;
    let mut group_haplotypes = contigs[cluster_index]
        .iter()
        .enumerate()
        .map(|(assembly_index, x)| SVGroupHaplotype {
            hap_id: GroupHaplotypeId {
                sample_index: 0,
                cluster_index,
                assembly_index,
            },
            contig_info: x.clone(),
        })
        .collect::<Vec<_>>();

    let (max_haplotype_count, expected_cn_info) = get_max_haplotype_count_for_regions(
        treat_single_copy_as_haploid,
        expected_copy_number_regions,
        &group_regions,
    );

    let debug = false;
    if debug {
        eprintln!("max_haplotype_count {max_haplotype_count}");
        eprintln!("contig_map {:?}", &first_arsv.contig_map);
    }

    let sample_expected_cn_info = vec![expected_cn_info];
    let mut sample_haplotype_list = vec![
        first_arsv
            .contig_map
            .iter()
            .take(max_haplotype_count)
            .copied()
            .collect::<Vec<_>>(),
    ];

    let mut tmp_arsvs = Vec::new();
    std::mem::swap(group_arsvs, &mut tmp_arsvs);
    let mut refined_svs = tmp_arsvs
        .into_iter()
        .map(|x| x.refined_sv)
        .collect::<Vec<_>>();
    let mut sv_haplotype_map = refined_svs
        .iter()
        .map(|x| x.id.assembly_index)
        .collect::<Vec<_>>();

    // Remove unused haplotypes -- the only haplotypes we keep should be those in the sample
    // haplotype list. Remap sample and sv haplotype lists to the new index values.
    //
    if sample_haplotype_list[0].len() < group_haplotypes.len() {
        if debug {
            eprintln!("BEFORE");
            eprintln!("sample_haplotype_list[0] {:?}", sample_haplotype_list[0]);
            eprintln!(
                "group_haplotypes {:?}",
                group_haplotypes
                    .iter()
                    .map(|x| &x.hap_id)
                    .collect::<Vec<_>>()
            );
            eprintln!("sv_haplotype_map {sv_haplotype_map:?}");
            eprintln!(
                "refined_svs {:?}",
                refined_svs.iter().map(|x| &x.id).collect::<Vec<_>>()
            );
        }

        // 1. Get mapping from old to new haplotype indexes:
        let haplotype_index_remap = {
            let mut x = vec![None; group_haplotypes.len()];
            for (new_haplotype_index, haplotype_index) in
                sample_haplotype_list[0].iter().enumerate()
            {
                x[*haplotype_index] = Some(new_haplotype_index);
            }
            x
        };

        // 2. Get mapping from old to new SV indexes
        let sv_index_remap = {
            let mut x = vec![None; sv_haplotype_map.len()];
            for (new_sv_index, (sv_index, _)) in sv_haplotype_map
                .iter()
                .enumerate()
                .filter(|&(_, &hap_index)| haplotype_index_remap[hap_index].is_some())
                .enumerate()
            {
                x[sv_index] = Some(new_sv_index);
            }
            x
        };

        // 3. Update sample_haplotype_list to reflect the new haplotype index values:
        for haplotype_index in sample_haplotype_list[0].iter_mut() {
            // The unwrap here is an assertion that all sample_haplotype_list values have been
            // remapped:
            let new_haplotype_index = haplotype_index_remap[*haplotype_index].unwrap();
            *haplotype_index = new_haplotype_index;
        }

        // 4. Update sv_haplotype_map to reflect new sv and haplotype index values:
        sv_haplotype_map = sv_haplotype_map
            .into_iter()
            .filter_map(|x| haplotype_index_remap[x])
            .collect();

        // 5. Shrink group haplotypes to remove filtered entries
        group_haplotypes = group_haplotypes
            .into_iter()
            .enumerate()
            .filter_map(|(i, x)| haplotype_index_remap[i].map(|_| x))
            .collect();

        // 6. Shrink refined SVs to remove filtered entries
        refined_svs = refined_svs
            .into_iter()
            .enumerate()
            .filter_map(|(i, x)| sv_index_remap[i].map(|_| x))
            .collect();

        if debug {
            eprintln!("AFTER");
            eprintln!("sample_haplotype_list[0] {:?}", sample_haplotype_list[0]);
            eprintln!(
                "group_haplotypes {:?}",
                group_haplotypes
                    .iter()
                    .map(|x| &x.hap_id)
                    .collect::<Vec<_>>()
            );
            eprintln!("sv_haplotype_map {sv_haplotype_map:?}");
            eprintln!(
                "refined_svs {:?}",
                refined_svs.iter().map(|x| &x.id).collect::<Vec<_>>()
            );
        }
    }

    if refined_svs.is_empty() {
        None
    } else {
        Some(SVGroup {
            group_regions,
            group_haplotypes,
            sample_haplotype_list,
            sample_expected_cn_info,
            sv_haplotype_map,
            refined_svs,
        })
    }
}

/// Take a set of annotated refined SVs from the same cluster, and try to convert these into an SVGroup
fn process_rsv_cluster(
    treat_single_copy_as_haploid: bool,
    contigs: &SampleAssemblies,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    group_arsvs: &mut Vec<AnnotatedRefinedSV>,
    sv_groups: &mut Vec<SVGroup>,
) {
    if !group_arsvs.is_empty()
        && let Some(sv_group) = process_rsv_cluster_to_sv_group(
            treat_single_copy_as_haploid,
            contigs,
            expected_copy_number_regions,
            group_arsvs,
        )
    {
        sv_groups.push(sv_group);
    }
}

/// Convert RefinedSV structures created from the candidate VCF, into the SVGroup structure used to
/// process overlapping SVs
///
fn group_single_sample_refined_svs(
    treat_single_copy_as_haploid: bool,
    contigs: &SampleAssemblies,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
    mut anno_refined_svs: Vec<AnnotatedRefinedSV>,
) -> Vec<SVGroup> {
    // Sort the refined svs to match the sorting process created during refinement. This also allows
    // the refined SVs to be joined together into the same single-sample SV groups formed during
    // refinement
    //
    anno_refined_svs.sort_by_key(|x| {
        let id = x.refined_sv.id.clone();
        (id.cluster_index, id.assembly_index, id.alignment_index)
    });

    // Convert refined sv list into sv groups. For single sample the grouping is simply all Svs
    // with the same cluster index
    //
    let mut sv_groups = Vec::new();

    // Cluster the annotated refined SVs into groups with matching cluster_index values, then
    // process each of these clusters into an sv_group
    let mut group_arsvs = Vec::new();
    let mut last_cluster_index = None;
    for arsv in anno_refined_svs.into_iter() {
        if last_cluster_index != Some(arsv.refined_sv.id.cluster_index) {
            process_rsv_cluster(
                treat_single_copy_as_haploid,
                contigs,
                expected_copy_number_regions,
                &mut group_arsvs,
                &mut sv_groups,
            );
            last_cluster_index = Some(arsv.refined_sv.id.cluster_index);
        }
        group_arsvs.push(arsv);
    }
    process_rsv_cluster(
        treat_single_copy_as_haploid,
        contigs,
        expected_copy_number_regions,
        &mut group_arsvs,
        &mut sv_groups,
    );

    // Check sv_group validity
    for sv_group in sv_groups.iter() {
        sv_group.assert_validity(treat_single_copy_as_haploid);
    }

    sv_groups
}

/// Walk through candidate SV records for sample and convert these into sv_groups for downstream processing
///
#[allow(clippy::too_many_arguments)]
fn process_sv_groups_from_candidate_sv_file(
    disable_small_indels: bool,
    treat_single_copy_as_haploid: bool,
    discover_dir: &Utf8Path,
    chrom_list: &ChromList,
    assembly_regions: &[Vec<GenomeSegment>],
    contigs: &SampleAssemblies,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
) -> Vec<SVGroup> {
    let candidate_sv_filename = discover_dir.join(CANDIDATE_SV_FILENAME);
    if !candidate_sv_filename.exists() {
        panic!("Can't find candidate SV file expected at {candidate_sv_filename}");
    }

    let mut reader = bcf::Reader::from_path(candidate_sv_filename).unwrap();
    let header = reader.header();

    // Get chrom names from header -- bcf chrom order might be different than our own chromosome
    // order:
    //
    let rid_to_chrom_index = {
        let mut rid_to_chrom_index = Vec::new();
        let chrom_count = header.contig_count();
        for rid in 0..chrom_count {
            let chrom_bytes = header.rid2name(rid).unwrap();
            let chrom_name = std::str::from_utf8(chrom_bytes).unwrap().to_string();
            let chrom_index = *chrom_list.label_to_index.get(&chrom_name).unwrap();
            rid_to_chrom_index.push(chrom_index);
        }
        rid_to_chrom_index
    };

    let mut anno_refined_svs = Vec::new();

    let mut rec = reader.empty_record();
    while let Some(r) = reader.read(&mut rec) {
        unwrap!(r, "Failed to parse variant record");

        let arsv = process_bcf_record_to_anno_refine_sv(
            chrom_list,
            assembly_regions,
            contigs,
            &rid_to_chrom_index,
            &rec,
        );
        if let Some(anno_refined_sv) = arsv {
            if disable_small_indels && anno_refined_sv.refined_sv.single_region_refinement {
                continue;
            }
            anno_refined_svs.push(anno_refined_sv);
        }
    }

    group_single_sample_refined_svs(
        treat_single_copy_as_haploid,
        contigs,
        expected_copy_number_regions,
        anno_refined_svs,
    )
}

/// Read SV Group list from candidate SV bcf
///
pub fn get_sample_sv_groups(
    disable_small_indels: bool,
    treat_single_copy_as_haploid: bool,
    discover_dir: &Utf8Path,
    chrom_list: &ChromList,
    target_regions: &GenomeRegions,
    expected_copy_number_regions: Option<&GenomeRegionsByChromIndex>,
) -> Vec<SVGroup> {
    // First read in assembly regions and contigs from discovery data, then we have enough
    // information to reconstruct the sv groups from the candidate SV file
    let assembly_regions = read_assembly_regions_from_bed(discover_dir, chrom_list);
    let contigs = get_candidate_sv_contig_info(discover_dir, chrom_list);
    let mut sv_groups = process_sv_groups_from_candidate_sv_file(
        disable_small_indels,
        treat_single_copy_as_haploid,
        discover_dir,
        chrom_list,
        &assembly_regions,
        &contigs,
        expected_copy_number_regions,
    );

    // Target regions is used for debugging only, so the intersection criteria here are very simple:
    // - The first begin pos of the first SV in a group needs to intersect the target regions
    if !target_regions.is_empty() {
        // Reformat target_regions to index on chrom id instead of name:
        let mut chrom_id_target_regions = Vec::new();
        for chrom_info in chrom_list.data.iter() {
            chrom_id_target_regions.push(target_regions.chroms.get(&chrom_info.label));
        }

        sv_groups.retain(|x| {
            let segment = &x.refined_svs[0].bp.breakend1.segment;
            if let Some(regions) = chrom_id_target_regions[segment.chrom_index] {
                regions.intersect_pos(segment.range.start)
            } else {
                false
            }
        });
    }

    sv_groups
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bnd_alt_allele() {
        let chrom_list = {
            let mut x = ChromList::default();
            x.add_chrom("chr1", 1000);
            x
        };

        let alt_allele = b"TCT]chr1:500]";

        let homlen = 10;
        let bnd_alt_info = parse_bnd_alt_allele(&chrom_list, alt_allele, homlen);

        assert_eq!(
            bnd_alt_info.breakend1_direction,
            BreakendDirection::LeftAnchor
        );
        assert_eq!(
            bnd_alt_info.breakend2_direction,
            BreakendDirection::LeftAnchor
        );
        assert_eq!(bnd_alt_info.breakend2_chrom_index, 0);
        assert_eq!(bnd_alt_info.breakend2_pos, 489);
        assert_eq!(bnd_alt_info.insert_seq, b"CT");
    }

    #[test]
    fn test_parse_bnd_alt_allele_hla() {
        let chrom_list = {
            let mut x = ChromList::default();
            x.add_chrom("HLA-DRB1*10:01:01", 1000);
            x
        };

        let alt_allele = b"TCT]HLA-DRB1*10:01:01:500]";

        let homlen = 10;
        let bnd_alt_info = parse_bnd_alt_allele(&chrom_list, alt_allele, homlen);

        assert_eq!(
            bnd_alt_info.breakend1_direction,
            BreakendDirection::LeftAnchor
        );
        assert_eq!(
            bnd_alt_info.breakend2_direction,
            BreakendDirection::LeftAnchor
        );
        assert_eq!(bnd_alt_info.breakend2_chrom_index, 0);
        assert_eq!(bnd_alt_info.breakend2_pos, 489);
        assert_eq!(bnd_alt_info.insert_seq, b"CT");
    }
}
