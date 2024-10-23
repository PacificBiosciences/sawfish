use std::path::{Path, PathBuf};

use log::info;
use rust_htslib::bcf;
use rust_vc_utils::{rev_comp_in_place, ChromList, GenomeRef};
use strum::EnumCount;

use crate::breakpoint::{
    get_breakpoint_vcf_sv_type, Breakend, BreakendDirection, Breakpoint, InsertInfo, VcfSVType,
};
use crate::cli;
use crate::expected_ploidy::SVLocusPloidy;
use crate::genome_segment::GenomeSegment;
use crate::refine_sv::{
    get_rsv_id_label, AlleleType, Genotype, RefinedSV, SVPhaseStatus, SVSampleScoreInfo,
    SVScoreInfo,
};
use crate::score_sv::QualityModelAlleles;
use crate::sv_group::SVGroup;
use crate::sv_id::get_bnd_sv_id_label;
use crate::vcf_utils;

pub const CONTIG_POS_INFO_KEY: &str = "CONTIG_POS";
pub const OVERLAP_ASM_INFO_KEY: &str = "OVERLAP_ASM";

fn get_sv_vcf_header(
    settings: &VcfSettings,
    chrom_list: &ChromList,
    sample_names: &[&str],
    candidate_mode: bool,
) -> bcf::Header {
    let mut header =
        vcf_utils::get_basic_vcf_header(&settings.ref_filename, chrom_list, sample_names);

    let contig_pos = format!("##INFO=<ID={},Number=1,Type=Integer,Description=\"SV position in the contig sequence before breakend1\">", CONTIG_POS_INFO_KEY);

    let overlap_asm = format!("##INFO=<ID={},Number=.,Type=Integer,Description=\"Assembly indices of overlapping haplotypes\">", OVERLAP_ASM_INFO_KEY);

    let min_qual_filter = format!(
        "##FILTER=<ID=MinQUAL,Description=\"QUAL score is less than {}\">",
        &settings.min_qual
    );

    // Note that the `PASS` and `.` FILTER records below are not typically included in the header,
    // but there's an oddity that rust-htslib interface that forces these to be present
    //
    let mut records: Vec<&[u8]> = vec![
        br#"##ALT=<ID=DEL,Description="Deletion">"#,
        br#"##ALT=<ID=INS,Description="Insertion">"#,
        br#"##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">"#,
        br#"##ALT=<ID=INV,Description="Inversion">"#,
        //br#"##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">"#,
        //br#"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">"#,
        br#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">"#,
        br#"##INFO=<ID=EVENT,Number=A,Type=String,Description="ID of associated event">"#,
        br#"##INFO=<ID=EVENTTYPE,Number=A,Type=String,Description="Type of associated event">"#,
        br#"##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">"#,
        br#"##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">"#,
        br#"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">"#,
        br#"##INFO=<ID=INSLEN,Number=.,Type=Integer,Description="Insertion length">"#,
        br#"##INFO=<ID=INSSEQ,Number=.,Type=String,Description="Insertion sequence (for symbolic alt alleles)">"#,
        br#"##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">"#,
        br#"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">"#,
        br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#,];

    if candidate_mode {
        records.push(contig_pos.as_bytes());
        records.push(overlap_asm.as_bytes());
    } else {
        let mut score_records: Vec<&[u8]> = vec![
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
            br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#,
            br#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">"#,
            br#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">"#,
            //br#"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Read depth for each allele on the forward strand">"#,
            //br#"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Read depth for each allele on the reverse strand">"#,
        ];
        if settings.enable_phasing {
            score_records.push(
                br#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#,
            );
        }
        records.append(&mut score_records);
    }

    let mut std_records: Vec<&[u8]> = vec![
        br#"##FILTER=<ID=.,Description="Unknown filtration status">"#,
        br#"##FILTER=<ID=PASS,Description="All filters passed">"#,
    ];
    records.append(&mut std_records);

    if !candidate_mode {
        records.push(br#"##FILTER=<ID=InvBreakpoint,Description="Breakpoint represented as part of an inversion record (same EVENT ID)">"#);
        records.push(min_qual_filter.as_bytes());
        records.push(br#"##FILTER=<ID=MaxScoringDepth,Description="SV candidate exceeds max scoring depth">"#);
        records.push(br#"##FILTER=<ID=ConflictingBreakpointGT,Description="Genotypes of breakpoints in a multi-breakpoint event conflict in the majority of cases">"#);
    }

    records.into_iter().for_each(|x| {
        header.push_record(x);
    });

    header
}

fn get_sv_type_vcf_label(sv_type: VcfSVType) -> &'static str {
    use VcfSVType::*;
    match sv_type {
        Deletion => "DEL",
        Inversion => "INV",
        Insertion => "INS",
        Duplication => "DUP",
        Breakpoint => "BND",
        SingleBreakend => "BND",
    }
}

/// A layer on top of the bcf Record struct which encodes the VCF record sorting logic
struct VcfRecord(bcf::Record);

impl PartialEq for VcfRecord {
    fn eq(&self, other: &Self) -> bool {
        self.0.rid() == other.0.rid() && self.0.pos() == other.0.pos()
    }
}

impl Eq for VcfRecord {}

impl PartialOrd for VcfRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VcfRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Sort first by rid and chrom to make tabix-able VCF output, then sort by other factors
        // for VCF records at the same position, to make the output deterministic
        self.0
            .rid()
            .cmp(&other.0.rid())
            .then(self.0.pos().cmp(&other.0.pos()))
            .then(self.0.id().cmp(&other.0.id()))
    }
}

/// Add the SVTYPE info tag to a bcf record
fn add_sv_type(sv_type: VcfSVType, record: &mut bcf::Record) {
    let sv_type_label = get_sv_type_vcf_label(sv_type);
    record
        .push_info_string("SVTYPE".as_bytes(), &[sv_type_label.as_bytes()])
        .unwrap();
}

/// Set HOMLEN and HOMSEQ
fn add_homology_tags(seg: &GenomeSegment, hom_seq: &[u8], record: &mut bcf::Record) {
    let hom_len = seg.range.size() - 1;
    assert_eq!(hom_len as usize, hom_seq.len());
    if hom_len > 0 {
        record
            .push_info_integer("HOMLEN".as_bytes(), &[hom_len as i32])
            .unwrap();

        record
            .push_info_string("HOMSEQ".as_bytes(), &[hom_seq])
            .unwrap();
    }
}

/// Set QUAL
fn add_qual(score: &SVScoreInfo, record: &mut bcf::Record) {
    record.set_qual(match score.alt_score {
        Some(x) => x.round(),
        None => bcf::record::Numeric::missing(),
    });
}

/// candidate_mode - True if writing out the candidate VCF in discovery mode, prior to genotyping
/// is_inversion_filter - True for breakends represented as a higher-level INV record, which will be filtered out
fn add_vcf_filters(
    settings: &VcfSettings,
    score: &SVScoreInfo,
    record: &mut bcf::Record,
    candidate_mode: bool,
    is_inversion_filter: bool,
    is_conflicting_gt: bool,
) {
    if candidate_mode {
        record.push_filter(".".as_bytes()).unwrap();
    } else {
        //TODO change this to a SAMPLE filter
        let is_max_scoring_depth = score.samples.iter().any(|x| x.is_max_scoring_depth);
        if is_max_scoring_depth {
            record.push_filter("MaxScoringDepth".as_bytes()).unwrap();
        }
        if let Some(alt_score) = score.alt_score {
            if alt_score < settings.min_qual as f32 {
                record.push_filter("MinQUAL".as_bytes()).unwrap();
            }
        }
        if is_inversion_filter {
            record.push_filter("InvBreakpoint".as_bytes()).unwrap();
        }
        if is_conflicting_gt {
            record
                .push_filter("ConflictingBreakpointGT".as_bytes())
                .unwrap();
        }
        if record.filters().peekable().peek().is_none() {
            record.push_filter("PASS".as_bytes()).unwrap();
        }
    }
}

fn add_inslen(bp: &Breakpoint, record: &mut bcf::Record) {
    let inslen = bp.insert_info.size() as i32;
    if inslen > 0 {
        record
            .push_info_integer("INSLEN".as_bytes(), &[inslen])
            .unwrap();
    }
}

fn add_insseq(bp: &Breakpoint, record: &mut bcf::Record) {
    if let InsertInfo::Seq(seq) = &bp.insert_info {
        if !seq.is_empty() {
            record
                .push_info_string("INSSEQ".as_bytes(), &[seq.as_slice()])
                .unwrap();
        }
    }
}

fn add_event(rsv: &RefinedSV, record: &mut bcf::Record) {
    if let Some(inversion_id) = rsv.ext.inversion_id {
        let event_tag = format!("INV{inversion_id}");
        record
            .push_info_string("EVENT".as_bytes(), &[event_tag.as_bytes()])
            .unwrap();
        record
            .push_info_string("EVENTTYPE".as_bytes(), &["INV".as_bytes()])
            .unwrap();
    }
}

fn add_contig_pos(contig_pos: usize, record: &mut bcf::Record) {
    record
        .push_info_integer(CONTIG_POS_INFO_KEY.as_bytes(), &[contig_pos as i32])
        .unwrap();
}

fn add_overlap_indexes(group_contig_indexes: &[usize], record: &mut bcf::Record) {
    if group_contig_indexes.len() > 1 {
        let group_contig_indexes = group_contig_indexes
            .iter()
            .map(|&x| x as i32)
            .collect::<Vec<_>>();
        record
            .push_info_integer(
                OVERLAP_ASM_INFO_KEY.as_bytes(),
                group_contig_indexes.as_slice(),
            )
            .unwrap();
    }
}

/// Get GT for sample
///
/// Return a 2-tuple of (1) htslib-encoded GT value (2) bool indicating if the GT is phased
///
fn get_genotype(
    sample_score: &SVSampleScoreInfo,
    format_haploid_samples_as_diploid: bool,
    enable_phasing: bool,
) -> (Vec<i32>, bool) {
    use bcf::record::GenotypeAllele::*;
    match sample_score.ploidy {
        SVLocusPloidy::Diploid => match &sample_score.gt {
            Some(gt) => {
                let (a0, a1) = match gt {
                    Genotype::Ref => (0, 0),
                    Genotype::Het => (0, 1),
                    Genotype::Hom => (1, 1),
                };
                let phase = if enable_phasing && (*gt == Genotype::Het) {
                    sample_score.phase
                } else {
                    SVPhaseStatus::Unknown
                };
                match phase {
                    SVPhaseStatus::Primary => {
                        (vec![i32::from(Phased(a0)), i32::from(Phased(a1))], true)
                    }
                    SVPhaseStatus::Secondary => {
                        (vec![i32::from(Phased(a1)), i32::from(Phased(a0))], true)
                    }
                    SVPhaseStatus::Unknown => (
                        vec![i32::from(Unphased(a0)), i32::from(Unphased(a1))],
                        false,
                    ),
                }
            }
            None => (
                vec![i32::from(UnphasedMissing), i32::from(UnphasedMissing)],
                false,
            ),
        },
        SVLocusPloidy::Haploid => match &sample_score.gt {
            Some(gt) => match gt {
                Genotype::Ref => {
                    if format_haploid_samples_as_diploid {
                        (vec![i32::from(Unphased(0)), i32::from(Unphased(0))], false)
                    } else {
                        (
                            vec![i32::from(Unphased(0)), vcf_utils::VECTOR_END_INTEGER],
                            false,
                        )
                    }
                }
                Genotype::Hom => {
                    if format_haploid_samples_as_diploid {
                        (vec![i32::from(Unphased(1)), i32::from(Unphased(1))], false)
                    } else {
                        (
                            vec![i32::from(Unphased(1)), vcf_utils::VECTOR_END_INTEGER],
                            false,
                        )
                    }
                }
                Genotype::Het => {
                    panic!(
                        "Unexpected haploid GT value. sample_score: {:?}",
                        sample_score
                    );
                }
            },
            None => (
                vec![i32::from(UnphasedMissing), vcf_utils::VECTOR_END_INTEGER],
                false,
            ),
        },
    }
}

fn get_genotype_quality(sample_score: &SVSampleScoreInfo) -> i32 {
    match &sample_score.gt {
        Some(_) => sample_score.gt_qscore,
        None => bcf::record::Numeric::missing(),
    }
}

fn get_genotype_lhoods(
    sample_score: &SVSampleScoreInfo,
    format_haploid_samples_as_diploid: bool,
) -> Vec<i32> {
    // This is defined even when GT is unknown, because in the absense of data, likelihoods are known even if uninformative
    let gt_lhoods = sample_score.gt_lhood_qscore.clone();

    if (!format_haploid_samples_as_diploid) && sample_score.ploidy == SVLocusPloidy::Haploid {
        vec![
            gt_lhoods[Genotype::Ref as usize],
            gt_lhoods[Genotype::Hom as usize],
            vcf_utils::VECTOR_END_INTEGER,
        ]
    } else {
        gt_lhoods
    }
}

fn get_sample_allele_depths(sample_score: &SVSampleScoreInfo) -> Vec<i32> {
    /*
    let adf = ad.iter().map(|x| x.fwd_strand as i32).collect::<Vec<_>>();
    record.push_format_integer("ADF".as_bytes(), &adf).unwrap();
    let adr = ad.iter().map(|x| x.rev_strand as i32).collect::<Vec<_>>();
    record.push_format_integer("ADR".as_bytes(), &adr).unwrap();
     */

    if sample_score.is_max_scoring_depth || sample_score.allele_depth.is_empty() {
        let mut x = vec![vcf_utils::VECTOR_END_INTEGER; QualityModelAlleles::COUNT];
        x[0] = bcf::record::Numeric::missing();
        x
    } else {
        // Compress Ref and Overlap allele states down to just Ref by convention:
        let ref_ad = (sample_score.allele_depth[AlleleType::Ref as usize].any_strand
            + sample_score.allele_depth[AlleleType::Overlap as usize].any_strand)
            as i32;
        let alt_ad = (sample_score.allele_depth[AlleleType::Alt as usize].any_strand) as i32;
        vec![ref_ad, alt_ad]
    }
}

fn add_sample_info(sv_score: &SVScoreInfo, record: &mut bcf::Record, enable_phasing: bool) {
    // Some tools don't want to accept haploid GT records, so this flag will setup VCF output to
    // write them as diploid instead if true.
    //
    let format_haploid_samples_as_diploid = false;

    let mut gts = Vec::new();
    let mut gqs = Vec::new();
    let mut pls = Vec::new();
    let mut ads = Vec::new();
    let mut pss = Vec::new();
    let mut is_any_phased = false;
    for sample_score in sv_score.samples.iter() {
        let (gt, is_sample_phased) = get_genotype(
            sample_score,
            format_haploid_samples_as_diploid,
            enable_phasing,
        );
        if is_sample_phased {
            is_any_phased = true;
        }

        gts.extend(gt);
        if enable_phasing {
            if let Some(phase_set) = sv_score.phase_set {
                pss.push(if is_sample_phased {
                    phase_set
                } else {
                    bcf::record::Numeric::missing()
                });
            }
        }

        gqs.push(get_genotype_quality(sample_score));
        pls.extend(get_genotype_lhoods(
            sample_score,
            format_haploid_samples_as_diploid,
        ));
        ads.extend(get_sample_allele_depths(sample_score));
    }
    record.push_format_integer(b"GT", &gts).unwrap();
    record.push_format_integer(b"GQ", &gqs).unwrap();
    record.push_format_integer(b"PL", &pls).unwrap();
    record.push_format_integer(b"AD", &ads).unwrap();
    if !pss.is_empty() && is_any_phased {
        record.push_format_integer(b"PS", &pss).unwrap();
    }
}

/// Add a scoring and sample fields
///
/// candidate_mode - True if writing out the candidate VCF in discovery mode, prior to genotyping
/// is_bnd - True if the record is being output in BND format
///
fn add_vcf_record_sample_scoring(
    settings: &VcfSettings,
    rsv: &RefinedSV,
    group_contig_indexes: &[usize],
    candidate_mode: bool,
    is_bnd: bool,
    mut record: bcf::Record,
) -> bcf::Record {
    add_qual(&rsv.score, &mut record);

    let is_inversion_filter = is_bnd && rsv.ext.inversion_id.is_some();
    let is_conflicting_gt = rsv.ext.is_inversion && rsv.ext.is_conflicting_gt;
    add_vcf_filters(
        settings,
        &rsv.score,
        &mut record,
        candidate_mode,
        is_inversion_filter,
        is_conflicting_gt,
    );

    if candidate_mode {
        add_contig_pos(rsv.contig_pos_before_breakend1, &mut record);
        add_overlap_indexes(group_contig_indexes, &mut record);
    } else {
        add_sample_info(&rsv.score, &mut record, settings.enable_phasing);
    }

    record
}

/// SV type classifier which moves up from just the bnd level to classify INV as well
///
fn get_rsv_vcf_sv_type(rsv: &RefinedSV) -> VcfSVType {
    let sv_type = get_breakpoint_vcf_sv_type(&rsv.bp);
    if sv_type == VcfSVType::Breakpoint && rsv.ext.is_inversion {
        VcfSVType::Inversion
    } else {
        sv_type
    }
}

fn convert_refined_sv_to_non_bnd_vcf_record(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    vcf: &bcf::Writer,
    rsv: &RefinedSV,
    group_contig_indexes: &[usize],
    candidate_mode: bool,
) -> Option<bcf::Record> {
    let sv_type = get_rsv_vcf_sv_type(rsv);
    let is_allowed_type = sv_type == VcfSVType::Deletion
        || sv_type == VcfSVType::Insertion
        || sv_type == VcfSVType::Duplication
        || sv_type == VcfSVType::Inversion;
    assert!(is_allowed_type);

    let seg1 = &rsv.bp.breakend1.segment;
    let start = seg1.range.start;

    // Bail out of difficult cases that can occasionally occur in circular chromosomes:
    // TODO: Better handle origin-spanning SVs in mitochondria
    if start < 0 {
        return None;
    }

    let mut record = vcf.empty_record();
    record.set_rid(Some(seg1.chrom_index as u32));
    record.set_pos(start);

    // Set ID
    {
        let id_str = get_rsv_id_label(rsv);
        record.set_id(id_str.as_bytes()).unwrap();
    }

    let end = rsv.bp.breakend2.as_ref().unwrap().segment.range.start;
    assert!(end >= start);

    let del_size = end - start;
    let ins_size = rsv.bp.insert_info.size() as i64;

    let is_symbolic_del = (sv_type == VcfSVType::Deletion)
        && (del_size >= settings.min_symbolic_deletion_size as i64);
    let is_symbolic_dup = sv_type == VcfSVType::Duplication;
    let is_symbolic_inv = sv_type == VcfSVType::Inversion;
    let is_symbolic_alt = is_symbolic_del || is_symbolic_dup || is_symbolic_inv;

    // Set REF and ALT
    let chrom_info = &chrom_list.data[seg1.chrom_index];
    let chrom_seq = genome_ref.chroms.get(chrom_info.label.as_str()).unwrap();
    if let InsertInfo::SizeRange(_) = &rsv.bp.insert_info {
        assert_ne!(
            sv_type,
            VcfSVType::Insertion,
            "Unexpected CandidateSV format"
        );
        let ref_seq = &chrom_seq[start as usize..(start + 1) as usize];
        let alt_seq = "<INS>".as_bytes();
        record.set_alleles(&[ref_seq, alt_seq]).unwrap();
    } else {
        let ref_seq_len = if is_symbolic_alt { 0 } else { del_size };
        let ref_seq = &chrom_seq[start as usize..(start + 1 + ref_seq_len) as usize];
        let alt_seq = if is_symbolic_del {
            b"<DEL>".to_vec()
        } else if is_symbolic_dup {
            b"<DUP:TANDEM>".to_vec()
        } else if is_symbolic_inv {
            b"<INV>".to_vec()
        } else {
            let mut x = ref_seq[0..1].to_vec();
            if let InsertInfo::Seq(seq) = &rsv.bp.insert_info {
                x.extend(seq);
            };
            x
        };
        record.set_alleles(&[ref_seq, alt_seq.as_slice()]).unwrap();
    }

    add_sv_type(sv_type, &mut record);

    // Set END - (use value+1 to convert to 1-indexed position)
    record
        .push_info_integer("END".as_bytes(), &[(end + 1) as i32])
        .unwrap();

    // Set SVLEN
    {
        // Get a positive SVLEN value for DUP and INV to be forward-compatible with VCF 4.4:
        let sv_len = match sv_type {
            VcfSVType::Deletion | VcfSVType::Insertion => ins_size - del_size,
            _ => del_size,
        };
        record
            .push_info_integer("SVLEN".as_bytes(), &[sv_len as i32])
            .unwrap();
    }

    add_homology_tags(seg1, &rsv.breakend1_homology_seq, &mut record);

    if !rsv.bp.is_precise {
        record.push_info_flag("IMPRECISE".as_bytes()).unwrap();
    }

    add_inslen(&rsv.bp, &mut record);

    if is_symbolic_alt {
        add_insseq(&rsv.bp, &mut record);
    }

    add_event(rsv, &mut record);

    let is_bnd = false;
    let record = add_vcf_record_sample_scoring(
        settings,
        rsv,
        group_contig_indexes,
        candidate_mode,
        is_bnd,
        record,
    );

    Some(record)
}

/// Get breakend from a breakpoint where breakend2 is already known to exist
///
fn get_bp_breakend(bp: &Breakpoint, is_breakend1: bool) -> &Breakend {
    if is_breakend1 {
        &bp.breakend1
    } else {
        bp.breakend2.as_ref().unwrap()
    }
}

/// Right-anchored BNDs are represented in VCF one position AFTER the breakend
fn get_bnd_vcf_start_pos(be: &Breakend) -> i64 {
    be.segment.range.start
        + if be.dir == BreakendDirection::RightAnchor {
            1
        } else {
            0
        }
}

/// Build the infamous BND alt sequence
///
fn get_bnd_alt_seq(
    chrom_list: &ChromList,
    bp: &Breakpoint,
    is_breakend1: bool,
    ref_seq: &[u8],
) -> Option<String> {
    let mate_breakend = get_bp_breakend(bp, !is_breakend1);
    let mate_segment = &mate_breakend.segment;
    let is_same_breakend_orientation = bp.same_orientation();

    // Getting mate_pos right is tricky so break down the steps here:
    let mate_pos = {
        // First, adjust which end of the homology range we pull from depending on whether this has
        // an inverted orientation:
        let mut pos = if is_same_breakend_orientation {
            mate_segment.range.end - 1
        } else {
            mate_segment.range.start
        };

        // Next, all right-anchored BND records are positioned on the base AFTER the break
        if mate_breakend.dir == BreakendDirection::RightAnchor {
            pos += 1;
        }

        // Finally, adjust this value to be 1-indexed since it's going into a string rather than a
        // (0-indexed) BCF position value:
        pos + 1
    };

    // Bail out of difficult cases that can occasionally occur in circular chromosomes:
    // TODO: Better handle origin-spanning SVs in mitochondria
    if mate_pos < 1 {
        return None;
    }

    // BND alt format includes sequence, including any breakpoint insertion sequence, as either
    // a prefix or suffix of the expression:
    let (prefix_seq, suffix_seq) = {
        let mut insert_seq = match &bp.insert_info {
            InsertInfo::Seq(s) => {
                let reverse_insert = !is_breakend1 && is_same_breakend_orientation;
                if reverse_insert {
                    let mut s2 = s.clone();
                    rev_comp_in_place(&mut s2);
                    s2
                } else {
                    s.clone()
                }
            }
            InsertInfo::NoInsert => Vec::new(),
            InsertInfo::SizeRange(_) => {
                panic!("Invalid breakpoint format");
            }
        };
        let breakend = get_bp_breakend(bp, is_breakend1);
        match &breakend.dir {
            BreakendDirection::LeftAnchor => {
                let mut alt_part = ref_seq.to_vec();
                alt_part.extend(insert_seq);
                (
                    std::str::from_utf8(alt_part.as_slice())
                        .unwrap()
                        .to_string(),
                    String::new(),
                )
            }
            BreakendDirection::RightAnchor => {
                insert_seq.extend(ref_seq);
                (
                    String::new(),
                    std::str::from_utf8(insert_seq.as_slice())
                        .unwrap()
                        .to_string(),
                )
            }
        }
    };

    let bracket = match &mate_breakend.dir {
        BreakendDirection::LeftAnchor => ']',
        BreakendDirection::RightAnchor => '[',
    };

    let mate_chrom = &chrom_list.data[mate_segment.chrom_index].label;

    Some(format!(
        "{prefix_seq}{bracket}{mate_chrom}:{mate_pos}{bracket}{suffix_seq}"
    ))
}

/// * vcf - Used to link the corresponding bam header to the vcf record
///
#[allow(clippy::too_many_arguments)]
fn convert_refined_sv_to_bnd_vcf_record(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    vcf: &bcf::Writer,
    rsv: &RefinedSV,
    group_contig_indexes: &[usize],
    is_breakend1: bool,
    candidate_mode: bool,
) -> Option<bcf::Record> {
    assert!(rsv.bp.breakend2.is_some());

    let breakend = get_bp_breakend(&rsv.bp, is_breakend1);
    let seg = &breakend.segment;

    let start = get_bnd_vcf_start_pos(breakend);

    // Bail out of difficult cases that can occasionally occur in circular chromosomes:
    // TODO: Better handle origin-spanning SVs in mitochondria
    if start < 0 {
        return None;
    }

    let mut record = vcf.empty_record();
    record.set_rid(Some(seg.chrom_index as u32));
    record.set_pos(start);

    // Set ID
    {
        let id_str = get_bnd_sv_id_label(&rsv.id, is_breakend1);
        record.set_id(id_str.as_bytes()).unwrap();
    }

    //set REF and ALT
    let chrom_info = &chrom_list.data[seg.chrom_index];
    let chrom_seq = genome_ref.chroms.get(chrom_info.label.as_str()).unwrap();
    {
        let ref_seq = &chrom_seq[start as usize..(start + 1) as usize];
        let alt_seq = get_bnd_alt_seq(chrom_list, &rsv.bp, is_breakend1, ref_seq);

        if let Some(alt_seq) = alt_seq {
            record.set_alleles(&[ref_seq, alt_seq.as_bytes()]).unwrap();
        } else {
            return None;
        }
    }

    add_sv_type(VcfSVType::Breakpoint, &mut record);

    // Set HOMLEN and HOMSEQ
    {
        let mut hom_seq = rsv.breakend1_homology_seq.clone();
        let reverse_hom_seq = !is_breakend1 && rsv.bp.same_orientation();
        if reverse_hom_seq {
            rev_comp_in_place(&mut hom_seq);
        }
        add_homology_tags(seg, &hom_seq, &mut record);
    }

    // Set MATEID
    {
        let mate_id_str = get_bnd_sv_id_label(&rsv.id, !is_breakend1);
        record
            .push_info_string("MATEID".as_bytes(), &[mate_id_str.as_bytes()])
            .unwrap();
    }

    if !rsv.bp.is_precise {
        record.push_info_flag("IMPRECISE".as_bytes()).unwrap();
    }

    add_inslen(&rsv.bp, &mut record);

    add_event(rsv, &mut record);

    let is_bnd = true;
    let record = add_vcf_record_sample_scoring(
        settings,
        rsv,
        group_contig_indexes,
        candidate_mode,
        is_bnd,
        record,
    );

    Some(record)
}

fn is_same_pos_record(r1: &bcf::Record, r2: &bcf::Record) -> bool {
    r1.rid() == r2.rid() && r1.pos() == r2.pos()
}

fn get_bcf_record_homlen(rec: &bcf::Record) -> i64 {
    let x = rec.info(b"HOMLEN").integer().unwrap();
    if let Some(x) = x {
        x.first().cloned().unwrap() as i64
    } else {
        0
    }
}

fn is_duplicate_record(r1: &bcf::Record, r2: &bcf::Record) -> bool {
    is_same_pos_record(r1, r2)
        && (r1.alleles() == r2.alleles())
        && (get_bcf_record_homlen(r1) == get_bcf_record_homlen(r2))
}

/// Filter out duplicate vcf records
///
/// Returns the de-duplicated records and the count of filtered duplicate records
///
/// Duplicated variants are supposed to be prevented by the pipeline through good calling logic, but
/// in case any make it this far we remove them here as a final filtration step. Note that because
/// this is just an 'emergency backup filter', it is conservative.
///
fn dedup_records(records: Vec<VcfRecord>) -> (Vec<VcfRecord>, usize) {
    let mut new_records: Vec<VcfRecord> = Vec::new();
    let mut first_matching_pos_index = 0;
    let mut duplicate_record_count = 0;
    for record in records {
        // Update first_matching_pos_index
        if !(new_records.is_empty()
            || is_same_pos_record(&record.0, &new_records[first_matching_pos_index].0))
        {
            first_matching_pos_index = new_records.len();
        }

        // Check this record against all accepted records at the same position to decide if it will
        // be filtered.
        //
        let mut filter_duplicate = false;
        for accepted_record in new_records[first_matching_pos_index..new_records.len()].iter() {
            if is_duplicate_record(&record.0, &accepted_record.0) {
                filter_duplicate = true;
                duplicate_record_count += 1;
                break;
            }
        }
        if !filter_duplicate {
            new_records.push(record);
        }
    }
    (new_records, duplicate_record_count)
}

#[allow(clippy::too_many_arguments)]
fn convert_sv_group_to_vcf_records(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_count: usize,
    vcf: &bcf::Writer,
    candidate_mode: bool,
    sv_group: &SVGroup,
    records: &mut Vec<VcfRecord>,
) {
    let group_contig_indexes = if candidate_mode {
        sv_group.sample_haplotype_list[0].clone()
    } else {
        Vec::new()
    };

    for refined_sv in sv_group.refined_svs.iter().filter(|x| !x.filter_sv()) {
        assert_eq!(refined_sv.score.samples.len(), sample_count);

        // Force the SV type to breakpoint for various cases
        let sv_type = {
            let x = get_rsv_vcf_sv_type(refined_sv);
            let force_breakpoint = if x == VcfSVType::SingleBreakend {
                false
            } else if refined_sv.ext.force_breakpoint_representation {
                true
            } else if candidate_mode && x == VcfSVType::Deletion {
                let be1 = &refined_sv.bp.breakend1;
                let be2 = refined_sv.bp.breakend2.as_ref().unwrap();
                std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) > 1_000
            } else {
                candidate_mode && x == VcfSVType::Duplication
            };

            if force_breakpoint {
                VcfSVType::Breakpoint
            } else {
                x
            }
        };

        match sv_type {
            VcfSVType::Deletion
            | VcfSVType::Insertion
            | VcfSVType::Duplication
            | VcfSVType::Inversion => {
                let record = convert_refined_sv_to_non_bnd_vcf_record(
                    settings,
                    genome_ref,
                    chrom_list,
                    vcf,
                    refined_sv,
                    &group_contig_indexes,
                    candidate_mode,
                );
                if let Some(record) = record {
                    records.push(VcfRecord(record));
                }
            }
            VcfSVType::Breakpoint => {
                let record0 = convert_refined_sv_to_bnd_vcf_record(
                    settings,
                    genome_ref,
                    chrom_list,
                    vcf,
                    refined_sv,
                    &group_contig_indexes,
                    true,
                    candidate_mode,
                );
                let record1 = convert_refined_sv_to_bnd_vcf_record(
                    settings,
                    genome_ref,
                    chrom_list,
                    vcf,
                    refined_sv,
                    &group_contig_indexes,
                    false,
                    candidate_mode,
                );
                if let (Some(record0), Some(record1)) = (record0, record1) {
                    records.push(VcfRecord(record0));
                    records.push(VcfRecord(record1));
                }
            }
            VcfSVType::SingleBreakend => {
                panic!("No support for single-ended breakend output to VCF");
            }
        }
    }
}

/// Returns the de-duplicated records and the count of filtered duplicate records
///
fn get_vcf_records(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_count: usize,
    sv_groups: &[SVGroup],
    vcf: &bcf::Writer,
    candidate_mode: bool,
) -> (Vec<VcfRecord>, usize) {
    let mut records = Vec::new();

    for sv_group in sv_groups {
        convert_sv_group_to_vcf_records(
            settings,
            genome_ref,
            chrom_list,
            sample_count,
            vcf,
            candidate_mode,
            sv_group,
            &mut records,
        );
    }

    records.sort();

    let enable_dedup_records = !(candidate_mode || settings.no_vcf_dedup);

    if enable_dedup_records {
        dedup_records(records)
    } else {
        (records, 0)
    }
}

/// This file outputs just the SV information, consistent with a traditional SV caller
///
/// A separate 'large-variant' VCF will output the joint information from both SV and CNV
///
/// # Arguments
/// * candidate_mode - If true, don't write out scoring information, and use bcf format
///
#[allow(clippy::too_many_arguments)]
fn write_sv_vcf_file(
    settings: &VcfSettings,
    filename: &Path,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_names: &[&str],
    sv_groups: &[SVGroup],
    candidate_mode: bool,
) -> VcfStats {
    let format = if candidate_mode {
        bcf::Format::Bcf
    } else {
        bcf::Format::Vcf
    };

    let header = get_sv_vcf_header(settings, chrom_list, sample_names, candidate_mode);
    let mut vcf = bcf::Writer::from_path(filename, &header, false, format).unwrap();

    let sample_count = sample_names.len();
    let (vcf_records, vcf_duplicate_record_count) = get_vcf_records(
        settings,
        genome_ref,
        chrom_list,
        sample_count,
        sv_groups,
        &vcf,
        candidate_mode,
    );

    for record in vcf_records.iter() {
        vcf.write(&record.0).unwrap();
    }

    VcfStats {
        output_record_count: vcf_records.len(),
        duplicate_record_count: vcf_duplicate_record_count,
    }
}

pub struct VcfSettings {
    pub ref_filename: String,
    pub output_dir: PathBuf,
    pub min_qual: i32,

    /// Disable end-stage VCF duplicate record filter (if not in candidate_mode)
    pub no_vcf_dedup: bool,

    /// All deletions larger than this size will be represented with a symbolic allele
    pub min_symbolic_deletion_size: usize,

    /// Enable phased genotype output if true
    pub enable_phasing: bool,
}

impl VcfSettings {
    pub fn new(
        ref_filename: &str,
        output_dir: &Path,
        min_qual: i32,
        no_vcf_dedup: bool,
        enable_phasing: bool,
    ) -> Self {
        Self {
            ref_filename: ref_filename.to_string(),
            output_dir: output_dir.to_owned(),
            min_qual,
            no_vcf_dedup,
            min_symbolic_deletion_size: 100_001,
            enable_phasing,
        }
    }
}

pub struct VcfStats {
    pub output_record_count: usize,
    pub duplicate_record_count: usize,
}

/// # Arguments
/// * candidate_mode - If true, don't write out scoring information, and use bcf format
///
#[allow(clippy::too_many_arguments)]
pub fn write_indexed_sv_vcf_file(
    shared_settings: &cli::SharedSettings,
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_names: &[&str],
    sv_groups: &[SVGroup],
    candidate_mode: bool,
) -> VcfStats {
    assert!(sample_names.is_empty() || (!candidate_mode));

    let filename = settings.output_dir.join(if candidate_mode {
        "candidate.sv.bcf"
    } else {
        "genotyped.sv.vcf.gz"
    });

    let label = if candidate_mode {
        "candidate"
    } else {
        "genotyped"
    };

    info!(
        "Writing {} structural variants to file: '{}'",
        label,
        filename.display()
    );

    let vcf_stats = write_sv_vcf_file(
        settings,
        &filename,
        genome_ref,
        chrom_list,
        sample_names,
        sv_groups,
        candidate_mode,
    );

    let build_tbi = !candidate_mode;
    vcf_utils::build_bcf_index(&filename, shared_settings.thread_count, build_tbi).unwrap();

    vcf_stats
}
