use camino::{Utf8Path, Utf8PathBuf};
use log::info;
use rust_htslib::bcf;
use rust_vc_utils::{ChromList, GenomeRef, rev_comp_in_place};
use strum::EnumCount;

use crate::breakpoint::{
    Breakend, BreakendDirection, BreakendNeighbor, Breakpoint, InsertInfo, VcfSVType,
    get_breakpoint_vcf_sv_type,
};
use crate::cli;
use crate::expected_ploidy::SVLocusPloidy;
use crate::genome_segment::GenomeSegment;
use crate::refine_sv::{
    AlleleType, Genotype, RefinedSV, SVPhaseStatus, SVSampleScoreInfo, SVScoreInfo,
    get_rsv_id_label,
};
use crate::refined_cnv::{CNVSampleScoreInfo, CNVScoreInfo, RefinedCNV, SharedSampleScoreInfo};
use crate::score_sv::QualityModelAlleles;
use crate::sv_group::SVGroup;
use crate::sv_id::get_bnd_id_label;
use crate::vcf_utils;

pub const CONTIG_POS_INFO_KEY: &str = "CONTIG_POS";
pub const OVERLAP_ASM_INFO_KEY: &str = "OVERLAP_ASM";
pub const BND0_NEIGHBOR_INFO_KEY: &str = "BND0_NEIGHBOR";
pub const BND1_NEIGHBOR_INFO_KEY: &str = "BND1_NEIGHBOR";

fn get_sv_vcf_header(
    settings: &VcfSettings,
    chrom_list: &ChromList,
    sample_names: &[&str],
    candidate_mode: bool,
) -> bcf::Header {
    let mut header =
        vcf_utils::get_basic_vcf_header(&settings.ref_filename, chrom_list, sample_names);

    let mut records = Vec::<&[u8]>::new();

    records.append(&mut vec![
        br#"##ALT=<ID=DEL,Description="Deletion">"#,
        br#"##ALT=<ID=INS,Description="Insertion">"#,
        br#"##ALT=<ID=DUP,Description="Duplication">"#,
        br#"##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">"#,
        br#"##ALT=<ID=INV,Description="Inversion">"#,
        br#"##ALT=<ID=CNV,Description="Copy number variable region">"#,
    ]);

    // Note that the `PASS` and `.` FILTER records below are not typically included in the header,
    // but there's an oddity in programatically creating vcfs in htslib that forces these to be present,
    // and specifically ordered before INFO/FORMAT.
    //
    records.push(br#"##FILTER=<ID=PASS,Description="All filters passed">"#);
    records.push(br#"##FILTER=<ID=.,Description="Unknown filtration status">"#);

    let min_qual_filter = format!(
        "##FILTER=<ID=MinQUAL,Description=\"QUAL score is less than {}\">",
        &settings.min_qual
    );

    if !candidate_mode {
        records.push(br#"##FILTER=<ID=InvBreakpoint,Description="Breakpoint represented as part of an inversion record (same EVENT ID)">"#);
        records.push(min_qual_filter.as_bytes());
        records.push(br#"##FILTER=<ID=MaxScoringDepth,Description="SV candidate exceeds max scoring depth in at least one sample">"#);
        records.push(br#"##FILTER=<ID=ConflictingBreakpointGT,Description="Genotypes of breakpoints in a multi-breakpoint event conflict in the majority of cases">"#);
    }

    records.append(&mut vec![
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
        br#"##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively">"#,
        br#"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">"#,
        br#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">"#,]);

    let contig_pos = format!(
        "##INFO=<ID={CONTIG_POS_INFO_KEY},Number=1,Type=Integer,Description=\"SV position in the contig sequence before breakend1\">"
    );
    let overlap_asm = format!(
        "##INFO=<ID={OVERLAP_ASM_INFO_KEY},Number=.,Type=Integer,Description=\"Assembly indices of overlapping haplotypes\">"
    );
    let bnd0_neighbor = format!(
        "##INFO=<ID={BND0_NEIGHBOR_INFO_KEY},Number=1,Type=String,Description=\"Breakend0 neighbor cluster and breakend indexes\">"
    );
    let bnd1_neighbor = format!(
        "##INFO=<ID={BND1_NEIGHBOR_INFO_KEY},Number=1,Type=String,Description=\"Breakend1 neighbor cluster and breakend indexes\">"
    );

    if candidate_mode {
        records.push(contig_pos.as_bytes());
        records.push(overlap_asm.as_bytes());
        records.push(bnd0_neighbor.as_bytes());
        records.push(bnd1_neighbor.as_bytes());
    } else {
        let mut score_records: Vec<&[u8]> = vec![
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
            br#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#,
            br#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">"#,
            br#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">"#,
            //br#"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Read depth for each allele on the forward strand">"#,
            //br#"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Read depth for each allele on the reverse strand">"#,
            br#"##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number">"#,
            br#"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality">"#,
        ];
        if settings.enable_phasing {
            score_records.push(
                br#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#,
            );
        }
        records.append(&mut score_records);
    }

    for x in records.into_iter() {
        header.push_record(x);
    }

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
        CopyNumberVariant => "CNV",
    }
}

/// A layer on top of the bcf Record struct which encodes the VCF record sorting logic
struct VcfRecord {
    record: bcf::Record,
    is_inversion: bool,
}

impl PartialEq for VcfRecord {
    fn eq(&self, other: &Self) -> bool {
        self.record.rid() == other.record.rid() && self.record.pos() == other.record.pos()
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
        // Sort first by rid and chrom to make tabix-able VCF output, then sort by:
        // 1. is_inversion : This makes inversion records precede their corresponding BND at the same position
        // 2. id : This is the final sort term, to keep the output deterministic since all IDs should be unique
        self.record
            .rid()
            .cmp(&other.record.rid())
            .then(self.record.pos().cmp(&other.record.pos()))
            .then(other.is_inversion.cmp(&self.is_inversion))
            .then(self.record.id().cmp(&other.record.id()))
    }
}

/// Add the SVTYPE info tag to a bcf record
fn add_sv_type(sv_type: VcfSVType, record: &mut bcf::Record) {
    let sv_type_label = get_sv_type_vcf_label(sv_type);
    record
        .push_info_string(b"SVTYPE", &[sv_type_label.as_bytes()])
        .unwrap();
}

/// Set HOMLEN and HOMSEQ
fn add_homology_tags(seg: &GenomeSegment, hom_seq: &[u8], record: &mut bcf::Record) {
    let hom_len = seg.range.size() - 1;
    assert_eq!(hom_len as usize, hom_seq.len());
    if hom_len > 0 {
        record
            .push_info_integer(b"HOMLEN", &[hom_len as i32])
            .unwrap();

        record.push_info_string(b"HOMSEQ", &[hom_seq]).unwrap();
    }
}

/// Set QUAL
fn add_qual(score: &SVScoreInfo, record: &mut bcf::Record) {
    record.set_qual(match score.alt_score() {
        Some(x) => x.round(),
        None => bcf::record::Numeric::missing(),
    });
}

/// candidate_mode - True if writing out the candidate VCF in discovery mode, prior to genotyping
/// is_inversion_filter - True for breakends represented as a higher-level INV record, which will be filtered out
fn add_sv_vcf_filters(
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
        if let Some(alt_score) = score.alt_score()
            && alt_score < settings.min_qual as f32
        {
            record.push_filter("MinQUAL".as_bytes()).unwrap();
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
        record.push_info_integer(b"INSLEN", &[inslen]).unwrap();
    }
}

fn add_insseq(bp: &Breakpoint, record: &mut bcf::Record) {
    if let InsertInfo::Seq(seq) = &bp.insert_info
        && !seq.is_empty()
    {
        record
            .push_info_string(b"INSSEQ", &[seq.as_slice()])
            .unwrap();
    }
}

fn add_event(rsv: &RefinedSV, record: &mut bcf::Record) {
    if let Some(inversion_id) = &rsv.ext.inversion_id {
        record
            .push_info_string(b"EVENT", &[inversion_id.as_bytes()])
            .unwrap();
        record.push_info_string(b"EVENTTYPE", &[b"INV"]).unwrap();
    }
}

fn add_svclaim(
    sv_type: VcfSVType,
    depth_support: bool,
    bp_support: bool,
    record: &mut bcf::Record,
) {
    if sv_type == VcfSVType::Deletion || sv_type == VcfSVType::Duplication {
        let mut svc = Vec::new();
        if depth_support {
            svc.push(b'D');
        }
        if bp_support {
            svc.push(b'J');
        }
        record
            .push_info_string(b"SVCLAIM", &[svc.as_slice()])
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

fn add_bnd_neighbor(
    breakend0_neighbor: Option<&BreakendNeighbor>,
    breakend1_neighbor: Option<&BreakendNeighbor>,
    record: &mut bcf::Record,
) {
    if let Some(x) = breakend0_neighbor {
        let val = format!("{}:{}", x.cluster_index, x.breakend_index);
        record
            .push_info_string(BND0_NEIGHBOR_INFO_KEY.as_bytes(), &[val.as_bytes()])
            .unwrap();
    }
    if let Some(x) = breakend1_neighbor {
        let val = format!("{}:{}", x.cluster_index, x.breakend_index);
        record
            .push_info_string(BND1_NEIGHBOR_INFO_KEY.as_bytes(), &[val.as_bytes()])
            .unwrap();
    }
}

/// Get SV htslib-encoded genotype for sample
///
/// # Arguments:
/// * `treat_single_copy_as_haploid` - If true, treat expected single-copy regions as haploid
/// * `format_haploid_samples_as_diploid` - If true format GT values as "1/1" and "0/0" instead of "1" and "0" in haploid regions
///
/// Return a 2-tuple of (1) htslib-encoded GT value (2) bool indicating if the GT is phased
///
fn get_sv_genotype(
    sample_score: &SVSampleScoreInfo,
    format_haploid_samples_as_diploid: bool,
    treat_single_copy_as_haploid: bool,
    enable_phasing: bool,
) -> (Vec<i32>, bool) {
    use bcf::record::GenotypeAllele::*;
    match sample_score
        .expected_cn_info
        .ploidy(treat_single_copy_as_haploid)
    {
        SVLocusPloidy::Diploid => match &sample_score.shared.gt {
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
        SVLocusPloidy::Haploid => match &sample_score.shared.gt {
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
                    panic!("Unexpected haploid GT value. sample_score: {sample_score:?}");
                }
            },
            None => (
                vec![i32::from(UnphasedMissing), vcf_utils::VECTOR_END_INTEGER],
                false,
            ),
        },
    }
}

fn get_genotype_quality(sample_score: &SharedSampleScoreInfo) -> i32 {
    match &sample_score.gt {
        Some(_) => sample_score.gt_qscore,
        None => bcf::record::Numeric::missing(),
    }
}

/// # Arguments
/// * `treat_single_copy_as_haploid` - It true, treat expected single-copy regions as haploid
///
fn get_genotype_lhoods(
    sample_score: &SVSampleScoreInfo,
    format_haploid_samples_as_diploid: bool,
    treat_single_copy_as_haploid: bool,
) -> Vec<i32> {
    // This is defined even when GT is unknown, because in the absense of data, likelihoods are known even if uninformative
    let gt_lhoods = sample_score.shared.gt_lhood_qscore.clone();

    if (!format_haploid_samples_as_diploid)
        && sample_score
            .expected_cn_info
            .ploidy(treat_single_copy_as_haploid)
            == SVLocusPloidy::Haploid
    {
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
    record.push_format_integer(b"ADF", &adf).unwrap();
    let adr = ad.iter().map(|x| x.rev_strand as i32).collect::<Vec<_>>();
    record.push_format_integer(b"ADR", &adr).unwrap();
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

/// # Arguments
/// * `treat_single_copy_as_haploid` - It true, treat expected single-copy regions as haploid
/// * `has_cnv_tags` - True if SV was matched to CNV and has additional CNV vcf tags
///
fn add_sv_sample_info(
    sv_score: &SVScoreInfo,
    record: &mut bcf::Record,
    enable_phasing: bool,
    treat_single_copy_as_haploid: bool,
    has_cnv_tags: bool,
) {
    // Some tools don't want to accept haploid GT records, so this flag will setup VCF output to
    // write them as diploid instead if true.
    //
    let format_haploid_samples_as_diploid = false;

    let mut gts = Vec::new();
    let mut gqs = Vec::new();
    let mut pls = Vec::new();
    let mut ads = Vec::new();
    let mut cns = Vec::new();
    let mut cnqs = Vec::new();
    let mut pss = Vec::new();
    let mut is_any_phased = false;
    for sample_score in sv_score.samples.iter() {
        let (gt, is_sample_phased) = get_sv_genotype(
            sample_score,
            format_haploid_samples_as_diploid,
            treat_single_copy_as_haploid,
            enable_phasing,
        );
        if is_sample_phased {
            is_any_phased = true;
        }

        gts.extend(gt);
        if enable_phasing && let Some(phase_set) = sv_score.phase_set {
            pss.push(if is_sample_phased {
                phase_set
            } else {
                bcf::record::Numeric::missing()
            });
        }

        gqs.push(get_genotype_quality(&sample_score.shared));
        pls.extend(get_genotype_lhoods(
            sample_score,
            format_haploid_samples_as_diploid,
            treat_single_copy_as_haploid,
        ));
        ads.extend(get_sample_allele_depths(sample_score));

        if has_cnv_tags {
            cns.push(get_copy_number(&sample_score.shared));
            cnqs.push(get_copy_number_quality(&sample_score.shared));
        }
    }
    record.push_format_integer(b"GT", &gts).unwrap();
    record.push_format_integer(b"GQ", &gqs).unwrap();
    record.push_format_integer(b"PL", &pls).unwrap();
    record.push_format_integer(b"AD", &ads).unwrap();
    if has_cnv_tags {
        record.push_format_integer(b"CN", &cns).unwrap();
        record.push_format_integer(b"CNQ", &cnqs).unwrap();
    }
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
    add_sv_vcf_filters(
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
        add_sv_sample_info(
            &rsv.score,
            &mut record,
            settings.enable_phasing,
            settings.treat_single_copy_as_haploid,
            rsv.ext.is_cnv_match,
        );
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
    //
    // Note that we keep this even though its deprecated in VCF 4.4, to reduce transitional confusion
    record
        .push_info_integer(b"END", &[(end + 1) as i32])
        .unwrap();

    // Set SVLEN
    //
    // SVLEN is different in the VCF 4.4 spec, always positive and follows the symbolic allele
    {
        // Ensure SVLEN is always positive to better align with the VCF 4.4 spec:
        let sv_len = match sv_type {
            VcfSVType::Insertion => ins_size,
            _ => del_size,
        };
        record
            .push_info_integer(b"SVLEN", &[sv_len as i32])
            .unwrap();
    }

    add_homology_tags(seg1, &rsv.breakend1_homology_seq, &mut record);

    if !rsv.bp.is_precise {
        record.push_info_flag(b"IMPRECISE").unwrap();
    }

    add_inslen(&rsv.bp, &mut record);

    if is_symbolic_alt {
        add_insseq(&rsv.bp, &mut record);
    }

    add_event(rsv, &mut record);

    add_svclaim(sv_type, rsv.ext.is_cnv_match, true, &mut record);

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
        let id_str = get_bnd_id_label(&rsv.id, is_breakend1);
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
        let mate_id_str = get_bnd_id_label(&rsv.id, !is_breakend1);
        record
            .push_info_string(b"MATEID", &[mate_id_str.as_bytes()])
            .unwrap();
    }

    if !rsv.bp.is_precise {
        record.push_info_flag(b"IMPRECISE").unwrap();
    }

    add_inslen(&rsv.bp, &mut record);

    add_event(rsv, &mut record);

    let is_bnd = true;
    let mut record = add_vcf_record_sample_scoring(
        settings,
        rsv,
        group_contig_indexes,
        candidate_mode,
        is_bnd,
        record,
    );

    if candidate_mode {
        add_bnd_neighbor(
            rsv.bp.breakend1_neighbor.as_ref(),
            rsv.bp.breakend2_neighbor.as_ref(),
            &mut record,
        );
    }

    Some(record)
}

fn is_same_pos_record(r1: &bcf::Record, r2: &bcf::Record) -> bool {
    r1.rid() == r2.rid() && r1.pos() == r2.pos()
}

fn get_bcf_record_homlen(rec: &bcf::Record) -> Option<i32> {
    let x = rec.info(b"HOMLEN").integer().unwrap();
    if let Some(x) = x {
        x.first().cloned()
    } else {
        None
    }
}

fn get_bcf_record_end(rec: &bcf::Record) -> Option<i32> {
    let x = rec.info(b"END").integer().unwrap();
    if let Some(x) = x {
        x.first().cloned()
    } else {
        None
    }
}

fn get_bcf_record_svlen(rec: &bcf::Record) -> Option<i32> {
    let x = rec.info(b"SVLEN").integer().unwrap();
    if let Some(x) = x {
        x.first().cloned()
    } else {
        None
    }
}

fn is_duplicate_record(r1: &bcf::Record, r2: &bcf::Record) -> bool {
    is_same_pos_record(r1, r2)
        && (r1.alleles() == r2.alleles())
        && (get_bcf_record_homlen(r1) == get_bcf_record_homlen(r2))
        && (get_bcf_record_svlen(r1) == get_bcf_record_svlen(r2))
        && (get_bcf_record_end(r1) == get_bcf_record_end(r2))
}

/// Filter out duplicate vcf records
///
/// Returns the de-duplicated records and the count of filtered duplicate records
///
/// Duplicated variants are supposed to be prevented by the pipeline through good calling logic, but
/// in case any make it this far we remove them here as a final filtration step. Note that because
/// this is just an 'emergency backup filter', it is conservative, and some duplications may make it
/// through.
///
fn dedup_records(records: Vec<VcfRecord>) -> (Vec<VcfRecord>, usize) {
    let mut deduplicated_records: Vec<VcfRecord> = Vec::new();
    let mut first_matching_pos_index = 0;
    let mut duplicate_record_count = 0;
    for record in records {
        // Update first_matching_pos_index
        if !(deduplicated_records.is_empty()
            || is_same_pos_record(
                &record.record,
                &deduplicated_records[first_matching_pos_index].record,
            ))
        {
            first_matching_pos_index = deduplicated_records.len();
        }

        // Check this record against all accepted records at the same position to decide if it will
        // be filtered.
        //
        let mut filter_duplicate = false;
        for accepted_record in
            deduplicated_records[first_matching_pos_index..deduplicated_records.len()].iter()
        {
            if is_duplicate_record(&record.record, &accepted_record.record) {
                filter_duplicate = true;
                duplicate_record_count += 1;
                break;
            }
        }
        if !filter_duplicate {
            deduplicated_records.push(record);
        }
    }
    (deduplicated_records, duplicate_record_count)
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
                    records.push(VcfRecord {
                        record,
                        is_inversion: refined_sv.ext.is_inversion,
                    });
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
                    let is_inversion = false;
                    records.push(VcfRecord {
                        record: record0,
                        is_inversion,
                    });
                    records.push(VcfRecord {
                        record: record1,
                        is_inversion,
                    });
                }
            }
            VcfSVType::SingleBreakend => {
                panic!("No support for single-ended breakend output to VCF");
            }
            VcfSVType::CopyNumberVariant => {
                panic!("Invalid type in sv_group: CopyNumberVariant");
            }
        }
    }
}

/// Get CNV htslib-encoded genotype for sample
///
/// Return htslib-encoded GT value
///
fn get_cnv_genotype(
    treat_single_copy_as_haploid: bool,
    cnv_sample_score: &CNVSampleScoreInfo,
) -> Vec<i32> {
    let is_haploid = treat_single_copy_as_haploid && cnv_sample_score.expected_copy_number == 1;
    let copy_number = cnv_sample_score
        .shared
        .copy_info
        .as_ref()
        .map(|x| x.copy_number);

    use bcf::record::GenotypeAllele::*;
    if !is_haploid {
        // diploid (default) case:
        match copy_number {
            Some(0) => vec![i32::from(Unphased(1)), i32::from(Unphased(1))],
            Some(1) => vec![i32::from(Unphased(0)), i32::from(Unphased(1))],
            Some(x) if x > 1 => vec![i32::from(UnphasedMissing), i32::from(Unphased(1))],
            _ => vec![i32::from(UnphasedMissing), i32::from(UnphasedMissing)],
        }
    } else {
        // haploid case:
        match copy_number {
            Some(x) if x == 0 || x > 1 => {
                vec![i32::from(Unphased(1)), vcf_utils::VECTOR_END_INTEGER]
            }
            _ => vec![i32::from(UnphasedMissing), vcf_utils::VECTOR_END_INTEGER],
        }
    }
}

fn get_copy_number(sample_score: &SharedSampleScoreInfo) -> i32 {
    match &sample_score.copy_info {
        Some(x) => x.copy_number as i32,
        None => bcf::record::Numeric::missing(),
    }
}

fn get_copy_number_quality(sample_score: &SharedSampleScoreInfo) -> i32 {
    match &sample_score.copy_info {
        Some(x) => x.copy_number_qscore,
        None => bcf::record::Numeric::missing(),
    }
}

fn add_cnv_sample_info(settings: &VcfSettings, cnv_score: &CNVScoreInfo, record: &mut bcf::Record) {
    let mut gts = Vec::new();
    let mut cns = Vec::new();
    let mut cnqs = Vec::new();

    for cnv_sample_score in cnv_score.samples.iter() {
        cns.push(get_copy_number(&cnv_sample_score.shared));
        cnqs.push(get_copy_number_quality(&cnv_sample_score.shared));

        let gt = get_cnv_genotype(settings.treat_single_copy_as_haploid, cnv_sample_score);
        gts.extend(gt);
    }
    record.push_format_integer(b"GT", &gts).unwrap();
    record.push_format_integer(b"CN", &cns).unwrap();
    record.push_format_integer(b"CNQ", &cnqs).unwrap();
}

/// Set CNV QUAL
fn add_cnv_qual(score: &CNVScoreInfo, record: &mut bcf::Record) {
    record.set_qual(match score.alt_score {
        Some(x) => x.round(),
        None => bcf::record::Numeric::missing(),
    });
}

fn add_cnv_vcf_filters(settings: &VcfSettings, score: &CNVScoreInfo, record: &mut bcf::Record) {
    if let Some(alt_score) = score.alt_score
        && alt_score < settings.min_qual as f32
    {
        record.push_filter("MinQUAL".as_bytes()).unwrap();
    }
    if record.filters().peekable().peek().is_none() {
        record.push_filter("PASS".as_bytes()).unwrap();
    }
}

fn add_vcf_record_cnv_sample_scoring(
    settings: &VcfSettings,
    score: &CNVScoreInfo,
    record: &mut bcf::Record,
) {
    add_cnv_qual(score, record);

    add_cnv_vcf_filters(settings, score, record);

    add_cnv_sample_info(settings, score, record);
}

fn get_cnv_id_label(cnv: &RefinedCNV) -> String {
    format!(
        "{}:CNV:{}:{}:{}",
        env!("CARGO_PKG_NAME"),
        cnv.source_sample_index,
        cnv.segment.chrom_index,
        cnv.segment.range.start
    )
}

fn convert_cnv_to_vcf_record(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_count: usize,
    vcf: &bcf::Writer,
    cnv: &RefinedCNV,
) -> Option<VcfRecord> {
    assert_eq!(sample_count, cnv.score.samples.len());

    let sv_type = match cnv.get_cnv_type() {
        Some(x) => x,
        None => {
            return None;
        }
    };

    let seg = &cnv.segment;

    // CNV segmentation format marks the start of the CNV range as zero-indexed inclusive. Here we translate to a zero-indexed version of the
    // standard VCF position, which starts one base before the CNV, per the spec for all symbolic alleles, so need to subtract 1.
    let start = seg.range.start - 1;
    let end = seg.range.end - 1;

    // Although a VCF stat POS of 0 (1-indexed) is VCF spec compliant, it is de-facto non-spec, b/c htslib and htsjdk don't really support it, so
    // this case gets special-cased up to 1:
    let start = std::cmp::max(start, 0);

    assert!(end > start);
    let cnv_size = end - start;

    // Bail out of difficult cases
    if start < -1 {
        return None;
    }

    let mut record = vcf.empty_record();
    record.set_rid(Some(seg.chrom_index as u32));
    record.set_pos(start);

    // Set ID
    {
        let id_str = get_cnv_id_label(cnv);
        record.set_id(id_str.as_bytes()).unwrap();
    }

    // Set REF and ALT
    let chrom_info = &chrom_list.data[seg.chrom_index];
    let chrom_seq = genome_ref.chroms.get(chrom_info.label.as_str()).unwrap();

    let ref_seq = if start >= 0 {
        &chrom_seq[start as usize..(start + 1) as usize]
    } else {
        b"N"
    };

    let alt_seq = match sv_type {
        VcfSVType::Deletion => b"<DEL>",
        VcfSVType::Duplication => b"<DUP>",
        _ => b"<CNV>",
    };
    record.set_alleles(&[ref_seq, alt_seq]).unwrap();

    record.push_info_flag(b"IMPRECISE").unwrap();

    add_sv_type(sv_type, &mut record);

    // Set END - (use value+1 to convert to 1-indexed position)
    record
        .push_info_integer(b"END", &[(end + 1) as i32])
        .unwrap();

    record
        .push_info_integer(b"SVLEN", &[cnv_size as i32])
        .unwrap();

    add_svclaim(sv_type, true, false, &mut record);

    add_vcf_record_cnv_sample_scoring(settings, &cnv.score, &mut record);

    let is_inversion = false;
    Some(VcfRecord {
        record,
        is_inversion,
    })
}

/// Convert `SVGroup`s into `VcfRecord`s
///
/// Returns the de-duplicated VCF records and the count of filtered duplicate VCF records
///
#[allow(clippy::too_many_arguments)]
fn get_vcf_records(
    settings: &VcfSettings,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_count: usize,
    sv_groups: &[SVGroup],
    cnvs: &[RefinedCNV],
    vcf: &bcf::Writer,
    candidate_mode: bool,
) -> (Vec<VcfRecord>, usize) {
    let mut records = Vec::new();

    for sv_group in sv_groups.iter() {
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
    let (mut deduplicated_records, duplicate_record_count) = if enable_dedup_records {
        dedup_records(records)
    } else {
        (records, 0)
    };

    // Add CNV records AFTER deduplication because the CNV segmentation process prevents the
    // style of FP duplicate we experience from breakpoint-based calling.
    //
    for cnv in cnvs {
        let cnv_vcf_record =
            convert_cnv_to_vcf_record(settings, genome_ref, chrom_list, sample_count, vcf, cnv);
        if let Some(x) = cnv_vcf_record {
            deduplicated_records.push(x);
        }
    }

    deduplicated_records.sort();

    (deduplicated_records, duplicate_record_count)
}

/// # Arguments
/// * `candidate_mode` - If true, don't write out scoring information, and use bcf format
///
#[allow(clippy::too_many_arguments)]
fn write_sv_vcf_file(
    settings: &VcfSettings,
    filename: &Utf8Path,
    genome_ref: &GenomeRef,
    chrom_list: &ChromList,
    sample_names: &[&str],
    sv_groups: &[SVGroup],
    cnvs: &[RefinedCNV],
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
        cnvs,
        &vcf,
        candidate_mode,
    );

    for record in vcf_records.iter() {
        vcf.write(&record.record).unwrap();
    }

    VcfStats {
        output_record_count: vcf_records.len(),
        duplicate_record_count: vcf_duplicate_record_count,
    }
}

pub struct VcfSettings {
    pub ref_filename: String,
    pub output_dir: Utf8PathBuf,
    pub min_qual: i32,

    /// Disable end-stage VCF duplicate record filter (if not in candidate_mode)
    pub no_vcf_dedup: bool,

    /// All deletions larger than this size will be represented with a symbolic allele
    pub min_symbolic_deletion_size: usize,

    /// Enable phased genotype output if true
    pub enable_phasing: bool,

    /// If true, treat SVs in expected single-copy regions as haploid
    pub treat_single_copy_as_haploid: bool,
}

impl VcfSettings {
    pub fn new(
        ref_filename: &str,
        output_dir: &Utf8Path,
        min_qual: i32,
        no_vcf_dedup: bool,
        enable_phasing: bool,
        treat_single_copy_as_haploid: bool,
    ) -> Self {
        Self {
            ref_filename: ref_filename.to_string(),
            output_dir: output_dir.to_owned(),
            min_qual,
            no_vcf_dedup,
            min_symbolic_deletion_size: 100_001,
            enable_phasing,
            treat_single_copy_as_haploid,
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
    cnvs: &[RefinedCNV],
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

    info!("Writing {label} structural variants to file: '{filename}'");

    let vcf_stats = write_sv_vcf_file(
        settings,
        &filename,
        genome_ref,
        chrom_list,
        sample_names,
        sv_groups,
        cnvs,
        candidate_mode,
    );

    let build_tbi = !candidate_mode;
    vcf_utils::build_bcf_index(&filename, shared_settings.thread_count, build_tbi).unwrap();

    vcf_stats
}
