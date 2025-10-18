use camino::Utf8Path;
use log::{debug, error, info};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
};
use rust_vc_utils::cigar::has_aligned_segments;
use rust_vc_utils::int_range::IntRange;
use rust_vc_utils::{ChromList, bam_reg2bin, get_alignment_end};

use crate::filenames::CONTIG_ALIGNMENT_FILENAME;
use crate::globals::PROGRAM_VERSION;

pub const SA_AUX_TAG: &[u8] = b"SA";
pub const CONTIG_AUX_TAG: &[u8] = b"sf";

#[derive(Clone)]
pub struct ContigAlignmentSegment {
    /// Reference chromosome ID
    pub tid: i32,

    /// Start position of alignment
    pub pos: i64,

    pub cigar: Vec<Cigar>,

    pub is_fwd_strand: bool,
}

impl std::fmt::Debug for ContigAlignmentSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "ContigAlignmentSegment: tid: {:?} pos: {:?} fwd?: {:?} cigar: {}",
            self.tid,
            self.pos,
            self.is_fwd_strand,
            CigarString(self.cigar.clone())
        )
    }
}

/// Assembly contig alignment information
///
#[derive(Clone)]
pub struct ContigAlignmentInfo {
    /// Contig source sample
    pub sample_index: usize,

    /// Index of SV candidate cluster
    pub cluster_index: usize,

    /// Index of all assemblies from one cluster
    pub assembly_index: usize,

    /// Sequence in fwd-strand orientation, this needs to be revcomped for any supplementary
    /// segments
    pub seq: Vec<u8>,

    /// Replicate this whole structure for each segment and add the zero-indexed segment id here
    pub segment_id: usize,

    /// Multiple alignment segments, mark the first one as primary and the others as supplemental
    pub segments: Vec<ContigAlignmentSegment>,

    /// Number of reads used to assemble the contig
    ///
    /// This isn't required for bam visualization, but used as part of joint genotyping
    ///
    pub supporting_read_count: usize,

    /// Region of the contig assembled into a consensus sequence by POA
    ///
    /// The contig may also include optional left and right flanks at lower consensus quality
    ///
    /// This isn't required for bam visualization, but used as part of joint genotyping.
    ///
    pub high_quality_contig_range: IntRange,
}

impl ContigAlignmentInfo {
    fn get_segment(&self) -> &ContigAlignmentSegment {
        &self.segments[self.segment_id]
    }
}

/// Produce one split read segment in the format required within the SA aux tag
///
fn get_sa_tag_segment(
    chrom_list: &ChromList,
    segment: &ContigAlignmentSegment,
    mapq: u8,
) -> String {
    let chrom = &chrom_list.data[segment.tid as usize].label;
    let schar = if segment.is_fwd_strand { '+' } else { '-' };
    format!(
        "{chrom},{},{schar},{},{mapq},0;",
        segment.pos + 1,
        CigarString(segment.cigar.clone()),
    )
}

fn convert_contig_alignment_to_bam_record(
    chrom_list: &ChromList,
    pkg_name: &str,
    contig_alignment: ContigAlignmentInfo,
) -> bam::Record {
    const CONTIG_MAPQ: u8 = 20;
    let qname = format!(
        "{}:{}:{}:{}",
        pkg_name,
        contig_alignment.sample_index,
        contig_alignment.cluster_index,
        contig_alignment.assembly_index
    );
    let mut record = bam::Record::new();
    let quals = vec![255u8; contig_alignment.seq.len()];
    let this_segment = contig_alignment.get_segment();
    record.set(
        qname.as_bytes(),
        Some(&CigarString(this_segment.cigar.clone())),
        contig_alignment.seq.as_slice(),
        &quals,
    );
    record.unset_unmapped();
    record.set_tid(this_segment.tid);
    record.set_pos(this_segment.pos);
    record.set_mapq(CONTIG_MAPQ);
    {
        let mut flag = 0u16;
        if !this_segment.is_fwd_strand {
            flag |= rust_htslib::htslib::BAM_FREVERSE as u16;
        }
        if contig_alignment.segment_id != 0 {
            flag |= rust_htslib::htslib::BAM_FSUPPLEMENTARY as u16;
        }
        if flag > 0 {
            record.set_flags(flag);
        }
    }
    let end = get_alignment_end(&record) as usize;
    record.set_bin(bam_reg2bin(this_segment.pos as usize, end));

    // Add SA tags:
    if contig_alignment.segments.len() > 1 {
        let mut aux_str = String::new();
        for (segment_index, segment) in contig_alignment.segments.iter().enumerate() {
            if segment_index == contig_alignment.segment_id {
                continue;
            }
            if !has_aligned_segments(&segment.cigar) {
                panic!("invalid CIGAR string in contig SA output for contig qname {qname}");
            }
            aux_str += get_sa_tag_segment(chrom_list, segment, CONTIG_MAPQ).as_str();
        }
        record.push_aux(SA_AUX_TAG, Aux::String(&aux_str)).unwrap();
    }

    // Add aux tag for contig metadata used during joint genotyping:
    {
        let hqrange = &contig_alignment.high_quality_contig_range;
        let aux_str = format!(
            "n_reads:{};hq_range:{}-{};",
            contig_alignment.supporting_read_count, hqrange.start, hqrange.end,
        );
        record
            .push_aux(CONTIG_AUX_TAG, Aux::String(&aux_str))
            .unwrap();
    }

    record
}

/// Create bam header for the sawfish contig alignments
///
/// This is a simple header containing just the contig info, and the sawfish commandline as a "PG" entry
///
fn get_contig_alignment_header(chrom_list: &ChromList, pkg_name: &str) -> bam::header::Header {
    let mut new_header = bam::header::Header::new();

    let mut hd_record = bam::header::HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", "1.6");
    hd_record.push_tag(b"SO", "coordinate");
    new_header.push_record(&hd_record);

    for chrom_info in chrom_list.data.iter() {
        let mut sq_record = bam::header::HeaderRecord::new(b"SQ");
        sq_record.push_tag(b"SN", &chrom_info.label);
        sq_record.push_tag(b"LN", chrom_info.length);
        new_header.push_record(&sq_record);
    }

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    let mut pg_record = bam::header::HeaderRecord::new(b"PG");
    pg_record.push_tag(b"PN", pkg_name);
    pg_record.push_tag(b"ID", format!("{pkg_name}-{PROGRAM_VERSION}"));
    pg_record.push_tag(b"VN", PROGRAM_VERSION);
    pg_record.push_tag(b"CL", &cmdline);

    new_header.push_record(&pg_record);
    new_header
}

/// Completes writing (and closing) the contig bam file itself
///
fn write_contig_alignments_bam(
    filename: &Utf8Path,
    thread_count: usize,
    chrom_list: &ChromList,
    contig_alignments: Vec<ContigAlignmentInfo>,
) {
    let pkg_name = env!("CARGO_PKG_NAME");

    // Setup new bam writer for debug output:
    let mut bam_writer = {
        let output_bam_header = get_contig_alignment_header(chrom_list, pkg_name);
        bam::Writer::from_path(filename, &output_bam_header, bam::Format::Bam).unwrap()
    };
    bam_writer.set_threads(thread_count).unwrap();

    for contig_alignment in contig_alignments.into_iter() {
        let record = convert_contig_alignment_to_bam_record(chrom_list, pkg_name, contig_alignment);
        bam_writer.write(&record).unwrap();
    }
}

/// The bam_header should not already have a sawfish PG entry, to avoid conflicting PG entires which cause problems in IGV
///
pub fn write_contig_alignments(
    output_dir: &Utf8Path,
    thread_count: usize,
    chrom_list: &ChromList,
    contig_alignment_type_label: &str,
    mut contig_alignments: Vec<ContigAlignmentInfo>,
) {
    // Sort contig alignments
    contig_alignments.sort_by(|a, b| {
        let aseg = a.get_segment();
        let bseg = b.get_segment();
        aseg.tid.cmp(&bseg.tid).then(aseg.pos.cmp(&bseg.pos))
    });

    let filename = output_dir.join(CONTIG_ALIGNMENT_FILENAME);

    info!("Writing {contig_alignment_type_label} contig alignments to bam file: '{filename}'");

    write_contig_alignments_bam(&filename, thread_count, chrom_list, contig_alignments);

    // Index bam file
    //

    // Samtools describes this value as "Set minimum interval size for CSI indices to 2^INT [14]""
    let csi_min_interval_pow2 = 14;
    match bam::index::build(
        &filename,
        None,
        bam::index::Type::Csi(csi_min_interval_pow2),
        thread_count as u32,
    ) {
        Ok(()) => {
            debug!("Finished building index for bam file {filename}.");
        }
        Err(e) => {
            error!("Error while building index for {filename}: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };
}
