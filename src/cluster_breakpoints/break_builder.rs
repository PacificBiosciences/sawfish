use std::collections::BTreeMap;

use rust_htslib::bam;
use rust_vc_utils::aux::is_aux_tag_found;
use rust_vc_utils::bam_utils::get_seq_order_read_position;
use rust_vc_utils::cigar::{
    get_cigar_ref_offset, get_hard_clipped_read_clip_positions,
    update_ref_and_hard_clipped_read_pos,
};
use rust_vc_utils::ChromList;

use crate::bam_sa_parser::{get_seq_order_read_split_segments, SeqOrderSplitReadSegment};
use crate::bam_utils::{test_read_for_large_insertion_soft_clip, LargeInsertionSoftClipState};
use crate::breakpoint::*;
use crate::genome_segment::{get_segment_distance, GenomeSegment};
use crate::int_range::IntRange;

/// Temporary representation used inside of break_builder to handle CIGAR candidates
#[derive(Debug)]
struct CandidateIndel {
    ref_pos: i64,
    del_size: u32,
    ins_size: u32,
}

impl CandidateIndel {
    fn new() -> Self {
        Self {
            ref_pos: 0,
            del_size: 0,
            ins_size: 0,
        }
    }

    fn is_empty(&self) -> bool {
        self.del_size == 0 && self.ins_size == 0
    }

    fn clear(&mut self) {
        self.ref_pos = 0;
        self.del_size = 0;
        self.ins_size = 0;
    }
}

pub struct BreakpointObservationMap {
    pub data: BTreeMap<GenomeSegment, Vec<BreakpointObservation>>,
}

impl BreakpointObservationMap {
    pub fn new() -> Self {
        Self {
            data: BTreeMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn len(&self) -> usize {
        self.data.values().flatten().count()
    }

    fn get_entry(&mut self, segment: &GenomeSegment) -> &mut Vec<BreakpointObservation> {
        self.data.entry(segment.clone()).or_default()
    }

    pub fn add(&mut self, bpo: BreakpointObservation) {
        assert!(bpo.breakpoint.is_standardized());
        self.get_entry(&bpo.breakpoint.breakend1.segment).push(bpo);
    }

    /// Merge other into self, destroying other
    ///
    pub fn merge(&mut self, other: &mut Self) {
        if self.is_empty() {
            std::mem::swap(self, other);
        } else {
            for (segment, obs) in other.data.iter_mut() {
                self.get_entry(segment).append(obs);
            }
        }
    }

    /// Sort the observations within each genome segment such that the downstream clustering
    /// method can produce stable results. This could later be a more limited version of sorting
    /// which only addresses what clustering needs.
    ///
    pub fn normalize(&mut self) {
        self.data.values_mut().for_each(|x| x.sort());
    }
}

/// Genome breakend observations
///
/// Raw breakend evidence prior to clustering and refinement. Exported by BreakBuider.
///
pub struct BreakObservations {
    /// The minimum observations to create a breakpoint cluster
    pub min_cluster_evidence_count: usize,

    /// Distance below which breakpoint observations will be clustered together
    ///
    /// Still defining precisely what this 'distance' between breakpoints means
    pub cluster_distance: usize,

    /// Individual observations that haven't been clustered yet
    ///
    /// We may not have observed enough of the genome yet to run clustering on some of these
    /// observations.
    pub singles: BreakpointObservationMap,
}

impl BreakObservations {
    pub fn new() -> Self {
        Self {
            min_cluster_evidence_count: 2,
            cluster_distance: 500,
            singles: BreakpointObservationMap::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.singles.is_empty()
    }

    /// Merge other into self, destroying other
    ///
    pub fn merge(&mut self, other: &mut Self) {
        if self.is_empty() {
            std::mem::swap(self, other);
        } else {
            self.singles.merge(&mut other.singles);
        }
    }
}

/// Accumulates breakpoint information in one segment of the genome
///
pub struct BreakBuilder<'a> {
    /// Minimum indel size of initial indel evidence
    ///
    /// Note this may be lower than the min_indel_size to account for noise in the size observations
    min_evidence_indel_size: u32,

    min_sv_mapq: u32,

    max_close_breakend_distance: usize,

    /// Soft-clipping on one side of the read must be at least this long to trigger as breakpoint
    /// evidence
    min_large_insertion_soft_clip_len: usize,

    worker_region: IntRange,

    observations: BreakObservations,

    chrom_list: &'a ChromList,

    is_targeted_scan: bool,
}

impl<'a> BreakBuilder<'a> {
    /// begin/end are the boundaries of the chromosome segment being analyzed
    pub fn new(
        min_evidence_indel_size: u32,
        min_sv_mapq: u32,
        max_close_breakend_distance: usize,
        begin: i64,
        end: i64,
        chrom_list: &'a ChromList,
        is_targeted_scan: bool,
    ) -> Self {
        Self {
            min_evidence_indel_size,
            min_sv_mapq,
            max_close_breakend_distance,
            min_large_insertion_soft_clip_len: 500,
            worker_region: IntRange::from_pair(begin, end),
            observations: BreakObservations::new(),
            chrom_list,
            is_targeted_scan,
        }
    }

    /// Check if candidate indel meets the minimum size, and if so return a breakpoint observation
    ///
    /// # Arguments
    /// * `record` - bam record for error/debug messages
    /// * `read_pos` - Position in read coordinates on the alignment strand for the position immediately after the indel
    ///
    fn process_cand_indel_impl(
        &self,
        record: &bam::Record,
        read_pos: usize,
        cand_indel: &CandidateIndel,
    ) -> Option<BreakpointObservation> {
        let indel_in_worker_region = self.worker_region.intersect_pos(cand_indel.ref_pos);
        let indel_large_enough = cand_indel.del_size >= self.min_evidence_indel_size
            || cand_indel.ins_size >= self.min_evidence_indel_size;

        if !(indel_in_worker_region && indel_large_enough) {
            return None;
        }

        let breakend1 = Breakend {
            segment: GenomeSegment {
                chrom_index: record.tid() as usize,
                range: IntRange::from_int(cand_indel.ref_pos),
            },
            dir: BreakendDirection::LeftAnchor,
        };

        let breakend2 = Some(Breakend {
            segment: GenomeSegment {
                chrom_index: record.tid() as usize,
                range: IntRange::from_int(cand_indel.ref_pos + cand_indel.del_size as i64),
            },
            dir: BreakendDirection::RightAnchor,
        });

        let breakpoint = Breakpoint {
            breakend1,
            breakend2,
            insert_info: InsertInfo::SizeRange(IntRange::from_int(cand_indel.ins_size as i64)),
            is_precise: false,
        };

        let breakend1_read_pos = match read_pos.checked_sub((cand_indel.ins_size + 1) as usize) {
            Some(x) => x,
            None => {
                let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                panic!(
                    "Unexpected insertion size {} and read_pos {read_pos} in read: {qname}",
                    cand_indel.ins_size
                );
            }
        };
        let breakend2_read_pos = read_pos;

        // Reset the read coordinates to sequencer order instead of alignment fwd-strand order
        let breakend1_seq_order_read_pos = get_seq_order_read_position(record, breakend1_read_pos);
        let breakend2_seq_order_read_pos = get_seq_order_read_position(record, breakend2_read_pos);

        Some(BreakpointObservation {
            breakpoint,
            evidence: BreakpointEvidenceType::ReadAlignment,
            qname: None,
            breakend1_seq_order_read_pos,
            breakend2_seq_order_read_pos,
            breakend1_neighbor: None,
            breakend2_neighbor: None,
        })
    }

    /// Check if candidate indel meets the minimum size, and if so push it into `all_bpo`
    ///
    /// # Arguments
    /// * `record` - bam record for error/debug messages
    /// * `read_pos` - Position in read coordinates on the alignment strand for the position immediately after the indel
    ///
    fn process_cand_indel(
        &self,
        record: &bam::Record,
        read_pos: usize,
        cand_indel: &mut CandidateIndel,
        all_bpo: &mut Vec<BreakpointObservation>,
    ) {
        if let Some(bpo) = self.process_cand_indel_impl(record, read_pos, cand_indel) {
            all_bpo.push(bpo);
        }
        cand_indel.clear();
    }

    /// Find large indels from the cigar string
    ///
    /// Fuse any adjacent deletion/insertion tags into one event and push qualifying events into `all_bpo`
    ///
    fn parse_breaks_from_cigar_alignment(
        &self,
        record: &bam::Record,
        all_bpo: &mut Vec<BreakpointObservation>,
    ) {
        use rust_htslib::bam::record::Cigar::*;

        let mut ref_pos = record.pos();
        let mut read_pos = 0;

        let mut cand_indel = CandidateIndel::new();

        for c in record.cigar().iter() {
            match c {
                Ins(len) => {
                    if cand_indel.is_empty() {
                        cand_indel.ref_pos = ref_pos - 1;
                    }
                    cand_indel.ins_size += *len;
                }
                Del(len) => {
                    if cand_indel.is_empty() {
                        cand_indel.ref_pos = ref_pos - 1;
                    }
                    cand_indel.del_size += *len;
                }
                _ => {
                    self.process_cand_indel(record, read_pos, &mut cand_indel, all_bpo);
                }
            }
            update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);

            // Bail out if we're past the end of the target region:
            if cand_indel.is_empty() && ref_pos >= self.worker_region.end {
                break;
            }
        }
        self.process_cand_indel(record, read_pos, &mut cand_indel, all_bpo);
    }

    /// Process 2 adjacent split read segments, ordered such that s1 comes before s2 in the read
    /// sequencing order.
    fn process_adjacent_split_segments(
        &self,
        s1: &SeqOrderSplitReadSegment,
        s2: &SeqOrderSplitReadSegment,
        all_bpo: &mut Vec<BreakpointObservation>,
    ) {
        if s1.mapq < self.min_sv_mapq as u8 || s2.mapq < self.min_sv_mapq as u8 {
            return;
        }

        fn get_breakend(s: &SeqOrderSplitReadSegment, is_s1: bool) -> Breakend {
            let right_anchor = s.is_fwd_strand ^ is_s1;
            let dir = if right_anchor {
                BreakendDirection::RightAnchor
            } else {
                BreakendDirection::LeftAnchor
            };
            let ref_pos = if right_anchor {
                s.pos
            } else {
                s.pos + get_cigar_ref_offset(&s.cigar)
            } - 1;
            Breakend {
                segment: GenomeSegment {
                    chrom_index: s.chrom_index,
                    range: IntRange::from_int(ref_pos),
                },
                dir,
            }
        }

        fn get_breakend_seq_order_read_pos(s: &SeqOrderSplitReadSegment, is_s1: bool) -> usize {
            let right_anchor = s.is_fwd_strand ^ is_s1;
            if right_anchor {
                s.seq_order_read_start
            } else {
                s.seq_order_read_end - 1
            }
        }

        let mut breakpoint = Breakpoint {
            breakend1: get_breakend(s1, true),
            breakend2: Some(get_breakend(s2, false)),
            insert_info: InsertInfo::SizeRange(IntRange::from_int(0)),
            is_precise: false,
        };
        breakpoint.standardize();
        all_bpo.push(BreakpointObservation {
            breakpoint,
            evidence: BreakpointEvidenceType::SplitRead,
            qname: None,
            breakend1_seq_order_read_pos: get_breakend_seq_order_read_pos(s1, true),
            breakend2_seq_order_read_pos: get_breakend_seq_order_read_pos(s2, false),
            breakend1_neighbor: None,
            breakend2_neighbor: None,
        });
    }

    /// Find general SV breakpoint signatures from any split reads annotated on the read
    ///
    fn parse_breaks_from_split_reads(
        &self,
        record: &bam::Record,
        all_bpo: &mut Vec<BreakpointObservation>,
    ) {
        // Only split reads are used in this method:
        const SA_AUX_TAG: &[u8] = b"SA";
        if !is_aux_tag_found(record, SA_AUX_TAG) {
            return;
        }

        // Procedure for parsing split read fragments:
        //
        // 1. First step is to put the split and primary alignment segments in order in read
        // coordinates and sequencer order
        //
        // 2. Then we can iterate through in order, and form a BreakpointObservation that can be
        // added to singles:
        //
        // self.observations.singles.add(bpo);
        //
        // How to determine 'ownership' when two regions both touch a split read?
        // - In non-targeted mode this is simple: we can just say that the bin containing the start
        // position of the primary read owns that whole read and enters all of the breakends.
        // - In targeted mode, perhaps every encounter with the read is entered, with qname attached
        // and the duplicate observations are removed later?
        //
        if self.is_targeted_scan {
            // Inclusion criteria are changed substantially in targeted mode. We'll defer this for
            // now and just understand that it will be broken for certain cases in targeted
            // analysis:
            //TODO: fix targeted split read search later to guarantee split reads are found in any
            // targeting arrangement
            if record.is_supplementary() {
                return;
            }
        } else {
            // First determine if this is the primary or supplemental alignment, we want to parse
            // all the breakpoints from a read once, so working exclusively with the primary
            // alignment should allow this.
            //
            // Requiring the read to start in the worker region prevents double entries (when
            // running in WGS mode).
            //
            let read_starts_in_worker_region = self.worker_region.intersect_pos(record.pos());
            if record.is_supplementary() || (!read_starts_in_worker_region) {
                return;
            }
        }

        let seq_order_read_split_segments =
            get_seq_order_read_split_segments(self.chrom_list, record);

        let debug = false;
        if debug {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            eprintln!("SA segments for read: {qname}");
            for seg in seq_order_read_split_segments.iter() {
                eprintln!("Sorted Segment: {:?}", seg);
            }
        }

        // Now iterate through the sorted list and convert each junction into a breakpoint observation:
        //
        let split_count = seq_order_read_split_segments.len();
        for split_index in 0..(split_count - 1) {
            self.process_adjacent_split_segments(
                &seq_order_read_split_segments[split_index],
                &seq_order_read_split_segments[split_index + 1],
                all_bpo,
            );
        }
    }

    /// Find general single-side SV breakpoints from read soft-clips
    ///
    fn parse_breaks_from_read_soft_clipping(
        &self,
        record: &bam::Record,
        all_bpo: &mut Vec<BreakpointObservation>,
    ) {
        // Criteria
        // - Read has at least N bases of soft clip
        // - Read is not soft-clipped on both ends.
        //
        // Originally, this method filtered out SA tag reads but that looks like it will remove most
        // insertions, so first see how many we can call, then go back and work on the coordination
        // vs. BND calls
        //
        let large_insertion_soft_clip_state =
            test_read_for_large_insertion_soft_clip(record, self.min_large_insertion_soft_clip_len);

        use LargeInsertionSoftClipState::*;
        let dir = match large_insertion_soft_clip_state {
            Null => {
                return;
            }
            Left => BreakendDirection::RightAnchor,
            Right => BreakendDirection::LeftAnchor,
        };

        let ref_pos = record.pos()
            + if large_insertion_soft_clip_state == Left {
                0
            } else {
                get_cigar_ref_offset(&record.cigar())
            };

        // Get the sequenced strand read position of the soft-clip breakpoint:
        let breakend1_seq_order_read_pos = {
            let (left_sclip_len, right_sclip_start, _read_size) =
                get_hard_clipped_read_clip_positions(&record.cigar());
            let read_pos = if large_insertion_soft_clip_state == Left {
                left_sclip_len
            } else {
                right_sclip_start - 1
            };
            get_seq_order_read_position(record, read_pos)
        };

        // Get a theoretical breakend2 position, even though we don't have a breakend location here:
        let breakend2_seq_order_read_pos =
            if (large_insertion_soft_clip_state == Left) ^ record.is_reverse() {
                breakend1_seq_order_read_pos - 1
            } else {
                breakend1_seq_order_read_pos + 1
            };

        let segment = GenomeSegment {
            chrom_index: record.tid() as usize,
            range: IntRange::from_int(ref_pos),
        };

        let breakend1 = Breakend { segment, dir };

        let breakpoint = Breakpoint {
            breakend1,
            breakend2: None,
            insert_info: InsertInfo::NoInsert,
            is_precise: false,
        };

        let debug = false;
        if debug {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            eprintln!("record {qname} accepted as soft-clip candidate. Breakpoint {breakpoint:?}");
        }

        all_bpo.push(BreakpointObservation {
            breakpoint,
            evidence: BreakpointEvidenceType::SoftClip,
            qname: None,
            breakend1_seq_order_read_pos,
            breakend2_seq_order_read_pos,
            breakend1_neighbor: None,
            breakend2_neighbor: None,
        });
    }

    /// Check for and annotate neighboring breakends
    fn annotate_neighboring_breakends(
        all_bpo: &mut [BreakpointObservation],
        max_close_breakend_distance: usize,
    ) {
        // Step 1: Come up with a filtered bpo list that removes cases we won't be considering for neighbor analysis:
        let min_neighbor_deletion_size = 100;
        let mut cand_bpo = all_bpo
            .iter_mut()
            .filter(|x| x.evidence != BreakpointEvidenceType::SoftClip)
            .filter(|x| match x.breakpoint.deletion_size() {
                Some(x) => x >= min_neighbor_deletion_size,
                None => true,
            })
            .collect::<Vec<_>>();

        // Step 2: Walk through filtered list in read pos order and test adjacent cases
        cand_bpo.sort_by_key(|x| x.breakend1_seq_order_read_pos);

        fn get_breakend_index(bpo: &BreakpointObservation, is_upstream_breakend: bool) -> usize {
            let is_be0_upstream =
                bpo.breakend1_seq_order_read_pos > bpo.breakend2_seq_order_read_pos;
            if is_be0_upstream ^ is_upstream_breakend {
                1
            } else {
                0
            }
        }

        fn get_breakend(bpo: &BreakpointObservation, is_upstream_breakend: bool) -> &Breakend {
            if get_breakend_index(bpo, is_upstream_breakend) > 0 {
                bpo.breakpoint.breakend2.as_ref().unwrap()
            } else {
                &bpo.breakpoint.breakend1
            }
        }

        let cand_bpo_count = cand_bpo.len();
        for index1 in 1..cand_bpo_count {
            let index0 = index1 - 1;
            {
                let be0 = get_breakend(cand_bpo[index0], true);
                let be1 = get_breakend(cand_bpo[index1], false);

                // Check if the breakends are close:
                match get_segment_distance(&be0.segment, &be1.segment) {
                    Some(x) => {
                        if x > max_close_breakend_distance {
                            continue;
                        };
                    }
                    None => continue,
                };

                // Check for expected breakend arrangement
                let (downstream_be, upstream_be) =
                    if be0.segment.range.start < be1.segment.range.start {
                        (be0, be1)
                    } else {
                        (be1, be0)
                    };
                if downstream_be.dir != BreakendDirection::RightAnchor
                    || upstream_be.dir != BreakendDirection::LeftAnchor
                {
                    continue;
                }
            }

            // Connect the breakends
            let bpo0_be_index = get_breakend_index(cand_bpo[index0], true);
            let bpo1_be_index = get_breakend_index(cand_bpo[index1], false);

            fn update_target(
                cand_bpo: &mut [&mut BreakpointObservation],
                target_index: usize,
                neighbor_index: usize,
                target_be_index: usize,
                neighbor_be_index: usize,
            ) {
                let target_neighbor = {
                    let bp = &cand_bpo[neighbor_index].breakpoint;
                    if neighbor_be_index == 0 {
                        Some(bp.breakend1.clone())
                    } else {
                        bp.breakend2.clone()
                    }
                };

                if target_be_index == 0 {
                    cand_bpo[target_index].breakend1_neighbor = target_neighbor;
                } else {
                    cand_bpo[target_index].breakend2_neighbor = target_neighbor;
                }
            }

            update_target(&mut cand_bpo, index0, index1, bpo0_be_index, bpo1_be_index);
            update_target(&mut cand_bpo, index1, index0, bpo1_be_index, bpo0_be_index);
        }
    }

    fn process_read_bpos(
        record: &bam::Record,
        all_bpo: &mut [BreakpointObservation],
        max_close_breakend_distance: usize,
    ) {
        let std_bp_count = all_bpo
            .iter()
            .filter(|x| x.evidence != BreakpointEvidenceType::SoftClip)
            .count();

        if std_bp_count <= 1 {
            return;
        }

        // qname label all bpos to prevent them from being counted as independent evidence in downstream steps
        for bpo in all_bpo.iter_mut() {
            bpo.qname = Some(record.qname().to_vec());
        }

        Self::annotate_neighboring_breakends(all_bpo, max_close_breakend_distance);
    }

    pub fn process_bam_record(&mut self, record: &bam::Record) {
        /*
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        eprintln!("Break builder processing read {}", qname);
         */
        if record.mapq() < self.min_sv_mapq as u8 {
            return;
        }

        // Get all breakpoint observations for this read:
        let mut all_bpo = Vec::new();
        self.parse_breaks_from_cigar_alignment(record, &mut all_bpo);
        self.parse_breaks_from_split_reads(record, &mut all_bpo);
        self.parse_breaks_from_read_soft_clipping(record, &mut all_bpo);

        // If multiple indel/split candidates are found from one read, then update the bpos with
        // more information for downstream phasing and dependent observation analysis:
        Self::process_read_bpos(record, &mut all_bpo, self.max_close_breakend_distance);

        for bpo in all_bpo {
            self.observations.singles.add(bpo);
        }
    }

    /// Process any buffered data remaining in the object, then return the summarized results
    ///
    pub fn complete_processing(&mut self) -> BreakObservations {
        std::mem::replace(&mut self.observations, BreakObservations::new())
    }
}
