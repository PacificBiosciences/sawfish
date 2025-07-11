use std::collections::BTreeSet;
use std::fmt;

use rust_vc_utils::{ChromList, rev_comp_in_place};
pub use strum::EnumCount;

use crate::genome_segment::{GenomeSegment, get_segment_dir_distance, get_segment_distance};
use crate::int_range::IntRange;

/// Direction of a breakend
///
/// 'LeftAnchor' means that the read to the left side of the breakend is locally mapped.
/// A 'LeftAnchor' breakend would correspond to the left side of a simple deletion
///
#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, strum::EnumCount)]
pub enum BreakendDirection {
    LeftAnchor,
    RightAnchor,
}

/// A single breakend, half of a breakpoint observation
///
#[derive(Clone, Eq, PartialEq, PartialOrd, Ord)]
pub struct Breakend {
    /// The segment range can represent the uncertainty in the breakend location (if this is a
    /// non-precise breakend candidate), or breakend homology in a refined candidate.
    ///
    /// By convention, any given position in the range represents the position immediately before
    /// the breakend occurs.
    ///
    /// For example, using the convention of zero-indexing all internal coordinates: a range of
    /// [99,100) would represent a breakend occurring between internal positions 99 and 100, or
    /// positions 100 and 101 in a 1-indexed output like VCF.
    ///
    /// One important implication of this convention is that a right-anchored breakend can have a
    /// start position of -1. This represents a BND ligating the contig end to another location starting
    /// from its first base. This case is actually somewhat common to find in the decoy sequences.
    ///
    pub segment: GenomeSegment,
    pub dir: BreakendDirection,
}

impl Breakend {
    pub fn new() -> Self {
        Self {
            segment: GenomeSegment::new(),
            dir: BreakendDirection::LeftAnchor,
        }
    }

    /// Get homology length
    ///
    /// This is pretty simple but standardizing to a function should help to reduce off by one
    /// errors, since this is used so frequently.
    ///
    pub fn get_homology_len(&self) -> i64 {
        let hom_len = self.segment.range.size() - 1;
        assert!(hom_len >= 0);
        hom_len
    }

    /// Return false if the breakend is improperly arranged off the edge of the contig
    ///
    /// See documentation of Breakend type above for notes on the -1 position adjustments for
    /// right-anchored breakends.
    ///
    pub fn is_valid_range(&self, chrom_list: &ChromList) -> bool {
        let chrom_size = chrom_list.data[self.segment.chrom_index].length as i64;
        let (min_start, max_end) = if self.dir == BreakendDirection::RightAnchor {
            (-1, chrom_size - 1)
        } else {
            (0, chrom_size)
        };
        let range = &self.segment.range;
        range.start < range.end && range.start >= min_start && range.end <= max_end
    }

    pub fn merge(&mut self, other: &Breakend) {
        assert_eq!(self.dir, other.dir);
        assert_eq!(self.segment.chrom_index, other.segment.chrom_index);
        self.segment.range.merge(&other.segment.range);
    }
}

impl fmt::Debug for Breakend {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Breakend: {:?} dir: {:?}", self.segment, self.dir)
    }
}

/// Get distance between two breakends, the distance is the number of bases separating the two
/// ranges.
///
/// Returns None if the breakends are on different chroms or have different directions. Returns a
/// distance of zero if they are adjacent or intersect.
///
pub fn get_breakend_distance(be1: &Breakend, be2: &Breakend) -> Option<usize> {
    if be1.dir != be2.dir {
        None
    } else {
        get_segment_distance(&be1.segment, &be2.segment)
    }
}

/// Get (direction, distance) between two breakends, the distance is the number of bases separating
/// the two ranges, and the direction is true if be2 is 'ahead' of be1.
///
/// Returns None if the breakends are on different chroms or have different directions. Returns a
/// distance of zero if they are adjacent or intersect.
///
pub fn get_breakend_dir_distance(be1: &Breakend, be2: &Breakend) -> Option<(bool, usize)> {
    if be1.dir != be2.dir {
        None
    } else {
        get_segment_dir_distance(&be1.segment, &be2.segment)
    }
}

#[derive(Clone, Copy, Debug, EnumCount, Eq, PartialEq, PartialOrd, Ord)]
pub enum BreakpointEvidenceType {
    ReadAlignment,
    SplitRead,
    SoftClip,
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd)]
pub enum InsertInfo {
    /// Represents an insertion when the size is uncertain and the sequence content is not known
    SizeRange(IntRange),

    /// A fully assembled insertion haplotype
    Seq(Vec<u8>),

    /// Use this instead of SizeRange(0,0) to indicate that the SV is known to have no insertion
    NoInsert,
}

impl InsertInfo {
    pub fn size(&self) -> usize {
        match self {
            InsertInfo::SizeRange(r) => ((r.start + r.end) / 2) as usize,
            InsertInfo::Seq(s) => s.len(),
            InsertInfo::NoInsert => 0,
        }
    }

    pub fn max_size(&self) -> usize {
        match self {
            InsertInfo::SizeRange(r) => r.end as usize,
            InsertInfo::Seq(s) => s.len(),
            InsertInfo::NoInsert => 0,
        }
    }

    /// Merge insertion size ranges
    pub fn merge(&mut self, other: &InsertInfo) {
        if matches!(self, InsertInfo::Seq(_)) || matches!(other, InsertInfo::Seq(_)) {
            panic!("Unexpected insert_info format for merge");
        }

        if let InsertInfo::SizeRange(osr) = other {
            if let InsertInfo::SizeRange(ssr) = self {
                ssr.merge(osr);
            } else if let InsertInfo::NoInsert = self {
                *self = other.clone();
            }
        }
    }
}

impl fmt::Debug for InsertInfo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InsertInfo::SizeRange(x) => {
                write!(f, "SizeRange({x:?})")
            }
            InsertInfo::Seq(x) => {
                write!(f, "Seq({:?})", std::str::from_utf8(x).unwrap())
            }
            InsertInfo::NoInsert => {
                write!(f, "NoInsert")
            }
        }
    }
}

#[derive(Clone, Eq, PartialEq, PartialOrd, Ord)]
pub struct Breakpoint {
    pub breakend1: Breakend,
    pub breakend2: Option<Breakend>,

    /// Information on inserted sequence between the breakends, either sequence or just length
    ///
    /// In case both breakends are anchored in the same direction, the inserted sequence should be
    /// oriented to match the direction of breakend1
    ///
    pub insert_info: InsertInfo,

    /// Once a breakpoint is marked as precise, its location should be estimated to single base
    /// accuracy, and the location range is interpreted as breakend homology
    pub is_precise: bool,
}

impl Breakpoint {
    pub fn new() -> Self {
        Self {
            breakend1: Breakend::new(),
            breakend2: None,
            insert_info: InsertInfo::SizeRange(IntRange::from_int(0)),
            is_precise: false,
        }
    }

    pub fn get_breakend(&self, breakend_index: usize) -> &Breakend {
        if breakend_index == 0 {
            &self.breakend1
        } else {
            self.breakend2.as_ref().unwrap()
        }
    }

    /// If true, the SV has an inverted breakend pattern. This could be due to a literal inversion,
    /// an inverted breakend pattern part of a complex local SV, or an inverted pattern in a
    /// translocation.
    ///
    /// Extra considerations for the same orientation state include (1) as breakend location is
    /// adjusted, the two breakends move in opposite directions (2) any inserted sequence has to
    /// be represented wrt the strand of either one of the breakends, but can't be represented for
    /// both.
    ///
    pub fn same_orientation(&self) -> bool {
        if self.breakend2.is_none() {
            return false;
        }
        self.breakend1.dir == self.breakend2.as_ref().unwrap().dir
    }

    /// If true, both breakends are defined, on the same chromosome and with a breakend direction
    /// pattern consistent with an indel
    pub fn is_indel(&self) -> bool {
        assert!(self.is_standardized());
        if self.breakend2.is_none() {
            return false;
        }
        let bnd1 = &self.breakend1;
        let bnd2 = self.breakend2.as_ref().unwrap();
        if bnd1.segment.chrom_index != bnd2.segment.chrom_index {
            return false;
        }
        bnd1.dir == BreakendDirection::LeftAnchor && bnd2.dir == BreakendDirection::RightAnchor
    }

    pub fn deletion_size(&self) -> Option<usize> {
        if !self.is_indel() {
            None
        } else {
            let be1 = &self.breakend1;
            let be2 = self.breakend2.as_ref().unwrap();
            let size = std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) as usize;
            Some(size)
        }
    }

    /// If true, both breakends are defined, on the same chromosome and with a breakend direction
    /// pattern consistent with an indel or duplication
    pub fn is_indel_or_dup(&self) -> bool {
        assert!(self.is_standardized());
        if self.breakend2.is_none() {
            return false;
        }
        let bnd1 = &self.breakend1;
        let bnd2 = self.breakend2.as_ref().unwrap();
        if bnd1.segment.chrom_index != bnd2.segment.chrom_index {
            return false;
        }
        bnd1.dir != bnd2.dir
    }

    pub fn deldup_size(&self) -> Option<usize> {
        if !self.is_indel_or_dup() {
            None
        } else {
            let be1 = &self.breakend1;
            let be2 = self.breakend2.as_ref().unwrap();
            let size = std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) as usize;
            Some(size)
        }
    }

    /// Test for standardized order of breakends 1 and 2
    ///
    /// If breakend2 is defined,
    /// then breakend2 must have a lower segment order,
    /// or they must be ordered Left,Right if segments are the same and in opposite directions
    ///
    /// The last condition is important to normalize insertions
    ///
    pub fn is_standardized(&self) -> bool {
        self.breakend2.is_none() || {
            let breakend2 = self.breakend2.as_ref().unwrap();
            self.breakend1.segment < breakend2.segment
                || self.breakend1.segment == breakend2.segment
                    && !(self.breakend1.dir == BreakendDirection::RightAnchor
                        && breakend2.dir == BreakendDirection::LeftAnchor)
        }
    }

    /// Standardize order of breakends 1 and 2
    ///
    pub fn standardize(&mut self) {
        if !self.is_standardized() {
            // flip breakends:
            std::mem::swap(&mut self.breakend1, self.breakend2.as_mut().unwrap());

            // revcomp the inserted sequence if required:
            if self.same_orientation() {
                if let InsertInfo::Seq(seq) = &mut self.insert_info {
                    rev_comp_in_place(seq);
                }
            }
        }
    }

    pub const DIR_TYPE_COUNT: usize =
        BreakendDirection::COUNT + BreakendDirection::COUNT * BreakendDirection::COUNT;

    pub fn get_dir_type_index(&self) -> usize {
        let i1 = self.breakend1.dir as usize;
        if self.breakend2.is_none() {
            i1
        } else {
            let i2 = self.breakend2.as_ref().unwrap().dir as usize;
            BreakendDirection::COUNT + i1 * BreakendDirection::COUNT + i2
        }
    }

    pub fn merge(&mut self, other: &Breakpoint) {
        // Assume from breakpoint standardization that we always match breakend 1 to 1 and 2 to 2
        //
        // TODO, should we assert on standardization here?
        //
        assert!(!(self.is_precise || other.is_precise));

        self.breakend1.merge(&other.breakend1);
        assert_eq!(self.breakend2.is_none(), other.breakend2.is_none());
        if let (Some(this_bp2), Some(other_bp2)) =
            (self.breakend2.as_mut(), other.breakend2.as_ref())
        {
            this_bp2.merge(other_bp2);
        }

        self.insert_info.merge(&other.insert_info);
    }
}

impl fmt::Debug for Breakpoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Breakpoint:\n be1: {:?}\n be2: {:?}\n insert: {:?}\n insert_size: {}\nis_precise: {}",
            self.breakend1,
            self.breakend2,
            self.insert_info,
            self.insert_info.size(),
            self.is_precise
        )
    }
}

/*
/// This is a more general version of the distance function which cross checks breakend2 matching
/// breakend1 of the other breakpoint. With breakpoint standardization its possible we can assume
/// this case never needs checking.
///
pub fn get_breakpoint_manhattan_distance(bp1: &Breakpoint, bp2: &Breakpoint) -> Option<usize> {
    // Handle single breakend breakpoints first:
    if bp1.breakend2.is_none() && bp2.breakend2.is_none() {
        get_breakend_distance(&bp1.breakend1, &bp2.breakend1)
    } else if bp1.breakend2.is_none() || bp2.breakend2.is_none() {
        None
    } else {
        // Normal breakpoints
        let bp1_be2 = bp1.breakend2.as_ref().unwrap();
        let bp2_be2 = bp2.breakend2.as_ref().unwrap();
        let match_dist = {
            let d11 = get_breakend_distance(&bp1.breakend1, &bp2.breakend1);
            let d22 = get_breakend_distance(bp1_be2, bp2_be2);
            if let (Some(d11), Some(d22)) = (d11, d22) {
                Some(d11 + d22)
            } else {
                None
            }
        };

        let switch_dist = {
            let d12 = get_breakend_distance(&bp1.breakend1, bp2_be2);
            let d21 = get_breakend_distance(bp1_be2, &bp2.breakend1);
            if let (Some(d12), Some(d21)) = (d12, d21) {
                Some(d12 + d21)
            } else {
                None
            }
        };

        if let (Some(match_dist), Some(switch_dist)) = (match_dist, switch_dist) {
            Some(std::cmp::min(match_dist, switch_dist))
        } else if match_dist.is_some() {
            match_dist
        } else if switch_dist.is_some() {
            switch_dist
        } else {
            None
        }
    }
}
 */

/// Breakpoint 'manhattan' distance finds the actual manhattan distance of breakends 1 and 2 only in
/// the worst case. Any movement of the second breakend in concert with the first, in the expected
/// direction for the direction pattern, is free.
///
pub fn get_breakpoint_manhattan_distance(bp1: &Breakpoint, bp2: &Breakpoint) -> Option<usize> {
    assert!(bp1.is_standardized() && bp2.is_standardized());
    assert_eq!(
        bp1.get_dir_type_index(),
        bp2.get_dir_type_index(),
        "No distance metric for different breakpoint dir types"
    );

    if bp1.breakend2.is_none() && bp2.breakend2.is_none() {
        // Handle single breakend breakpoints
        get_breakend_distance(&bp1.breakend1, &bp2.breakend1)
    } else {
        // Handle normal breakpoints:
        let bp1_be2 = bp1.breakend2.as_ref().unwrap();
        let bp2_be2 = bp2.breakend2.as_ref().unwrap();
        let d11 = get_breakend_dir_distance(&bp1.breakend1, &bp2.breakend1);
        let d22 = get_breakend_dir_distance(bp1_be2, bp2_be2);
        if let (Some((dir11, dist11)), Some((dir22, dist22))) = (d11, d22) {
            // When the two breakends of the breakpoint shift together in the expected pattern
            // relative to another breakpoint, the distance between those representations is
            // reduced, as it is more likely they correspond to the same event.
            //
            // For example, if bp1 represents a 100 base deletion starting at position 100, and bp2 represents
            // a 100 base deletion starting at position 110., the two breakends both shift in the expected direction
            // together, so the distance between the breakpoints is 10.
            //
            // Extending this example, if bp1 represents a 100 base deletion starting at position 100, and bp2 represents
            // a 120 base deletion start at position 90, now the breakends are each 10 bases appart as in the previous
            // example, but they've moved in opposite directions, so are further penalized with manahatan distances across
            // the two breakends, to define a distance between the breakpoints of 20.
            //
            let is_expected_shift_pattern = bp1.same_orientation() ^ (dir11 == dir22);
            if is_expected_shift_pattern {
                Some(std::cmp::max(dist11, dist22))
            } else {
                Some(dist11 + dist22)
            }
        } else {
            None
        }
    }
}

/// A single breakpoint observation
///
/// A breakend 'observation' is a breakpoint proposal from a single piece of evidence
#[derive(Eq, PartialEq, PartialOrd, Ord)]
pub struct BreakpointObservation {
    pub breakpoint: Breakpoint,
    pub evidence: BreakpointEvidenceType,
    pub qname: Option<Vec<u8>>,

    /// Sequencer-order read position of the last mapped base on the anchor side of breakend1
    pub breakend1_seq_order_read_pos: usize,

    /// Sequencer-order read position of the last mapped base on the anchor side of breakend2
    pub breakend2_seq_order_read_pos: usize,

    /// If a close breakend was found on the same haplotype in this read, record it here for the breakend 1 or 2 side of the breakpoint
    /// so that we can account for these later.
    ///
    pub breakend1_neighbor: Option<Breakend>,
    pub breakend2_neighbor: Option<Breakend>,
}

impl fmt::Debug for BreakpointObservation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "bp: {:?} evidence_type: {:?}",
            self.breakpoint, self.evidence
        )
    }
}

#[derive(Clone)]
pub struct ConsolidatedAssemblySegmentInfo {
    pub segment: GenomeSegment,
    pub insert_info: InsertInfo,
    pub cluster_ids: Vec<usize>,
}

#[derive(Clone)]
pub struct LargeInsertionInfo {
    /// Segment spanning between the two soft-clip clusters composing the insertion candidate
    ///
    pub segment: GenomeSegment,

    /// First soft-clip cluster index of the insertion
    ///
    /// ie. The first cluster encountered during a left to right scan
    ///
    pub root_cluster_index: usize,

    /// Second soft-clip cluster index of the insertion
    ///
    /// ie. The second cluster encountered during a left to right scan
    ///
    #[allow(unused)]
    pub paired_cluster_index: usize,
}

#[derive(Clone)]
pub struct NeighborBreakend {
    pub cluster_index: usize,
    pub breakend_index: usize,
    pub breakpoint: Breakpoint,
}

/// A cluster of breakpoint observations
///
/// A cluster aggregates multiple similar breakpoint observations.
///
#[derive(Clone)]
pub struct BreakpointCluster {
    pub breakpoint: Breakpoint,

    /// Count the types of direct read-based breakpoint observations supporting this breakpoint candidate
    pub evidence: [usize; BreakpointEvidenceType::COUNT],

    /// qnames for a subset of the evidence reads
    ///
    /// qnames are only stored for reads which have multiple breakpoint observations (currently ignoring softclip breakends)
    ///
    /// These are used:
    /// (1) to prevent double-counting of these observations as independent support of the same haplotype
    /// (2) for some forms of early hapltoype phasing, such as modifying the alt allele model based on a close
    /// neighboring breakend
    ///
    pub evidence_qnames: BTreeSet<Vec<u8>>,

    pub assembly_segments: Vec<GenomeSegment>,

    /// In some cases indel Breakpoint clusters are so close together that these are
    /// consolidated into one assembly. If so the assembly region is provided here for one
    /// representative breakpoint cluster.
    ///
    /// This consolidation does not include assembly flank windows, so it can be treated similarly
    /// to the assembly segments above.
    ///
    pub consolidated_assembly_segment: Option<ConsolidatedAssemblySegmentInfo>,

    /// If this cluster is part of a consolidation, point to the consolidated cluster id here:
    pub consolidated_cluster_index: Option<usize>,

    pub large_insertion_info: Option<LargeInsertionInfo>,

    pub breakend1_neighbor_observations: Vec<Breakend>,
    pub breakend2_neighbor_observations: Vec<Breakend>,

    pub breakend1_neighbor: Option<NeighborBreakend>,
    pub breakend2_neighbor: Option<NeighborBreakend>,
}

impl BreakpointCluster {
    pub fn new() -> Self {
        Self {
            breakpoint: Breakpoint::new(),
            evidence: [0; BreakpointEvidenceType::COUNT],
            evidence_qnames: BTreeSet::new(),
            assembly_segments: Vec::new(),
            consolidated_assembly_segment: None,
            consolidated_cluster_index: None,
            large_insertion_info: None,
            breakend1_neighbor_observations: Vec::new(),
            breakend2_neighbor_observations: Vec::new(),
            breakend1_neighbor: None,
            breakend2_neighbor: None,
        }
    }

    pub fn from_breakpoint_observation(bpo: &BreakpointObservation) -> Self {
        let mut bpc = BreakpointCluster::new();
        bpc.merge_breakpoint_observation(bpo);
        bpc
    }

    pub fn is_empty(&self) -> bool {
        self.evidence_count() == 0
    }

    pub fn evidence_count(&self) -> usize {
        self.evidence.iter().sum()
    }

    pub fn merge_breakpoint_observation(&mut self, bpo: &BreakpointObservation) {
        let is_new_evidence = if let Some(qname) = &bpo.qname {
            self.evidence_qnames.insert(qname.clone())
        } else {
            true
        };

        if self.is_empty() {
            self.breakpoint = bpo.breakpoint.clone();
        } else {
            self.breakpoint.merge(&bpo.breakpoint);
        }
        if is_new_evidence {
            self.evidence[bpo.evidence as usize] += 1;
        }

        if let Some(x) = &bpo.breakend1_neighbor {
            self.breakend1_neighbor_observations.push(x.clone());
        }
        if let Some(x) = &bpo.breakend2_neighbor {
            self.breakend2_neighbor_observations.push(x.clone());
        }
    }
}

impl fmt::Debug for BreakpointCluster {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "bp: {:?} evidence: {:?}", self.breakpoint, self.evidence)
    }
}

/// Higher level SV classifications that can be inferred from a single breakpoint
///
/// These are the classifications that will be used for VCF output
///
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VcfSVType {
    Deletion,
    Insertion,
    Duplication,
    Inversion,
    Breakpoint,
    SingleBreakend,

    /// This type is only used for a variant that can't be described as a deletion or duplication
    CopyNumberVariant,
}

/// Get chromosome span from breakend 1 to 2, if it can be defined for this breakpoint
///
pub fn get_breakpoint_span(bp: &Breakpoint) -> Option<usize> {
    let be2 = bp.breakend2.as_ref()?;
    let be1 = &bp.breakend1;
    if be1.segment.chrom_index != be2.segment.chrom_index {
        None
    } else {
        Some(std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) as usize)
    }
}

pub fn get_breakpoint_vcf_sv_type(bp: &Breakpoint) -> VcfSVType {
    use VcfSVType::*;
    if let Some(be2) = &bp.breakend2 {
        let be1 = &bp.breakend1;
        let chrom_diff = be1.segment.chrom_index != be2.segment.chrom_index;
        let same_be_dirs = be1.dir == be2.dir;
        if chrom_diff || same_be_dirs {
            Breakpoint
        } else if be1.dir == BreakendDirection::RightAnchor {
            Duplication
        } else {
            let ins_size = bp.insert_info.size();
            let del_size =
                std::cmp::max(be2.segment.range.start - be1.segment.range.start, 0) as usize;
            if ins_size > del_size {
                Insertion
            } else {
                Deletion
            }
        }
    } else {
        SingleBreakend
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_breakpoint_manhattan_distance() {
        let bp1 = Breakpoint {
            breakend1: Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(100),
                },
                dir: BreakendDirection::LeftAnchor,
            },
            breakend2: Some(Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(200),
                },
                dir: BreakendDirection::RightAnchor,
            }),
            insert_info: InsertInfo::NoInsert,
            is_precise: false,
        };

        let bp2 = Breakpoint {
            breakend1: Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(110),
                },
                dir: BreakendDirection::LeftAnchor,
            },
            breakend2: Some(Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(210),
                },
                dir: BreakendDirection::RightAnchor,
            }),
            insert_info: InsertInfo::NoInsert,
            is_precise: false,
        };
        let bp_dist = get_breakpoint_manhattan_distance(&bp1, &bp2);
        assert_eq!(bp_dist, Some(9));

        let bp2 = Breakpoint {
            breakend1: Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(90),
                },
                dir: BreakendDirection::LeftAnchor,
            },
            breakend2: Some(Breakend {
                segment: GenomeSegment {
                    chrom_index: 0,
                    range: IntRange::from_int(210),
                },
                dir: BreakendDirection::RightAnchor,
            }),
            insert_info: InsertInfo::NoInsert,
            is_precise: false,
        };
        let bp_dist = get_breakpoint_manhattan_distance(&bp1, &bp2);
        assert_eq!(bp_dist, Some(18));
    }
}
