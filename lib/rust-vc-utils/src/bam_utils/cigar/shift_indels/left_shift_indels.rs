use rust_htslib::bam::record::Cigar;

use super::cigar_indel_shifter::{CigarShiftBuilder, ShiftDirection};
use crate::cigar::{clean_up_cigar_edge_indels, compress_cigar, update_ref_and_read_pos};

/// Left-shift all indels in the input cigar
///
/// Does not preserve X/= match states in left-shifted output
///
/// Indels may be fused into combined into insertion deletion events. All such events
/// will be output in "nImD" format (ie. insertion before deletion state in cigar string)
///
/// Returns a 2-tuple of:
/// 1. Shifted ref pos
/// 2. Shifted alignment cigar
///
pub fn left_shift_indels(
    ref_pos: i64,
    cigar: &[Cigar],
    ref_seq: &[u8],
    read_seq: &[u8],
) -> (i64, Vec<Cigar>) {
    let ignore_hard_clip = false;

    let mut ref_head_pos = ref_pos;
    let mut read_head_pos = 0;

    let mut cigar_block_info = CigarShiftBuilder::new(ShiftDirection::Left, ref_seq, read_seq);

    for c in cigar.iter() {
        cigar_block_info.add_element(c, ref_head_pos, read_head_pos);
        update_ref_and_read_pos(c, &mut ref_head_pos, &mut read_head_pos, ignore_hard_clip);
    }
    let mut shift_cigar = cigar_block_info.get_cigar();

    let ref_pos_shift = clean_up_cigar_edge_indels(&mut shift_cigar);
    let shift_cigar = compress_cigar(&shift_cigar);
    (ref_pos + ref_pos_shift as i64, shift_cigar)
}
