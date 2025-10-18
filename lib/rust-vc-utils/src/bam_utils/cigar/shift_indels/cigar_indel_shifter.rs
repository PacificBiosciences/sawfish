use crate::{IntRange, get_indel_breakend_homology_info};
use rust_htslib::bam::record::Cigar;

#[derive(PartialEq)]
pub enum ShiftDirection {
    Left,
    Right,
}

pub struct CigarShiftBuilder<'a> {
    dir: ShiftDirection,

    ref_seq: &'a [u8],
    read_seq: &'a [u8],

    match_block_size: u32,

    is_in_indel_block: bool,

    indel_block_ref_start: i64,
    indel_block_read_start: usize,

    indel_block_del_size: u32,
    indel_block_ins_size: u32,

    shift_cigar: Vec<Cigar>,
}

impl<'a> CigarShiftBuilder<'a> {
    pub fn new(dir: ShiftDirection, ref_seq: &'a [u8], read_seq: &'a [u8]) -> Self {
        Self {
            dir,
            ref_seq,
            read_seq,
            match_block_size: 0,
            is_in_indel_block: false,
            indel_block_ref_start: 0,
            indel_block_read_start: 0,
            indel_block_del_size: 0,
            indel_block_ins_size: 0,
            shift_cigar: Vec::new(),
        }
    }

    pub fn add_element(&mut self, cigar_seg: &Cigar, ref_pos: i64, read_pos: usize) {
        use Cigar::*;
        match cigar_seg {
            Del(len) => self.add_del(*len, ref_pos, read_pos),
            Ins(len) => self.add_ins(*len, ref_pos, read_pos),
            Match(len) | Equal(len) | Diff(len) => self.add_match(*len),
            _ => self.add_other(Some(cigar_seg)),
        }
    }

    pub fn get_cigar(&mut self) -> Vec<Cigar> {
        self.add_other(None);
        if ShiftDirection::Right == self.dir {
            self.shift_cigar.reverse();
        }
        std::mem::take(&mut self.shift_cigar)
    }

    fn add_indel(&mut self, ref_pos: i64, read_pos: usize) {
        if ShiftDirection::Right == self.dir || !self.is_in_indel_block {
            self.indel_block_ref_start = ref_pos;
            self.indel_block_read_start = read_pos;
            if !self.is_in_indel_block {
                self.is_in_indel_block = true;
            }
        }
    }

    fn add_del(&mut self, len: u32, ref_pos: i64, read_pos: usize) {
        if len > 0 {
            self.add_indel(ref_pos, read_pos);
            self.indel_block_del_size += len;
        }
    }

    fn add_ins(&mut self, len: u32, ref_pos: i64, read_pos: usize) {
        if len > 0 {
            self.add_indel(ref_pos, read_pos);
            self.indel_block_ins_size += len;
        }
    }

    fn push_del_segment(&mut self) {
        if self.indel_block_del_size > 0 {
            self.shift_cigar.push(Cigar::Del(self.indel_block_del_size));
            self.indel_block_del_size = 0;
        }
    }

    fn push_ins_segment(&mut self) {
        if self.indel_block_ins_size > 0 {
            self.shift_cigar.push(Cigar::Ins(self.indel_block_ins_size));
            self.indel_block_ins_size = 0;
        }
    }

    fn end_indel(&mut self) {
        if !self.is_in_indel_block {
            return;
        }

        self.is_in_indel_block = false;

        let shift_len = {
            let ref_range = IntRange::from_pair(
                self.indel_block_ref_start,
                self.indel_block_ref_start + self.indel_block_del_size as i64,
            );
            let read_range = IntRange::from_pair(
                self.indel_block_read_start as i64,
                self.indel_block_read_start as i64 + self.indel_block_ins_size as i64,
            );
            let (range, _) = get_indel_breakend_homology_info(
                self.ref_seq,
                &ref_range,
                self.read_seq,
                &read_range,
            );
            std::cmp::max(
                0,
                match self.dir {
                    ShiftDirection::Left => -range.start,
                    ShiftDirection::Right => range.end,
                },
            ) as u32
        };

        let actual_shift_len = std::cmp::min(self.match_block_size, shift_len);
        let shifted_match_block_size = self.match_block_size - actual_shift_len;
        if shifted_match_block_size > 0 {
            self.shift_cigar
                .push(Cigar::Match(shifted_match_block_size));
        }
        self.match_block_size = actual_shift_len;

        // Arrange order so that combined insertion/deletion events will always be output in "nImD" format
        if self.dir == ShiftDirection::Left {
            self.push_ins_segment();
        }
        self.push_del_segment();
        if self.dir == ShiftDirection::Right {
            self.push_ins_segment();
        }
    }

    fn add_match(&mut self, len: u32) {
        self.end_indel();
        self.match_block_size += len;
    }

    fn add_other(&mut self, cigar_seg: Option<&Cigar>) {
        self.end_indel();
        if self.match_block_size > 0 {
            self.shift_cigar.push(Cigar::Match(self.match_block_size));
            self.match_block_size = 0;
        }
        if let Some(x) = cigar_seg {
            self.shift_cigar.push(*x);
        }
    }
}
