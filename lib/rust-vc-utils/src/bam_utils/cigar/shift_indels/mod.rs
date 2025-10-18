mod cigar_indel_shifter;
mod left_shift_indels;
mod right_shift_indels;

pub use left_shift_indels::*;
pub use right_shift_indels::*;

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::Cigar;

    use super::*;
    use Cigar::*;

    #[test]
    fn test_shift_alignment_match() {
        let ref_pos = 2;
        let cigar = vec![Match(6)];
        let ref_seq = b"XXABCCDEXX";
        let read_seq = b"ABCCDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);

        let (shift_ref_pos, shift_cigar) = right_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);
    }

    #[test]
    fn test_shift_alignment_softclip() {
        let ref_pos = 4;
        let cigar = vec![SoftClip(2), Match(2), SoftClip(2)];
        let ref_seq = b"XXABCCDEXX";
        let read_seq = b"ABCCDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);

        let (shift_ref_pos, shift_cigar) = right_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);
    }

    #[test]
    fn test_shift_alignment_ins() {
        let ref_pos = 2;
        let cigar = vec![Match(3), Ins(1), Match(2)];
        let ref_seq = b"XXABCDEXX";
        let read_seq = b"ABCCDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(2), Ins(1), Match(3)]);

        let (shift_ref_pos, shift_cigar) =
            right_shift_indels(shift_ref_pos, &shift_cigar, ref_seq, read_seq);

        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);
    }

    #[test]
    fn test_shift_alignment_ins_to_edge() {
        let ref_pos = 4;
        let cigar = vec![Match(1), Ins(1), Match(2)];
        let ref_seq = b"XXABCDEXX";
        let read_seq = b"CCDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![SoftClip(1), Match(3)]);

        let cigar = vec![Match(2), Ins(1), Match(1)];
        let read_seq = b"CDEE";

        let (shift_ref_pos, shift_cigar) = right_shift_indels(ref_pos, &cigar, ref_seq, read_seq);

        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(3), SoftClip(1)]);
    }

    #[test]
    fn test_shift_alignment_del() {
        let ref_pos = 2;
        let cigar = vec![Match(3), Del(1), Match(2)];
        let ref_seq = b"XXABCCDEXX";
        let read_seq = b"ABCDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(2), Del(1), Match(3)]);

        let (shift_ref_pos, shift_cigar) =
            right_shift_indels(shift_ref_pos, &shift_cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, shift_cigar);
    }

    #[test]
    fn test_shift_alignment_del_on_interrupted_hpol() {
        let ref_pos = 2;
        let cigar = vec![Match(3), Del(3), Match(2)];
        let ref_seq = b"XXABBCBBBAXX";
        let read_seq = b"ABBBA";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(1), Del(3), Match(4)]);

        let (shift_ref_pos, shift_cigar) =
            right_shift_indels(shift_ref_pos, &shift_cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar);
    }

    #[test]
    fn test_shift_alignment_del_to_edge() {
        let ref_pos = 4;
        let cigar = vec![Match(1), Del(1), Match(2)];
        let ref_seq = b"XXABCCDEXX";
        let read_seq = b"CDE";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, 5);
        assert_eq!(shift_cigar, vec![Match(3)]);

        let cigar = vec![Match(2), Del(1), Match(1)];
        let ref_seq = b"XXABCDEEXX";

        let (shift_ref_pos, shift_cigar) = right_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(3)]);
    }

    #[test]
    fn test_shift_alignment_multi_indel() {
        let ref_pos = 2;
        let cigar = vec![Match(3), Ins(1), Match(2), Del(1), Match(1)];
        let ref_seq = b"XXABCDEEFXX";
        let read_seq = b"ABCCDEF";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(
            shift_cigar,
            vec![Match(2), Ins(1), Match(2), Del(1), Match(2)]
        );

        let (shift_ref_pos, shift_cigar) =
            right_shift_indels(shift_ref_pos, &shift_cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, cigar,);
    }

    #[test]
    fn test_shift_alignment_indel_cluster() {
        let ref_pos = 2;
        let cigar = vec![Match(4), Del(2), Ins(2), Match(1)];
        let ref_seq = b"XXABBBABFXX";
        let read_seq = b"ABBBBBF";

        let (shift_ref_pos, shift_cigar) = left_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(3), Ins(2), Del(2), Match(2)]);

        let cigar = vec![Match(3), Del(2), Ins(2), Match(2)];
        let (shift_ref_pos, shift_cigar) = right_shift_indels(ref_pos, &cigar, ref_seq, read_seq);
        assert_eq!(shift_ref_pos, ref_pos);
        assert_eq!(shift_cigar, vec![Match(4), Ins(2), Del(2), Match(1)]);
    }
}
