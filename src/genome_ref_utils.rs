use rust_vc_utils::{ChromList, GenomeRef, GenomeSegment};

/// Get a slice of the reference sequence corresponding to genome_segment
///
pub fn get_ref_segment_seq<'a>(
    chrom_list: &ChromList,
    genome_ref: &'a GenomeRef,
    genome_segment: &GenomeSegment,
) -> &'a [u8] {
    let chrom_info = &chrom_list.data[genome_segment.chrom_index];
    let segment_chrom_name = chrom_info.label.as_str();
    let range = &genome_segment.range;
    assert!(range.start >= 0);
    assert!(range.end > range.start);

    let chrom_ref = genome_ref.chroms.get(segment_chrom_name).unwrap();

    assert!(range.end <= chrom_ref.len() as i64);

    &chrom_ref[range.start as usize..range.end as usize]
}
