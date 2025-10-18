use rust_htslib::bcf::header::HeaderView;
use rust_vc_utils::ChromList;

/// Get chrom names from bcf header and provide a map to switch the chrom indexes
/// back to those in chrom_list
///
pub fn bcf_chrom_index_map(chrom_list: &ChromList, header: &HeaderView) -> Vec<usize> {
    let mut rid_to_chrom_index = Vec::new();
    let chrom_count = header.contig_count();
    for rid in 0..chrom_count {
        let chrom_bytes = header.rid2name(rid).unwrap();
        let chrom_name = std::str::from_utf8(chrom_bytes).unwrap().to_string();
        let chrom_index = *chrom_list.label_to_index.get(&chrom_name).unwrap();
        rid_to_chrom_index.push(chrom_index);
    }
    rid_to_chrom_index
}
