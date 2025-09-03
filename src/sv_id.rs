/// Details used to generate a unique ID for every SV in the sample
#[derive(Clone, Eq, PartialEq, PartialOrd, Ord)]
pub struct SVUniqueIdData {
    /// Index of sample in this analysis
    ///
    /// This will always be zero in discovery, but allows the discovery candidates to be merged
    /// at the joint call step while retaining a globally unique id on VCF output
    pub sample_index: usize,

    /// Index of SV candidate cluster
    pub cluster_index: usize,

    /// Index of all assemblies from one cluster
    pub assembly_index: usize,

    /// Index of all SVs parsed out of one assembly alignment
    pub alignment_index: usize,
}

impl std::fmt::Debug for SVUniqueIdData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", get_sv_id_label(self))
    }
}

/// Get core component of ID label basd on indexes only
fn get_sv_id_core_label(id: &SVUniqueIdData) -> String {
    format!(
        "{}:{}:{}:{}",
        id.sample_index, id.cluster_index, id.assembly_index, id.alignment_index
    )
}

/// The default ID label for 'standard' SVs
pub fn get_sv_id_label(id: &SVUniqueIdData) -> String {
    let core_id = get_sv_id_core_label(id);
    format!("{}:{core_id}", env!("CARGO_PKG_NAME"),)
}

/// A variant of the default SV label used for inversions corresponding to the given id:
pub fn get_inv_id_label(id: &SVUniqueIdData) -> String {
    let core_id = get_sv_id_core_label(id);
    format!("{}:INV:{core_id}", env!("CARGO_PKG_NAME"),)
}

/// A variant of the default SV label used for breakends
pub fn get_bnd_id_label(id: &SVUniqueIdData, is_breakend1: bool) -> String {
    let bnd_index = if is_breakend1 { 0 } else { 1 };
    format!("{}:{}", get_sv_id_label(id), bnd_index)
}

pub fn get_sv_id_from_label(id: &str) -> (SVUniqueIdData, Option<bool>) {
    let parts = id.split(':').collect::<Vec<_>>();
    let plen = parts.len();
    assert!((5..=6).contains(&plen));
    let sample_index = parts[1].parse::<usize>().unwrap();
    let cluster_index = parts[2].parse::<usize>().unwrap();
    let assembly_index = parts[3].parse::<usize>().unwrap();
    let alignment_index = parts[4].parse::<usize>().unwrap();

    let sv_id = SVUniqueIdData {
        sample_index,
        cluster_index,
        assembly_index,
        alignment_index,
    };

    let is_breakend1 = if plen > 5 {
        let breakend_index = parts[5].parse::<usize>().unwrap();
        Some(breakend_index == 0)
    } else {
        None
    };

    (sv_id, is_breakend1)
}
