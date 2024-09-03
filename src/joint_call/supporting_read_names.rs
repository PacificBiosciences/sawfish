use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use flate2::{write::GzEncoder, Compression};
use log::info;
use unwrap::unwrap;

use crate::breakpoint::{get_breakpoint_vcf_sv_type, VcfSVType};
use crate::refine_sv::get_rsv_id_label;
use crate::sv_group::SVGroup;
use crate::sv_id::get_bnd_sv_id_label;

/// Key on sample name
type VariantSupportingReadNames = BTreeMap<String, Vec<String>>;

/// Key on the variant id field
type SupportingReadNames = BTreeMap<String, VariantSupportingReadNames>;

fn get_supporting_read_names_from_sv_groups(
    sample_names: &[&str],
    sv_groups: &[SVGroup],
) -> SupportingReadNames {
    let sample_count = sample_names.len();

    let mut supporting_read_names = SupportingReadNames::default();

    for sv_group in sv_groups {
        for refined_sv in sv_group.refined_svs.iter().filter(|x| !x.filter_sv()) {
            assert_eq!(refined_sv.score.samples.len(), sample_count);

            // Simplified Refined SV type interpretation that skips inversions:
            if refined_sv.ext.is_inversion {
                continue;
            }

            // Get the VCF label for this SV:
            let sv_type = get_breakpoint_vcf_sv_type(&refined_sv.bp);
            let is_breakpoint = (sv_type == VcfSVType::Breakpoint)
                || refined_sv.ext.force_breakpoint_representation;

            let labels = if is_breakpoint {
                vec![
                    get_bnd_sv_id_label(&refined_sv.id, true),
                    get_bnd_sv_id_label(&refined_sv.id, false),
                ]
            } else {
                vec![get_rsv_id_label(refined_sv)]
            };

            let mut variant_read_names = VariantSupportingReadNames::default();
            for (sample_index, sample_score_info) in refined_sv.score.samples.iter().enumerate() {
                variant_read_names.insert(
                    sample_names[sample_index].to_string(),
                    sample_score_info.supporting_read_names.clone(),
                );
            }

            for variant_label in labels.iter() {
                supporting_read_names.insert(variant_label.clone(), variant_read_names.clone());
            }
        }
    }
    supporting_read_names
}

fn write_supporting_read_names_file(output_dir: &Path, supporting_read_names: SupportingReadNames) {
    let filename = output_dir.join("supporting_reads.json.gz");

    info!(
        "Writing supporting read names to file: '{}'",
        filename.display()
    );

    let fp = unwrap::unwrap!(
        std::fs::File::create(&filename),
        "Unable to create supporting read names json file: '{}'",
        filename.display()
    );

    let srn_string = unwrap!(
        serde_json::to_string_pretty(&supporting_read_names),
        "Failed to serialize supporting read names"
    );

    GzEncoder::new(fp, Compression::default())
        .write_all(srn_string.as_bytes())
        .unwrap();
}

/// Write supporting read list for each variant
pub fn write_supporting_read_names(
    output_dir: &Path,
    sample_names: &[&str],
    sv_groups: &[SVGroup],
) {
    // First step is the build the supportingReadNames structure:
    let supporting_read_names = get_supporting_read_names_from_sv_groups(sample_names, sv_groups);

    write_supporting_read_names_file(output_dir, supporting_read_names);
}
