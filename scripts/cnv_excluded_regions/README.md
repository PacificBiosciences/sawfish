# CNV excluded region generation

This directory contains scripts and instructions for generating regions where CNV calling should be excluded in sawfish.
For users simply interested in using regions that have already been pre=computed, see the user guide [cnv excluded
regions section](../../docs/user_guide.md#cnv-excluded-regions).

The information in this directory provides both documentation for how sawfish's pre-computed excluded region tracks were
generated, as well as a template for users interested in extending these methods to new references sequences or to
different exclusion stringency levels, etc.

## Excluded region overview

The sawfish pre-computed cnv excluded regions are currently generated from two information sources. The first is
annotation tracks for the reference genome for regions such as centromeres and assembly gaps. The second is regions
where CNVs are found to be common across a diverse sample cohort, indicating an enrichment for false positives.

Ideally the excluded regions used in an analysis should include data from both sources, but results from a cohort
aligned to the same reference genome are not always available.

In procedure outlined below, we describe how the primary exclusion track for GRCh38 was generated and combined from both
annotation tracks and common CNV calls in a cohort.

## Conda environment setup

Required tools for all excluded region generation scripts in this directory are provided in the accompanying conda
environment specification. This conda environment can be setup as follows:

```
conda env create --file excluded_region_tools.yaml
```

## Annotation

The scripts used to produce sawfish's pre-computed 'annotation_only' exclusion regions are provided in the [annotation
directory](annotation). The following can be used to produce the excluded region tracks hg38:

```
conda activate excluded_region_tools

annotation/get_excluded_regions_from_ref_annotation.bash
```

This script includes a simple `ref` value which can be modified to produce the excluded region track for hg19 as well,
and provides a template which could be adapted to other genomes supported on the UCSC browser or similar annotation
sources.

The script `convert_hg19_regions_to_hs37d5.bash` can be used to further convert the hg19 excluded regions to hs37d5.

## Common CNV

Generating the common CNV excluded regions track requires sawfish results from a diverse sample cohort. It has been
found to be more effective to run sawfish is run on these samples without any `--cnv-excluded-regions` option, so that
all difficult CNV regions are attempted for the purpose of finding the common CNV regions. Note that only larger scale
events are being used for this operation so the `--disable-small-indels` flag can be used to accelerate the sawfish runs
used for this purpose.

All sawfish cohort VCFs can be input into the following script found in the `common_cnv/` subdirectory as follows:

```
conda activate excluded_region_tools

./common_cnv/convert_sawfish_vcfs_to_common_cnv_exclusion_bed.bash ${COHORT_SAMPLE_PATH}/*/sawfish_joint-call_output/genotyped.sv.vcf.gz > common_cnv.bed
```

The cohort samples can be distributed in any number of VCFs: one sample per VCFs, all samples in one VCF or any other
combination so long as each sample is input only once.

The resulting `common_cnv.bed` file written to stdout by this script identifies regions where 50% or more of cohort
samples contained the same type of CNV (deletion or duplication). These regions can be combined with annotation-based
excluded regions as described in the final step below.

For the pre-computed GRCh38 excluded regions track provided with sawfish, the common CNV mask was generated from the [47
HPRC year1 cohort samples](https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt).

## Combining results to create a full excluded region set

The final integration step combines the annotation and common CNV excluded regions:

1. Combine annotation and common CNV regions
2. Fill in any gaps of less than 20kb between these regions to make a fully expanded exclusion track
3. Combine the annotation and common CNV regions with the expanded exclusion track, to retain labels which explain each
   excluded region.

The script segment below describes this operation, assuming the annotation based excluded regions are from the
pre-computed track provided with sawfish, and the common CNV excluded regions come from the step described above.

```
conda activate excluded_region_tools

gzip -dc ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_only.hg38.bed.gz >|\
cat - common_cnv.bed |\
bedtools sort -i - >|\
combined.bed

bedtools merge -d 20000 -i combined.bed |\
bedtools subtract -a - -b combined.bed |\
sed "s/$/\texcluded_gap_fill/" |
cat - combined.bed |\
bedtools sort -i - | bgzip -c >|\
annoation_and_common_cnv.hg38.bed.gz

tabix -p bed annoation_and_common_cnv.hg38.bed.gz
```
