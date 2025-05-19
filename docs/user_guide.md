# Sawfish User Guide

## Table of Contents

* [Overview](#overview)
* [Accuracy evaluation](#accuracy-evaluation)
* [Getting started](#getting-started)
* [Outputs](#outputs)
* [Other usage details](#other-usage-details)

## Overview

Sawfish is a joint structural variant (SV) and copy number variant (CNV) caller for mapped HiFi sequencing reads. It
discovers germline structural variants from local sequence assembly and jointly genotypes these variants across multiple
samples. Sawfish additionally applies copy number segmentation on each sample's sequencing coverage levels,
synchronizing structural variant breakpoints with copy number change boundaries in the process to improve classification
of breakpoint-based structural variant calls, in addition to calling copy number variants.

Key features:
- Combined assessment of all large variants in each sample.
  - Sawfish provides a unified view of both SVs and CNVs in each sample, with redundant calls merged into single
    variants describing both breakpoint and copy number detail.
- High SV discovery and genotyping accuracy
  - All breakpoint-based structural variants are modeled and genotyped as local haplotypes, yielding substantial
    accuracy gains on modern SV truth sets such as the GIAB HG002 T2T SVs. 
- High resolution
  - All breakpoint-based structural variants are assembled to basepair resolution and reported with breakpoint homology
    and insertion details.
- Integrated copy number segmentation
  - Integrated copy number segmentation with GC-bias correction is used to: (1) call CNVs independent of any breakpoint
    support, and (2) improve the classification of large structural variant deletion and duplication calls, any such
    calls lacking consistent depth support are reclassified as breakends.
- Simple multi-threaded workflow
  - A single command-line is used for each of the discover and joint-call steps

Breakpoint-based SVs are reported as deletions, insertions, duplications and inversions when supported by the
corresponding breakpoint and depth pattern, otherwise the breakpoint itself is reported. Copy number variants are
reported as deletions and duplications. The minimum variant size is 35 bases (configurable). A maximum size is only
applied to inversions (100kb).

## Accuracy evaluation

Recommended methods for breakpoint-based SV accuracy assessment and benchmarking results are described in the [sawfish
App Note in *Bioinformatics*](https://doi.org/10.1093/bioinformatics/btaf136). A step-by-step overview of these
benchmarking methods is provided on the [accuracy assessment page](accuracy.md).

## Getting started

### Installation

Sawfish binaries are available for 64-bit Linux platforms. These can be installed either directly from the GitHub
release tarball, or via conda as described below.

#### Install from GitHub

To install sawfish from github, download the latest release tarball compiled for 64-bit Linux on the [github release
channel](https://github.com/PacificBiosciences/sawfish/releases/latest), then unpack the tar file. Using v0.12.1 as an
example, the tar file can be downloaded and unpacked as follows:

    wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.1/sawfish-v0.12.1-x86_64-unknown-linux-gnu.tar.gz
    tar -xzf sawfish-v0.12.1-x86_64-unknown-linux-gnu.tar.gz

The sawfish binary is found in the `bin/` directory of the unpacked file distribution. This can be run with the help
option to test the binary and review latest usage details:

    sawfish-v0.12.1-x86_64-unknown-linux-gnu/bin/sawfish --help

#### Install from conda

For [conda](https://github.com/conda/conda) users, installing sawfish on conda may be a more convenient option. Sawfish
is available for conda on Linux from the `bioconda` channel. A new conda environment with the latest sawfish release can
be created as follows:

    conda create -n sawfish -c bioconda sawfish

### Analysis steps

Sawfish analyzes samples in 2 steps:

1. `discover` - The discover step needs to be run once on each sample prior to entering that sample in the the
   joint-call step. The discover step completes several tasks:
   - Identifies candidate structural variant (SV) regions
    - Assembles each candidate SV into a local SV haplotype.
    - Builds a binned track of sequencing coverage over the whole genome
    - Runs initial copy number segmentation on the depth track, and iterates the segmentation to estimate and adjust for
      sample-specific GC-bias levels.
2. `joint-call` - The joint call step takes the output of the sawfish 'discover' step for one-to-many samples and
   provides jointly genotyped SV calls over the sample set. Joint calling includes the following operations:
    - Merge duplicate SV haplotypes
    - Associate deduplicated SV haplotypes with samples
    - Evaluate SV read support in each sample
    - Genotype quality assessment
    - Identification of breakends used to reward additional copy-number change boundary points.
    - Re-segmentation of depth in all samples.
    - Merging of SV and CNV signatures corresponding to the same variant.
    - Merged SV and CNV results written out to VCF.

### Calling SVs and CNVs on one sample

To call SVs and CNVs in one sample, run `discover` on the mapped sample bam, and then run `joint-call` on the output
directory of the discover step.

The following example shows how this is done for a mapped sample bam named `HG002.GRCh38.bam`, using 16 threads for both
the `discover` and the `joint-call` steps.

    sawfish discover \
      --threads 16  \
      --ref GRCh38.fa \
      --bam HG002.GRCh38.bam \
      --expected-cn expected_cn.hg38.XY.bed \
      --cnv-excluded-regions cnv.excluded_regions_inc_sawfish_common_cnv.hg38.bed.gz \
      --output-dir HG002_discover_dir
    
    sawfish joint-call \
      --threads 16 \
      --sample HG002_discover_dir \
      --output-dir HG002_joint_call_dir

The primary output of the `joint-call` step can be found in `HG002_joint_call_dir/genotyped.sv.vcf.gz`. See the
[outputs](#outputs) section below for discussion of this and all other output files.

Note that the reference fasta and sample bam path specified in the `discover` step are stored and reused in the
subsequent `joint-call` step.

### Joint calling SVs and CNVs across a set of samples

To call SVs and CNVs on a set of samples, run `discover` separately on each mapped sample bam, and then run `joint-call`
on all `discover` step output directories.

The following example shows how this is done for mapped sequences from the HG002 trio, given the following bam files:
`HG004.GRCh38.bam`, `HG003.GRCh38.bam`, `HG002.GRCh38.bam`.

As a first step, `discover` needs to be run on all 3 samples. In the example below 16 threads are used to process each
sample. Note that these 3 commands are independent and could be run in parallel.

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG004.GRCh38.bam \
      --expected-cn expected_cn.hg38.XY.bed \
      --cnv-excluded-regions cnv.excluded_regions_inc_sawfish_common_cnv.hg38.bed.gz \
      --output-dir HG004_discover_dir

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG003.GRCh38.bam \
      --expected-cn expected_cn.hg38.XY.bed \
      --cnv-excluded-regions cnv.excluded_regions_inc_sawfish_common_cnv.hg38.bed.gz \
      --output-dir HG003_discover_dir

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG002.GRCh38.bam \
      --expected-cn expected_cn.hg38.XY.bed \
      --cnv-excluded-regions cnv.excluded_regions_inc_sawfish_common_cnv.hg38.bed.gz \
      --output-dir HG002_discover_dir

After all discover steps have completed, joint calling can be run over all 3 samples using the following command:

    sawfish joint-call \
      --threads 16 \
      --sample HG004_discover_dir \
      --sample HG003_discover_dir \
      --sample HG002_discover_dir \
      --output-dir HG002_trio_joint_call_dir

The primary output of the `joint-call` step can be found in `HG002_trio_joint_call_dir/genotyped.sv.vcf.gz`. See the
[outputs](#outputs) section below for discussion of this and all other output files.

Just as in the single-sample case, note that the reference fasta and all 3 sample bam paths specified in the `discover`
steps are stored and reused in the subsequent `joint-call` step.

## Inputs

### Input read alignments

Read alignments for the query sample must be supplied in BAM or CRAM format as an argument in the discover step.

Sawfish has been tested with sequencing reads mapped by [pbmm2](https://github.com/PacificBiosciences/pbmm2). In general
it is designed to work on supplementary alignments without hard-clipping. If this requirement is fulfilled it may work
with other mappers, but no others are tested or supported.

When joint-calling over multiple samples, all input alignment files must have been mapped to the same reference genome.

### Reference fasta

A genome reference sequence file in fasta format is required as input for every run at the discovery step. Every
chromosome name in the input read alignment file must be be present in the reference sequence file. Their is no
reciprocal requirement, the reference fasta may contain chromosome names not present in the input bam file.

### Expected copy number

An expected copy number file can be provided for each sample during the discover step, to set expected copy number per
region of the genome. Any regions not specified will have a default expected copy number of 2. If no file is specified
the default expected copy number of 2 will apply to the whole genome.

The expected copy number is important in determining which copy number segments will be output as CNV deletions or
duplications, and is especially useful to indicate the expected copy number for mammalian sex chromosomes.

During the joint-call step, the expected copy number inputs from the discover phase are retained per sample, so for
instance in a human pedigree analysis the expected copy number on the chrX nPAR region can vary by sample sex chromosome
complement.

The expected copy number file must be in BED format, with the first 3 columns used to specify regions following standard
BED format, and with expected copy number provided in column 5. Column 4 is ignored and can be used as a region label.

Example expected copy number BED files are provided for some common human reference genomes and sex chromosome
complements in the sawfish [expected_cn](../data/expected_cn/) directory. These can be used directly or serve as
templates for other sample configurations.

Note that in earlier versions of sawfish, before CNV output was introduced, the expected copy number input would change
the ploidy used by the genotyper for breakpoint-based calls. This ploidy change is no longer made by default, and for
general purpose calling this change is no longer recommended. If the previous ploidy-change behavior is preferred, the
`--treat_single_copy_as_haploid` option can be provided in the joint-call step to cause any region with copy number 1 to
be treated as haploid (all other cases will continue to be treated as diploid).

### CNV excluded regions

Certain regions of each reference genome may present inherent difficulties to the prediction of meaningful CNV calls.
Such regions can be marked as excluded for the purpose of CNV calling by providing the regions in BED file format using
the discover step argument `--cnv-excluded-regions`. The way that excluded regions impact CNV calling is summarized in
the section further below, but note that these regions do not change the behavior of breakpoint-based SV calling.

Pre-computed CNV excluded regions are provided in the sawfish [cnv_excluded_regions](../data/cnv_excluded_regions/)
directory for some common human reference genomes. The recommended exclusion track for GRCh38 is
`cnv.excluded_regions_inc_sawfish_common_cnv.hg38.bed.gz`. This is the only reference genome for which complete excluded
regions are provided. The complete excluded regions are a combination of (1) regions annotated as assembly gaps,
centromeres, and alpha satellite sequences (2) regions where sawfish calls the same CNV type in at least 50% of the 47
HPRC year1 cohort samples
(https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt).

For other genomes, 'base'/minimal excluded region tracks are provided for the first excluded region type only (regions
annotated as assembly gaps, centromeres, and alpha satellite sequences). These files are named with the pattern
`cnv.excluded_regions.${genome_tag}.bed.gz`, where the genome_tag value may be hg38, hg19 or hs37d5. For reference
genomes other than GRCh38, these files provide some level of exclusion to reduce false positives, and should produce a
better result than not excluding any regions. For other reference genomes, it is recommended to develop a similarly
expanded exclusion region track including common sawfish CNV calls from a background cohort.

The scripts used to produce the 'base' exclusion regions are provided in the [cnv_excluded_regions script
directory](../scripts/cnv_excluded_regions). Here `get_cnv_exclusion_regions.bash` is used to produce the exclusion
tracks for hg19 and hg38, and provides a template which could be adapted to other genomes supported on the UCSC browser.
The script `convert_hg19_regions_to_hs37d5.bash` is used to convert the hg19 regions to hs37d5.

#### How excluded regions influence copy number calling

Excluded regions are designed to prevent any copy number segments from changing within the excluded region. While copy
number variants spanning excluded regions are penalized, there an allowance for a longer copy number segment to span
through a relatively small excluded region without interruption. This allows a megabase-scale copy number gain to be
represented as a continuous CNV over a reference assembly gap, for example.

A summary of excluded region behavior is as follows:
- All depth bins intersecting an excluded region are removed from the depth bins track. 
- All minor allele frequency evidence intersecting an excluded region are removed from the minor allele frequency track.
- Segmentation will treat any depth bins intersecting an excluded region as having a small bias in favor of a special
unknown copy-number state -- the probability of all other copy number states are equal, but lower than the unknown
state. This means that a copy number change can span through a short excluded region if there is sufficient evidence on
the left or right flank, but longer excluded regions should be segmented into an unknown state.

### Minor allele frequency

On the discover step command line for each sample, a small-variant VCF file for that sample can be provided as an
argument to `--maf`. The given VCF will be parsed to create a minor allele frequency track for the genome. This
information is written to an IGV visualization track for assessment and interpretation of the CNV output. This
information may also be used for improved segmentation and CNV calling in the future, although it is not used in the
current release.

Any small-variant caller VCF with an `AD` entry for the small-variant calls could work for this purpose, but the feature
is tested and best supported for the output from [DeepVariant](https://github.com/google/deepvariant).

## Outputs

### Joint call step

The primary output of the joint-calling step are the SV and CNV calls for all samples in VCF 4.4 format, written to
`${OUTPUT_DIR}/genotyped.sv.vcf.gz`. Details of the SV and representation in this file are provided below.

#### Quality scores

The primary quality metrics for each SV call are:
1. `QUAL` - This is the phred-scaled confidence that the given SV allele exists in the set of genotyped samples.
2. `GQ` - This value is provided once for each sample. It is the phred-scaled confidence that the given sample genotype
   is correct.

The primary quality metrics for each CNV call are:
1. `QUAL` - This is the phred-scaled confidence that the segment copy number is not the expected copy number in at least
   one sample.
2. `CNQ` - This value is provided once for each sample. It is the phred-scaled confidence that the given copy number
   (`CN`) value is correct for this segment.

All phred-scaled quality scores in the VCF output have a maximum value of 999.

#### Filters

The following filters may be applied to each VCF record:

- `MinQUAL` - The variant quality score (`QUAL`) is less than 10
- `MaxScoringDepth` - Read depth at an SV locus exceeds 1000x, so all scoring and genotyping steps were disabled.
- `InvBreakpoint` - This breakpoint is represented as part of a separate VCF inversion record (the inversion record
  shares the same `EVENT` ID)
- `ConflictingBreakpointGT` - Genotypes of breakpoints in a multi-breakpoint event conflict in the majority of cases
  (This filter is only relevant to inversions at present)

#### SV and CNV types

Notes on formatting and representation of SVs and CNVs are listed below for each major type.

##### Deletions

Deletion records include both breakpoint-based SV deletions and copy-number loss CNVs, these all have an INFO entry of
`SVTYPE=DEL`. Per the VCF 4.4 spec, the `SVCLAIM` value is used to distinguish the type of support for each call, where
"D" indicates depth support from copy-number segmentation, "J" indicates breakpoint-based support from read assembly,
and "DJ" indicates that the call is supported by both.

All breakpoint-based deletions of 100kb or smaller are represented by directly writing the deleted sequence in the VCF
`REF` field and any breakpoint insertion sequence in `ALT`. Deletions larger than 100kb, or lacking breakpiont support
are written as symbolic alleles using the ALT value of `<DEL>`.

All candidate breakpoint-based deletions at least 50kb in length without depth support from copy-number segmentation
will be reported in the VCF output as a pair of breakend (`BND`) records instead.

###### Insertions

Any indel-like SVs where the length of sequence inserted at the breakpoint exceeds the length of deleted sequence will
be formatted as an insertion in the VCF output if it is possible to fully assemble the inserted sequence, and will be
formatted as a duplication otherwise. If represented as an insertion the full inserted sequence assembly will be written
to the VCF `ALT` field.

##### Duplications

Duplication records include both breakpoint-based SV duplications and copy-number gain CNVs, these all have an INFO
entry of `SVTYPE=DUP`. Per the VCF 4.4 spec, the `SVCLAIM` value is used to distinguish the type of support for each
call, where "D" indicates depth support from copy-number segmentation, "J" indicates breakpoint-based support from read
assembly, and "DJ" indicates that the call is supported by both.

Very large insertions with long breakpoint homology will be represented as duplications in the VCF output only if they
cannot be output as insertions. These will be written to the VCF output using the symbolic ALT value of `<DUP:TANDEM>`.
Copy-number gain CNV records without breakpoint-based support use a symbolic ALT value of `<DUP>`.

 All candidate breakpoint-based duplications at least 50kb in length without depth support from copy-number segmentation
 will be reported in the VCF output as a pair of breakend (`BND`) records instead.

##### Breakpoints

All SV breakpoints which can't be modeled as one of the simple SV types above will be output as a pair of breakend
(`BND`) records.

##### Inversions

Sawfish will currently identify one type of multi-breakpoint complex SV signature, corresponding to that of a simple (or
balanced) inversion. Inversions are identified when two intra-chromosomal inverted breakpoints of opposite orientation
have overlapping spans with at least a 60% reciprocal overlap. The longer span must not be greater than 100kb.

When an inversion is found, a VCF record will be output using the `<INV>` symbolic allele summarizing the inversion in
as much detail as possible. It is not possible to retain the details of all 4 breakends in this format such as all
breakend positions and breakpoint insertion sequences. For this reason the corresponding breakend records are retained
in the VCF output but marked as filtered, such that full breakend details remain available in the output. The inversion
record and the filtered breakend records are given a shared VCF `EVENT` label so that their relationship can be
identified.

#### Overlapping SV formatting

All sawfish SVs are output so that only one allele is described in each VCF record, even if an overlapping SV allele is
output at the same locus. The internal SV calling model accounts for up to 2 overlapping alleles per sample during
genotyping and quality scoring. Reads which support a 2nd alternate allele at any given locus will be counted as
supporting the reference in output fields such as allele depth (`AD`). This protocol matches standard SV caller
formatting conventions. Users interested in a more detailed output format, such as representing overlapping read support
on the VCF `<*>` allele can request this for prioritization.

#### Phasing

Sawfish adds short-range phasing information to clarify the relationship of heterozygous SVs called from the same or
overlapping SV haplotypes. This does not have the range of general read-backed phasing and will only result in phased
genotype output for smaller insertions and deletions. Each local cluster of phased genotypes corresponds to a phase set
as annotated using the VCF `PS` tag. The phase set ID is the `POS` value of the first SV called from the SV haplotype
cluster.

#### Per-sample output

In addition to the final merged SV and CNV VCF file output, the `joint-call` step also provides certain output files for
each sample. These files are primarily associated with copy number segmentation or visualizing/interpreting the CNV
caller output. The per-sample output is written to the directory `${OUTPUT_DIR}/samples`. Within this directory there is
one subdirectory per sample following the pattern `sample{sample_index}_{sample_name}`, where sample index reflects the
order that samples are listed on command-line for the `joint-call` step.

##### Copy-number segmentation track

The final copy-number segmentation result for the given sample is provided in `copynum.bedgraph`, where the copy number
value is listed in column 4. Not that any region segmented into the 'excluded' state will be represented as an uncovered
gap in the region coverage.

##### Depth track

The depth bin values used as input to the segmentation process for each sample are provided in bigwig format in the file
`depth.bw`. This track can be useful to visualize and interpret the CNV output.

##### Minor allele frequency track

When a minor allele frequency input file is provided for the sample, the corresponding minor allele frequency track will
be output in bigwig format in the file `maf.bw`.This track can be useful to visualize and interpret the CNV output.

#### Optional variant read support output

To help show which reads support each SV allele, the optional `--report-supporting-reads` argument can be added to the
joint-call command line. When this is used a compressed json output file is provided in
`${OUTPUT_DIR}/supporting_reads.json.gz`.

In this json output file, the top-level objects are variant IDs matching those provided in the ID field of the VCF
output. Nested under each variant ID are sample IDs. For each sample ID associated with a variant, the array of
supporting read QNAME values are provided. A simplified example output is shown below for two variants:

```
{
  "sawfish:0:1041:0:0": {
    "HG002": [
      "m84005_220919_232112_s2/22021538/ccs",
      "m84005_220919_232112_s2/108659098/ccs",
      "m84005_220919_232112_s2/166989308/ccs"
    ]
  },
  "sawfish:0:1051:0:0": {
    "HG002": [
      "m84005_220919_232112_s2/130223022/ccs",
      "m84005_220919_232112_s2/9113818/ccs",
      "m84005_220919_232112_s2/84214835/ccs",
      "m84005_220919_232112_s2/116654499/ccs"
    ]
  }
}
```

Note that the number of read QNAME entries should often match the supporting AD count for the alternate allele from the
same variant/sample entry in the VCF, but this is not always an exact match. Also to keep a consistent relationship
between supporting reads and variants, no output is provided for VCF records with the inversion (`<INV>`) allele type,
but the supporting reads for the breakends comprising each inversion are provided.

### Discover step

The discover step produces a number of output files in the discover output directory used by sawfish during the
subsequent joint calling step. Although these are not documented or intended for end users, some of the more important
files are noted below:

- `assembly.regions.bed` - Describes each region of the genome targeted for assembly.
- `candidate.sv.bcf` - These are the candidate SVs expressed in a simplified format for each sample. These are used as
  input for joint genotyping together with the aligned candidate contigs.
- `discover.settings.json` - Various parameters from the discover step (either user input or default) are recorded in
  this file. Some of the paths to files like the sample bam and reference fasta will be reused in the joint call step.

### Debug outputs from either step

In either run step, the following files are produced to help debug problematic runs or SV calls:

- `${OUTPUT_DIR}/sawfish.log` - High level logging output
- `${OUTPUT_DIR}/run_stats.json` - Run statistics and component timings
- `${OUTPUT_DIR}/contig.alignment.bam` - Contigs for assembled SV haplotypes aligned back to the reference genome. For
  the joint-call output this file shows the contigs used for the final VCF output, after all haplotype merging across
  samples has been completed.

#### Contig alignment format

SV haplotype contig alignments are output to `${OUTPUT_DIR}/contig.alignment.bam` in either the discover or joint-call
steps, and can be useful for reviewing SV calls. For instance, this file can be viewed in alignment browsers such as
IGV.

Aligned contigs are provided for all single-breakpoint SV calls. To find the contig for a given SV, locate the SV's VCF
ID field, such as `sawfish:0:2803:1:2`, and take the prefix from this ID that includes the first three digits, in this
case `sawfish:0:2803:1`. This is the `QNAME` value of the corresponding SV haplotype alignment(s) in the contig
alignment BAM file.

Contigs are not available for CNVs or multi-breakpoint events such as inversions. For the latter case, contigs are
available for each individual breakpoint comprising the event.

In addition to standard sequence and alignment information, each contig BAM record includes a custom aux field called
`sf` which provides a list of key/value properties associated with the contig, for instance:

    sf:Z:n_reads:15;hq_range:1500-2311;

The properties are:

- `n_reads` - The number or reads used to assemble the contig
- `hq_range` - The high-quality assembled region of the contig, prior to appending any flanking read sequence anchors.

An example contig alignment bam record is:

```
sawfish:0:92:0   0       chr1    1649635 20      211=1X108=2I38=1X201=1X86=20I89=1I12=1X281=6D17=1X65=1D22=1X359=1X129=1X154=53D99=1X210=1X140=1X134=1X73=1X55=1X400=1X186=1X133=1X233=1X44=1D74=1635D9= *       0       0       TCCCTAATGAGAAATAAAGTGTCATGCAAAGAAACCTCACTTCAAAAATTTCACATGAAGCCGGGCACGGAGGCTTATGCCTGTAATCCTAGCACTTTGGGAGGCTGAGGCGGGCGGATCACCTGAGGTCAGGAGTTCAAGGCCATCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCTGGGCGTGGTGGCAGACACCTGTAATCCCAGCTACTAAGCGAGGCTGAGGCAGAAGAATTGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCACGCCACTGCACTACAGCCTGGGCAAAAAAAAAAAAAAAAAACCCACGTGAAACTGAAATTAAGGCCGGGCGCGGTGGCTCACGCCTGTAATTCCAGCACTCTGGGAGGCCGAGGTGGGCGGATCACAAGGTCAGATCGGGACCATCCTGGCTAACACGGTGAAACCCCATCTCTACTAAAAATACAAAAAATTAGCTGGGTGTGGTGGCGGGCACCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCGGGAGAATAGCGTGAACCCGGGAGATGGAATTTGCAGTGAGCTGAGATTGCGCCACTGTACTCCAGCCTGGGTGACAAGCAAGACTCCGTCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAGAAATTAAATCAAGAACAGTAAATATTTAAATAAATATTTAAATAATGATGTTAACGTTAAGTAATCCTAATTTTTCTTTTTTTTCTTTTTTTTTTTTTTGAGATGGAGTCTTGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGCAATCTCGGCTCACTGCAAGCTCCGCCTCCCGTGTTCACACCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGTACTACAGGCGCCTGCCACCACGCCCGGCTAATTTTTTGTCTTTTTAGTAGAGACGGGGTTCCACCATGTTAGCCAGGATGGTCTTGATGTTCTGACCTCGTGATCTGCCGGCCTCGGCCTCCCAAAGTGTGGGGATTACAGTTGTGAGCCACCGCGCCCGGCCTTTTTTTTTTTTTTTTTAAAAGAGACAGGGTCTCGCTATATTGGCCAGGCTGGTCTTGAACTCCCGGACTCAATTGATCCTCCAAGTGCTGGGATTACAGGCCTGAGCCACTGCACCCAGCCGAATAATCATGATTTTATGTTAAATAAAAAACTTTGAAAATAAAAAACTATCTGCAGTAAGCGTCTAATTATGAAGAAAGAAGAAAAAAGAAAAACAATTCTGCTATCACAGAAGAGCAGAATTGTAATATTCGTCTTTTAAAACTTTTCCATACTGAATAAACTATAATTATCAGTTTTATAATACAAAAATCACTCTTCACAAAGACTACAGAACAAAGCTTTGCTATCAGTGGGCTTCTCCACTGTGCAATGAAGCCACATTAATTAATCAAGTGTATTTATAATCATGACATTTCAATCGGGCTCCAGGTCCAATTTTCCTAACACCCGTAAGAACCTCTTGATGTTGGTACGAGATCAAACTGCTCAAGCCAAACCCATTCTTTGGACTTGAGCAAATACCCATTTTGGGGTGTGTTTTTCTCCTATACTTGTTGAATTCAGGTCATTTTAAATGTAAACAAACTGCTCCCAAACAATATAATGGGGGAGAGAAAACCCCAGAGGAAAAATGGACTAGCCATTCTGAATGGTCTGTGACACACGCACGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTCGGTCTGTGATACACGCATGCTTTCAGCTGGAGTTTGCTCTCTGTAGCCCCTCTGAATGGTCTGTGACACATGCACGCTTTCAGCTAGAGTACTCTCTCTATAGCCCTTCTGAATGGTCTGTGACACACGCATGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTCGGTCTGGGACACATGCATGCTTTTAGCTAGAGTTTGCTCTGTATAGCCCTTCTGAACGGTCTGTGACACACGCATGCTTTCAGCTGGAGTTTGCTCTCTATAGCCCCTCTGAATGGTCTGTGACACACGCATGCTTTCAGATAGAGTATTCTCTCTATAGCCCTTCTGAATGGTCTGTAACACACGCAAGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTTGGTCTGTGACACACGCATGCTTTTAGCTAGAGTTTGCTCTGTATAGCCCTTCTGAATGGTCTGTGACACATGCATGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTCAGTCTGTGACACACACATGCTTTTAGCTAGAGTTTGCTCTGTATAGCCCTTCTGAATGGTCTGTGACACACGCGTGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTTGGTCTGTGACACACGCATGCTTTTAGCTAGTTTGCTCTCATAGCCCTTCTGAACGGTCTGTGACACATGCATGCTTTCAGCTATTCTCTCTATAGCCATTGTGAATGGTCTGTGACACACGCACGCTTTCAGCTAGAGTTTGCTCTTTCTGGTTTTTGGTCTGTGACACACGCATGCTTTCAGCTAGAGTTTGCTCTCTCTGGTTTTCGGTCTGTGACGCACGCATGCTTTTAGCTAGAGTATTCTCTCTATAGCCATTCTGAACGGTCTGTGACACACGTATGCTTTCAGCTAGAGTTTGCTTTCTCTGGTTTTTCAGTGGTGCTCTGGGGAAGGCAGAAGAGTAGGAACAGGAAAGAAACCACACTTGAACATGATGTCAAAGAAAGTAAATGCTTCTGTACCCCCTTCTGCTGAATGGCTACGATGCCTACGTTTCTCTTTTCTCTTTTCATCTTTTCTGTGATGAGCTTTTTCTTTCCGAGACATTTGCTGGGGTGGTTTGATGGCCAAAGAATCATCTTCTTCTCCTCTGAAATAAAACACAACAGCACTGCGTCATGCTTGAGAAAGTGCAAAGAAGCATCAGGCTATTATAAGGTTTCTTCAACCCAGAAAAATGCATGATTCAGACAGGAACAAAGCTGAAACATCATTTAAAAAATTACATTAATTCTCCAACTTTAGGCATCTTTTTTTTCTTTTTTTCTTTTTTTTAGACAGTCTCGCTCTGTTGCCCGGGCTGTAGTGGCACGATCTCGGCTCACTGCAATCTCCACCCTCCGGGTTCATGCCATTCTCTTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCGCCACGCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTTACCATGTTAGCCAGGATGGTCTTGGTCTCCTGACCTCATGATCCGCCCACCTCGGTCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGCGCCCGGCCTGTATTTATTTTTTTGAGACGGAGTCTCGCTCTGTTGCCCAGGCTGGAATGCAGTGGTACGATCTCGGCTCATTGCAACCTCCCCTTCCAGTCCCAGGTTCAAGCAATTCTCCTGCCTCTGCCTCAGGAGTAGCTGGGATTACAGGCATGCGCCACCACACCCGGCTAATTTTTTATTTTTAGTAGAGACGGGGTTTCACCATATTGGTCAGGCTGGTCTCAAACTTGTGACATCATGATCCACCCACCTCG     *       sf:Z:n_reads:12;hq_range:1500-2103;
```

In this case, the VCF record for the corresponding SV call generated from the contig is:

```
chr1    1651421 sawfish:0:92:0:0        GTAGCCCCTCTGAACGGTCTGTGACACACGCATGCTTTCAGCTAGAGTACTCTA  G       504     PASS    SVTYPE=DEL;END=1651474;SVLEN=53;HOMLEN=13;HOMSEQ=TAGCCCCTCTGAA;SVCLAIM=J        GT:GQ:PL:AD     0/1:387:537,0,387:9,12
```

## Other usage details

### Analyzing CNV/Large variants only

Sawfish has an option to run selecting for large-scale breakpoint events and CNVs only by adding the
`--disable-small-indels` flag in the discover step of every sample. In this mode all CNVs and non-indel SVs will be
retained as well as large scale breakpoint-based deletions and duplications. Not that this mode is based on calling
mechanism rather than a strict size cutoff, but in general all content above ~10kb should still be present in this mode.
This is an experimental feature, meaning that the exact behavior of this mode may change in the future based on user
feedback, but in general when running with this option sawfish can be thought of as a breakpoint-enhanced CNV caller.

### Determinism

Sawfish should always produce the same output from a given command-line and input file set (allowing for expected
changes in timestamps, benchmark timers and similar metadata).

### Output directory / clobber

Each step of the pipeline accepts the argument `--output-dir` where all files from the step will be written. If not
specified the default of either `sawfish_discover_output` or `sawfish_joint-call_output` will be used. Sawfish will not
proceed if the output directory already exists, unless the `--clobber` argument is given as well.

### VCF ID field

The entries in the output VCF ID field (such as `sawfish:0:2803:1:2` and `sawfish:INV4`) are designed to guarantee a
unique identifier for each record in the VCF output. This identifier isn't meant to convey useful details about the call
and may be reformatted in future releases.

### Expected compute requirements

#### Runtime

For a typical ~30x HiFi sample analyzed on 16 threads, the `discover` step should complete in about 30-40 minutes and
the `joint-call` step should complete in about 5 minutes.

Running the `joint-call` step on 10 samples at 30-100x depth completes in about 1 hour on 64 threads.

In general, runtime response to thread count is expected to be nearly linear. The current joint calling scheme has been
designed with pedigree-scale analysis in mind. Sawfish joint calling has completed on 47 HPRC samples in testing, but
substantially larger cohorts would be difficult without further changes to the joint-calling design.

#### Memory

The `discover` step should typically require less than 8Gb/thread so long as at least several threads are selected. The
`joint-call` step should require substantially less memory but hasn't been tested at scale with less than 1Gb/thread.

### Discover step input file path storage

Sawfish accesses several files associated with each sample during joint-genotyping in the `joint-call` step. For
instance, this is done to test read support for each allele by accessing the sample alignment file.

To find these files for each sample, input file paths are stored from the `discover` step in a configuration file
written to the output directory here:

    ${OUTPUT_DIR}/discover.settings.json

These input file paths are normally canonicalized, so that relative paths can be reliably reused after any change to the
working directory. In some cases it may be more convenient to store relative file paths. To do so the `discover` step
option `--disable-path-canonicalization` can be used to store all input paths as-is. This may be useful if e.g., the
`discover` and `joint-call` steps are being run in different directory structures.

Note that for even more complex situations, the paths in the above discover settings json file can be manually edited
before running the `joint-call` step.
