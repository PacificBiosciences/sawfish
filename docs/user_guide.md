# Sawfish User Guide

## Table of Contents

* [Overview](#overview)
* [Accuracy evaluation](#accuracy-evaluation)
* [Getting started](#getting-started)
* [Inputs](#inputs)
* [Outputs](#outputs)
* [Other usage details](#other-usage-details)

## Overview

Sawfish is a joint structural variant (SV) and copy number variant (CNV) caller for mapped HiFi sequencing reads. It
discovers germline structural variants from local sequence assembly and jointly genotypes these variants across multiple
samples. Sawfish additionally applies copy number segmentation on each sample's sequencing coverage levels,
synchronizing structural variant breakpoints with copy number change boundaries in the process to improve both
structural variant and CNV calling accuracy.

Key features:
- Combined assessment of all large variants in each sample.
  - Sawfish provides a unified view of SVs and CNVs in each sample, with each jointly-supported variant merged into one
    record detailing both breakpoint and copy number support.
- High SV discovery and genotyping accuracy
  - All breakpoint-based structural variants are modeled and genotyped as local haplotypes, yielding substantial
    accuracy gains on modern SV truth sets such as the GIAB HG002 T2T SVs.
- High resolution
  - All breakpoint-based structural variants are assembled to basepair resolution and reported with breakpoint homology
    and insertion details.
- Integrated copy number segmentation
  - Integrated copy number segmentation with GC-bias correction is used to: (1) independently call CNVs and (2) improve
    the classification of large SV deletion and duplication calls, any such calls lacking consistent depth support are
    reclassified as breakends.
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
channel](https://github.com/PacificBiosciences/sawfish/releases/latest), then unpack the tar file. Using v2.0.4 as an
example, the tar file can be downloaded and unpacked as follows:

    wget https://github.com/PacificBiosciences/sawfish/releases/download/v2.0.0/sawfish-v2.0.0-x86_64-unknown-linux-gnu.tar.gz
    tar -xzf sawfish-v2.0.4-x86_64-unknown-linux-gnu.tar.gz

The sawfish binary is found in the `bin/` directory of the unpacked file distribution. This can be run with the help
option to test the binary and review latest usage details:

    sawfish-v2.0.4-x86_64-unknown-linux-gnu/bin/sawfish --help

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
    - Runs initial 'draft' copy number segmentation on the depth track, and iterates the segmentation to estimate and
      adjust for sample-specific GC-bias levels.
2. `joint-call` - Given the 'discover' command results from each sample, the joint call step merges and genotypes SVs,
   calls CNVs and unifies SV/CNV calls from one or more samples. Joint calling includes the following operations:
    - Merge duplicate SV haplotypes
    - Associate deduplicated SV haplotypes with samples
    - Evaluate SV read support in each sample
    - Genotype quality assessment
    - Identification of breakends used to reward additional copy-number change boundary points.
    - Re-segmentation of depth in all samples.
    - Merging of SV and CNV signatures corresponding to the same variant.
    - Merged SV and CNV results written out to VCF.

### Quickstart Example 1: Calling SVs and CNVs on one sample

To call SVs and CNVs in one sample, run `discover` on the mapped sample bam, and then run `joint-call` on the output
directory of the discover step.

The following example shows how this is done for a mapped sample bam named `HG002.GRCh38.bam`, using 16 threads for both
the `discover` and the `joint-call` steps.

#### Single sample discover step

In this example the discover step is run on the given bam file input. The additional tracks for expected copy number and
regions excluded from CNV calling are optional but substantially improve the value of sawfish's CNV output, these tracks are
described in detail in the [inputs](#inputs) section of the user guide further below.

    sawfish discover \
      --threads 16  \
      --ref GRCh38.fa \
      --bam HG002.GRCh38.bam \
      --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XY.bed \
      --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
      --output-dir HG002_discover_dir

#### Single sample joint-call step

The `joint-call` command can be specified in one of two styles, the first is to provide a path to each input discover
directory using one or more `--sample` command-line entries. For the example case the joint-call command in this style
would be as follows:

    sawfish joint-call \
      --threads 16 \
      --sample HG002_discover_dir \
      --output-dir HG002_joint_call_dir

Note that when the `joint-call` step is called using the commands demonstrated in the above example, the reference fasta
and sample bam path specified in the `discover` step are stored and reused in the subsequent `joint-call` step. This can
be convenient for quick analysis of a small number of files on a persistent filesystem.

The `joint-call` step can also be configured in a second style, by directly specifying the reference and specifying the
bam path for each sample in a sample CSV file. In this style the joint-call command would be run as follows:

    cat << END > all_samples.csv
    HG002_discover_dir, HG002.GRCh38.bam
    END

    sawfish joint-call \
      --threads 16 \
      --ref GRCh38.fa \
      --sample-csv all_samples.csv \
      --output-dir HG002_joint_call_dir

Here the `--sample-csv` argument is used instead of the `--sample` argument. Besides allowing the specification of
new bam file paths, this can be convenient as a way to specify larger numbers of input samples, see the trio example
below for an example of how this file is input for a multi-sample analysis.

Whichever style is used to run the `joint-call` command, the primary output of this step can be found in
`HG002_joint_call_dir/genotyped.sv.vcf.gz`. See the [outputs](#outputs) section below for discussion of this and all
other output files.

### Quickstart Example 2: Joint calling SVs and CNVs on a trio

To call SVs and CNVs jointly on multiple samples, run `discover` separately on each mapped sample bam, and then run `joint-call`
on all `discover` step output directories.

The following example shows how this is done for mapped sequences from the HG002 trio, given the following bam files:
`HG004.GRCh38.bam`, `HG003.GRCh38.bam`, `HG002.GRCh38.bam`.

#### Trio discover step

As a first step, `discover` needs to be run on all 3 samples. In the example below 16 threads are used to process each
sample. Note that these 3 commands are independent and could be run in parallel.

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG004.GRCh38.bam \
      --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XX.bed \
      --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
      --output-dir HG004_discover_dir

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG003.GRCh38.bam \
      --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XY.bed \
      --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
      --output-dir HG003_discover_dir

    sawfish discover \
      --threads 16 \
      --ref GRCh38.fa \
      --bam HG002.GRCh38.bam \
      --expected-cn ${DISTRO_ROOT_DIR}/data/expected_cn/expected_cn.hg38.XY.bed \
      --cnv-excluded-regions ${DISTRO_ROOT_DIR}/data/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz \
      --output-dir HG002_discover_dir

#### Trio joint-call step

After all discover steps have completed, everything is ready to run the joint-call step. As discussed for the
single-sample example, there are two command-line styles that can be used to provide the sample inputs.

The first style is shown in the command below, where the `--sample` option is provided multiple times to specify
the 3 discover step results. When using this approach the reference and per-sample bam paths provided in the discover
steps above will be re-used for joint-calling:

    sawfish joint-call \
      --threads 16 \
      --sample HG004_discover_dir \
      --sample HG003_discover_dir \
      --sample HG002_discover_dir \
      --output-dir HG002_trio_joint_call_dir

Just as in the single-sample case, note that the reference fasta and all 3 sample bam paths specified in the `discover`
steps are stored and reused in the subsequent `joint-call` step.

The `joint-call` step can also be configured in a second style, by directly specifying the reference and specifying the
bam path for each sample in a sample CSV file. In this style the joint-call command would be run as follows:

    cat << END > all_samples.csv
    HG002_discover_dir, HG002.GRCh38.bam
    HG003_discover_dir, HG003.GRCh38.bam
    HG004_discover_dir, HG004.GRCh38.bam
    END

    sawfish joint-call \
      --threads 16 \
      --ref GRCh38.fa \
      --sample-csv all_samples.csv \
      --output-dir HG002_trio_joint_call_dir

Here, as in the single-sample example above, the `--sample-csv` argument is used instead of the `--sample` argument,
allowing direct specification of all bam sample paths and unifying all input sample information to one file.

Whichever style is used to run the `joint-call` command, the primary output of this step can be found in
`HG002_trio_joint_call_dir/genotyped.sv.vcf.gz`. See the [outputs](#outputs) section below for discussion of this and
all other output files.

### Upgrading to sawfish v2 from previous sawfish versions

The sawfish version 2 release adds a substantial new CNV calling and integration feature. For users switching from
previous sawfish versions this can largely be treated as a gradual change, in that the accuracy of smaller SVs remains
just as high and computational resource demands are similar. However, the following differences should be considered:

1. The final output is now VCF v4.4, which introduces some subtle changes, notably `SVLEN` is now the absolute value of
the SV size.

2. Specifying expected copy number regions with the  `--expected-cn` argument is now an important configuration input
for getting meaningful SV/CNV results from the sex chromosomes. Sawfish's previous behavior to change SV ploidy as a
function of expected copy number is no longer the default behavior, see full details in the [expected copy number
section](#expected-copy-number).

3. Specifying excluded CNV regions with the `--cnv-excluded-regions` argument is another important new configuration
input for improving CNV precision. See the [CNV excluded regions section](#cnv-excluded-regions) for full details.

## Inputs

### Discover step inputs

#### Input read alignments

HiFi read alignments for the query sample must be supplied in BAM or CRAM format as an argument in the discover step.

Sawfish has been tested with HiFi sequencing reads mapped by [pbmm2](https://github.com/PacificBiosciences/pbmm2). In
general it is designed to work on supplementary alignments without hard-clipping. If this requirement is fulfilled it
may work with other mappers, but no others are tested or supported.

When joint-calling over multiple samples, all input alignment files must have been mapped to the same reference genome.

#### Reference fasta

A genome reference sequence file in fasta format is required as input for every run at the discover step, as specified
by the `--ref` argument. Every chromosome name in the input read alignment file must be be present in the reference
sequence file. Their is no reciprocal requirement, the reference fasta may contain chromosome names not present in the
input bam file.

#### Expected copy number

An BED file can be provided for each sample during the discover step to set expected copy number per region of the
genome, by using the `--expected-cn` option. Any regions not specified will have a default expected copy number of 2. If
no file is specified the default expected copy number of 2 will apply to the whole genome.

The expected copy number is important in determining which copy number segments will be output as CNV deletions or
duplications, and is especially useful to indicate the expected copy number for mammalian sex chromosomes.

Pre-generated expected copy number BED files are provided for some common human reference genomes and sex chromosome
complements in the sawfish [expected_cn](../data/expected_cn/) directory. These can be used directly or serve as
templates for other sample configurations. As an example, the BED file for hg38 and karyotype XY is:

```
chrX	0	2781479	chrX_PAR_1	2
chrX	2781479	155701382	chrX_uniq_1	1
chrX	155701382	156040895	chrX_PAR_2	2
chrY	0	2781479	chrY_PAR_1	0
chrY	2781479	56887902	chrY_uniq_1	1
chrY	56887902	57227415	chrY_PAR_2	0
```

As demonstrated in this example, the expected copy number file must be in BED format, with the first 3 columns used to
specify regions following standard BED format. Expected copy number must be provided in column 5. Column 4 is ignored
and can be used as a region label.

During the joint-call step, the expected copy number inputs from the discover phase are retained per sample allowing,
e.g. the expected copy number on the chrX nPAR region to vary by sample sex chromosome complement in a human pedigree
analysis.
ß
Note that in earlier versions of sawfish before CNV output was introduced, the expected copy number input would change
the ploidy used by the genotyper for breakpoint-based calls. Since sawfish v2, this ploidy change is no longer made by
default, and for general purpose calling this change is no longer recommended. If the previous ploidy-change behavior is
preferred, the `--treat-single-copy-as-haploid` option can be provided in the joint-call step to cause any region with
copy number 1 to be treated as haploid (all other cases will continue to be treated as diploid).

#### CNV excluded regions

Certain regions of each reference genome may present inherent difficulties to the prediction of meaningful CNV calls.
Such regions can be marked as excluded for the purpose of CNV calling by providing the regions in BED file format using
the discover step argument `--cnv-excluded-regions`. The way that excluded regions impact CNV calling is summarized in
the section further below, but note that these regions do not change the behavior of breakpoint-based SV calling.

Pre-computed CNV excluded regions are provided in the sawfish [cnv_excluded_regions](../data/cnv_excluded_regions/)
directory for some common human reference genomes. The recommended exclusion track for GRCh38 is
`annotation_and_common_cnv.hg38.bed.gz`. This is the only reference genome for which a combined annotation and common CNV based
exclude region set is provided. Annotation-based excluded regions include assembly gaps,centromeres, and alpha satellite
sequences. Common CNV excluded regions specify where sawfish calls the same CNV type in a high fraction of samples within a
diverse sample cohort. In this case the common CNV regions indicate that the CNV type is present in at least 50% of the 47
[HPRC year1 cohort samples](https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/sample_metadata/hprc_year1_sample_metadata.txt).

For other genomes, 'annotation_only' excluded region tracks are provided. These files are named with the pattern
`annotation_only.${genome_tag}.bed.gz`, where the genome_tag value may be `hg38`, `hg19` or `hs37d5`. For reference
genomes other than GRCh38, these files provide some level of exclusion to reduce false positives, and should produce a
better result than not excluding any regions. For other reference genomes, it is recommended to develop a similarly
expanded exclusion region track including common sawfish CNV calls from a background cohort. The methods used to
produce both annotation and common-cnv based excluded regions are described in the [cnv_excluded_regions script
directory](../scripts/cnv_excluded_regions), which can be used to extend support to additional reference genomes.

##### How excluded regions influence copy number calling

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

#### Minor allele frequency

On the discover step command line for each sample, a small-variant VCF file for that sample can be provided as an
argument to `--maf`. The given VCF will be parsed to create a minor allele frequency track for the genome. This
information is written to an IGV visualization track for assessment and interpretation of the CNV output. This
information may also be used for improved segmentation and CNV calling in the future, although it is not used in the
current release.

Any small-variant caller VCF with an `AD` entry for the small-variant calls could work for this purpose, but the feature
is tested and best supported for the output from [DeepVariant](https://github.com/google/deepvariant).

### Joint-call step inputs

#### Reference fasta

A reference fasta file path can optionally be specified by the `--ref` argument for the `joint-call` step. If not specified,
the reference fasta file path will be taken from the first sample discover step data.

Directly specifying the reference path can be helpful when the filesystem context of the `discover` and `joint-call` steps is
not the same and therefore the file is located on a different path. Note that the reference sequence contained in the file
should be the same as that specified to all input sample discover-step runs.

#### Sample input

One or more samples can be specified to the joint-call step. There are two ways to specify the input samples, only one
of these two approaches can be used at a time.

##### Specify sample input directly on the command-line

In this approach, each sample's sawfish discover output directory is specified using the `--sample` command-line
argument. This argument can be specified once for each input sample. When using this approach, the bam path used in the
sawfish discover step will be re-used for the joint-call step, so it is assumed the bam is still available in the same
location.

##### Specify sample input using a CSV file

Using this approach, all samples inputs are listed in a CSV file, which is provided using the `--sample-csv`
command-line argument. This approach also (optionally) allows a new bam file path to be provided for each sample, which
can be useful if the input bam files are in a different location compared to when the sawfish discover step was run.

The sample CSV file format is designed to be relatively flexible and has the following requirements:
- The first column of each record must describe each sample's sawfish discover output directory
- The second column is optional, and can be used to specify the sample's bam file path
    - If the column is not present or blank, the bam file path will be extracted from the discover output directory (ie. it will be the same path provided on
the discover step commmand-line)
- Records can have differing number of columns, and for any record with more than 2 columns, the extra content will be
ignored
- The file can contain comments starting with the `#` character

## Outputs

### Joint-call step outputs

The primary output of the joint-calling step are the SV and CNV calls for all samples in VCF 4.4 format, written to
`${OUTPUT_DIR}/genotyped.sv.vcf.gz`. Details of the SV and representation in this file are provided below.

#### Quality scores

The primary quality metrics for each variant call are:
- `QUAL` - This is the phred-scaled confidence that the given alternate allele exists in the set of analyzed samples
    - For SV calls supported only by breakpoint evidence, it reflects the probability of any non-reference genotype in
      any sample.
    - For CNV calls supported only by depth evidence, it reflects the probability of that the segment copy number is not
      the expected copy number in any sample.
    - For merged SV/CNV calls, QUAL reflects the maximum of the breakpoint and depth-based QUAL values from the merged
      components
- `GQ` - This value is provided once for each sample in a breakpont-based call. It is the phred-scaled confidence that
   the given sample genotype is correct in this sample based on supporting read evidence at the breakpoint.
- `CNQ` - This value is provided once for each sample in a depth-based call. It is the phred-scaled confidence that the
   given copy number (`CN`) value is correct for this segment in the given sample.

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
`REF` field and any breakpoint insertion sequence in `ALT`. Deletions larger than 100kb, or lacking breakpoint support
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

##### CNVs

Per the above sections, most CNVs will be described as deletions `SVTYPE=DEL` or duplications `SVTYPE=DUP` whether or
not they are merged with a breakpoint-based SV call.

For CNVs that are not merged to a breakpoint-based SV in a multi-sample analysis, it is possible for some samples to show a copy
number gain and other samples to show a copy number loss of the same genomic interval. Where such cases occur the output record
will be given `SVTYPE=CNV`, with a `<CNV>` symbolic alt allele.

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

#### Genotype formatting

The sawfish genotype output is designed to follow the VCF 4.4 spec wherever possible, but the following notes should
supplement the spec to help interpret these results.

##### Overlapping SV formatting

All sawfish SVs are output so that only one allele is described in each VCF record, even if an overlapping SV allele is
output at the same locus. The internal SV calling model accounts for up to 2 overlapping alleles per sample during
genotyping and quality scoring. Reads which support a 2nd alternate allele at any given locus will be counted as
supporting the reference in output fields such as allele depth (`AD`). This protocol matches standard SV caller
formatting conventions. Users interested in a more detailed output format, such as representing overlapping read support
on the VCF `<*>` allele can request this for prioritization.

##### CNV genotypes

CNVs that have not been merged to a breakpoint-based SV call follow a slightly different genotype formatting convention
compared to other SVs in the sawfish output. For a diploid region of the genome, all copy number 0 calls will have
genotype `0/0`, and copy number 1 calls will have genotype `0/1`. For any copy number gain (copy number 3 or higher),
the genotype will be `./1`, reflecting that sawfish has only analyzed the aggregate sample copy number without any
allele-specific copy number estimate. The `./1` genotype reflects that one allele is duplicated, and the other allele's
copy number status is unknown, it may be lost, unchanged, or, for copy number or 4 and up, duplicated as well.

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
value is listed in column 4. Note that any region segmented into the 'excluded' state will be represented as an uncovered
gap in the region coverage.

This file will not appear in the output for any sample run with the `--disable-cnv` discover step option.

##### GC-bias corrected depth track

The bigwig file `gc_bias_corrected_depth.bw` provides binned depth values enumerated from the sample alignment file and
rescaled to correct for the GC-bias pattern inferred from the sample.

This track can be especially useful to visualize and interpret CNV calls. Note that the copy number segmentation model
does not directly operate on the values in this track -- instead it uses the original depth values together with the
local GC-bias estimate for each bin. During segmentation, the GC-bias estimate is used to modify the expected depth
rather than scaling the observed depth.

This file will not appear in the output for any sample run with the `--disable-cnv` discover step option.

##### Depth track

The binned depth values enumerated from the sample alignment file and used as input to the segmentation process are
provided in bigwig format in the file `depth.bw`.

##### Minor allele frequency track

When a minor allele frequency input file is provided for the sample, the corresponding minor allele frequency track will
be output in bigwig format in the file `maf.bw`.This track can be useful to visualize and interpret the CNV output.

#### Variant read support output

To show which reads support each SV allele, the optional `--report-supporting-reads` argument can be added to the
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

### Discover step outputs

The discover step produces a number of output files in the discover output directory used by sawfish during the
subsequent joint calling step. Although these are not fully documented or intended for end users, some of the more
important files are noted below:

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

### Disabling CNV / non-WGS input

Sawfish CNV calling was designed with an assumption of HiFi WGS input. On various types of targeted data the GC-bias
estimation and depth segmentation routines may not be able to complete, or may produce unhelpful results. All CNV
processing can be disabled for given a sample in these cases by specifying the `--disable-cnv` option on the sawfish
discover step.

### Faster analysis for CNV/Large variants only

Sawfish has a faster CNV-focused mode which can be enabled by using the `--fast-cnv-mode` flag in the discover step of
every sample. With this setting, sawfish analyzes only the larger-scale SV breakpoint evidence that could be useful to
improve CNV/large-variant accuracy, together with depth-based CNV analysis. Smaller assembly regions are skipped, which
will remove all insertions and most breakpoint-based deletions below about 1kb.

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

In general, runtime response to thread count is expected to be nearly linear for both sawfish discover and joint-call steps.

For a typical ~30x HiFi sample analyzed on 16 threads, the `discover` step should complete in about 30-40 minutes and
the `joint-call` step should complete in about 5 minutes.

The current joint calling scheme has been designed with pedigree-scale analysis in mind, so runtimes for typical small
pedigrees should be practical. However the runtime is super-linear with sample count, so the method is not practical for
larger cohorts at this time. The following examples should give an idea of what runtimes to expect for different
joint-calling scenarios:

|                | samples | sample type | threads | wall-time | core-hours |
|----------------|:-------:|:-----------:|:-------:|:---------:|:----------:|
| HG002          |    1    |  ~30x human |    16   |   ~5min   |   ~1.3     |
| Plat Ped g2+g3 |   10    | >=30x human |    64   |   ~29min  |   ~31      |
| HPRC Year1     |   47    |  ~30x human |    64   |   ~3hr    |   ~192     |

##### Factors that could lead to longer runtimes

If a given case shows runtime scaling that is considerably longer than the above guidelines for either sawfish step, the
points below may be helpful. One of the factors that could extend runtime is sawfish's alignment file access pattern. In
both the discover and joint-call steps, sawfish will randomly access segments of the alignment file containing reads
associated with candidate SV call breakends. This random access pattern relies on both good file I/O and reasonably fast
decompression of alignment file segments. The following two cases should be considered in this context:

- CRAM input files

CRAM files, and especially 'archival' CRAM with higher compression levels may lead to considerably slower runtimes due
to the burden of random access decompression. It may be worth making a temporary BAM copy of the sample
CRAM file to use during the sawfish analysis in these cases.

- Network storage

Various types of network/cloud file storage system may have poor I/O or poor random-access I/O, even if they perform
well in the context of a caller which reads the alignment file end-to-end. In these cases it may be worth copying the
alignment files onto the local compute node storage (such as a /scratch drive) during sawfish analysis.

#### Memory

The `discover` step should typically require less than 8Gb/thread so long as at least several threads are selected. The
`joint-call` step should require substantially less memory but hasn't been tested at scale with less than 1Gb/thread.

### Discover step input file path storage

In the `joint-call` step, sawfish primarily relies on the files it has written to the discover step output directory for
each sample. For two of the file paths provided as input to the `discover` step, sawfish may rely on being able to
access the original file path provided during the discover step. These two files are the input alignment file (specified
with `--bam`), and the reference fasta file (specified with `--ref`).

Note that the sawfish joint-call step can be configured to eliminate any such original path reuse as follows:
- If the `--ref` argument is provided in the `joint-call` command-line, then the reference file path will not be reused
  from the discover step.
- If the input samples are specified using the `--sample-csv` option, and a bam file path is provided in column 2 for
  every sample, then no bam paths will be reused from the input discover steps.

In the event that the above conditions do not apply, then the following details on how paths are reused may be helpful.

The original file paths used in the discover step for each sample are stored in a configuration file written to the
discover step output directory here:

    ${DISCOVER_STEP_OUTPUT_DIR}/discover.settings.json

These input file paths are normally canonicalized, so that relative paths can be reliably reused after any change to the
working directory. In some cases it may be more convenient to store relative file paths. To do so the `discover` step
option `--disable-path-canonicalization` can be used to store all input paths as-is. This may be useful if e.g., the
`discover` and `joint-call` steps are being run in different directory structures.

Note that for even more complex situations, the paths in the above discover settings json file can be manually edited
before running the `joint-call` step, but in general the above mentioned `--sample-csv` sample specification option
should provide a simpler path customization option.
