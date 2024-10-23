# Sawfish User Guide

## Table of Contents

* [Overview](#overview)
* [Accuracy Evaluation](#accuracy-evaluation)
* [Getting Started](#getting-started)
* [Outputs](#outputs)
* [Usage Details](#usage-details)

## Overview

Sawfish calls structural variants from mapped HiFi sequencing reads. It discovers germline variants from local sequence assembly and jointly genotypes variants across multiple samples.

Key features:
- High SV discovery and genotyping accuracy
  - All variants are modeled and genotyped as local haplotypes, yielding substantial accuracy gains on modern SV truth sets such as the GIAB HG002 T2T SVs. 
- High resolution
  - All structural variants are assembled to basepair resolution and reported with breakpoint homology and insertion details.
- Integrated depth assessment
  - Integrated depth estimation with GC-bias correction is used to classify large deletion and duplication calls for higher precision.
- Simple multi-threaded workflow
  - A single command-line is used for each of the discover and joint-genotyping steps

All SVs are modeled internally as breakpoints, but will be reported as deletions, insertions, duplications and inversions when supported by the corresponding breakpoint and depth pattern, otherwise the breakpoint itself is reported. The minimum variant size is 35 bases (configurable). A maximum size is only applied to inversions (100kb).

## Accuracy Evaluation

Recommended methods for SV accuracy assessment and benchmarking results are described in the [sawfish preprint](https://doi.org/10.1101/2024.08.19.608674). An earlier simplified assessment and benchmarking approach is described in [accuracy](accuracy.md).

## Getting Started

### Installation

Sawfish binaries are available for 64-bit Linux platforms. These can be installed either directly from the GitHub release tarball, or via conda as described below.

#### Install from GitHub

To install sawfish from github, download the latest release tarball compiled for 64-bit Linux on the
[github release channel](https://github.com/PacificBiosciences/sawfish/releases/latest), then unpack the tar file.
Using v0.12.1 as an example, the tar file can be downloaded and unpacked as follows:

    wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.1/sawfish-v0.12.1-x86_64-unknown-linux-gnu.tar.gz
    tar -xzf sawfish-v0.12.1-x86_64-unknown-linux-gnu.tar.gz

The sawfish binary is found in the `bin/` directory of the unpacked file distribution.
This can be run with the help option to test the binary and review latest usage details:

    sawfish-v0.12.1-x86_64-unknown-linux-gnu/bin/sawfish --help

#### Install from conda

For [conda](https://github.com/conda/conda) users, installing sawfish on conda may be a more convenient option. Sawfish is available for conda on Linux from the `bioconda` channel. A new conda environment with the latest sawfish release can be created as follows:

    conda create -n sawfish -c bioconda sawfish

### Analysis Steps

Sawfish analyzes samples in 2 steps:

1. `discover` - The discover step identifies candidate structural variant (SV) regions and assembles each local SV haplotype.
2. `joint-call` - The joint call step takes the output of the sawfish 'discover' step for one to many samples and provides jointly genotyped SV calls over the sample set. Joint calling includes the following operations:
    - Merge duplicate SV haplotypes
    - Associate deduplicated SV haplotypes with samples
    - Evaluate SV read support in each sample
    - Genotype quality assessment and VCF output

### Calling SVs on one sample

To call SVs in one sample, run `discover` on the mapped sample bam, and then run joint call on the output directory of the discover step.

The following example shows how this is done for a mapped sample bam named `HG002.GRCh38.bam`, using 16 threads for both the `discover` and the `joint-call` steps.

    sawfish discover --threads 16 --ref GRCh38.fa --bam HG002.GRCh38.bam --output-dir HG002_discover_dir
    sawfish joint-call --threads 16 --sample HG002_discover_dir --output-dir HG002_joint_call_dir

The final joint calling output can be found in `HG002_joint_call_dir/genotyped.sv.vcf.gz`. See the [outputs](#outputs) section below for discussion of the VCF contents.

Note that the reference fasta and sample bam specified in the `discover` step are still used in the subsequent `joint-call` step, they simply don't need to be specified on the command-line.

### Joint Calling SVs across a set of samples

To call SVs on a set of samples, run `discover` separately on each mapped sample bam, and then run joint call on all `discover` step output directories.

The following example shows how this is done for mapped sequences from the HG002 trio, given the following bam files: `HG004.GRCh38.bam`, `HG003.GRCh38.bam`, `HG002.GRCh38.bam`.

As a first step, `discover` needs to be run on all 3 samples. In the example below 16 threads are used to process each sample. Note that these 3 command-lines could be run in parallel.

    sawfish discover --threads 16 --ref GRCh38.fa --bam HG004.GRCh38.bam --output-dir HG004_discover_dir
    sawfish discover --threads 16 --ref GRCh38.fa --bam HG003.GRCh38.bam --output-dir HG003_discover_dir
    sawfish discover --threads 16 --ref GRCh38.fa --bam HG002.GRCh38.bam --output-dir HG002_discover_dir

After all discover steps have completed, joint calling can be run over all 3 samples using the following command:

    sawfish joint-call --threads 16 --sample HG004_discover_dir --sample HG003_discover_dir --sample HG002_discover_dir --output-dir HG002_trio_joint_call_dir

The final joint calling output can be found in `HG002_trio_joint_call_dir/genotyped.sv.vcf.gz`. See the (outputs)[#outputs] section below for detailed discussion of the output VCF contents.

Just as in the single-sample case, note that the reference fasta and all 3 sample bams specified in the `discover` steps are still used in the subsequent `joint-call` step, but they don't need to be specified on the command-line, since their paths are recorded in the metadata of each discover output data.

## Outputs

### Joint call step

The primary user output of the sawfish SV caller is the SV VCF produced by the joint-call step. This file lists all SVs in VCF 4.2 format. Details of the SV representation in this file are provided below.

#### Quality scores

The primary quality metrics for each SV call are:
1. `QUAL` - This is the phred-scaled confidence that the given SV allele exists in the set of genotyped samples.
2. `GQ` - This value is provided once for each sample. It is the phred-scaled confidence that the given sample genotype is correct.

All phred-scaled quality scores in the VCF output have a maximum value of 999.

#### Filters

The following filters may be applied to each VCF record:

- `ConflictingBreakpointGT` - Genotypes of breakpoints in a multi-breakpoint event conflict in the majority of cases (This filter is only relevant to inversions at present)
- `MinQUAL` - The SV allele quality score (`QUAL`) is less than 10
- `MaxScoringDepth` - Read depth at the SV locus exceeds 1000x, so all scoring and genotyping steps were disabled.
- `InvBreakpoint` - This breakpoint is represented as part of a separate VCF inversion record (the inversion record shares the same EVENT ID)

#### SV Types

Notes on formatting and representation of SVs are listed below for each major type.

##### Deletions

All deletions of 100kb or smaller are represented by directly writing the deleted sequence in the VCF `REF` field and any breakpoint insertion sequence in `ALT`. Deletions larger than 100kb are written as symbolic alleles using the ALT value of `<DEL>`. All candidate deletions at least 50kb in length will be checked for a supporting depth signature, if this support is not found the candidate deletion will be reported in the VCF output as a breakend (`BND`) pair instead.

###### Insertions

Any indel-like SVs where the length of sequence inserted at the breakpoint exceeds the length of deleted sequence will be formatted as an insertion in the VCF output if it is possible to fully assemble the inserted sequence, and will be formatted as a duplication otherwise. If represented as an insertion the full inserted sequence assembly will be written to the VCF `ALT` field.

##### Duplications

Very large insertions with long breakpoint homology will be represented as duplications in the VCF output only if they cannot be output as insertions. These will be written to the VCF output using the symbolic ALT value of `<DUP:TANDEM>`. All candidate duplications at least 50kb in length will be checked for a supporting depth signature, if this support is not found the candidate duplication will be reported in the VCF output as a breakend (`BND`) pair instead.

##### Breakpoints

All SV breakpoints which can't be modeled as one of the simple SV types above will be output as a pair of breakend (`BND`) records.

##### Inversions

Sawfish will currently identify one type of multi-breakpoint complex SV signature, corresponding to that of a simple inversion. Inversions are identified when two intra-chromosomal inverted breakpoints of opposite orientation have overlapping spans with at least an 80% reciprocal overlap. The longer span must not be greater than 100kb.

When an inversion is found, a VCF record will be output using the `<INV>` symbolic allele summarizing the inversion in as much detail as possible. It is not possible to retain the details of all 4 breakends in this format such as all breakend positions and breakpoint insertion sequences. For this reason the corresponding breakend records are retained in the VCF output but marked as filtered, such that full breakend details remain available in the output. The inversion record and the filtered breakend records are given a shared VCF `EVENT` label so that their relationship can be identified.

#### Overlapping SV formatting

All sawfish SVs are output so that only one allele is described in each VCF record, even if an overlapping SV allele is output at the same locus. The internal SV calling model accounts for up to 2 overlapping alleles per sample during genotyping and quality scoring, reads supporting a 2nd alternate allele at any given locus will be counted as support the reference in output fields such as allele depth (`AD`). This protocol matches standard SV caller formatting conventions. Users interested in a more detailed output format, such as representing overlapping read support on the VCF `<*>` allele can request this for prioritization.

#### Phasing

Sawfish adds short-range phasing information to clarify the relationship of heterozygous SVs called from the same or overlapping SV haplotypes. This does not have the range of general read-backed phasing and will only result in phased genotype output for smaller insertions and deletions. Each local cluster of phased genotypes corresponds to a phase set as annotated using the VCF `PS` tag. The phase set ID is the `POS` value of the first SV called from the SV haplotype cluster.

#### Optional variant read support output

To help show which reads support each SV allele, the optional `--report-supporting-reads` argument can be added to the joint-call command line. When this is used a compressed json output file is provided in `${OUTPUT_DIR}/supporting_reads.json.gz`.

In this json output file, the top-level objects are variant IDs matching those provided in the ID field of the VCF output. Nested under each variant ID are sample IDs. For each sample ID associated with a variant, the array of supporting read QNAME values are provided. A simplified example output is shown below for two variants:

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

Note that the number of read QNAME entries should often match the supporting AD count for the alternate allele from the same variant/sample entry in the VCF, but this is not always an exact match. Also to keep a consistent relationship between supporting reads and variants, no output is provided for VCF records with the inversion (`<INV>`) allele type, but the supporting reads for the breakends comprising each inversion are provided.

### Discover step

The discover step produces a number of output files in the discover output directory used by sawfish during the subsequent joint calling step. Although these are not intended for direct use, some of the important files are described here:

- `assembly.regions.bed` - Describes each region of the genome targeted for assembly.
- `contig.alignment.bam` - This is a BAM file containing the SV locus contigs aligned back to the genome to create candidate SVs for each sample.
- `candidate.sv.bcf` - These are the candidate SVs expressed in a simplified format for each sample. These are used as input for joint genotyping together with the aligned candidate contigs.
- `discover.settings.json` - Various parameters from the discover step (either user input or default) are recorded in this file. Some of the paths to files like the sample bam and reference fasta will be reused in the joint call step.

### Debug outputs

In either run step, the following files are produced to help debug problematic runs:
1. `${OUTPUT_DIR}/sawfish.log` - High level logging output
2. `${OUTPUT_DIR}/run_stats.json` - Run statistics and component timings

## Usage Details

### Read mapper

Sawfish has been tested with sequencing reads mapped by [pbmm2](https://github.com/PacificBiosciences/pbmm2). In general it is designed to work on supplementary alignments without hard-clipping. If this requirement is fulfilled it may work with other mappers, but no others are tested or supported.

### Determinism

Sawfish should always produce the same output from a given command-line and input file set (allowing for expected changes in timestamps, benchmark timers and similar metadata).

### Output directory / clobber

Each step of the pipeline accepts the argument `--output-dir` where all files from the step will be written. If not specified the default of either `sawfish_discover_output` or `sawfish_joint-call_output` will be used. Sawfish will not proceed if the output directory already exists, unless the `--clobber` argument is given as well.

### VCF ID Field

The entries in the output VCF ID field (such as `sawfish:0:2803:1:2` and `sawfish:INV4`) are designed to guarantee a unique identifier for each record in the VCF output. This identifier isn't meant to convey useful details about the call and may be reformatted in future releases.

### Expected compute requirements

#### Runtime

For a typical ~30x HiFi sample analyzed on 16 threads, the `discover` step should complete in about 30-40 minutes and the `joint-call` step should complete in about 5 minutes.

Running the `joint-call` step on 10 samples at 30-100x depth completes in about 1 hour on 64 threads.

In general, runtime response to thread count is expected to be nearly linear. The current joint calling scheme has been designed with pedigree-scale analysis in mind. Sawfish joint calling has completed on 47 HPRC samples in testing, but substantially larger cohorts would be difficult without further changes to the joint-calling design.

#### Memory

The `discover` step should typically require less than 8Gb/thread so long as at least several threads are selected. The `joint-call` step should require substantially less memory but hasn't been tested at scale with less than 1Gb/thread.

### Haploid regions

The SV caller `discover` step accepts a specially formatted BED file format which specifies expected copy number/ploidy by genome region. By default all regions of the genome are treated as diploid, so these files only need to specify non-diploid regions.  

In the copy number BED file, the first 3 columns are used to specify regions following standard BED format, and expected copy number will be read from column 5 of the input BED file. Column 4 is ignored and can be used as a region label. For the purpose of SV calling, any region with copy number 1 will be treated as haploid and all other values will be treated as diploid.

Expected copy number BED files are typically used to specify ploidy in the non-PAR regions of the sex chromosomes. For example, in the example discover step for HG002, we can additionally specify an `--expected-cn` argument as follows:

    sawfish discover --threads 16 --ref GRCh38.fa --bam HG002.GRCh38.bam --output-dir HG002_discover_dir --expected-cn ${SAWFISH_DIR}/data/expected_cn/expected_cn.hg38.XY.bed

The file `expected_cn.hg38.XY.bed` contains:
```
chrX    0       2781479 chrX_PAR_1      2
chrX    2781479 155701382       chrX_uniq_1     1
chrX    155701382       156040895       chrX_PAR_2      2
chrY    0       2781479 chrY_PAR_1      0
chrY    2781479 56887902        chrY_uniq_1     1
chrY    56887902        57227415        chrY_PAR_2      0
```

...expected sex chromosome copy number files for this and other references can be found in the [expected_cn](../data/expected_cn) directory.

All expected copy number files submitted for each sample at the discover phase are saved in the discover directory and used to select per-sample ploidy in the specified regions during the joint-calling step.

### Discover step input file path storage

Sawfish accesses several files associated with each sample during joint-genotyping in the `joint-call` step. For instance, this is done to test read support for each allele by accessing the sample alignment file.

To find these files for each sample, input file paths are stored from the `discover` step in a configuration file written to the output directory here:

    ${OUTPUT_DIR}/discover.settings.json

These input file paths are normally canonicalized, so that relative paths can be reliably reused after any change to the working directory. In some cases it may be more convenient to store relative file paths. To do so the `discover` step option `--disable-path-canonicalization` can be used to store all input paths as-is. This may be useful if e.g., the `discover` and `joint-call` steps are being run in different directory structures.

Note that for even more complex situations, the paths in the above discover settings json file can be manually edited before running the `joint-call` step.
