# Sawfish accuracy assessment

Sawfish tends to assemble and describe multiple overlapping SV alleles at loci where older SV benchmarks tend to compress overlapping variants down to a single allele, which can complicate interpretation of results. For this reason we describe some effective sawfish accuracy evaluation approaches here.

## Assessing SV accuracy on HG002

Single-sample SV calling performance on HG002 can be assessed with multiple GIAB SV benchmark sets. We demonstrate results for CMRG v1.0 and a recent draft T2T SV benchmark set below.

### Results

#### HG002 GIAB CMRG assessment

SV accuracy in medically relevant genes can be assessed by comparing calls against the GIAB CMRG SV benchmark set. Sawfish SV calls are generated and then compared to the benchmark using truvari, as detailed below in [methods](#methods).

This accuracy assessment shows a reduction in false SV calls from sawfish:

| Method  | Recall | #FN | Precision | #FP | F1    |
|:-------:|:------:|:---:|:---------:|:---:|:-----:|
| sawfish | 0.991  | 2   | 0.995     | 1   | 0.993 |
| pbsv    | 0.958  | 9   | 0.985     | 3   | 0.971 |


#### HG002 GIAB T2T SV assessment using hap_eval

We assess general SV accuracy against the GIAB T2T SV benchmark set for HG002. Both sawfish and the GIAB T2T benchmark set include more detailed overlapping SV alleles compared to previous SV benchmark sets, complicating assessment of accuracy between the SV caller output and the benchmark set. Optimizing assessment methods for this case is an active area of community methods development. For the time being, we have found `hap-eval` to handle this assessment problem in a fairly effective manner without an SV phasing requirement. We additionally restrict evaluation to autosomes to simplify assessment. Details of SV call generation and assessment methods are provided in [methods](#methods) below.

The resulting accuracy assessment shows that sawfish yields a substantial improvement in both recall and precision:

| Method  | Recall | Precision | F1    |
|:-------:|:------:|:---------:|:-----:|
| sawfish | 0.962  | 0.981     | 0.971 |
| pbsv    | 0.920  | 0.940     | 0.930 |


### Methods

The above analyses use the sawfish v0.10.0 release binary. HG002 HiFi data and reference files were downloaded and mapped with pbmm2 as described below. The following commands were then used to produce the benchmarked sawfish SV call output:

    sawfish discover --threads 16 --ref human_GRCh38_no_alt_analysis_set.fasta --bam m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam
    sawfish joint-call --threads 16 --sample sawfish_discover_output

After running these commands the final sawfish SV output will be found in `sawfish_joint-call_output/genotyped.sv.vcf.gz`.

For CMRG assessment, the SVs were assessed against the GIAB CMRG v1.0 benchmark set using truvari with following command:

    truvari bench -r 1000 --passonly --pick ac -f human_GRCh38_no_alt_analysis_set.fasta -b HG002_GRCh38_CMRG_SV_v1.00.vcf.gz --includebed HG002_GRCh38_CMRG_SV_v1.00.bed -c sawfish_joint-call_output/genotyped.sv.vcf.gz -o out

Benchmark set links and truvari installation details are given below.

For general accuracy assessment, the SVs are assessed against a recent GIAB T2T benchmark set. The benchmark set download link and sex chromosome exclusion details are described below. Instructions to download the version of `hap-eval` used here together with an example installation procedure are also described below.

Given the sawfish SV calls, GIAB T2T benchmark set files and a version of `hap-eval` activated in the current environment, the benchmarking command is:

    hap_eval --reference human_GRCh38_no_alt_analysis_set.fasta --base GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz --interval GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.only_autosomes.bed --comp sawfish_joint-call_output/genotyped.sv.vcf.gz

SV calls from pbsv are assessed on CMRG and T2T benchmark sets using the same methods.

#### Downloading the GRCh38 reference

All analyses use the reference `human_GRCh38_no_alt_analysis_set.fasta`, [described here](https://github.com/PacificBiosciences/reference_genomes/tree/main/reference_genomes/human_GRCh38_no_alt_analysis_set).

The reference fasta file can be downloaded and unpacked as follows:

    wget https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
    tar -xzf human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz

After unpacking, the reference fasta path will be `human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta`.

#### Creating conda environment for other PacBio tools

The analysis uses pbmm2 1.13.1 for read mapping and pbsv 2.9.0 for SV comparison. The recommended way to recreate this analysis is to create a new conda environment for these
dependencies:

    conda create -y -n sawfish_accuracy_test_v1_test -c bioconda pbmm2=1.13.1 pbsv=2.9.0

#### Downloading and mapping HiFi Data

Download the example HiFi Revio sequencing reads for HG002 here:

    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/m84011_220902_175841_s1.hifi_reads.bam

These reads are mapped using pbmm2 1.13.1, see the above section on conda environment creation for a recommended installation scheme.

In the `sawfish_accuracy_test_v1` conda environment, the example sequencing reads can be mapped as follows

    pbmm2 --log-level INFO align -j 16 --preset CCS --sample HG002 --sort --unmapped human_GRCh38_no_alt_analysis_set.fasta m84011_220902_175841_s1.hifi_reads.bam m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam

The resulting output bam `m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam` is used for all downstream SV calling steps.

#### Generating pbsv calls

This analysis includes SV calls generated using pbsv 2.9.0, see the above section for a recommended installation scheme.

The pbsv calls were generated by first downloading the GRCh38 tandem repeat track:

    wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/077573d4cacf5632d3ecf125220526e12ef61309/annotations/human_GRCh38_no_alt_analysis_set.trf.bed

From the `sawfish_accuracy_test_v1` conda environment, pbsv calling was run as follows:

    pbsv discover --hifi -b human_GRCh38_no_alt_analysis_set.trf.bed m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.svsig.gz
    pbsv call --hifi -j 8 -t INS,DEL human_GRCh38_no_alt_analysis_set.fasta m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.svsig.gz m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.pbsv.vcf

#### Downloading GIAB CMRG SV benchmark

The CMRG SV benchmark was obtained from the following URL:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00

The specific CMRG benchmark files for GRCh38 SVs are:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed

#### Downloading and modifying the GIAB draft T2T SV benchmark

We assess general SV accuracy compared to the GIAB draft SV benchmark (V0.012-20231107) based on the T2T-HG002-Q100v1.0 diploid assembly. These benchmark files were obtained from the following URL:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107

The specific T2T benchmark files for GRCh38 SVs are:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed

For this analysis the confident regions are additionally restricted to autosomes as follows:

    grep -v chr[XY] GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed > GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.only_autosomes.bed

#### Installing truvari

Truvari 4.1.0 is used for benchmarking calls against the GIAB CMRG truth set. Truvari is run through a conda `truvari_v4.1.0` environment configured as follows:

    conda create -y -n truvari_v4.1.0 -c conda-forge -c bioconda truvari=4.1.0

#### Installing hap-eval

hap-eval is used for benchmarking calls against the GIAB T2T truth set. The latest version of hap-eval can be obtained from Sentieon here: https://github.com/Sentieon/hap-eval. To match the hap-eval version from this analysis, use the following branch:

    https://github.com/ctsa/hap-eval/tree/sawfish-t2t-example-stable

Follow the hap-eval README for full installation instructions. An example of the hap-eval version we used into a new conda `hap_eval_0327909` environment is:

    conda create -y -n hap_eval_0327909 python=3.10
    conda activate hap_eval_0327909
    git clone --branch sawfish-t2t-example-stable --recurse-submodules https://github.com/ctsa/hap-eval.git
    pip install ./hap-eval
