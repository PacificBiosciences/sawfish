# Sawfish accuracy assessment

Sawfish tends to assemble and describe multiple overlapping SV alleles at loci where older truth sets may compress overlapping variants down to a single allele, which can complicate interpretation of benchmarking results. For this reason we describe some improved sawfish accuracy evaluation approaches here.

## Assessing genome-wide SV accuracy on HG002

We assess general SV accuracy against the new GIAB T2T SV truth set for HG002. Both sawfish and the GIAB T2T truth set include more detailed overlapping SV alleles compared to previous
truth sets, complicating benchmarking. Optimizing benchmarking methods for this case is an active area of community methods development. For the time being, we have found `hap-eval` to handle this benchmarking problem in a fairly effective manner without an SV phasing requirement. We additionally restrict evaluation to autosomes to simplify assessment.

Result using this benchmarking approach are:

#### HG002 GIAB T2T 1.0 SV assessment using hap_eval

| Method  | Recall | Precision | F1    |
|:-------:|:------:|:---------:|:-----:|
| sawfish | 0.958  | 0.981     | 0.969 |
| pbsv    | 0.922  | 0.943     | 0.932 |


### Details to reproduce HG002 accuracy assessment

The above analysis used the sawfish v0.10.0 release binary. HG002 HiFi data and reference files were downloaded as described below. The following commands were used to produce the benchmarked SV call output:

    sawfish discover --threads 16 --ref human_GRCh38_no_alt_analysis_set.fasta --bam HG002.m84011_220902_175841_s1.GRCh38.bam
    sawfish joint-call --threads 16 --sample sawfish_discover_output

After running these commands the final sawfish SV output will be found in `sawfish_joint-call_output/genotyped.sv.vcf.gz`.

The SVs are benchmarked against the GIAB T2T v1.0 truth set. Download and sex chromosome exclusion details are described below.
Instructions to download the version of `hap-eval` used here together with an example installation procedure are also described below.

Given the sawfish SV calls, GIAB T2T truth set files and a version of `hap-eval` activated in the current environment, the benchmarking command is:

    hap_eval --reference human_GRCh38_no_alt_analysis_set.fasta --base GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz --interval GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.only_autosomes.bed --comp sawfish_joint-call_output/genotyped.sv.vcf.gz

SV calls from pbsv (see below for download paths) are benchmarked using the same method.

#### HiFi Data

Download the example HiFi Revio data for HG002 from the following paths:

    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam.bai

#### GRCh38 reference

The linked HiFi BAMs and all analyses use the reference `human_GRCh38_no_alt_analysis_set.fasta`, [described here](https://github.com/PacificBiosciences/reference_genomes/tree/main/reference_genomes/human_GRCh38_no_alt_analysis_set).

The reference fasta file can be downloaded and unpacked as follows:

    wget https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
    tar -xzf human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz

After unpacking, the reference fasta path will be `human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta`.

#### pbsv calls

pbsv calls for the same HG002 sample are available from the following paths:

    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.pbsv.vcf.gz
    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.pbsv.vcf.gz.tbi

#### GIAB T2T SV truth set

The SV truth set is obtained from the following GIAB directory:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/

The specific T2T truth set version for GRCh38 SVs is:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed

For this analysis the confident regions are additionally restricted to autosomes only as follows:

    grep -v chr[XY] GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed > GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.only_autosomes.bed

#### hap-eval

The latest version of `hap-eval` can be obtained from Sentieon here: https://github.com/Sentieon/hap-eval. To match the exact version from this analysis, use
the following branch:

    https://github.com/ctsa/hap-eval/tree/sawfish-t2t-example-stable

Follow the `hap-eval` README for full installation instructions. An example installation procedure using conda is:

    conda create -y -n hap_eval_0327909 python=3.10
    conda activate hap_eval_0327909
    git clone --branch sawfish-t2t-example-stable --recurse-submodules https://github.com/ctsa/hap-eval.git
    pip install ./hap-eval
