# Sawfish SV accuracy assessment on HG002

Sawfish assembles and describes multiple overlapping SV alleles at loci where older SV benchmarks tend to compress overlapping variants down to a single allele, which can complicate interpretation of benchmarking results. We describe best benchmarking practices to mitigate these interpretation challenges here. These methods are also detailed in the [sawfish App Note in *Bioinformatics*](https://doi.org/10.1093/bioinformatics/btaf136). The exact SV VCF files, assessment results, and more information on assessment script details from the App Note are also [available on Zenodo](https://doi.org/10.5281/zenodo.14898462).

## Results

Single-sample SV calling performance on HG002 can be assessed with multiple GIAB SV benchmark sets. We demonstrate results for CMRG v1.0 and a recent draft T2T SV benchmark set below to
demonstrate performance on both medically-relevant and genome-wide SV calling problems. The expected outcome of the methods below should match those reported for sawfish in Tables S1
and S4 from the App Note as follows:

### CMRG

| Method  | F1     | Recall | #FN | Precision | #FP |
|:-------:|:------:|:------:|:---:|:---------:|:---:|
| sawfish | 0.9929 | 0.9907 | 2   | 0.9951    | 1   |

### T2T

| Method  | F1     | Recall | Precision |
|:-------:|:------:|:------:|:---------:|
| sawfish | 0.9763 | 0.9625 | 0.9905    |


## Methods

Below we describe several methods required to setup the analysis for both GIAB benchmark sets on HG002, then detail the separate assessment procedures used for the CMRG and T2T benchmark sets.

### Downloading the GRCh38 reference

All analyses use the reference `human_GRCh38_no_alt_analysis_set.fasta`, [described here](https://github.com/PacificBiosciences/reference_genomes/tree/main/reference_genomes/human_GRCh38_no_alt_analysis_set).

The reference fasta file can be downloaded and unpacked as follows:

    wget https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
    tar -xzf human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz

After unpacking, the reference fasta path will be `human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta`.

### Creating conda environment for supporting analysis tools

The analysis uses pbmm2 1.13.1 for read mapping, and sawfish 0.12.9 for SV calling. These software versions match those used in the sawfish App Note. Assessment
is run in a truvari singularity image, so singularity will be required as well. The recommended way to recreate this analysis is to make a new conda environment for
these dependencies as follows.

    conda create -y -n sawfish_accuracy_test_v2 python=3.9 bioconda::pbmm2=1.13.1 bioconda::sawfish=0.12.9 conda-forge::singularity

### Downloading truvari singularity image

All assessment below uses truvari 4.2.2 with MAFFT for refinement. The CMRG assessment uses the truvari bench command alone, but for the T2T assessment it is important
to use the bench and refine command, with the bench/refine results consolidated as described further below.

We have found that the most reliable way to run this version of truvari with MAFFT is from a singularity image. The singularity image file used in the App Note can be
pulled from `cloud.sylabs.io` as follows:

    conda activate sawfish_accuracy_test_v2
    singularity pull --arch amd64 library://ctsa/truvari/truvari:4.2.2

Store the downloaded singularity image file `truvari_4.2.2.sif` in a location it can be referenced in the subsequent steps to this guide.

### Downloading and mapping HiFi Data

Download the example HiFi Revio sequencing reads for HG002 here:

    https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/m84011_220902_175841_s1.hifi_reads.bam

These reads are mapped using pbmm2 1.13.1, see the above section on conda environment creation for a recommended installation scheme.

In the `sawfish_accuracy_test_v2` conda environment, the example sequencing reads can be mapped as follows:

    conda activate sawfish_accuracy_test_v2
    pbmm2 --log-level INFO align -j 16 --preset CCS --sample HG002 --sort --unmapped human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta m84011_220902_175841_s1.hifi_reads.bam m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam

Thread count can be adjusted to fit local compute resources. The resulting output bam `m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam` is used for all downstream SV calling steps.

### Calling SVs with sawfish

SVs are called from mapped reads using sawfish 0.12.9, see the above section on conda environment creation for a recommended installation scheme. The following
commands can be used to replicate the analysis from the App Note:

    conda activate sawfish_accuracy_test_v2
    sawfish discover --threads 16 --ref human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta --bam m84011_220902_175841_s1.pbmm2-1.13.1.GRCh38.bam
    sawfish joint-call --threads 16 --sample sawfish_discover_output

Thread count can be adjusted to fit local compute resources. After running these commands the final sawfish SV output will be found in `sawfish_joint-call_output/genotyped.sv.vcf.gz`.


### CMRG Assessment methods

#### Downloading GIAB CMRG SV benchmark

The CMRG SV benchmark was obtained from the following URL:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00

The specific CMRG benchmark files for GRCh38 SVs are:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed


#### Running CMRG assessment with truvari bench

After completing all of the download and configuration steps above, truvari can be run to assess the sawfish HG002 results using the script segment described below.
This script also demonstrates how truvari can be run through a singularity image. Note that we illustrate singularity path bindings that will work for the specified
input file paths, which assume all above download and configuration operations have been run from a common working directory. The binding scheme will likely require
adaption for each user's file organization scheme:

```
# Activate conda environment to use singularity:
#
conda activate sawfish_accuracy_test_v2

sif=truvari_4.2.2.sif

ref=human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta
bench_vcf=HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
bench_bed=HG002_GRCh38_CMRG_SV_v1.00.bed
query_vcf=sawfish_joint-call_output/genotyped.sv.vcf.gz

output_dir=truvari_cmrg_result

singularity exec --bind $working_dir:/data  --cleanenv --pwd /data $sif \
truvari bench \
    --reference $ref \
    --includebed $bench_bed \
    --base $bench_vcf \
    --comp $query_vcf \
    --output $output_dir \
    --passonly \
    --pick ac \
    --dup-to-ins
```

On completion, the analysis results can be found in `truvari_cmrg_result/summary.json`.


### T2T Assessment methods

#### Downloading and modifying the GIAB draft T2T SV benchmark

We assess genome-wide SV accuracy using the GIAB draft SV benchmark (V0.019-20241113) based on the T2T-HG002-Q100v1.1 diploid assembly. These benchmark files were obtained from the following URL:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113

The specific T2T benchmark files for GRCh38 SVs are:

    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz.tbi
    https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed

These files are modified for the app note analysis as follows.

1. The benchmark files are modified to select for structural variants only as follows:

    bcftools view -i 'INFO/SVTYPE!="."' -Oz -o GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
    bcftools index -t GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz

2. The confident regions are restricted to autosomes as follows:

    grep -v chr[XY] GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed > GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed


#### Running T2T assessment with truvari bench/refine

After completing all of the download and configuration steps above, truvari can be run to assess the sawfish HG002 results using the script segment described below.
This script also demonstrates how truvari can be run through a singularity image. Note that we illustrate singularity path bindings that will work for the specified
input file paths, which assume all above download and configuration operations have been run from a common working directory. The binding scheme will likely require
adaption for each user's file organization scheme:

```
# Activate conda environment to use singularity:
#
conda activate sawfish_accuracy_test_v2

sif=truvari_4.2.2.sif

ref=human_GRCh38_no_alt_analysis_set/human_GRCh38_no_alt_analysis_set.fasta
bench_vcf=GRCh38_HG2-T2TQ100-V1.1_stvar.svtype.vcf.gz
bench_bed=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.only_autosomes.bed
query_vcf=sawfish_joint-call_output/genotyped.sv.vcf.gz

# Truvari bench requires that the output directory does not already exist:
output_dir=truvari_t2t_result
rm -rf $output_dir

# Copy this truvari post-processing script from the sawfish repo to the local directory to simplify usage inside of the singularity container:
cp -f ${SAWFISH_REPO_ROOT_DIR}/scripts/truvari_utils/process_truvari_ga4gh_vcfs.py .

working_dir=$(pwd)

# MAFFT (used in the refine step) can be sensitive to the space available in the default tmp directory, so create one on the regular filesystem:
tmp_dir=$(mktemp -d)
trap 'rm -rf $tmp_dir' EXIT

# Step 1: Run bench
singularity exec --bind $working_dir:/data --cleanenv --pwd /data $sif \
truvari bench \
    --reference $ref \
    --includebed $bench_bed \
    --base $bench_vcf \
    --comp $query_vcf \
    --output $output_dir \
    --passonly \
    --pick ac \
    --dup-to-ins

# Step 2: Run refine
singularity exec --bind $working_dir:/data --bind /tmp:$tmp_dir --cleanenv --pwd /data $sif \
truvari refine \
    --reference $ref \
    --regions $output_dir/candidate.refine.bed \
    --recount \
    --use-region-coords \
    --use-original-vcfs \
    --threads 16 \
    --align mafft \
    $output_dir

# Step 3: Produce final assessment VCF outputs from the combined bench/refine assessment
ga4gh_prefix=${output_dir}/ga4gh_with_refine
singularity exec --bind $working_dir:/data --cleanenv --pwd /data $sif \
truvari ga4gh \
    --input $output_dir \
    --output $ga4gh_prefix \
    --with-refine

# Step 4: Summarize performance from final VCF outputs
#
# Note this script is copied from the sawfish repository scripts directory above.
#
./process_truvari_ga4gh_vcfs.py --truth-vcf ${ga4gh_prefix}_truth.vcf.gz --query-vcf ${ga4gh_prefix}_query.vcf.gz >| ${ga4gh_prefix}.size_stratified.accuracy.stats.txt
```

On completion, the analysis results can be found in `truvari_t2t_result/ga4gh_with_refine.size_stratified.accuracy.stats.txt`.
