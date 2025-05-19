#!/usr/bin/env bash

set -o errexit
set -o nounset

sawfish_vcf_to_cnv_bedgraphs() {
  if [ -z ${vcf_index+x} ]; then
    declare -g vcf_index=-1
  fi
  vcf_index=$((vcf_index + 1))
  zvi=$(printf "%04d" $vcf_index)

  # 1. Select for all pass calls with depth support
  # 2. Count the number of samples where the variant is found
  # 3. Output variant span and counts to bedgraph separately for DEL and DUP
  bcftools view -i 'FILTER="PASS" && SVCLAIM~"D"' $1 |\
  bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\t%INFO/SVTYPE\t%N_PASS(GT~"1")\n' |\
  awk -v OFS='\t' \
   -v del_out=${vcf_work_dir}/del_cnv_count.${zvi}.bedgraph \
   -v dup_out=${vcf_work_dir}/dup_cnv_count.${zvi}.bedgraph \
  '{
    if ($4 == "DEL") {
      print $1,$2,$2+$3,$5 > del_out;
    } else if ($4 == "DUP") {
      print $1,$2,$2+$3,$5 > dup_out;
    }
  }'
}

consolidate_cnv_bedgraphs() {
  bedtools unionbedg -i $@ |\
  awk -v OFS='\t' '{
    cnv=0;
    for (i = 4; i <= NF; i++) { cnv += $i }
    print $1,$2,$3,cnv;
  }' |\
  bedtools sort -i -
}

#
#
#

if [ $# -eq 0 ]; then
  cat << EOF
usage $0 [list/glob of sawfish cohort VCFs] > common_cnv.bed

Compute common CNV mask from sawfish VCFs representing a sample cohort. Input VCFs can be multi-sample, so long as samples are
not repeated in the input VCF set.

Assumes recent bcftools and bedtools are in the path, via the accompanying 'excluded_region_tools' conda environment
EOF
  exit 1
fi

# Create a temporary working file location:
work_dir=$(mktemp -d)
trap 'rm -rf $work_dir' EXIT

vcf_subdir=per_vcf
vcf_work_dir=${work_dir}/${vcf_subdir}
mkdir -p $vcf_work_dir

# Cycle through all input VCFs to make del/dup beds
sample_count=0
for arg in "$@"; do
  sawfish_vcf_to_cnv_bedgraphs $arg
  sample_count=$((sample_count + $(bcftools query -l $arg | wc -l)))
done

cd $work_dir

consolidate_cnv_bedgraphs ${vcf_subdir}/del_cnv_count.*.bedgraph >| del_cnv_count.bedgraph
consolidate_cnv_bedgraphs ${vcf_subdir}/dup_cnv_count.*.bedgraph >| dup_cnv_count.bedgraph

# Convert common CNV frequency into an exclusion bed track
cutoff_freq=0.5
min_count=$(echo "$cutoff_freq * $sample_count" | bc)
cat del_cnv_count.bedgraph | awk -v mc=$min_count '$4>=mc' | bedtools merge -i - | sed "s/$/\tDEL_ge_$cutoff_freq/" >| common_del.bed
cat dup_cnv_count.bedgraph | awk -v mc=$min_count '$4>=mc' | bedtools merge -i - | sed "s/$/\tDUP_ge_$cutoff_freq/" >| common_dup.bed

# Combine without merging to retain separate DEL/DUP labels, and write final common cnv bed to stdout
cat common_del.bed common_dup.bed | bedtools sort -i -
