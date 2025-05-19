#!/usr/bin/env bash

set -o errexit
set -o nounset

# This script yields a recommended set of exclusion regions for CNV calling on hs37d5, by
# converting excluded regions from hg19.
#
# Assumes tabix and bgzip are in the path, via the accompanying 'excluded_region_tools' conda environment
#

hg19_excluded_regions=annotation_only.hg19.bed.gz

hg19_renamed() {
  gzip -dc $hg19_excluded_regions |\
  sed s/^chr// |\
  awk '$1~/^([1-2]?[0-9]|[XY]|MT)$/'
}

other() {
  wget -O - http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai |\
  awk '$1!~/^([1-2]?[0-9]|[XY]|MT)$/ {printf "%s\t0\t%s\tother\n",$1,$2}'
}


label=annotation_only.hs37d5

cat <(hg19_renamed) <(other) | bgzip -c >| $label.bed.gz
tabix -p bed $label.bed.gz

