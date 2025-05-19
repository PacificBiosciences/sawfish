#!/usr/bin/env bash

set -o errexit
set -o nounset

# This script builds excluded regions for CNV calling based on UCSC annotation tracks for a given reference genome. It
# is recommended to supplement these regions with empirical regions based on common CNV calls in a background cohort.
#
# The reference defaults to human hg38, but can be adapted for other references depending on UCSC track availability.
#
# The script assumes tabix and bgzip are in the path, via the accompanying 'excluded_region_tools' conda environment.
#

# The reference tag for the genome version to use. This script has only been tested for ref values 'hg19' and 'hg38'. It
# may work with other reference versions in the UCSC genome browser database, but this hasn't been tested.
ref=hg38


outfile=annotation_only.${ref}.bed.gz

base_url=http://hgdownload.cse.ucsc.edu/goldenPath/$ref

# Get UCSC Gaps track, simplify labels and convert to bed format:
get_gaps() {
  wget -O - $base_url/database/gap.txt.gz |\
  gzip -dc |\
  awk -v OFS='\t' '{
    if ($8 ~ /^(clone|contig|scaffold|short_arm)$/) {
       $8 = "gap";
    }
    print $2,$3,$4,$8;
  }'
}

# Get UCSC Centromere track if it exists, simplify labels and convert to bed format
#
# This file is optional because in at least some older genomes, it may not exist. In such cases the centromere regions
# are annotated in the gap track instead.
#
get_centromeres() {
  url=$base_url/database/centromeres.txt.gz
  if wget --quiet --spider $url; then
    wget -O - $url |\
    gzip -dc |\
    awk -v OFS='\t' '{print $2,$3,$4,"centromere";}'
  fi
}

# Get alpha-satellite regions from the UCSC repeatmasker track
get_alpha() {
  wget -O - $base_url/database/rmsk.txt.gz |\
  gzip -dc |\
  awk -v OFS='\t' '$11=="ALR/Alpha" {
    print $6,$7,$8,$11;
  }'
}

cat <(get_gaps) <(get_centromeres) |\
cat - <(get_alpha) |\
sort -k1,1 -k2,2g | bgzip -c >|\
$outfile

tabix -p bed $outfile
