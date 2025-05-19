<h1 align="center"><img width="250px" src="img/logo.svg"/></h1>

<h1 align="center">Sawfish</h1>

<h4 align="center">Joint structural variant and copy number variant caller for HiFi sequencing data</h3>

Sawfish is a joint structural variant (SV) and copy number variant (CNV) caller for mapped HiFi sequencing reads. It discovers germline structural variants from local sequence assembly and jointly genotypes these variants across multiple samples. Sawfish additionally applies copy number segmentation on each sample's sequencing coverage levels, synchronizing structural variant breakpoints with copy number change boundaries in the process to improve classification of breakpoint-based structural variant calls, in addition to calling copy number variants.

Key features:
- Combined assessment of all large variants in each sample.
  - Sawfish provides a unified view of both SVs and CNVs in each sample, with redundant calls merged into single variants describing both breakpoint and copy number detail.
- High SV discovery and genotyping accuracy
  - All breakpoint-based structural variants are modeled and genotyped as local haplotypes, yielding substantial accuracy gains on modern SV truth sets such as the GIAB HG002 T2T SVs.
- High resolution
  - All breakpoint-based structural variants are assembled to basepair resolution and reported with breakpoint homology and insertion details.
- Integrated copy number segmentation
  - Integrated copy number segmentation with GC-bias correction is used to: (1) call CNVs independent of any breakpoint support, and (2) improve the classification of large structural variant deletion and duplication calls, any such calls lacking consistent depth support are reclassified as breakends.
- Simple multi-threaded workflow
  - A single command-line is used for each of the discover and joint-call steps

Breakpoint-based SVs are reported as deletions, insertions, duplications and inversions when supported by the corresponding breakpoint and depth pattern, otherwise the breakpoint itself is reported. Copy number variants are reported as deletions and duplications. The minimum variant size is 35 bases (configurable). A maximum size is only applied to inversions (100kb).

## Getting started

See the [getting started](docs/user_guide.md#getting-started) section in the [user guide](docs/user_guide.md) to start using sawfish.

## Citation

Sawfish is published in *Bioinformatics*, this can be cited as:

[Saunders, C. T., Holt, J. M., Baker, D. N., Lake, J. A., Belyeu, J. R. Kronenberg, Z., Rowell, W. J., & Eberle, M. (2025). Sawfish: Improving long-read structural variant discovery and genotyping with local haplotype modeling. *Bioinformatics*](https://doi.org/10.1093/bioinformatics/btaf136)

## Support

Create a new [issue ticket](https://github.com/PacificBiosciences/sawfish/issues) on this repo for support with current capabilities or new feature requests.

## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
