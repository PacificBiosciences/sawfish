<h1 align="center"><img width="250px" src="img/logo.svg"/></h1>

<h1 align="center">Sawfish</h1>

<h4 align="center">Structural variant caller for HiFi sequencing data based on local haplotype assembly</h3>

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

## Getting started

See the [getting started](docs/user_guide.md#getting-started) section in the [user guide](docs/user_guide.md) to start using sawfish.

## Citation

Sawfish is detailed in the following pre-print:

[Saunders, C. T., Holt, J. M., Baker, D. N., Lake, J. A., Belyeu, J. R. Kronenberg, Z., Rowell, W. J., & Eberle, M. (2024). Sawfish: Improving long-read structural variant discovery and genotyping with local haplotype modeling. *bioRxiv*](https://doi.org/10.1101/2024.08.19.608674)

## Support

Create a new [issue ticket](https://github.com/PacificBiosciences/sawfish/issues) on this repo for support with current capabilities or new feature requests.

Note that this is a research method in early development. The interface and behaviors are still being refined and may change in subsequent feature releases.

## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
