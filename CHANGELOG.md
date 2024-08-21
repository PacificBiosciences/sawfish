# Change Log

## v0.12.2 - 2024-08-21

### Fixed
-CR-377 Remove `rq` tag requirement in input alignment records

## v0.12.1 - 2024-08-08

### Fixed
- CR-375 Fix infrequent discover mode failures due to invalid breakpoint ranges.
  - Issue seems to have started with v0.12.0 via CR-340

## v0.12.0 - 2024-05-23

### Added
- CR-327 Add local SV phasing
  - Short-range phasing provided in VCF records wherever multiple hets are genotyped on one or more overlapping SV haplotypes

### Changed
- CR-340 Improve handling of multi-breakpoint haplotypes
  - Improves detection of complex SV breakpoints, especially for small inversions
- CR-333 Adjust alignment for long breakpoint homology
  - Improves detection of high homology inversions
- CR-330 Expand trimmed read search region for small SV regions
  - Improves detection of soft-clipped evidence for low homology insertions

## v0.11.0 - 2024-04-18

Improve a number of discovery and scoring features related to inversions, inverted breakpoints and large deletions.

### Changed
- CR-323 Improve scoring for a number of cases relevant to inversions and inverted breakpoints
- CR-321 Change contig flank size handling to improve large inversion calling
- CR-318 Expand trimmed read search region for large SV candidates
  - Improves recall for inverted breakpoints and duplications
- CR-317 Standardize on csi indexing for bam output

### Fixed
- CR-324 Fix non-deterministic inversion output

## v0.10.0 - 2024-03-05

Initial github release
