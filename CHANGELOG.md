# Change Log

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
