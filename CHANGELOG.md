# Change Log

## v2.0.2 - 2025-07-11

### Added
- Add support for cram file inputs which include bzip2 and lzma codec blocks [PacificBiosciences/sawfish#19]

## v2.0.1 - 2025-06-18

### Added
- CR-497 Add GC-bias corrected depth track to the per-sample output from joint-call

## v2.0.0 - 2025-05-19

This is the first official release of SV/CNV integration in sawfish. This version adds full depth segmentation into the
SV calling process to provide several benefits:

1. Consistent view of all large variants in a single output, including calls based on breakpoint evidence, calls from
depth evidence and calls jointly supported by both.

2. A breakpoint-enhanced CNV caller: Sawfish provides standard depth-based CNV calls without requiring breakpoint
evidence, but when available, breakpoints can be used to provide more accurate segmentation and higher sensitivity at
lower sequencing depths.

3. An improved SV caller, via more precise calling of large, unbalanced SVs: All large deletions and duplications must
have some level of depth segmentation support, enabling more precise large variant output.

## v1.0.1 - 2025-04-13

### Fixed
- CR-468 Fix high memory usage when joint-genotyping from CRAM [PacificBiosciences/sawfish#15]
  - This issue only impacts analysis from CRAM, all BAM-input analysis is unaffected.
  - The excess memory usage had only a minor effect on single-sample analysis, but could become very large at higher sample and thread counts.
- CR-469 Fix very low-frequency non-deterministic genotyping result
  - This also fixes non-deterministic read order in the optional supporting reads output

## v1.0.0 - 2025-04-10
This release marks the initial stabilization of sawfish with its publication in [Bioinformatics](https://doi.org/10.1093/bioinformatics/btaf136).
It is functionally identical to v0.12.10 but includes updated in-source documentation for both accuracy
assessment and methods.

### Added
- CR-467 Add in-source pandoc sawfish methods document.

### Changed
- CR-462 Updated accuracy assessment page from the original sawfish release to match all updated methods
reported in the sawfish paper.

## v0.12.10 - 2025-02-24

### Fixed
- Fix handling of chromosome names with colons, eg. 'HLA-DRB1*10:01:01' [PacificBiosciences/sawfish#11]

## v0.12.9 - 2025-01-10

### Added
- CR-418 Provide BAM output of assembled SV contig alignments in the joint-call step
  - BAM output reflects all assembled SV haplotypes used during genotyping, which can be useful for reviewing SV calls

## v0.12.8 - 2024-12-13

### Fixed
- Increase system open file limit [PacificBiosciences/sawfish#9]
  - May simplify joint-call for large pedigrees at high thread counts
- Improve error message when split reads map to an unknown chromosome [PacificBiosciences/sawfish#8]

## v0.12.7 - 2024-10-23

### Fixed
- Fix "Illegal Instruction" error reported for some use cases [PacificBiosciences/sawfish#3]
  - Removed gcc native cpu optimization in WFA2-lib believed to be causing this issue
  - Added additional `--debug` output to discover mode contig alignment logic

## v0.12.6 - 2024-10-15

### Fixed
- Fix discover mode input path canonicalization (reversed flag logic)

## v0.12.5 - 2024-10-14

### Added
- Add new joint-call `--report-supporting-reads` option to report read names supporting each variant

### Changed
- CR-390 Canonicalize all discover mode input paths
  - Also provide new `--disable-path-canonicalization` discover step option to store input paths as-is
- CR-391 Don't create output directory until command line is validated

## v0.12.4 - 2024-09-13

### Added
- CR-384 Add debug logging option
  - High detail level intended to improve crash reports from external users
  - Initially populated for breakpoint refinement only, debug log coverage will be expanded as required
- Expose clobber option to overwrite existing output directory

### Fixed
- CR-385 Improve error message for unexpected alignment patterns from VACmap
  - This also adjusts some off-by-one errors in neighbor extension handling, which could cause infrequent changes to inversion output
  - Note VACmap is still unsupported; working towards clear error messages for problematic alignments from any source

## v0.12.3 - 2024-09-03

### Fixed
- CR-378 Improve error message for hard-clipped split read input
- CR-379 Improve error message when discover directory is missing

## v0.12.2 - 2024-08-21

### Fixed
- CR-377 Remove `rq` tag requirement in input alignment records

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
