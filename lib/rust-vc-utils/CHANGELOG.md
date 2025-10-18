## v0.33.0 - 2025-XX-XX

### Added
- Add RegionMap
- Add drop_true from sawfish
- New methods to remove clipped regions from cigars

### Changed

- In all cigar shift routines, reverse cigar indel cluster ordering such that insertion state is always ordered first
- Add method to restrict GenomeRef sequence characters
- Make reverse comp methods tolerant to unexpected characters
  - Still not supporting IUPAC ambiguity codes

## v0.32.0 - 2025-08-07

### Added
- Add alignment edit distance and reorg gap-compressed identity
- Add alignment clipping routines from sawfish
- Add indel left and right-shift methods from portello
- Add sawfish `get_indel_breakend_homology_info` method

## v0.31.0 - 2025-07-09

This update includes numerous additions to support the initial release of portello.

### Added
- Add `assert_bam_eof` method to assert that bam/cram files are not truncated
- Add `remove_aux_if_found` method for removing aux records from bam that may not exist
- Add `clean_up_cigar_edge_indels` method to standardize all edge indels in a cigar string
- Add `cigar_edge_insertion_to_softclip` method to standardize edge insertions in a cigar string
- Add subset of sawfish GenomeSegment methods
- Add `get_read_segment_to_ref_pos_map` method to generate read-to-ref position maps from split reads.
- Add sawfish int range methods
- Add sawfish split read parse and ordering methods
  - Update to support degenerate case

### Changed
- Breaking interface change -- all cigar functions processing read lengths have changed names and arguments.
  - Cigar read-length functions have previously been supported as two functions, one to include hard-clips and one to
    ignore hard-clips.
  - All methods are consolidated to one function taking an `ignore_hard_clip` bool argument instead.
  - Update guide:
    - Replace any "complete" methods with `ignore_hard_clip` set to `false`
```
let offset = get_cigar_complete_read_offset(cigar); // Before
let offset = get_cigar_read_offset(cigar, false); // After
```
    - Replace any "hard_clipped" methods with `ignore_hard_clip` set to `true`
```
let offset = get_cigar_hard_clipped_read_offset(cigar); // Before
let offset = get_cigar_read_offset(cigar, true); // After
```

## v0.30.0 - 2025-06-09

### Changed
- Add immutable revcomp and add lower-case to revcomp to fit vcf DNA string spec.

## v0.29.0 - 2024-10-31

### Added
- Add simple mean tracker
- Add new force-periodic option to `ProgressReporter`, breaks api
- Add new bam record utils

## v0.28.0 - 2024-09-04

### Changed
- Add changes to bam aux tag parsing
    - CR-380 Reduce the chance of a 'panic-in-panic' leading to a SIGILL
    - Add new parser for float aux tags

## v0.27.0 - 2024-08-29

### Changed
- Changed to github version of rust-libbigwig

## v0.26.0 - 2024-03-21

### Fixed
- Fix fasta reader error in formatted message

## v0.25.0 - 2024-03-15

### Added
- Add progress reporter

## v0.24.1 - 2023-10-05

### Changed
- Updated build to remove unused rust-htslib dependencies

## v0.24.0 - 2023-08-21

### Added
- Add new sparse window sum container

### Fixed
- Fix basemod parsing error messages

## v0.23.0 - 2023-07-26

### Added
- Add new cigar methods

### Changed
- Changed Cargo.toml to more flexibly allow dependency changes from client code, principally so that this library no
  longer imposes such an exact rust-htslib requirement.
- Multiple API breaking changes:
  - Refactored parameterless new methods into default
  - Refactored bam_util submodule names

## v0.22.0 - 2023-06-21

### Added
- Add new cigar convenience functions

### Changed
- Updated to rust-htslib 0.44.1

## v0.21.0 - 2023-06-06

### Added
- Added deterministic vector downsample
- Added new cigar processing methods.
    - Hard-clip test
    - Individual cigar segment offset functions

### Changed
- Updated cigar interfaces
  - Used `hard_clip` more consistently in all functions with an underscore. This update **intentionally breaks api**.
  - Accept `&[Cigar]` instead of `CigarString` to generalize input cigar requirements to more data types.

## v0.20.0 - 2023-05-25

### Added
- Added several new cigar processing methods.
  - This update **intentionally breaks api** for all cigar processing functions where hard-clip handling was previously
    unspecified, and any downstream bam processing functions in the library which used the previously unspecified
    behavior. In most of these cases they have been replaced with hard-clip and non-hard-clip versions. Updating code will
    have to explicitly consider which version to update to.

## v0.19.0 - 2023-05-25

### Added
- Added new cigar clipping method

### Fixed
- Reorganized internal bam utility structure

## v0.18.0 - 2023-05-25

### Added
- Added bam aux string method
- Added bam reg2bin

## v0.17.0 - 2023-05-18

### Added
- New ChromList ctor and test

## v0.16.0 - 2023-05-10

### Added
- Added new aux parse options

## v0.15.0 - 2023-05-01

### Added
- Added bam record read to ref map

## v0.14.0 - 2023-04-13

### Added
- Added new bam aux tag utils
- Updated to 2021 rest edition

## v0.13.0 - 2023-04-07

### Added
- rev comp

## v0.12.0 - 2023-04-07

### Added
- New bam cigar tracking util

## v0.11.0 - 2023-04-04

### Fixed
- Fixed bug in `decode_cpg_meth_info` that would produce an invalid read_pos for a final C in the read with methylation when the read was reverse-mapped.

