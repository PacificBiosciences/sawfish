# Change Log

## v0.10 - 2024-03-05

Changed:
- CR-302 Change zero-depth sample output
  -  When a sample has no supporting reads for any allele, define GT and GQ as unknown in the output VCF
- CR-298 Improve haplotype merge
  - Remove filters preventing unlikely haplotype merge testing. Improved merging removes ~25% of FP trio denovo calls
  - Parallelize haplotype merge testing to compensate for the extra tests resulting from the above change, reducing total joint-call time for typical thread counts
- Change VCF output max Q-Score from 99 to 999

Fixed:
- Improve error messages for missing, unmapped and unindexed alignment file inputs
- CR-299 Fix non-deterministic haplotype selection due to unstable floating point comparison

## v0.9 - 2024-02-02

Added:
- CR-273 Add inversion formatting to VCF output
- CR-267 Add tandem duplication formatting to VCF output
- CR-265 Add periodic status update for long running tasks in non-tty contexts

Changed:
- CR-266 Reduce peak memory demand for large insertions
   - Stabilizes memory demand to no more than 8Gb/thread in all human test sets available.
- CR-269 Report very large deletions as symbolic alleles

Fixed:
- Fix inversion output to be deterministic

## v0.8 - 2023-12-06

This update substantially improves genotype concordance on family-scale pedigrees.

Added:
- CR-262 Add support for per-sample haploid genotyping regions
- Add parallel sample data processing during joint genotyping

Changed:
- CR-263 Improve discover runtime for high depth samples

Fixed:
- CR-259 Adjust overlapping contig registration during allele scoring
- CR-261 Fix breakend homology range handling in allele scoring
- CR-260 Generalize replicate SV handling for multi-sample
- CR-258 Improve multi-sample handling of allele candidates from other samples.
  - This substantially improves genotype consistency across larger pedigrees.
- Fixed blocked allele ref AD counts

## v0.7 - 2023-11-20

This update adds an initial version of multi-sample joint genotyping.

Added:
- CR-252 Enable multi-sample calling at trio scale, with partial accuracy optimization.
   - This covers most merging and multi-sample genotyping functionality, but accuracy is not fully optimized to the 
   single sample workflow level.
- CR-250 Break single-sample SV calling workflow into two separate workflows for candidate discovery and genotyping

Fixed:
- CR-257 Fix read allele alignment issue causing multisample false positives
- CR-255 Fix multi-allelic insertions occurring in multisample context

## v0.6.1 - 2023-10-23

Added:
- CR-249 Classify breakpoints consistent with large deletions based on observed change in sequencing depth.

## v0.6 - 2023-10-04

This update focuses on substantial improvements to genotype accuracy.

Added:
- Add progress bar to track breakpoint cluster refinement
- CR-243 Automate BAM and VCF output indexing

Fixed:
- CR-247 Remove large insertion duplicates
- CR-246 Improve genotype overcall pattern occurring on certain deletions
- CR-245 Fix spurious ref AD counts attributed to undercalling GT of homozygous insertions
- Fix read offset error in scoring routine
- CR-244 Fix AD counts for spurious overlapping insertion pattern
- CR-242 Fix AD counts in overlapping deletion regions

## v0.5 - 2023-09-22

Added:
- CR-235 Fix AD counts in complex haplotypes
- Add supporting read counts to contig debug bam
- CR-234 Add diploid quality scores
- CR-209 Add AD counts to all SV output

Changed:
- Adjust assembly tier2 parameters for noisy contexts
- Adjust large insertion trigger to allow for broader breakpoint arrangement
- Add soft-clip reads to standard insertion assemblies

## v0.4 - 2023-08-11

Added:
- Add second backup cluster stage with higher noise tolerance when initial clustering fails
- Add local check backup alignment for POA assembly to reduce false positive contig assemblies

Changed:
- Remove a false positive breakend cluster pattern triggered by noisy reads
- Remove duplicated large insertion assemblies where insertions can also be called by standard assembly to improve precision
- Reorder large insertion jobs and limit their inputs to reduce memory demand and improve thread utilization
- Adjust contig to ref alignment to increase SV clustering
- Remove indel calls made from fake flank regions to improve precision.
- Relax min reported indel size and min evidence indel size
- Relax criteria for 'fake flank' read selection to improve recall of a few difficult cases
- Adjust contig flank sizes and wfa2 clipping to improve recall at high homology breakends
- Filter non-variant reads from local indel assemblies to improve recall

Fixed:
- Fix bug in run stats output
- Fix minor bug in 'fake flank' read selection and improve other selection criteria, improving recall.

## v0.3 - 2023-07-26

Added
- Add option to reduce overlapping SV alleles so that comparisons can be made to simplified SV truth sets.
- Add new `*.run_stats.json` output to track various run metrics.

Changed
- Adjust contig-to-alt alignment parameters to increase convex gap open cost, improving GIAB precision and recall.
- Reverse contig-to-alt alignments to left-shift raw alignment results, creating more conventional representation of
  several SVs.

Fixed
- Fix incorrect breakpoint insertion sequences in all BND records.
- Fix splitting of breakpoint alignments to account for breakpoint insertions on either breakend.
- Fix multiple correctness and stability issues with BND records:
  - Fix start position of all right-anchored BND records to be one base higher, corresponding to the base after the
  break. This representation matches the spec.
  - Fix breakend homology sequence handling to sync both sides of all BND records and prevent infrequent out of range
  panic at end of chromosome sequence.
  - Fix breakend position and homology sequence for second BND record of all inverted breakpoints.
- Fix large deletion false positives due to flawed candidate read search
- Make SV VCF output deterministic

## v0.2 - 2023-06-27

Recall improvement update

- Added numerous features focused on improving recall of the SV candidate set prior to scoring. Changes include:
  - Added a flanking sequence structure to properly align variants with high breakend homology
  - Improve alignment scoring and assembly details
  - Added large insertion assembly. Largest called insertions in HG002 no go up to ~28kb
- Added left-shift and exact breakend homology annotations to all calls (using HOMLEN/HOMSEQ)

## v0.1 - 2023-06-13

Initial sawfish release. Current capabilities:

- VCF output covering all SV types, but mostly expressed in BND format.
  - SV VCF output is currently unscored and not integrated with CNV.
  - No breakend homology/left-shift

- VCF output currently provides all indels in DEL/INS format, and all other types in BND format
