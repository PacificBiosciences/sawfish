# Rust bindings for WFA2-Lib

Rust language bindings for the excellent
[WFA2-Lib](https://github.com/smarco/WFA2-lib) library.

Work in progress. Tests and features are not yet complete.

## Autovectorization
Remember to specify the correct C compiler! For me it is `CC=/usr/local/opt/llvm/bin/clang`.

## Usage

```rust
use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, MemoryModel, WFAlignerGapAffine};

let mut aligner =
    WFAlignerGapAffine::new(6, 4, 2, AlignmentScope::Alignment, MemoryModel::MemoryLow);

let pattern = b"TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
let text = b"TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
let status = aligner.align_end_to_end(pattern, text);
assert_eq!(status, AlignmentStatus::StatusSuccessful);
assert_eq!(aligner.score(), -24);
assert_eq!(aligner.cigar(), "MMMXMMMMDMMMMMMMIMMMMMMMMMXMMMMMM");
let (a, b, c) = aligner.matching(pattern, text);
assert_eq!(
    format!("{}\n{}\n{}", a, b, c),
    "TCTTTACTCGCGCGTT-GGAGAAATACAATAGT\n|||||||| ||||||| ||||||||||||||||\nTCTATACT-GCGCGTTTGGAGAAATAAAATAGT"
);
```
