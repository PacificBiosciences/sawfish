# sawfish patches

This fork carries the changes sawfish needs on top of upstream. Branch `sawfish-performance`:

| Commit | Content |
|---|---|
| `9a978f9` (tag `sawfish_v1`) | Base: the `ctsa/spoa-rs` API changes for sawfish |
| `d8df2ee` | Flatten the `spoa-sys/spoa` submodule (`rvaser/spoa` at `981ad59`) into this repository. No content change. |
| `0540f91` | Performance patches (below) |

All patches return the same alignments and scores as the unpatched code. They only change time and memory.

## Patched files

- `spoa-sys/spoa/src/sisd_alignment_engine.cpp`, `.hpp`: low-memory linear-gap kSW/kOV alignment
  (score rows kept only until their successors are computed, tiles recomputed on backtrack);
  pruned alignment that skips cells that cannot reach a required score, with a fallback to the full
  alignment; AVX2 row kernels (disable with `SPOA_DISABLE_AVX2=1`); 16-bit scores for sequences that
  fit; score buffers returned to the OS when freed. `SPOA_VERIFY_LOW_MEMORY_LINEAR=1` checks the
  low-memory kernel against the full-matrix kernel.
- `spoa-sys/spoa/src/alignment_engine.cpp`, `spoa-sys/spoa/include/spoa/alignment_engine.hpp`:
  `AlignPruned`, `AlignPrunedAuto` and `MaxPossibleScore` interface, with full-alignment defaults.
- `spoa-sys/include/bindings.hpp`, `spoa-sys/src/bindings.cpp`, `spoa-sys/src/lib.rs`, `src/lib.rs`:
  expose the above to Rust.

## Updating from upstream

Diff `spoa-sys/spoa` against the new `rvaser/spoa` revision, re-apply the changes to the files above,
and run the sawfish `test_data` comparison with `SPOA_VERIFY_LOW_MEMORY_LINEAR=1`.
