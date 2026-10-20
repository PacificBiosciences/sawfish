# sawfish patches

This fork carries the changes sawfish needs on top of upstream. Branch `sawfish-performance`:

| Commit | Content |
|---|---|
| `4342b3b` | Base: `ctsa/rust-wfa2` |
| `ce25305` | Flatten the `wfa2-sys/WFA2-lib` submodule (`ctsa/WFA2-lib` at `700a93e`) into this repository. No content change. |
| `e8ddc72` | Performance patches (below) |

Alignment results are identical to unpatched WFA2. The patches only change time and memory.

## Patched files

- `src/affine2p.rs` (new): exact low-memory reimplementation of WFA2 gap-affine 2-piece ends-free
  alignment, with a checkpointed backtrace for very long alignments. It follows WFA2 as of
  `4342b3b`, including tie-breaking, trimming, termination order and backtrace, so a WFA2 update
  can change which of several equally good alignments WFA2 returns, and this must be checked.
- `src/linear_dp.rs` (new): exact gap-linear ends-free score by anti-diagonal DP, for scoring
  dissimilar read flanks without WFA2.
- `src/aligner.rs`, `src/lib.rs`, `wfa2-sys/src/lib.rs`: expose the above and the memory budget.
- `wfa2-sys/WFA2-lib/system/mm_allocator.c`, `.h`: process-wide memory budget shared by all
  threads, with an oldest-claim rule so that some thread can always make progress. Disabled
  (zero) by default. It can deadlock if a long-lived aligner holds budget memory across
  waiting calls, so only enable it where aligners are short-lived.
- `wfa2-sys/WFA2-lib/utils/commons.h`, `wavefront/wavefront_compute_affine2p.c`,
  `wavefront/wavefront_extend_kernels.c`: `target_clones` for AVX2 on the affine2p compute and match
  extension kernels, selected at runtime.

## Updating from upstream

Diff `wfa2-sys/WFA2-lib` against the new upstream revision, re-apply the changes above, and run the
sawfish `test_data` comparison with `SAWFISH_VERIFY_FAST_WFA=1`.
