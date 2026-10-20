# sawfish change notes

This branch of `spoa-rs` includes some API customizations for the sawfish SV caller (https://github.com/PacificBiosciences/sawfish).

All edits are potential candidates for upstream if cleaned up a bit, but the branch is maintained primarily for quick API change trials.

Change summary:
- Update wrapper to include the score of each alignment
- Additional convenience functions:
    - `get_alignment_overlap_size`
    - `get_alignment_clip_size`
    - `add_alignment_noqual`
    - `get_sequence_count`

----

# spoa

A rust wrapper around the C++ [SPOA](https://github.com/rvaser/spoa) partial order alignment library.

Relevant DOIs of interest:
- SPOA C++ implementation: [10.1101/gr.214270.116](https://doi.org/10.1101/gr.214270.116)
- partial order alignment: [10.1093/bioinformatics/18.3.452](https://doi.org/10.1093/bioinformatics/18.3.452)
- consensus from partial order alignment: [10.1093/bioinformatics/btg109](https://doi.org/10.1093/bioinformatics/btg109)
