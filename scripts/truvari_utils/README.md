# Truvari utils

`process_truvari_ga4gh_vcfs.py` - This script processes the VCF files output by the truvari 'ga4gh' command, to produce standard accuracy metrics (precision/recall/f1) for the given sample. Metrics are generated for all SVs, as well as different SV size and type categories. This script provides a useful way to summarize the total performance of all SVs from a combined truvari bench and refine analysis. See the sawfish [accuracy page](../../docs/accuracy.md) for more information about the GIAB T2T assessment pipeline where this is used.

