"""Compare amino-acid group motif frequencies across taxonomy analyses."""
from __future__ import annotations

import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent
PROJECT_ROOT = BASE.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from libraries.compare_aa_group_motifs_lib import run_compare


# ---- Paths (edit here) ----
# Notes:
# - To add another taxonomy, append another dict to COMPARISONS.
# - "outputs_root" should point to the <order>/ncbi_fibroin_sequences folder that
#   contains data_analysis/aa_motifs/...
# - "taxonomy_terms" must match the taxonomy folder name under aa_motifs
#   (case-sensitive on some systems).
# - "label" is just the display name used in tables/plots.
# - Optional per-entry overrides:
#   * "protein_types": list like ["heavy chain"] (overrides PROTEIN_TYPES for that entry)
#   * "partial_full", "length_range", "length_threshold", "length_mode",
#     "longest_factor", "min_motif_length", "max_motif_length" (override filters per entry)
COMPARISONS = [
  {
        "label": "trichoptera",
        "outputs_root": Path(r"bioinformatic_silk\caddisfly\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Trichoptera"],
    },
   {
        "label": "Integripalpia",
        "outputs_root": Path(r"bioinformatic_silk\caddisfly\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Integripalpia"],
    },
    {
        "label": "Annulipalpia",
        "outputs_root": Path(r"bioinformatic_silk\caddisfly\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Annulipalpia"],
    },
]



DEFAULT_OUT_DIR = Path(r"bioinformatic_silk\caddisfly\ncbi_fibroin_sequences\data_analysis\aa_motifs_comparison")
PLOTS_OUT_DIR = Path(r"bioinformatic_silk\caddisfly\ncbi_fibroin_sequences\plots\aa_motifs_comparison")
# ---- End paths ----

# ---- Filters (edit after paths) ----
PARTIAL_FULL = None  # "partial", "full", or None for both
PROTEIN_TYPES = ["heavy chain"]  # restrict to specific protein types
LENGTH_RANGE: tuple[int, int] | None = None  # exact length range filter (min, max)
LENGTH_THRESHOLD = 1000  # threshold length value
LENGTH_MODE = "ge"  # "ge" (>=) or "le" (<=) for threshold filter
LONGEST_FACTOR = None  # optional factor for longest selection logic
MIN_MOTIF_LENGTH = 5  # shortest motif length to include
MAX_MOTIF_LENGTH = 50  # longest motif length to include
MIN_AVG_RELATIVE = 1  # minimum average relative percent to keep a motif in the comparison
USE_UNIQUE_MOTIFS = False  # use unique motifs instead of overlapping motif counts
# Window normalization uses: count / (L - k + 1) * 100, which is <= 100%.
USE_WINDOW_NORMALIZED = False  # normalize counts by window count per sequence

# Filter logic:
# - If MIN_AVG_RELATIVE_TAXA is set: keep motifs where >= MIN_AVG_RELATIVE_MIN_TAXA taxonomies meet the threshold
#   AND at least one of the listed taxa meets the threshold.
# - If MIN_AVG_RELATIVE_TAXA is empty: keep motifs where >= MIN_AVG_RELATIVE_MIN_TAXA taxonomies meet the threshold.
MIN_AVG_RELATIVE_TAXA =  ["Integripalpia", "Annulipalpia"]  # taxa that must meet MIN_AVG_RELATIVE
MIN_AVG_RELATIVE_MIN_TAXA = 2  # number of taxa that must meet MIN_AVG_RELATIVE
REQUIRE_NON_FILTER_TAXA = True  # require at least one non-filter taxon to meet MIN_AVG_RELATIVE
# ---- End filters ----


def main() -> int:
    filters = {
        "partial_full": PARTIAL_FULL,
        "protein_types": PROTEIN_TYPES,
        "length_range": LENGTH_RANGE,
        "length_threshold": LENGTH_THRESHOLD,
        "length_mode": LENGTH_MODE,
        "longest_factor": LONGEST_FACTOR,
        "min_motif_length": MIN_MOTIF_LENGTH,
        "max_motif_length": MAX_MOTIF_LENGTH,
        "min_avg_relative": MIN_AVG_RELATIVE,
        "use_unique_motifs": USE_UNIQUE_MOTIFS,
        "use_window_normalized": USE_WINDOW_NORMALIZED,
        "min_avg_relative_taxa": MIN_AVG_RELATIVE_TAXA,
        "min_avg_relative_min_taxa": MIN_AVG_RELATIVE_MIN_TAXA,
        "require_non_filter_taxa": REQUIRE_NON_FILTER_TAXA,
    }
    return run_compare(COMPARISONS, DEFAULT_OUT_DIR, PLOTS_OUT_DIR, filters)


if __name__ == "__main__":
    raise SystemExit(main())
