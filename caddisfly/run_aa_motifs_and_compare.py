"""Run aa_motifs_analysis for multiple taxonomies, then compare their motifs."""
from __future__ import annotations

import sys
import time
import json
from pathlib import Path
from typing import Dict, Iterable, Tuple

BASE = Path(__file__).resolve().parent
PROJECT_ROOT = BASE.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from libraries.amino_acid_group_motifs import (
    DEFAULT_GROUPS,
    analyze_group_motifs,
    dedupe_records_by_accession,
    resolve_aa_motifs_output_dir,
    select_longest_records,
    write_json,
    write_md,
    write_reflected_json,
    write_reflected_md,
)
from libraries.aminoacids_composition_analysis_lib import filter_records, load_records_from_root
from libraries.compare_aa_group_motifs_lib import run_compare


# ---- Groups (edit here) ----
# If GROUPS is None, defaults are used.

# GROUPS: Dict[str, Iterable[str]] | None = None     #   unnote it to have default groups

# Optional group legend override (used in compare plots/markdown).
GROUPS = {
     "π": ["A", "V", "L", "I", "M", "P"],
     "α": ["F", "W", "Y"],
     "θ": ["C", "N", "Q", "G","S", "T"],
     "ψ": ["D", "E"],
     "φ": ["K", "R", "H"],}
# Example:

# GROUPS = {
#     "π": ["A", "V", "L", "I", "M", "P"],
#     "α": ["F", "W", "Y"],
#     "η": ["S", "T"],
#     "θ": ["C", "N", "Q", "G"],
#     "ψ": ["D", "E"],
#     "φ": ["K", "R", "H"],
# }

# GROUPS = {
#     "H": ["A", "V", "L", "I", "M", "P", "F", "W", "Y"],
#     "P": ["S", "T", "C", "N", "Q", "G", "D", "E", "K", "R", "H"],
# }

# GROUP_LEGEND: Dict[str, str] | None = None          #  unnote it to have no legend override
# Optional group label to create an extra nesting folder (e.g., "default" or "nonpolar_polar_aromatics").

GROUP_LEGEND = {
     "H": "Hydrophobic",
     "α": "aromatic",
     "θ": "polar",
     "q": "charged",
 }

# Example:
# GROUP_LEGEND = {
#     "π": "nonpolar",
#     "α": "aromatic",
#     "η": "hydroxyl",
#     "θ": "polar",
#     "ψ": "negative charged",
#     "φ": "positive charged",
# }

#GROUP_LEGEND = {
#     "H": "hydrophobic + aromatic",
#     "P": "polar/charged",
# }

GROUP_LABEL ="Hydrophobic_aromatic_polar_charged"
# GROUP_LABEL = "default"           # unnote it to have default label
# GROUP_LABEL = "H_P"


# ---- End groups ----

# ---- Sync compare config (optional) ----
# If True, reuse settings from caddisfly/compare_aa_group_motifs.py so results match.
USE_COMPARE_CONFIG = False
# If True, skip recomputing aa_motifs outputs and only run the comparison.
REUSE_EXISTING_OUTPUTS = False
# ---- End sync ----

# ---- Analysis targets (edit here) ----
ANALYSES = [

    {
        "label": "Integripalpia",
        "outputs_root": Path(r"D:\python\Lab\caddisfly\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Integripalpia"],
    },
    {
        "label": "Annulipalpia",
        "outputs_root": Path(r"D:\python\Lab\caddisfly\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Annulipalpia"],
    },
    {
        "label": "Bombyx",
        "outputs_root": Path(r"D:\python\LHuge_data_for_bioinformatic_project\moths_and_butterflies\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Bombyx"],
    },
    {
        "label": "saturniidae",
        "outputs_root": Path(r"D:\python\LHuge_data_for_bioinformatic_project\moths_and_butterflies\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["saturniidae"],
    },
    {
        "label": "Nymphulinae",
        "outputs_root": Path(r"D:\python\LHuge_data_for_bioinformatic_project\moths_and_butterflies\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Nymphulinae"],
    },
    {
        "label": "Eumeta japonica",
        "outputs_root": Path(r"D:\python\LHuge_data_for_bioinformatic_project\moths_and_butterflies\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Eumeta japonica"],
    },
    {
        "label": "Chironomoidea",
        "outputs_root": Path(r"D:\python\Lab\Diptera\ncbi_fibroin_sequences"),
        "taxonomy_terms": ["Chironomoidea"],
        "length_threshold": None,
    },
]






COMPARE_OUT_DIR = Path(r"D:\python\Lab\caddisfly\ncbi_fibroin_sequences\data_analysis\aa_motifs_comparison")
COMPARE_PLOTS_DIR = Path(r"D:\python\Lab\caddisfly\ncbi_fibroin_sequences\plots\aa_motifs_comparison")
# ---- End analysis targets ----

# ---- Filters (edit after targets) ----
PARTIAL_FULL = None  # "partial", "full", or None for both
PROTEIN_TYPES = ["heavy chain"]  # restrict to specific protein types
LENGTH_RANGE: Tuple[int, int] | None = None  # exact length range filter (min, max)
LENGTH_THRESHOLD = 1000  # threshold length value
LENGTH_MODE = "ge"  # "ge" (>=) or "le" (<=) for threshold filter
LONGEST_FACTOR = None  # optional factor for longest selection logic
LONGEST_FACTOR_SCOPE = "species"  # scope for longest selection
USE_LONGEST_ONLY = False  # keep only the longest record per species if True

# Output controls.
INCLUDE_MOTIF_STRING = False  # include full motif strings per record (large output)

# Motif scan settings.
INCLUDE_MOTIF_COUNTS = True  # include motif counts/frequencies
MIN_MOTIF_LENGTH = 5  # shortest motif length to include
MAX_MOTIF_LENGTH = 50  # longest motif length to include
SKIP_UNKNOWN_MOTIFS = True  # skip motifs that contain "unknown"
INCLUDE_UNIQUE_MOTIFS = True  # include unique motif counts/frequencies
UNIQUE_MOTIF_LENGTH = 5  # fixed length for unique motifs
SKIP_COMBINED_MOTIFS = True  # drop motifs that are repeats/concats
INCLUDE_REFLECTED_SEQUENCE = True  # include reflected sequence output
INCLUDE_SIMILAR_MOTIFS = True  # include per-species shared motifs
INCLUDE_TAXONOMY_SHARED_MOTIFS = True  # include taxonomy shared motifs
PHYLO_TREE_PATH = None  # optional phylo tree path override

LONGEST_SCOPE = "species"

# Compare filters.
MIN_AVG_RELATIVE = 3  # minimum average relative percent to keep a motif
USE_UNIQUE_MOTIFS = False  # use unique motifs instead of overlapping motif counts
# Window normalization uses: count / (L - k + 1) * 100, which is <= 100%.
USE_WINDOW_NORMALIZED = False  # normalize counts by window count per sequence
MIN_AVG_RELATIVE_TAXA: list[str] = ["Integripalpia", "Annulipalpia"]  # taxa that must meet MIN_AVG_RELATIVE
MIN_AVG_RELATIVE_MIN_TAXA = 2  # number of taxa that must meet MIN_AVG_RELATIVE
REQUIRE_NON_FILTER_TAXA = True  # require at least one non-filter taxon to meet MIN_AVG_RELATIVE
# ---- End filters ----


def _resolve_root(base: Path, override: Path | None) -> Path:
    if override:
        return override if override.is_absolute() else (base / override)
    return base / "ncbi_fibroin_sequences"


def main() -> int:
    tic = time.perf_counter()
    groups = GROUPS or DEFAULT_GROUPS

    analyses = ANALYSES
    compare_out_dir = COMPARE_OUT_DIR
    compare_plots_dir = COMPARE_PLOTS_DIR
    compare_filters = {
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
        "group_legend_override": GROUP_LEGEND,
        "group_map_override": GROUPS,
        "group_label": GROUP_LABEL,
    }

    if USE_COMPARE_CONFIG:
        from caddisfly import compare_aa_group_motifs as compare_cfg

        analyses = compare_cfg.COMPARISONS
        compare_out_dir = compare_cfg.DEFAULT_OUT_DIR
        compare_plots_dir = compare_cfg.PLOTS_OUT_DIR
        compare_filters.update(
            {
                "partial_full": compare_cfg.PARTIAL_FULL,
                "protein_types": compare_cfg.PROTEIN_TYPES,
                "length_range": compare_cfg.LENGTH_RANGE,
                "length_threshold": compare_cfg.LENGTH_THRESHOLD,
                "length_mode": compare_cfg.LENGTH_MODE,
                "longest_factor": compare_cfg.LONGEST_FACTOR,
                "min_motif_length": compare_cfg.MIN_MOTIF_LENGTH,
                "max_motif_length": compare_cfg.MAX_MOTIF_LENGTH,
                "min_avg_relative": compare_cfg.MIN_AVG_RELATIVE,
                "use_unique_motifs": getattr(compare_cfg, "USE_UNIQUE_MOTIFS", False),
                "use_window_normalized": getattr(compare_cfg, "USE_WINDOW_NORMALIZED", False),
                "min_avg_relative_taxa": compare_cfg.MIN_AVG_RELATIVE_TAXA,
                "min_avg_relative_min_taxa": compare_cfg.MIN_AVG_RELATIVE_MIN_TAXA,
                "require_non_filter_taxa": compare_cfg.REQUIRE_NON_FILTER_TAXA,
            }
        )

    # Write effective config for debugging/repro.
    compare_out_dir_path = Path(compare_out_dir).expanduser().resolve()
    compare_out_dir_path.mkdir(parents=True, exist_ok=True)
    effective_config = {
        "use_compare_config": USE_COMPARE_CONFIG,
        "reuse_existing_outputs": REUSE_EXISTING_OUTPUTS,
        "analyses": analyses,
        "compare_out_dir": str(compare_out_dir_path),
        "compare_plots_dir": str(Path(compare_plots_dir).expanduser().resolve()),
        "compare_filters": compare_filters,
    }
    (compare_out_dir_path / "run_aa_motifs_and_compare_config.json").write_text(
        json.dumps(effective_config, indent=2, default=str),
        encoding="utf-8",
    )

    comparisons: list[dict[str, object]] = []
    for entry in analyses:
        data_root = Path(entry.get("data_root") or entry.get("outputs_root") or BASE).expanduser().resolve()
        outputs_root = Path(entry.get("outputs_root") or data_root).expanduser().resolve()
        taxonomy_terms = entry.get("taxonomy_terms")

        if not REUSE_EXISTING_OUTPUTS:
            resolved_root = _resolve_root(BASE, data_root)
            root = resolved_root.resolve()

            records = load_records_from_root(root, skip_tree=True)
            filtered_records = filter_records(
                records,
                taxonomy_terms=taxonomy_terms,
                protein_types=entry.get("protein_types", PROTEIN_TYPES),
                partial_full=entry.get("partial_full", PARTIAL_FULL),
                length_range=entry.get("length_range", LENGTH_RANGE),
                length_threshold=entry.get("length_threshold", LENGTH_THRESHOLD),
                length_mode=entry.get("length_mode", LENGTH_MODE),
                longest_factor=entry.get("longest_factor", LONGEST_FACTOR),
                longest_factor_scope=LONGEST_FACTOR_SCOPE,
            )
            filtered_records = dedupe_records_by_accession(filtered_records)
            if entry.get("use_longest_only", USE_LONGEST_ONLY):
                filtered_records = select_longest_records(filtered_records, scope=LONGEST_SCOPE)

            summary = analyze_group_motifs(
                filtered_records,
                groups=groups,
                include_motif_string=INCLUDE_MOTIF_STRING,
                include_motif_counts=INCLUDE_MOTIF_COUNTS,
                min_motif_length=entry.get("min_motif_length", MIN_MOTIF_LENGTH),
                max_motif_length=entry.get("max_motif_length", MAX_MOTIF_LENGTH),
                skip_unknown_motifs=SKIP_UNKNOWN_MOTIFS,
                include_unique_motifs=INCLUDE_UNIQUE_MOTIFS,
                unique_motif_length=UNIQUE_MOTIF_LENGTH,
                skip_combined_motifs=SKIP_COMBINED_MOTIFS,
                include_reflected_sequence=INCLUDE_REFLECTED_SEQUENCE,
                include_similar_motifs=INCLUDE_SIMILAR_MOTIFS,
                include_taxonomy_shared_motifs=INCLUDE_TAXONOMY_SHARED_MOTIFS,
                taxonomy_terms=taxonomy_terms,
                phylo_tree_path=entry.get("phylo_tree_path", PHYLO_TREE_PATH),
                root=root,
            )

            output_root = resolve_aa_motifs_output_dir(
                outputs_root,
                taxonomy_terms=taxonomy_terms,
                partial_full=entry.get("partial_full", PARTIAL_FULL),
                protein_types=entry.get("protein_types", PROTEIN_TYPES),
                length_range=entry.get("length_range", LENGTH_RANGE),
                length_threshold=entry.get("length_threshold", LENGTH_THRESHOLD),
                length_mode=entry.get("length_mode", LENGTH_MODE),
                longest_factor=entry.get("longest_factor", LONGEST_FACTOR),
                min_motif_length=MIN_MOTIF_LENGTH if INCLUDE_MOTIF_COUNTS else None,
                max_motif_length=MAX_MOTIF_LENGTH if INCLUDE_MOTIF_COUNTS else None,
            )
            write_json(output_root, summary, "aa_group_motifs_analysis.json")
            write_md(output_root, summary, "aa_group_motifs_analysis.md")
            if INCLUDE_REFLECTED_SEQUENCE:
                write_reflected_json(output_root, summary, "reflected_sequence.json")
                write_reflected_md(output_root, summary, "reflected_sequence.md")

        comparisons.append(
            {
                "label": entry.get("label") or str(taxonomy_terms[0]) if taxonomy_terms else "Unknown",
                "outputs_root": outputs_root,
                "taxonomy_terms": taxonomy_terms,
                "protein_types": entry.get("protein_types", PROTEIN_TYPES),
                "partial_full": entry.get("partial_full", PARTIAL_FULL),
                "length_range": entry.get("length_range", LENGTH_RANGE),
                "length_threshold": entry.get("length_threshold", LENGTH_THRESHOLD),
                "length_mode": entry.get("length_mode", LENGTH_MODE),
                "longest_factor": entry.get("longest_factor", LONGEST_FACTOR),
                "min_motif_length": entry.get("min_motif_length", MIN_MOTIF_LENGTH),
                "max_motif_length": entry.get("max_motif_length", MAX_MOTIF_LENGTH),
            }
        )

    run_compare(comparisons, compare_out_dir, compare_plots_dir, compare_filters, argv=None)

    toc = time.perf_counter()
    print(f"Elapsed: {toc - tic:0.3f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
