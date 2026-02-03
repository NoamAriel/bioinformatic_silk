"""Amino-acid group motif analysis helper script."""
from __future__ import annotations

from pathlib import Path
import sys
import time
from typing import Tuple


BASE = Path(__file__).resolve().parent
PROJECT_ROOT = BASE.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from libraries.amino_acid_group_motifs import (
    DEFAULT_GROUPS,
    analyze_group_motifs,
    dedupe_records_by_accession,
    select_longest_records,
    resolve_aa_motifs_output_dir,
    write_json,
    write_md,
    write_reflected_json,
    write_reflected_md,
)
from libraries.aminoacids_composition_analysis_lib import filter_records, load_records_from_root

# Optional overrides (absolute or relative to this file).
data_root: Path | None = None
# Optional overrides for where outputs should be written.
outputs_root: Path | None = None

# Filters (match aa_composition_analysis.py behavior).
taxonomy_terms = ["Chironomoidea"]  # or [] / None  Chironomoidea, Culicoidea, Tephritoidea, bactrocera
protein_types = ["heavy chain"]  # or None
partial_full = None
length_range: Tuple[int, int] | None = None
length_threshold = None
length_mode = "ge"
longest_factor = None
longest_factor_scope = "species"
use_longest_only = False  # set True to keep only the longest record per species. set false to keep all records.

# Output controls.
include_motif_string = False  # set True for small datasets

# Motif scan settings (group-label motifs like "ππφηπ").
include_motif_counts = True
min_motif_length = 5
max_motif_length = 50
skip_unknown_motifs = True  # skip windows containing "unknown"
include_unique_motifs = True
unique_motif_length = 5  # set to min_motif_length or any fixed length
skip_combined_motifs = True  # drop motifs that are repeats or concatenations
include_reflected_sequence = True
include_similar_motifs = True
include_taxonomy_shared_motifs = True
phylo_tree_path = Path("D:/python/bioinformatic_silk/Diptera/ncbi_fibroin_sequences/phylo_tree.json")

longest_scope = "species"  # "species" or "global"


tic = time.perf_counter()

# Root directory containing the sequences to analyze.
resolved_root = data_root if data_root and data_root.is_absolute() else (BASE / (data_root or "ncbi_fibroin_sequences"))
root = resolved_root.resolve()

records = load_records_from_root(root, skip_tree=True)
filtered_records = filter_records(
    records,
    taxonomy_terms=taxonomy_terms,
    protein_types=protein_types,
    partial_full=partial_full,
    length_range=length_range,
    length_threshold=length_threshold,
    length_mode=length_mode,
    longest_factor=longest_factor,
    longest_factor_scope=longest_factor_scope,
)
filtered_records = dedupe_records_by_accession(filtered_records)
if use_longest_only:
    filtered_records = select_longest_records(filtered_records, scope=longest_scope)

summary = analyze_group_motifs(
    filtered_records,
    groups=DEFAULT_GROUPS,
    include_motif_string=include_motif_string,
    include_motif_counts=include_motif_counts,
    min_motif_length=min_motif_length,
    max_motif_length=max_motif_length,
    skip_unknown_motifs=skip_unknown_motifs,
    include_unique_motifs=include_unique_motifs,
    unique_motif_length=unique_motif_length,
    skip_combined_motifs=skip_combined_motifs,
    include_reflected_sequence=include_reflected_sequence,
    include_similar_motifs=include_similar_motifs,
    include_taxonomy_shared_motifs=include_taxonomy_shared_motifs,
    taxonomy_terms=taxonomy_terms,
    phylo_tree_path=phylo_tree_path,
    root=root,
)

output_root = resolve_aa_motifs_output_dir(
    root,
    taxonomy_terms=taxonomy_terms,
    partial_full=partial_full,
    protein_types=protein_types,
    length_range=length_range,
    length_threshold=length_threshold,
    length_mode=length_mode,
    longest_factor=longest_factor,
    min_motif_length=min_motif_length if include_motif_counts else None,
    max_motif_length=max_motif_length if include_motif_counts else None,
)
json_path = write_json(output_root, summary, "aa_group_motifs_analysis.json")
md_path = write_md(output_root, summary, "aa_group_motifs_analysis.md")
if include_reflected_sequence:
    write_reflected_json(output_root, summary, "reflected_sequence.json")
    write_reflected_md(output_root, summary, "reflected_sequence.md")

print(json_path, md_path)

toc = time.perf_counter()
print(f"Elapsed: {toc - tic:0.3f} seconds")
