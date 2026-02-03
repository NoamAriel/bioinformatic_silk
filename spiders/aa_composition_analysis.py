# Amino acid composition analysis helper script.
# format_aa_composition_analysis.py

# This script performs amino acid composition analysis on protein sequences
# from a specified root directory. It uses the analyze_from_root function from the aminoacids_composition_analysis_lib library.
# Adjust the filter parameters as needed.
# Import necessary modules and set up the library path for custom libraries.
# Make sure to have the required libraries in the 'libraries' directory.
# You can modify the main function as needed.
# The script performs the following task:
# 1) Analyze amino acid composition in protein sequences from a specified root directory.   
# 2) Print the paths to the generated JSON and Markdown output files.
# 3) You can adjust the filter parameters such as letters, taxonomy_terms, protein_types,
#  partial_full, length_range, length_threshold, length_mode, and longest_factor as needed.
# 4) Finally, it prints the paths to the generated JSON and Markdown output files.  


from pathlib import Path
import sys
import time

from matplotlib.pylab import full

BASE = Path(__file__).resolve().parent
PROJECT_ROOT = BASE.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from libraries.aminoacids_composition_analysis_lib import analyze_from_root

tic = time.perf_counter()

# Optional overrides (absolute or relative to this file)
# - data_root: where the ncbi_fibroin_sequences live
# - plots_root: where plots should be written (defaults to <data_root>/plots)
# - tables_root: where tables should be written (defaults to <data_root>/tables)
data_root: Path | None = Path(r"D:\python\bioinformatic_silk\spiders\ncbi_fibroin_sequences")
plots_root: Path | None = Path(r"D:\python\bioinformatic_silk\spiders\ncbi_fibroin_sequences\plots")
phylo_tree_json = Path(r"D:\python\bioinformatic_silk\spiders\ncbi_fibroin_sequences\phylo_tree.json")
tables_root: Path | None = Path(r"D:\python\bioinformatic_silk\spiders\ncbi_fibroin_sequences\tables")

# Root directory containing the sequences to analyze.
resolved_root = data_root if data_root and data_root.is_absolute() else (BASE / (data_root or "ncbi_fibroin_sequences"))
root = resolved_root.resolve()

resolved_plots = None
if plots_root:
    resolved_plots = plots_root if plots_root.is_absolute() else (BASE / plots_root)
    resolved_plots = resolved_plots.resolve()

resolved_tables = None
if tables_root:
    resolved_tables = tables_root if tables_root.is_absolute() else (BASE / tables_root)
    resolved_tables = resolved_tables.resolve()


json_path, md_path = analyze_from_root(
    root=root,
    letters="wfy",       #stfwhydenqrk avli gp cm   krh nq ed ednqkrh nqkrh wpfh sty   wfy    # amino acids you care about. for example, "sae" or "c" or "gyplsrh"
    taxonomy_terms=["araneae"],           # or [] / None
    protein_types=["CrSp"],        #  MaSp1 MaSp2 MaSp3 MiSp Flag  AgSp     PySp      or [] / None
    partial_full=None,                  # "full", "partial", or None
length_range=None,                    # e.g., (100, 2450) or None
length_threshold=1000,                # e.g., 1500 or None
length_mode="ge",                     # optional, default is "ge". ge: greater equal, le: less equal
longest_factor=None,                   # optional default is 2.0
longest_factor_scope="species",        # "species" (per organism) or "global" (all records)
tables_dir=resolved_tables,
plots_dir=resolved_plots,
phylo_tree_json=phylo_tree_json,
)

print(json_path, md_path)

toc = time.perf_counter()
print(f"Elapsed: {toc - tic:0.3f} seconds")
