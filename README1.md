# bioinformatic_silk

Bioinformatic Silk studies adhesive aquatic silk, especially from caddisfly (Trichoptera), and compares it to other silk-producing animals to connect sequence features, amino-acid composition, and material properties.

**Using the GUI App**

Quick start:
1. Ensure you have Python 3 with Tkinter available.
2. Install Python packages: `requests`, `beautifulsoup4`, `pandas`, `numpy`, `matplotlib`.
3. From the termianl, run:

uv run ./gui_bioinformatic_silk.py


The GUI app opens a multi-tab window and writes logs to the bottom panel. Many fields include a `?` help button with short explanations.

**GUI Tabs**
- `NCBI Scraper`: Download protein records from NCBI for a taxonomic order (and optional families). Set `Order name`, optional `Family names`, `Protein terms`, `Expected types JSON`, and output settings, then click `Run NCBI Scraper`.
- `Taxonomy Tree`: Render a phylogenetic tree using `phylo_tree.json` and `Species_Index.md`. Choose output path/format and optional rank range, set filters, then click `Generate Taxonomy Tree`.
- `AA Composition`: Analyze amino-acid composition for selected letters. Choose data root, letters, optional plots/tables directories, set filters, then click `Run AA Composition`.
- `SXn Analysis`: Analyze serine S-Xn motifs. Choose data root, min/max motif length, optional plots/tables directories, set filters, then click `Run SXn Analysis + Plots`.
- `AA Motifs Compare`: Compare amino-acid group motifs across analyses. Provide `Analyses JSON` (list of dicts), optional group map/legend overrides, filters, motif lengths, toggles, and output directories, then click `Run AA Motifs Compare`.
- `Help`: Built-in reference for filters, tabs, and equations used by the analyses.

**Filters (used across analysis tabs)**
- `Taxonomy/species`: Comma-separated taxonomy terms or species names to include.
- `Protein types`: Comma-separated substrings to match protein types (e.g., heavy chain).
- `Partial/full`: Choose full, partial, or all records.
- `Length range`: Inclusive min-max sequence length.
- `Length threshold + mode`: Keep lengths >= threshold (`ge`) or <= threshold (`le`).
- `Longest factor` + `Longest scope`: Keep sequences near the longest length (per species or global scope).

**Outputs**
Results are written as JSON/MD summaries and plots under the selected data root (or the optional output directories you provide). The GUI also opens new plot windows for recent images.
