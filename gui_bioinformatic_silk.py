import json
import threading
import queue
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import tkinter as tk
from tkinter import ttk, filedialog, messagebox


PROJECT_ROOT = Path(__file__).resolve().parent

HELP_TEXTS: dict[str, str] = {
    "Order name": "Taxonomic order to search in NCBI (e.g., Trichoptera).",
    "Family names": "Optional comma-separated family names to also search.",
    "Protein terms": "Comma-separated protein terms for the NCBI query (e.g., fibroin).",
    "Expected types JSON": "JSON list or dict used to classify protein types from NCBI titles.",
    "NCBI API key": "Optional NCBI API key to increase request limits.",
    "Output root": "Root folder where NCBI sequences will be written.",
    "Sleep time (s)": "Delay between requests to NCBI (seconds). \nleave blank for default (0.2s).",
    "Path max length": "Optional max path length. \nleave blank for default (200).",
    "Allow long paths": "On Windows, allow extended-length paths (\\\\?\\ prefix).",
    "Data root": "Root folder containing ncbi_fibroin_sequences for the selected group.",
    "Phylo tree JSON": "Path to phylo_tree.json used for taxonomy ordering/plotting.",
    "Species Index": "Species_Index.md used to annotate taxonomy tree.",
    "Output path": "Output file path without extension (e.g., taxonomy_tree).",
    "Output format": "Image format for the taxonomy tree output.",
    "Rank from": "Start rank or rank name to display in the taxonomy tree.",
    "Rank to": "End rank or rank name to display in the taxonomy tree.",
    "Letters": "Amino acids to analyze, e.g. nqkrh or sty.",
    "Plots dir (opt)": "Optional override folder for plots output.",
    "Tables dir (opt)": "Optional override folder for tables output.",
    "Min motif length": "Minimum motif length (n). For SXn, this sets the smallest n.",
    "Max motif length": "Maximum motif length (n). For SXn, this sets the largest n.",
    "Taxonomy/species": "Comma-separated taxonomy terms or species names for filtering.",
    "Protein types": "Comma-separated substrings to match protein type (e.g., heavy chain).",
    "Partial/full": "full = only full sequences, partial = only partial, all = no filter.",
    "Length range": "Inclusive min-max sequence length filter.",
    "Length threshold": "Single length threshold used with Length mode.",
    "Length mode": "ge keeps lengths >= threshold, le keeps lengths <= threshold.",
    "Longest factor": "Keeps sequences near the longest length (per species or global). \nFor example: if L_max=3000 and f=2, keep lengths >= 1500; if f=4, keep lengths >= 2250.",
    "Longest scope": "Apply longest-factor per species or across all records.",
    "Compare out dir": "Where comparison tables/markdown are written.",
    "Compare plots dir": "Where comparison plots are written.",
    "Analyses JSON (list of dicts)": "List of analysis targets with label/outputs_root/taxonomy_terms.",
    "Group map JSON": "Override group symbol -> amino acids mapping.",
    "Group legend JSON": "Override group symbol -> label mapping for plots.",
    "Group label": "Optional folder label for grouping outputs.",
    "Unique motif length": "Fixed k for unique motif counting.",
    "Min avg relative": "Minimum average relative percent to keep a motif.",
    "Min avg relative taxa": "Taxa that must meet the min average relative threshold.",
    "Min avg min taxa": "Number of taxa that must meet the threshold.",
    "Include motif counts": "Include overlapping motif counts/frequencies.",
    "Include unique motifs": "Include fixed-length unique motif counts.",
    "Include reflected sequence": "Write reflected sequence outputs.",
    "Include similar motifs": "Include per-species shared motifs.",
    "Include taxonomy shared motifs": "Include motifs shared across taxa.",
    "Include motif strings": "Store full motif strings (large output).",
    "Skip unknown motifs": "Remove motifs containing 'unknown'.",
    "Skip combined motifs": "Remove repeated/concatenated motifs.",
    "Use longest only per species": "Keep only the longest record per species.",
    "Compare using unique motifs": "Compare unique motifs instead of overlapping motifs.",
    "Use window normalized": "Normalize motif counts by window count (L - k + 1).",
    "Require non-filter taxa": "Require at least one non-filter taxon to meet the threshold.",
    "Reuse existing outputs": "Skip recomputing motif outputs if they already exist.",
}


def _to_list(text: str) -> list[str]:
    if not text:
        return []
    parts = []
    for chunk in text.replace(";", ",").split(","):
        cleaned = chunk.strip()
        if cleaned:
            parts.append(cleaned)
    return parts


def _opt_int(text: str) -> int | None:
    text = text.strip()
    return int(text) if text else None


def _opt_float(text: str) -> float | None:
    text = text.strip()
    return float(text) if text else None


def _opt_path(text: str) -> Path | None:
    text = text.strip()
    return Path(text) if text else None


def _length_range(min_text: str, max_text: str) -> tuple[int, int] | None:
    if not min_text.strip() and not max_text.strip():
        return None
    if not min_text.strip() or not max_text.strip():
        raise ValueError("Length range needs both min and max.")
    return int(min_text.strip()), int(max_text.strip())


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Bioinformatic Silk GUI")
        self.geometry("1100x760")
        self.minsize(1000, 680)
        self._log_queue: queue.Queue[str] = queue.Queue()
        self._image_refs: list[tk.PhotoImage] = []
        self._build_ui()
        self.after(150, self._drain_log)

    def _build_ui(self) -> None:
        root = ttk.Frame(self, padding=10)
        root.pack(fill=tk.BOTH, expand=True)

        ttk.Label(root, text=f"Project root: {PROJECT_ROOT}").pack(anchor=tk.W)

        notebook = ttk.Notebook(root)
        notebook.pack(fill=tk.BOTH, expand=True, pady=(10, 6))

        self.tab_ncbi = ttk.Frame(notebook)
        self.tab_tax = ttk.Frame(notebook)
        self.tab_aa = ttk.Frame(notebook)
        self.tab_sxn = ttk.Frame(notebook)
        self.tab_motifs = ttk.Frame(notebook)
        self.tab_help = ttk.Frame(notebook)

        notebook.add(self.tab_ncbi, text="NCBI Scraper")
        notebook.add(self.tab_tax, text="Taxonomy Tree")
        notebook.add(self.tab_aa, text="AA Composition")
        notebook.add(self.tab_sxn, text="SXn Analysis")
        notebook.add(self.tab_motifs, text="AA Motifs Compare")
        notebook.add(self.tab_help, text="Help")

        self._build_ncbi_tab(self.tab_ncbi)
        self._build_tax_tab(self.tab_tax)
        self._build_aa_tab(self.tab_aa)
        self._build_sxn_tab(self.tab_sxn)
        self._build_motifs_tab(self.tab_motifs)
        self._build_help_tab(self.tab_help)

        log_frame = ttk.LabelFrame(root, text="Log")
        log_frame.pack(fill=tk.BOTH, expand=False)
        self.log_text = tk.Text(log_frame, height=8, wrap="word")
        self.log_text.pack(fill=tk.BOTH, expand=True)

    def _row(self, parent: ttk.Frame, label: str) -> ttk.Frame:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=2)
        ttk.Label(row, text=label, width=22).pack(side=tk.LEFT)
        help_text = HELP_TEXTS.get(label)
        if help_text:
            ttk.Button(row, text="?", width=2, command=lambda: messagebox.showinfo(label, help_text)).pack(
                side=tk.LEFT, padx=(0, 6)
            )
        return row

    def _entry(self, parent: ttk.Frame, label: str, default: str = "") -> ttk.Entry:
        row = self._row(parent, label)
        entry = ttk.Entry(row)
        entry.insert(0, default)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        return entry

    def _combo(self, parent: ttk.Frame, label: str, values: list[str], default: str) -> ttk.Combobox:
        row = self._row(parent, label)
        combo = ttk.Combobox(row, values=values, state="readonly")
        combo.set(default)
        combo.pack(side=tk.LEFT)
        return combo

    def _path(self, parent: ttk.Frame, label: str, default: str = "") -> ttk.Entry:
        row = self._row(parent, label)
        entry = ttk.Entry(row)
        entry.insert(0, default)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        def _browse() -> None:
            path = filedialog.askdirectory()
            if path:
                entry.delete(0, tk.END)
                entry.insert(0, path)

        ttk.Button(row, text="Browse", command=_browse).pack(side=tk.LEFT, padx=6)
        return entry

    def _file(self, parent: ttk.Frame, label: str, default: str = "") -> ttk.Entry:
        row = self._row(parent, label)
        entry = ttk.Entry(row)
        entry.insert(0, default)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

        def _browse() -> None:
            path = filedialog.askopenfilename()
            if path:
                entry.delete(0, tk.END)
                entry.insert(0, path)

        ttk.Button(row, text="Browse", command=_browse).pack(side=tk.LEFT, padx=6)
        return entry

    def _range(self, parent: ttk.Frame, label: str) -> tuple[ttk.Entry, ttk.Entry]:
        row = self._row(parent, label)
        min_entry = ttk.Entry(row, width=10)
        min_entry.pack(side=tk.LEFT)
        ttk.Label(row, text="to").pack(side=tk.LEFT, padx=4)
        max_entry = ttk.Entry(row, width=10)
        max_entry.pack(side=tk.LEFT)
        return min_entry, max_entry

    def _check(self, parent: ttk.Frame, label: str, default: bool) -> tk.BooleanVar:
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=2)
        var = tk.BooleanVar(value=default)
        ttk.Checkbutton(row, text=label, variable=var).pack(side=tk.LEFT)
        help_text = HELP_TEXTS.get(label)
        if help_text:
            ttk.Button(row, text="?", width=2, command=lambda: messagebox.showinfo(label, help_text)).pack(
                side=tk.LEFT, padx=(6, 0)
            )
        return var

    def _log(self, msg: str) -> None:
        self._log_queue.put(msg)

    def _drain_log(self) -> None:
        while not self._log_queue.empty():
            msg = self._log_queue.get_nowait()
            self.log_text.insert(tk.END, msg + "\n")
            self.log_text.see(tk.END)
        self.after(150, self._drain_log)

    def _show_image_popup(self, image_path: Path, title: str) -> None:
        try:
            if not image_path.exists():
                self._log(f"Image not found: {image_path}")
                return
            win = tk.Toplevel(self)
            win.title(title)
            container = ttk.Frame(win)
            container.pack(fill=tk.BOTH, expand=True)

            toolbar = ttk.Frame(container)
            toolbar.pack(fill=tk.X, padx=6, pady=4)

            canvas = tk.Canvas(container, background="white", highlightthickness=0)
            canvas.pack(fill=tk.BOTH, expand=True)

            orig = tk.PhotoImage(file=str(image_path))
            zoom_factor = {"value": 1.0}
            pan = {"x": 0, "y": 0, "start_x": 0, "start_y": 0}

            def render() -> None:
                w = canvas.winfo_width()
                h = canvas.winfo_height()
                if w <= 1 or h <= 1:
                    return
                base_scale = min(w / orig.width(), h / orig.height())
                scale = max(0.05, base_scale * zoom_factor["value"])
                frac = Fraction(scale).limit_denominator(10)
                img = orig.zoom(frac.numerator, frac.numerator).subsample(frac.denominator, frac.denominator)
                canvas.delete("all")
                canvas.create_image(w // 2 + pan["x"], h // 2 + pan["y"], image=img, anchor="center")
                win._img = img  # keep reference

            def zoom_in() -> None:
                zoom_factor["value"] *= 1.25
                render()

            def zoom_out() -> None:
                zoom_factor["value"] /= 1.25
                render()

            def zoom_reset() -> None:
                zoom_factor["value"] = 1.0
                pan["x"] = 0
                pan["y"] = 0
                render()

            ttk.Button(toolbar, text="+", width=3, command=zoom_in).pack(side=tk.LEFT)
            ttk.Button(toolbar, text="-", width=3, command=zoom_out).pack(side=tk.LEFT, padx=(4, 0))
            ttk.Button(toolbar, text="Reset", command=zoom_reset).pack(side=tk.LEFT, padx=(8, 0))

            def start_pan(event) -> None:
                pan["start_x"] = event.x
                pan["start_y"] = event.y

            def do_pan(event) -> None:
                dx = event.x - pan["start_x"]
                dy = event.y - pan["start_y"]
                pan["start_x"] = event.x
                pan["start_y"] = event.y
                pan["x"] += dx
                pan["y"] += dy
                render()

            def on_wheel(event) -> None:
                if event.state & 0x0004:  # Ctrl
                    if event.delta > 0:
                        zoom_in()
                    else:
                        zoom_out()

            canvas.bind("<Configure>", lambda _e: render())
            canvas.bind("<ButtonPress-1>", start_pan)
            canvas.bind("<B1-Motion>", do_pan)
            canvas.bind("<MouseWheel>", on_wheel)
            win.after(50, render)
        except Exception as exc:
            self._log(f"Failed to show image: {exc}")

    def _show_images_popups(self, image_paths: list[Path], title_prefix: str) -> None:
        for idx, path in enumerate(image_paths, start=1):
            self._show_image_popup(path, f"{title_prefix} ({idx}/{len(image_paths)})")

    def _latest_png(self, root: Path) -> Path | None:
        if not root.exists():
            return None
        pngs = list(root.rglob("*.png"))
        if not pngs:
            return None
        return max(pngs, key=lambda p: p.stat().st_mtime)

    def _all_pngs(self, root: Path) -> list[Path]:
        if not root.exists():
            return []
        pngs = list(root.rglob("*.png"))
        pngs.sort(key=lambda p: p.stat().st_mtime)
        return pngs

    def _new_pngs(self, root: Path, since_ts: float) -> list[Path]:
        if not root.exists():
            return []
        pngs = [p for p in root.rglob("*.png") if p.stat().st_mtime >= since_ts]
        pngs.sort(key=lambda p: p.stat().st_mtime)
        return pngs

    def _build_filters(self, parent: ttk.Frame, defaults: dict[str, str]) -> dict[str, Any]:
        fields: dict[str, Any] = {}
        fields["taxonomy_terms"] = self._entry(parent, "Taxonomy/species", defaults.get("taxonomy_terms", ""))
        fields["protein_types"] = self._entry(parent, "Protein types", defaults.get("protein_types", ""))
        fields["partial_full"] = self._combo(parent, "Partial/full", ["all", "full", "partial"], defaults.get("partial_full", "all"))
        fields["length_range"] = self._range(parent, "Length range")
        fields["length_threshold"] = self._entry(parent, "Length threshold", defaults.get("length_threshold", ""))
        fields["length_mode"] = self._combo(parent, "Length mode", ["ge", "le"], defaults.get("length_mode", "ge"))
        fields["longest_factor"] = self._entry(parent, "Longest factor", defaults.get("longest_factor", ""))
        fields["longest_scope"] = self._combo(parent, "Longest scope", ["species", "global"], defaults.get("longest_scope", "species"))
        return fields

    def _read_filters(self, fields: dict[str, Any]) -> dict[str, Any]:
        taxonomy_terms = _to_list(fields["taxonomy_terms"].get()) or None
        protein_types = _to_list(fields["protein_types"].get()) or None
        partial_full_val = fields["partial_full"].get().strip().lower()
        partial_full = None if partial_full_val == "all" else partial_full_val
        length_range = _length_range(fields["length_range"][0].get(), fields["length_range"][1].get())
        length_threshold = _opt_int(fields["length_threshold"].get())
        length_mode = fields["length_mode"].get().strip() or "ge"
        longest_factor = _opt_float(fields["longest_factor"].get())
        longest_scope = fields["longest_scope"].get().strip() or "species"
        return {
            "taxonomy_terms": taxonomy_terms,
            "protein_types": protein_types,
            "partial_full": partial_full,
            "length_range": length_range,
            "length_threshold": length_threshold,
            "length_mode": length_mode,
            "longest_factor": longest_factor,
            "longest_scope": longest_scope,
        }

    def _build_tax_tab(self, parent: ttk.Frame) -> None:
        frame = ttk.Frame(parent, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        self.tax_data_root = self._path(frame, "Data root", str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences"))
        self.tax_tree_path = self._file(frame, "Phylo tree JSON", str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences" / "phylo_tree.json"))
        self.tax_index_path = self._file(frame, "Species Index", str(PROJECT_ROOT / "caddisfly" / "Species_Index.md"))
        self.tax_out_path = self._path(frame, "Output path", str(PROJECT_ROOT / "caddisfly" / "taxonomy_tree"))
        self.tax_output_format = self._combo(frame, "Output format", ["png", "svg", "pdf"], "png")
        try:
            from libraries.generate_taxonomy_graph import RANK_ORDER
            rank_values = ["(any)"] + [r.title() for r in RANK_ORDER]
        except Exception:
            rank_values = [
                "(any)",
                "Kingdom",
                "Phylum",
                "Class",
                "Order",
                "Suborder",
                "Family",
                "Subfamily",
                "Genus",
                "Species",
                "Subspecies",
            ]
        self.tax_rank_from = self._combo(frame, "Rank from", rank_values, "(any)")
        self.tax_rank_to = self._combo(frame, "Rank to", rank_values, "(any)")
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        self.tax_filters = self._build_filters(frame, {"partial_full": "all"})
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        ttk.Button(frame, text="Generate Taxonomy Tree", command=self._run_taxonomy_tree).pack(anchor=tk.W)

    def _build_ncbi_tab(self, parent: ttk.Frame) -> None:
        frame = ttk.Frame(parent, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        self.ncbi_order = self._entry(frame, "Order name", "Trichoptera")
        self.ncbi_family = self._entry(frame, "Family names", "")
        self.ncbi_terms = self._entry(frame, "Protein terms", "fibroin")
        self.ncbi_api_key = self._entry(frame, "NCBI API key", "")
        self.ncbi_output_root = self._path(frame, "Output root", str(PROJECT_ROOT / "ncbi_sequences"))
        self.ncbi_sleep = self._entry(frame, "Sleep time (s)", "0.2")
        self.ncbi_path_max = self._entry(frame, "Path max length", "200")
        self.ncbi_allow_long = self._check(frame, "Allow long paths", False)

        row = self._row(frame, "Expected types JSON")
        default_expected = {
            "heavy chain": [
                "heavy chain",
                "fib-h",
                "h-fibroin",
                "h chain",
                "fibroin heavy chain",
            ],
            "light chain": [
                "light chain",
                "fib-l",
                "l-fibroin",
                "l chain",
                "fibroin light chain",
            ],
        }
        self.ncbi_expected = tk.Text(frame, height=6, wrap="word")
        self.ncbi_expected.insert("1.0", json.dumps(default_expected, indent=2))
        self.ncbi_expected.pack(fill=tk.X, expand=False, pady=(2, 6))

        ttk.Button(frame, text="Run NCBI Scraper", command=self._run_ncbi_scraper).pack(anchor=tk.W)

    def _build_aa_tab(self, parent: ttk.Frame) -> None:
        frame = ttk.Frame(parent, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        self.aa_data_root = self._path(frame, "Data root", str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences"))
        self.aa_letters = self._entry(frame, "Letters", "nqkrh")
        self.aa_phylo_tree = self._file(frame, "Phylo tree JSON", str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences" / "phylo_tree.json"))
        self.aa_plots_dir = self._path(frame, "Plots dir (opt)", "")
        self.aa_tables_dir = self._path(frame, "Tables dir (opt)", "")
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        self.aa_filters = self._build_filters(frame, {"taxonomy_terms": "trichoptera", "protein_types": "heavy chain", "partial_full": "full"})
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        ttk.Button(frame, text="Run AA Composition", command=self._run_aa_composition).pack(anchor=tk.W)

    def _build_sxn_tab(self, parent: ttk.Frame) -> None:
        frame = ttk.Frame(parent, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)
        self.sxn_data_root = self._path(frame, "Data root", str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences"))
        self.sxn_min_n = self._entry(frame, "Min motif length", "3")
        self.sxn_max_n = self._entry(frame, "Max motif length", "50")
        self.sxn_plots_dir = self._path(frame, "Plots dir (opt)", "")
        self.sxn_tables_dir = self._path(frame, "Tables dir (opt)", "")
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        self.sxn_filters = self._build_filters(frame, {"taxonomy_terms": "trichoptera", "protein_types": "heavy chain", "partial_full": "full"})
        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        ttk.Button(frame, text="Run SXn Analysis + Plots", command=self._run_sxn).pack(anchor=tk.W)

    def _build_motifs_tab(self, parent: ttk.Frame) -> None:
        outer = ttk.Frame(parent, padding=0)
        outer.pack(fill=tk.BOTH, expand=True)

        canvas = tk.Canvas(outer, highlightthickness=0)
        scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        frame = ttk.Frame(canvas, padding=10)
        window_id = canvas.create_window((0, 0), window=frame, anchor="nw")

        def _on_frame_configure(_event) -> None:
            canvas.configure(scrollregion=canvas.bbox("all"))

        def _on_canvas_configure(event) -> None:
            canvas.itemconfigure(window_id, width=event.width)

        frame.bind("<Configure>", _on_frame_configure)
        canvas.bind("<Configure>", _on_canvas_configure)
        canvas.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))

        defaults = {}
        try:
            from caddisfly import run_aa_motifs_and_compare as cfg

            defaults["compare_out_dir"] = str(Path(cfg.COMPARE_OUT_DIR))
            defaults["compare_plots_dir"] = str(Path(cfg.COMPARE_PLOTS_DIR))
            defaults["analyses_json"] = json.dumps(cfg.ANALYSES, indent=2, default=str)
            defaults["groups_json"] = json.dumps(cfg.GROUPS or {}, indent=2, ensure_ascii=True)
            defaults["group_legend_json"] = json.dumps(cfg.GROUP_LEGEND or {}, indent=2, ensure_ascii=True)
            defaults["group_label"] = str(cfg.GROUP_LABEL or "")
        except Exception:
            defaults["compare_out_dir"] = str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences" / "data_analysis" / "aa_motifs_comparison")
            defaults["compare_plots_dir"] = str(PROJECT_ROOT / "caddisfly" / "ncbi_fibroin_sequences" / "plots" / "aa_motifs_comparison")
            defaults["analyses_json"] = "[]"
            defaults["groups_json"] = "{}"
            defaults["group_legend_json"] = "{}"
            defaults["group_label"] = ""

        self.motif_compare_out = self._path(frame, "Compare out dir", defaults["compare_out_dir"])
        self.motif_compare_plots = self._path(frame, "Compare plots dir", defaults["compare_plots_dir"])

        ttk.Label(frame, text="Analyses JSON (list of dicts)").pack(anchor=tk.W)
        self.motif_analyses = tk.Text(frame, height=7, wrap="word")
        self.motif_analyses.insert("1.0", defaults["analyses_json"])
        self.motif_analyses.pack(fill=tk.X, expand=False)

        ttk.Label(frame, text="Group map JSON").pack(anchor=tk.W, pady=(6, 0))
        self.motif_groups = tk.Text(frame, height=4, wrap="word")
        self.motif_groups.insert("1.0", defaults["groups_json"])
        self.motif_groups.pack(fill=tk.X, expand=False)

        ttk.Label(frame, text="Group legend JSON").pack(anchor=tk.W, pady=(6, 0))
        self.motif_group_legend = tk.Text(frame, height=4, wrap="word")
        self.motif_group_legend.insert("1.0", defaults["group_legend_json"])
        self.motif_group_legend.pack(fill=tk.X, expand=False)

        self.motif_group_label = self._entry(frame, "Group label", defaults["group_label"])

        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        self.motif_filters = self._build_filters(frame, {"partial_full": "all", "protein_types": "heavy chain"})
        self.motif_min_len = self._entry(frame, "Min motif length", "5")
        self.motif_max_len = self._entry(frame, "Max motif length", "50")
        self.motif_unique_len = self._entry(frame, "Unique motif length", "5")
        checks = ttk.Frame(frame)
        checks.pack(fill=tk.X, pady=4)
        self.motif_include_counts = self._check(checks, "Include motif counts", True)
        self.motif_include_unique = self._check(checks, "Include unique motifs", True)
        self.motif_include_reflected = self._check(checks, "Include reflected sequence", True)
        self.motif_include_similar = self._check(checks, "Include similar motifs", True)
        self.motif_include_tax_shared = self._check(checks, "Include taxonomy shared motifs", True)
        self.motif_include_string = self._check(checks, "Include motif strings", False)
        self.motif_skip_unknown = self._check(checks, "Skip unknown motifs", True)
        self.motif_skip_combined = self._check(checks, "Skip combined motifs", True)
        self.motif_use_longest_only = self._check(checks, "Use longest only per species", False)

        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        self.motif_min_avg_relative = self._entry(frame, "Min avg relative", "3")
        self.motif_min_avg_taxa = self._entry(frame, "Min avg relative taxa", "Integripalpia, Annulipalpia")
        self.motif_min_avg_taxa_count = self._entry(frame, "Min avg min taxa", "2")
        self.motif_use_unique_compare = self._check(frame, "Compare using unique motifs", False)
        self.motif_use_window_norm = self._check(frame, "Use window normalized", False)
        self.motif_require_non_filter = self._check(frame, "Require non-filter taxa", True)
        self.motif_reuse_outputs = self._check(frame, "Reuse existing outputs", False)

        ttk.Separator(frame).pack(fill=tk.X, pady=6)
        ttk.Button(frame, text="Run AA Motifs Compare", command=self._run_motifs_compare).pack(anchor=tk.W)

    def _build_help_tab(self, parent: ttk.Frame) -> None:
        frame = ttk.Frame(parent, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        help_text = (
            "Overview\n"
            "This GUI wraps the analysis scripts in this repository. Choose a tab, set filters, and run.\n\n"
            "Filters (applies to AA Composition, SXn, and AA Motifs Compare)\n"
            "- Taxonomy/species: Comma-separated names. Matches taxonomy lineage or organism name.\n"
            "- Protein types: Comma-separated substrings (e.g., heavy chain, light chain).\n"
            "- Partial/full: full = only full sequences, partial = only partial, all = no filter.\n"
            "- Length range: Inclusive min-max sequence length.\n"
            "- Length threshold + mode: ge keeps >= threshold, le keeps <= threshold.\n"
            "- Longest factor: Keeps sequences near the longest within a species (or global if selected).\n\n"
            "NCBI Scraper tab\n"
            "- Downloads protein records from NCBI Protein for an order (and optional families).\n"
            "- Protein terms are combined with organism filters to find relevant sequences.\n"
            "- Expected types JSON helps classify heavy/light chain names from titles.\n"
            "- Output root is the folder where sequences and taxonomy trees are saved.\n\n"
            "Taxonomy Tree tab\n"
            "- Uses phylo_tree.json and Species_Index.md to render a tree image.\n"
            "- Rank From/To limits the displayed taxonomic range.\n"
            "- Filters apply before rendering to keep only species that match your criteria.\n\n"
            "AA Composition tab\n"
            "- Runs amino-acid composition analysis and writes JSON/MD plus plots/tables.\n"
            "- Letters: amino acids to analyze, e.g. nqkrh or sty.\n"
            "- Plots dir / Tables dir are optional overrides.\n\n"
            "SXn Analysis tab\n"
            "- Runs serine and S-Xn motif analysis and generates plots.\n"
            "- Min/Max motif length controls the n range.\n"
            "- Plots/tables directories are optional overrides.\n\n"
            "AA Motifs Compare tab\n"
            "- Analyses JSON: list of entries with label/outputs_root/taxonomy_terms, etc.\n"
            "- Group map/legend: override motif group definitions and labels.\n"
            "- Motif lengths and toggles control which motifs are counted and written.\n"
            "- Compare settings filter which motifs are shown in the comparison plots.\n"
            "- Reuse existing outputs skips recomputing per-taxonomy motif outputs.\n"
            "\nAA Motifs toggles\n"
            "- Include motif counts: write motif counts/frequencies for overlapping motifs.\n"
            "- Include unique motifs: include fixed-length unique motif counts.\n"
            "- Include reflected sequence: write reflected sequence outputs.\n"
            "- Include similar motifs: include per-species shared motifs.\n"
            "- Include taxonomy shared motifs: include motifs shared across taxa.\n"
            "- Include motif strings: store full motif strings (can be large).\n"
            "- Skip unknown motifs: remove motifs containing 'unknown'.\n"
            "- Skip combined motifs: remove repeated/concatenated motifs.\n"
            "- Use longest only per species: keep only the longest record per species before analysis.\n"
            "\nCompare toggles\n"
            "- Compare using unique motifs: compare unique motifs instead of overlapping motifs.\n"
            "- Use window normalized: normalize motif counts by window count per sequence.\n"
            "- Require non-filter taxa: require at least one non-filter taxon to pass the min threshold.\n"
            "\nEquations used\n"
            "Length threshold filter: keep L if (mode == ge and L >= T) or (mode == le and L <= T).\n"
            "Longest-factor filter (after other filters): keep L if L >= L_max - (L_max / f).\n"
            "Example: if L_max=3000 and f=2, keep lengths >= 1500. If f=4, keep lengths >= 2250.\n"
            "AA composition percent (per species, per letter a): percent_a = count_a / length * 100.\n"
            "Total AA percent (per species, letters set A): total_percent = sum(count_a for a in A) / length * 100.\n"
            "SXn serine percent: percent_S = count_S / length * 100.\n"
            "SXn motif window normalization (compare option): freq = count / (L - k + 1) * 100.\n"
            "Average relative percent across taxa: avg = sum(percent_taxon) / N.\n"
            "\nEquation symbols\n"
            "L: sequence length (amino acids).\n"
            "T: length threshold.\n"
            "f: longest_factor.\n"
            "L_max: maximum length in the filtered set (per species or global scope).\n"
            "a: a specific amino acid letter (e.g., S, T, Y).\n"
            "A: set of amino acids included in the analysis.\n"
            "count_a: number of occurrences of amino acid a in a sequence.\n"
            "count: motif count for a given motif in a sequence.\n"
            "k: motif length (e.g., SXn motif length or unique motif length).\n"
            "N: number of taxa used in the average.\n"
        )

        text = tk.Text(frame, wrap="word")
        text.insert("1.0", help_text)
        text.configure(state="disabled")
        text.pack(fill=tk.BOTH, expand=True)
        text.bind("<MouseWheel>", lambda event: text.yview_scroll(int(-1 * (event.delta / 120)), "units"))

    def _run_in_thread(self, name: str, fn) -> None:
        threading.Thread(target=fn, name=name, daemon=True).start()

    def _run_ncbi_scraper(self) -> None:
        def task() -> None:
            try:
                from libraries.ncbi_protein_scraper_lib import run_ncbi_protein_scraper
            except Exception as exc:
                self._log(f"Import error: {exc}")
                return
            try:
                order_name = self.ncbi_order.get().strip()
                if not order_name:
                    self._log("NCBI scraper: order name is required.")
                    return
                family_names = _to_list(self.ncbi_family.get()) or None
                protein_terms = _to_list(self.ncbi_terms.get()) or None
                expected_raw = self.ncbi_expected.get("1.0", tk.END).strip()
                expected_types = json.loads(expected_raw) if expected_raw else None
                api_key = self.ncbi_api_key.get().strip() or None
                output_root = self.ncbi_output_root.get().strip() or "ncbi_sequences"
                sleep_time = float(self.ncbi_sleep.get().strip() or "0.2")
                path_max_length = _opt_int(self.ncbi_path_max.get())
                allow_long_paths = bool(self.ncbi_allow_long.get())

                run_ncbi_protein_scraper(
                    order_name=order_name,
                    family_names=family_names,
                    protein_terms=protein_terms,
                    expected_types=expected_types,
                    api_key=api_key,
                    output_root=output_root,
                    sleep_time=sleep_time,
                    path_max_length=path_max_length,
                    allow_long_paths=allow_long_paths,
                )
                self._log("NCBI scraper completed.")
            except Exception as exc:
                self._log(f"NCBI scraper failed: {exc}")

        self._run_in_thread("ncbi_scraper", task)

    def _run_taxonomy_tree(self) -> None:
        def task() -> None:
            try:
                from libraries.aminoacids_composition_analysis_lib import filter_records, load_records_from_root
                from libraries.generate_taxonomy_graph import (
                    RANK_ORDER,
                    RANK_ORDER_MAP,
                    assign_keys,
                    parse_species_index,
                    parse_taxonomy_tree,
                    render_tree,
                )
            except Exception as exc:
                self._log(f"Import error: {exc}")
                return
            try:
                data_root = Path(self.tax_data_root.get()).expanduser().resolve()
                tree_path = Path(self.tax_tree_path.get()).expanduser().resolve()
                index_path = Path(self.tax_index_path.get()).expanduser().resolve()
                out_path = Path(self.tax_out_path.get()).expanduser().resolve()
                output_format = self.tax_output_format.get().strip().lower() or "png"
                rank_from = self.tax_rank_from.get().strip()
                rank_to = self.tax_rank_to.get().strip()
                rank_from = None if not rank_from or rank_from.lower() == "(any)" else rank_from
                rank_to = None if not rank_to or rank_to.lower() == "(any)" else rank_to
                filters = self._read_filters(self.tax_filters)

                tree_data = json.loads(tree_path.read_text(encoding="utf-8"))
                tree = parse_taxonomy_tree(tree_data)
                index = parse_species_index(index_path)
                assign_keys(tree, index)

                def record_species(rec: dict) -> str:
                    for key in ("organism_name", "organism", "species"):
                        val = rec.get(key)
                        if val:
                            return str(val)
                    return ""

                if any([filters["protein_types"], filters["partial_full"], filters["length_range"], filters["length_threshold"], filters["longest_factor"]]):
                    records = load_records_from_root(data_root, skip_tree=True)
                    filtered = filter_records(
                        records,
                        taxonomy_terms=None,
                        protein_types=filters["protein_types"],
                        partial_full=filters["partial_full"],
                        length_range=filters["length_range"],
                        length_threshold=filters["length_threshold"],
                        length_mode=filters["length_mode"],
                        longest_factor=filters["longest_factor"],
                        longest_factor_scope=filters["longest_scope"],
                    )
                    species_set = {record_species(r).lower() for r in filtered if record_species(r)}

                    def filter_tree(node: dict) -> dict | None:
                        children = [c for c in (filter_tree(ch) for ch in node["children"]) if c]
                        is_species = node.get("rank") == "species"
                        if is_species and node["name"].lower() in species_set:
                            return {"name": node["name"], "rank": node.get("rank"), "children": children, "key": ""}
                        if not is_species and children:
                            return {"name": node["name"], "rank": node.get("rank"), "children": children, "key": ""}
                        return None

                    tree = filter_tree(tree) or tree

                def find_rank(root: dict, name: str) -> str | None:
                    key = name.strip().lower()
                    nodes = [root]
                    ranks = set()
                    while nodes:
                        node = nodes.pop()
                        if node.get("name", "").strip().lower() == key and node.get("rank"):
                            ranks.add(node["rank"])
                        nodes.extend(node.get("children", []))
                    if len(ranks) > 1:
                        raise ValueError(f"Ambiguous rank name '{name}'.")
                    return next(iter(ranks)) if ranks else None

                def resolve_range() -> tuple[int, int] | None:
                    if not rank_from and not rank_to:
                        return None
                    start = 0
                    end = len(RANK_ORDER) - 1
                    if rank_from:
                        key = rank_from.strip().lower()
                        if key not in RANK_ORDER_MAP:
                            key = (find_rank(tree, rank_from) or key)
                        start = RANK_ORDER_MAP[key]
                    if rank_to:
                        key = rank_to.strip().lower()
                        if key not in RANK_ORDER_MAP:
                            key = (find_rank(tree, rank_to) or key)
                        end = RANK_ORDER_MAP[key]
                    return min(start, end), max(start, end)

                rank_range = resolve_range()
                if rank_range:
                    def filter_rank(node: dict, lo: int, hi: int) -> list[dict]:
                        kept = []
                        for child in node["children"]:
                            kept.extend(filter_rank(child, lo, hi))
                        rank = node.get("rank")
                        idx = RANK_ORDER_MAP.get(rank) if isinstance(rank, str) else None
                        if idx is None:
                            return kept
                        if lo <= idx <= hi:
                            return [{"name": node["name"], "rank": rank, "children": kept, "key": ""}]
                        return kept

                    kept_children = []
                    for child in tree["children"]:
                        kept_children.extend(filter_rank(child, rank_range[0], rank_range[1]))
                    tree = {"name": tree["name"], "rank": tree.get("rank"), "children": kept_children, "key": tree.get("key", "")}

                render_tree(tree, index, out_path, fmt=output_format)
                image_path = out_path.with_suffix(f".{output_format}")
                self._log(f"Taxonomy tree written: {image_path}")
                self._show_image_popup(image_path, "Taxonomy Tree")
            except Exception as exc:
                self._log(f"Taxonomy tree failed: {exc}")

        self._run_in_thread("taxonomy", task)

    def _run_aa_composition(self) -> None:
        def task() -> None:
            try:
                from libraries.aminoacids_composition_analysis_lib import analyze_from_root
            except Exception as exc:
                self._log(f"Import error: {exc}")
                return
            try:
                start_ts = time.time()
                root = Path(self.aa_data_root.get()).expanduser().resolve()
                filters = self._read_filters(self.aa_filters)
                letters = self.aa_letters.get().strip()
                plots_dir = _opt_path(self.aa_plots_dir.get())
                tables_dir = _opt_path(self.aa_tables_dir.get())
                phylo_tree = _opt_path(self.aa_phylo_tree.get())
                json_path, md_path = analyze_from_root(
                    root=root,
                    letters=letters,
                    taxonomy_terms=filters["taxonomy_terms"],
                    protein_types=filters["protein_types"],
                    partial_full=filters["partial_full"],
                    length_range=filters["length_range"],
                    length_threshold=filters["length_threshold"],
                    length_mode=filters["length_mode"],
                    longest_factor=filters["longest_factor"],
                    longest_factor_scope=filters["longest_scope"],
                    plots_dir=plots_dir,
                    tables_dir=tables_dir,
                    phylo_tree_json=phylo_tree,
                )
                self._log(f"AA composition JSON: {json_path}")
                self._log(f"AA composition MD  : {md_path}")
                plots_root = plots_dir or (root / "plots")
                images = self._new_pngs(plots_root, start_ts)
                if images:
                    self._show_images_popups(images, "AA Composition Plot")
            except Exception as exc:
                self._log(f"AA composition failed: {exc}")

        self._run_in_thread("aa_composition", task)

    def _run_sxn(self) -> None:
        def task() -> None:
            try:
                from libraries.serine_sxn_analysis_lib import analyze_from_root
                from libraries.serine_sxn_plot_lib import (
                    category_out_dir,
                    ensure_out_dir,
                    load_summary,
                    plot_motif_counts_and_fraction,
                    plot_phylo_types,
                    plot_serine,
                    plot_total_sxn,
                    plot_x_composition,
                    filter_records,
                )
            except Exception as exc:
                self._log(f"Import error: {exc}")
                return
            try:
                start_ts = time.time()
                root = Path(self.sxn_data_root.get()).expanduser().resolve()
                filters = self._read_filters(self.sxn_filters)
                min_n = int(self.sxn_min_n.get().strip())
                max_n = int(self.sxn_max_n.get().strip())
                plots_dir = _opt_path(self.sxn_plots_dir.get())
                tables_dir = _opt_path(self.sxn_tables_dir.get())

                json_path, md_path = analyze_from_root(
                    root,
                    max_n=max_n,
                    min_n=min_n,
                    taxonomy_terms=filters["taxonomy_terms"],
                    protein_types=filters["protein_types"],
                    partial_full=filters["partial_full"],
                    length_range=filters["length_range"],
                    length_threshold=filters["length_threshold"],
                    length_mode=filters["length_mode"],
                    longest_factor=filters["longest_factor"],
                    longest_factor_scope=filters["longest_scope"],
                )
                self._log(f"SXn analysis JSON: {json_path}")
                self._log(f"SXn analysis MD  : {md_path}")

                summary = load_summary(json_path)
                records = summary.get("analyzed_records", [])
                filtered = filter_records(
                    records,
                    taxonomy_terms=filters["taxonomy_terms"],
                    protein_types=filters["protein_types"],
                    partial_full=filters["partial_full"],
                    length_range=filters["length_range"],
                    length_threshold=filters["length_threshold"],
                    length_mode=filters["length_mode"],
                    longest_factor=filters["longest_factor"],
                    longest_factor_scope=filters["longest_scope"],
                )
                if not filtered:
                    self._log("No records remain after filtering; skipping plots.")
                    return

                out_dir = plots_dir or ensure_out_dir(json_path)
                categories = ["full", "partial", "unknown"]
                present = {r.get("partial_full", "unknown") for r in filtered}
                for cat in categories:
                    if cat not in present:
                        continue
                    subset = [r for r in filtered if r.get("partial_full", "unknown") == cat]
                    if not subset:
                        continue
                    cat_out = category_out_dir(out_dir, cat, min_n, max_n)
                    sub_summary = {**summary, "analyzed_records": subset}
                    plot_serine(sub_summary, cat_out, min_n, max_n, suffix=cat, tables_dir=tables_dir or (root / "tables"))
                    plot_total_sxn(sub_summary, cat_out, min_n, max_n, suffix=cat, tables_dir=tables_dir or (root / "tables"))
                    plot_motif_counts_and_fraction(sub_summary, cat_out, min_n, max_n, suffix=cat, tables_dir=tables_dir or (root / "tables"))
                    plot_x_composition(sub_summary, cat_out, min_n, max_n, suffix=cat, tables_dir=tables_dir or (root / "tables"))
                    plot_phylo_types(sub_summary, cat_out, min_n, max_n, suffix=cat)

                self._log(f"SXn plots written under: {out_dir}")
                images = [p for p in self._new_pngs(out_dir, start_ts) if "phylo_types" not in p.name.lower()]
                if images:
                    self._show_images_popups(images, "SXn Plot")
            except Exception as exc:
                self._log(f"SXn analysis failed: {exc}")

        self._run_in_thread("sxn", task)

    def _run_motifs_compare(self) -> None:
        def task() -> None:
            try:
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
            except Exception as exc:
                self._log(f"Import error: {exc}")
                return
            try:
                start_ts = time.time()
                analyses = json.loads(self.motif_analyses.get("1.0", tk.END) or "[]")
                if not isinstance(analyses, list):
                    raise ValueError("Analyses JSON must be a list of dicts.")
                group_map = json.loads(self.motif_groups.get("1.0", tk.END) or "{}")
                group_legend = json.loads(self.motif_group_legend.get("1.0", tk.END) or "{}")
                group_label = self.motif_group_label.get().strip() or None

                filters = self._read_filters(self.motif_filters)
                min_motif_length = int(self.motif_min_len.get().strip())
                max_motif_length = int(self.motif_max_len.get().strip())
                unique_motif_length = int(self.motif_unique_len.get().strip())
                include_counts = bool(self.motif_include_counts.get())
                include_unique = bool(self.motif_include_unique.get())
                include_reflected = bool(self.motif_include_reflected.get())
                include_similar = bool(self.motif_include_similar.get())
                include_tax_shared = bool(self.motif_include_tax_shared.get())
                include_string = bool(self.motif_include_string.get())
                skip_unknown = bool(self.motif_skip_unknown.get())
                skip_combined = bool(self.motif_skip_combined.get())
                use_longest_only = bool(self.motif_use_longest_only.get())

                min_avg_relative = float(self.motif_min_avg_relative.get().strip())
                use_unique_compare = bool(self.motif_use_unique_compare.get())
                use_window_norm = bool(self.motif_use_window_norm.get())
                min_avg_relative_taxa = _to_list(self.motif_min_avg_taxa.get())
                min_avg_relative_min_taxa = int(self.motif_min_avg_taxa_count.get().strip())
                require_non_filter = bool(self.motif_require_non_filter.get())
                reuse_outputs = bool(self.motif_reuse_outputs.get())

                compare_out_dir = Path(self.motif_compare_out.get()).expanduser().resolve()
                compare_plots_dir = Path(self.motif_compare_plots.get()).expanduser().resolve()
                groups = group_map or DEFAULT_GROUPS

                comparisons: list[dict[str, object]] = []
                for entry in analyses:
                    if not isinstance(entry, dict):
                        raise ValueError("Each analysis entry must be a dict.")
                    data_root = Path(entry.get("data_root") or entry.get("outputs_root") or PROJECT_ROOT)
                    outputs_root = Path(entry.get("outputs_root") or data_root)
                    taxonomy_terms = entry.get("taxonomy_terms")
                    entry_protein_types = entry.get("protein_types", filters["protein_types"])
                    entry_partial_full = entry.get("partial_full", filters["partial_full"])
                    entry_length_range = entry.get("length_range", filters["length_range"])
                    entry_length_threshold = entry.get("length_threshold", filters["length_threshold"])
                    entry_length_mode = entry.get("length_mode", filters["length_mode"])
                    entry_longest_factor = entry.get("longest_factor", filters["longest_factor"])
                    entry_use_longest_only = entry.get("use_longest_only", use_longest_only)

                    if not reuse_outputs:
                        root = (PROJECT_ROOT / data_root).resolve() if not data_root.is_absolute() else data_root.resolve()
                        records = load_records_from_root(root, skip_tree=True)
                        filtered = filter_records(
                            records,
                            taxonomy_terms=taxonomy_terms,
                            protein_types=entry_protein_types,
                            partial_full=entry_partial_full,
                            length_range=entry_length_range,
                            length_threshold=entry_length_threshold,
                            length_mode=entry_length_mode,
                            longest_factor=entry_longest_factor,
                            longest_factor_scope="species",
                        )
                        filtered = dedupe_records_by_accession(filtered)
                        if entry_use_longest_only:
                            filtered = select_longest_records(filtered, scope="species")

                        summary = analyze_group_motifs(
                            filtered,
                            groups=groups,
                            include_motif_string=include_string,
                            include_motif_counts=include_counts,
                            min_motif_length=entry.get("min_motif_length", min_motif_length),
                            max_motif_length=entry.get("max_motif_length", max_motif_length),
                            skip_unknown_motifs=skip_unknown,
                            include_unique_motifs=include_unique,
                            unique_motif_length=unique_motif_length,
                            skip_combined_motifs=skip_combined,
                            include_reflected_sequence=include_reflected,
                            include_similar_motifs=include_similar,
                            include_taxonomy_shared_motifs=include_tax_shared,
                            taxonomy_terms=taxonomy_terms,
                            phylo_tree_path=entry.get("phylo_tree_path"),
                            root=root,
                        )

                        output_root = resolve_aa_motifs_output_dir(
                            outputs_root,
                            taxonomy_terms=taxonomy_terms,
                            partial_full=entry_partial_full,
                            protein_types=entry_protein_types,
                            length_range=entry_length_range,
                            length_threshold=entry_length_threshold,
                            length_mode=entry_length_mode,
                            longest_factor=entry_longest_factor,
                            min_motif_length=min_motif_length if include_counts else None,
                            max_motif_length=max_motif_length if include_counts else None,
                        )
                        write_json(output_root, summary, "aa_group_motifs_analysis.json")
                        write_md(output_root, summary, "aa_group_motifs_analysis.md")
                        if include_reflected:
                            write_reflected_json(output_root, summary, "reflected_sequence.json")
                            write_reflected_md(output_root, summary, "reflected_sequence.md")

                    comparisons.append(
                        {
                            "label": entry.get("label") or (taxonomy_terms[0] if taxonomy_terms else "Unknown"),
                            "outputs_root": outputs_root,
                            "taxonomy_terms": taxonomy_terms,
                            "protein_types": entry_protein_types,
                            "partial_full": entry_partial_full,
                            "length_range": entry_length_range,
                            "length_threshold": entry_length_threshold,
                            "length_mode": entry_length_mode,
                            "longest_factor": entry_longest_factor,
                            "min_motif_length": entry.get("min_motif_length", min_motif_length),
                            "max_motif_length": entry.get("max_motif_length", max_motif_length),
                        }
                    )

                compare_filters = {
                    "partial_full": filters["partial_full"],
                    "protein_types": filters["protein_types"],
                    "length_range": filters["length_range"],
                    "length_threshold": filters["length_threshold"],
                    "length_mode": filters["length_mode"],
                    "longest_factor": filters["longest_factor"],
                    "min_motif_length": min_motif_length,
                    "max_motif_length": max_motif_length,
                    "min_avg_relative": min_avg_relative,
                    "use_unique_motifs": use_unique_compare,
                    "use_window_normalized": use_window_norm,
                    "min_avg_relative_taxa": min_avg_relative_taxa,
                    "min_avg_relative_min_taxa": min_avg_relative_min_taxa,
                    "require_non_filter_taxa": require_non_filter,
                    "group_legend_override": group_legend or None,
                    "group_map_override": group_map or None,
                    "group_label": group_label,
                }

                compare_out_dir.mkdir(parents=True, exist_ok=True)
                run_compare(comparisons, compare_out_dir, compare_plots_dir, compare_filters, argv=None)
                self._log(f"AA motifs comparison written to: {compare_out_dir}")
                self._log(f"AA motifs plots written to: {compare_plots_dir}")
                images = self._new_pngs(compare_plots_dir, start_ts)
                if images:
                    self._show_images_popups(images, "AA Motifs Compare Plot")
            except Exception as exc:
                self._log(f"AA motifs compare failed: {exc}")

        self._run_in_thread("motifs", task)


def main() -> None:
    app = App()
    app.mainloop()


if __name__ == "__main__":
    main()
