"""Library for comparing amino-acid group motif frequencies across taxonomies."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Iterable, List

from libraries.amino_acid_group_motifs import resolve_aa_motifs_output_dir


def _resolve_analysis_json(path: Path) -> Path:
    if path.is_file():
        return path
    if not path.exists():
        raise FileNotFoundError(f"Path not found: {path}")

    direct = path / "aa_group_motifs_analysis.json"
    if direct.exists():
        return direct

    matches = list(path.rglob("aa_group_motifs_analysis.json"))
    if not matches:
        raise FileNotFoundError(f"No aa_group_motifs_analysis.json found under {path}")
    if len(matches) > 1:
        preview = "\n".join(str(p) for p in matches[:10])
        raise ValueError(
            "Multiple aa_group_motifs_analysis.json files found. "
            "Please pass the specific file or a narrower directory.\n"
            f"{preview}"
        )
    return matches[0]


def _slugify_label(text: str) -> str:
    return "".join(ch if ch.isalnum() and ord(ch) < 128 else "_" for ch in text).strip("_") or "comparison"

def _min_avg_label(value: float) -> str:
    text = f"{value:g}"
    safe = "".join(ch if ch.isalnum() else "_" for ch in text).strip("_")
    return f"MIN_AVG_{safe or '0'}"


def _label_from_path(path: Path) -> str:
    parts = [p for p in path.parts]
    if "aa_motifs" in parts:
        idx = parts.index("aa_motifs")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    if path.is_file():
        return path.parent.name
    return path.name


def _load_summary(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve_from_outputs_root(
    outputs_root: Path,
    taxonomy_terms: Iterable[str] | None,
    partial_full: str | None,
    protein_types: Iterable[str] | None,
    length_range: tuple[int, int] | None,
    length_threshold: int | None,
    length_mode: str,
    longest_factor: float | None,
    min_motif_length: int | None,
    max_motif_length: int | None,
) -> Path:
    out_dir = resolve_aa_motifs_output_dir(
        outputs_root,
        taxonomy_terms=taxonomy_terms,
        partial_full=partial_full,
        protein_types=protein_types,
        length_range=length_range,
        length_threshold=length_threshold,
        length_mode=length_mode,
        longest_factor=longest_factor,
        min_motif_length=min_motif_length,
        max_motif_length=max_motif_length,
    )
    return out_dir / "aa_group_motifs_analysis.json"


def _taxonomy_shared_node(
    summary: Dict[str, object],
    taxonomy_terms: Iterable[str] | None,
) -> Dict[str, object] | None:
    if not taxonomy_terms:
        return None
    target = str(list(taxonomy_terms)[0]).lower()
    shared = (summary.get("taxonomy_shared_motifs") or {}).get("nodes") or []
    for node in shared:
        name = str(node.get("name", "") or "").lower()
        if name == target:
            return node
    return None


def _aggregate_motif_stats(
    records: List[Dict[str, object]],
    skip_missing: bool,
    use_unique_motifs: bool,
    use_window_normalized: bool,
    unique_motif_length: int | None,
) -> Dict[str, Dict[str, float | int]]:
    sums: Dict[str, float] = {}
    with_motif: Dict[str, int] = {}
    for rec in records:
        freqs: Dict[str, float] = {}
        seq_len = int(rec.get("length", 0) or 0)
        if use_window_normalized and seq_len > 0:
            if use_unique_motifs:
                counts = rec.get("unique_motifs") or {}
                motif_len = int(unique_motif_length or 0)
                denom = seq_len - motif_len + 1 if motif_len > 0 else 0
                if denom > 0:
                    freqs = {m: (c / denom * 100.0) for m, c in counts.items()}
            else:
                counts = rec.get("motif_counts") or {}
                lengths = rec.get("motif_lengths") or {}
                for motif, count in counts.items():
                    motif_len = int(lengths.get(motif, 0) or 0)
                    denom = seq_len - motif_len + 1 if motif_len > 0 else 0
                    if denom > 0:
                        freqs[motif] = count / denom * 100.0
        if not freqs:
            key = "unique_motif_frequencies" if use_unique_motifs else "motif_frequencies"
            freqs = rec.get(key) or {}
        for motif, val in freqs.items():
            sums[motif] = sums.get(motif, 0.0) + float(val or 0.0)
            with_motif[motif] = with_motif.get(motif, 0) + 1

    total_records = len(records)
    stats: Dict[str, Dict[str, float | int]] = {}
    for motif, total in sums.items():
        num_with = with_motif.get(motif, 0)
        denom = num_with if skip_missing else total_records
        avg = total / denom if denom else 0.0
        stats[motif] = {
            "avg_relative_percent": avg,
            "num_records": denom,
            "num_records_with_motif": num_with,
        }
    return stats


def _write_markdown(
    out_path: Path,
    motifs: List[str],
    rows: List[Dict[str, object]],
    skip_missing: bool,
    group_legend: Dict[str, str] | None,
    min_avg_relative: float,
    group_map: Dict[str, List[str]] | None,
) -> None:
    lines: List[str] = []
    lines.append("# Motif comparison")
    lines.append("")
    if skip_missing:
        lines.append("Averages exclude records where the motif is missing.")
    else:
        lines.append("Averages treat missing motifs as 0%.")
    if group_legend:
        lines.append("")
        lines.append("## Group legend")
        lines.append("")
        lines.append("| Group | Description |")
        lines.append("| --- | --- |")
        for key, desc in group_legend.items():
            lines.append(f"| {key} | {desc} |")
    if group_map:
        lines.append("")
        lines.append("## Group classification")
        lines.append("")
        lines.append("| Group | Amino acids |")
        lines.append("| --- | --- |")
        for key, acids in group_map.items():
            acids_str = ", ".join(acids)
            lines.append(f"| {key} | {acids_str} |")
    lines.append("")
    for motif in motifs:
        lines.append(f"## Motif: {motif}")
        lines.append("")
        lines.append("| Taxonomy | Avg Relative (%) | Records | Records with motif |")
        lines.append("| --- | --- | --- | --- |")
        for row in rows:
            m = row["motifs"][motif]
            if m["avg_relative_percent"] < min_avg_relative:
                continue
            lines.append(
                f"| {row['label']} | {m['avg_relative_percent']:.4f} | "
                f"{m['num_records']} | {m['num_records_with_motif']} |"
            )
        lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _plot_all_motifs(
    out_dir: Path,
    motifs: List[str],
    rows: List[Dict[str, object]],
    title: str | None,
    group_legend: Dict[str, str] | None,
    min_avg_relative: float,
    group_map: Dict[str, List[str]] | None,
) -> Path:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError("matplotlib is required for plotting.") from exc

    labels = [row["label"] for row in rows]
    num_groups = len(labels)
    num_motifs = len(motifs)
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * num_motifs), 4 + 0.3 * num_groups))

    bar_width = 0.8 / max(1, num_groups)
    x_positions = list(range(num_motifs))
    for i, row in enumerate(rows):
        offset = (i - (num_groups - 1) / 2) * bar_width
        values = [
            (
                row["motifs"][motif]["avg_relative_percent"]
                if row["motifs"][motif]["avg_relative_percent"] >= min_avg_relative
                else 0.0
            )
            for motif in motifs
        ]
        bars = ax.bar([x + offset for x in x_positions], values, width=bar_width, label=row["label"])
        for bar, val in zip(bars, values):
            if val <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(1, val * 0.01),
                f"{val:.2f}",
                ha="center",
                va="bottom",
                fontsize=6,
                rotation=90,
            )

    ax.set_ylabel("Avg Relative (%)")
    ax.set_title(title or "Motif comparison")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(motifs, rotation=90, fontsize=6)
    if num_groups > 1:
        ax.legend(fontsize=8, loc="upper right")
    if group_legend:
        legend_lines = [f"{k} = {v}" for k, v in group_legend.items()]
        ax.text(
            0.01,
            0.99,
            " | ".join(legend_lines),
            fontsize=7,
            va="top",
            ha="left",
            transform=ax.transAxes,
        )
    if group_map:
        class_lines = [f"{k}: {', '.join(v)}." for k, v in group_map.items()]
        ax.text(
            0.01,
            0.92,
            "\n".join(class_lines),
            fontsize=7,
            va="top",
            ha="left",
            transform=ax.transAxes,
        )
    fig.tight_layout()

    out_path = out_dir / "all_motifs_comparison.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def run_compare(
    comparisons: List[Dict[str, object]],
    default_out_dir: Path,
    plots_out_dir: Path,
    filters: Dict[str, object],
    argv: List[str] | None = None,
) -> int:
    tic = time.perf_counter()
    parser = argparse.ArgumentParser(
        description="Compare group motif frequencies across aa_group_motifs_analysis.json files."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="Paths to aa_group_motifs_analysis.json files or directories containing them.",
    )
    parser.add_argument(
        "--motifs",
        nargs="+",
        default=None,
        help="Motifs to compare (Greek group symbols allowed). If omitted, compare all motifs found.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(default_out_dir),
        help="Directory to write comparison outputs.",
    )
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Exclude records where the motif is missing from the average.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title prefix.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plotting (faster when comparing many motifs).",
    )
    args = parser.parse_args(argv)

    # Filters
    partial_full = filters.get("partial_full")
    protein_types = filters.get("protein_types")
    length_range = filters.get("length_range")
    length_threshold = filters.get("length_threshold")
    length_mode = filters.get("length_mode")
    longest_factor = filters.get("longest_factor")
    min_motif_length = filters.get("min_motif_length")
    max_motif_length = filters.get("max_motif_length")
    min_avg_relative = float(filters.get("min_avg_relative", 0.0))
    use_unique_motifs = bool(filters.get("use_unique_motifs", False))
    use_window_normalized = bool(filters.get("use_window_normalized", False))
    min_avg_relative_taxa = list(filters.get("min_avg_relative_taxa") or [])
    min_avg_relative_min_taxa = int(filters.get("min_avg_relative_min_taxa", 2))
    require_non_filter_taxa = bool(filters.get("require_non_filter_taxa", False))

    analysis_paths: List[Path] = []
    labels: List[str] = []
    if args.inputs:
        resolved_inputs = [Path(p).expanduser().resolve() for p in args.inputs]
        analysis_paths = [_resolve_analysis_json(p) for p in resolved_inputs]
        labels = [_label_from_path(p) for p in analysis_paths]
    else:
        for entry in comparisons:
            out_root = Path(entry["outputs_root"]).expanduser().resolve()
            json_path = _resolve_from_outputs_root(
                out_root,
                taxonomy_terms=entry.get("taxonomy_terms"),
                partial_full=entry.get("partial_full", partial_full),
                protein_types=entry.get("protein_types", protein_types),
                length_range=entry.get("length_range", length_range),
                length_threshold=entry.get("length_threshold", length_threshold),
                length_mode=entry.get("length_mode", length_mode),
                longest_factor=entry.get("longest_factor", longest_factor),
                min_motif_length=entry.get("min_motif_length", min_motif_length),
                max_motif_length=entry.get("max_motif_length", max_motif_length),
            )
            analysis_paths.append(json_path)
            labels.append(str(entry.get("label") or _label_from_path(json_path)))

    all_motifs: set[str] = set()
    rows: List[Dict[str, object]] = []
    group_legend: Dict[str, str] | None = None
    group_map: Dict[str, List[str]] | None = None
    for idx, path in enumerate(analysis_paths):
        if not path.exists():
            raise FileNotFoundError(f"Missing analysis json: {path}")
        summary = _load_summary(path)
        if group_legend is None:
            group_legend = summary.get("group_legend") or None
        if group_map is None:
            group_map_raw = summary.get("groups") or None
            if isinstance(group_map_raw, dict):
                group_map = {str(k): [str(x) for x in v] for k, v in group_map_raw.items()}
        records = summary.get("analyzed_records") or []
        unique_motif_length = summary.get("unique_motif_length")
        label = labels[idx] if idx < len(labels) else _label_from_path(path)

        taxonomy_terms = None
        if not args.inputs and idx < len(comparisons):
            taxonomy_terms = comparisons[idx].get("taxonomy_terms")
        shared_node = _taxonomy_shared_node(summary, taxonomy_terms)
        if shared_node:
            shared_stats = shared_node.get("shared_motifs") or {}
            stats = {
                motif: {
                    "avg_relative_percent": float(val or 0.0),
                    "num_records": int(shared_node.get("num_records", 0) or 0),
                    "num_records_with_motif": int(shared_node.get("num_records", 0) or 0),
                }
                for motif, val in shared_stats.items()
            }
            total_records = int(shared_node.get("num_records", 0) or 0)
        else:
            stats = _aggregate_motif_stats(
                records,
                args.skip_missing,
                use_unique_motifs,
                use_window_normalized,
                int(unique_motif_length) if isinstance(unique_motif_length, int) else None,
            )
            total_records = len(records)

        motif_keys = set(stats.keys())
        all_motifs.update(motif_keys)
        entry: Dict[str, object] = {
            "label": label,
            "source": str(path),
            "total_records": total_records,
            "motifs": stats,
        }
        rows.append(entry)

    if args.motifs:
        motifs = args.motifs
    else:
        motifs = sorted(all_motifs)
    if not motifs:
        raise ValueError("No motifs found to compare.")

    for row in rows:
        for motif in motifs:
            if motif in row["motifs"]:
                continue
            denom = 0 if args.skip_missing else int(row.get("total_records", 0))
            row["motifs"][motif] = {
                "avg_relative_percent": 0.0,
                "num_records": denom,
                "num_records_with_motif": 0,
            }

    for row in rows:
        row.pop("total_records", None)

    # Optional legend overrides from filters.
    legend_override = filters.get("group_legend_override")
    if isinstance(legend_override, dict):
        group_legend = {str(k): str(v) for k, v in legend_override.items()}
    group_map_override = filters.get("group_map_override")
    if isinstance(group_map_override, dict):
        group_map = {str(k): [str(x) for x in v] for k, v in group_map_override.items()}

    if min_avg_relative_taxa:
        filter_labels = {str(name) for name in min_avg_relative_taxa}
        motifs = [
            motif
            for motif in motifs
            if (
                sum(
                    1
                    for row in rows
                    if row["motifs"][motif]["avg_relative_percent"] >= min_avg_relative
                )
                >= min_avg_relative_min_taxa
            )
            and any(
                row["label"] in filter_labels
                and row["motifs"][motif]["avg_relative_percent"] >= min_avg_relative
                for row in rows
            )
            and (
                not require_non_filter_taxa
                or any(
                    row["label"] not in filter_labels
                    and row["motifs"][motif]["avg_relative_percent"] >= min_avg_relative
                    for row in rows
                )
            )
        ]
    else:
        motifs = [
            motif
            for motif in motifs
            if sum(
                1
                for row in rows
                if row["motifs"][motif]["avg_relative_percent"] >= min_avg_relative
            )
            >= min_avg_relative_min_taxa
        ]
    if not motifs:
        raise ValueError("No motifs remain after applying MIN_AVG_RELATIVE.")

    for row in rows:
        row["motifs"] = {motif: row["motifs"][motif] for motif in motifs}

    group_label = filters.get("group_label")
    group_part = _slugify_label(str(group_label)) if group_label else None
    comparison_label = "_".join(_slugify_label(label) for label in labels) or "comparison"
    min_avg_part = _min_avg_label(min_avg_relative)
    out_dir = Path(args.out_dir).expanduser().resolve()
    if group_part:
        out_dir = out_dir / group_part
    out_dir = out_dir / comparison_label
    out_dir = out_dir / min_avg_part
    out_dir.mkdir(parents=True, exist_ok=True)

    comparison = {
        "motifs": motifs,
        "entries": rows,
        "skip_missing": args.skip_missing,
    }
    json_path = out_dir / "aa_group_motif_comparison.json"
    json_path.write_text(json.dumps(comparison, indent=2), encoding="utf-8")

    md_path = out_dir / "aa_group_motif_comparison.md"
    _write_markdown(
        md_path,
        motifs,
        rows,
        args.skip_missing,
        group_legend,
        min_avg_relative,
        group_map,
    )

    plot_paths: List[str] = []
    if not args.no_plot:
        plots_dir = Path(plots_out_dir).expanduser().resolve()
        if group_part:
            plots_dir = plots_dir / group_part
        plots_dir = plots_dir / comparison_label
        plots_dir = plots_dir / min_avg_part
        plots_dir.mkdir(parents=True, exist_ok=True)
        plot_path = _plot_all_motifs(
            plots_dir,
            motifs,
            rows,
            args.title,
            group_legend,
            min_avg_relative,
            group_map,
        )
        plot_paths.append(str(plot_path))

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    for p in plot_paths:
        print(f"Wrote {p}")
    toc = time.perf_counter()
    print(f"Elapsed: {toc - tic:0.3f} seconds")
    return 0
