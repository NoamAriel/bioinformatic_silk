import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


def _slugify(text: str) -> str:
    """Make a filesystem-friendly slug."""
    return "".join(ch if ch.isalnum() else "_" for ch in text).strip("_") or "none"


def greedy_sxn_runs(seq: str, max_n: int, min_n: int) -> Dict[str, Any]:
    """Find non-overlapping [SX]_n runs using a longest-first strategy."""
    seq = seq.upper()
    working = seq
    motif_runs: Dict[int, int] = {}
    x_residue_counts: Dict[int, Dict[str, int]] = {}
    total_runs = 0
    total_motif_residues = 0

    for n in range(max_n, min_n - 1, -1):
        pattern = re.compile(rf"(S[A-Z]){{{n}}}")
        matches = list(pattern.finditer(working))
        if not matches:
            continue

        motif_runs[n] = len(matches)
        total_runs += len(matches)
        total_motif_residues += len(matches) * (n * 2)

        x_counts: Dict[str, int] = defaultdict(int)
        for m in matches:
            run = m.group(0)
            for i in range(1, len(run), 2):
                x_counts[run[i]] += 1
        x_residue_counts[n] = dict(x_counts)

        filler = "$" * (n * 2)
        temp = working
        for m in matches:
            s, e = m.span()
            temp = temp[:s] + filler + temp[e:]
        working = temp

    fraction = (total_motif_residues / len(seq) * 100) if seq else 0.0
    return {
        "motif_runs": motif_runs,
        "x_residue_counts": x_residue_counts,
        "total_runs": total_runs,
        "total_motif_residues": total_motif_residues,
        "fraction_motif_residues": fraction,
    }


def serine_stats(seq: str) -> Tuple[int, float]:
    seq_up = seq.upper()
    count = seq_up.count("S")
    frac = (count / len(seq_up) * 100) if seq_up else 0.0
    return count, frac


def analyze_records(
    records: List[Dict[str, Any]],
    max_n: int = 50,
    min_n: int = 2,
) -> Dict[str, Any]:
    """
    Analyze serine composition and [SX]_n motifs for each record.
    Expects each record to have: origin_sequence, accession, organism_name,
    taxonomy_from_araneae or taxonomy_full, partial_full, and type.
    """
    max_n = max(1, min(100, max_n))
    min_n = max(1, min(max_n, min_n))

    analyzed: List[Dict[str, Any]] = []
    type_counts: Dict[str, int] = defaultdict(int)
    type_pf_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    species_pf_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    species_info: Dict[str, Dict[str, Any]] = {}

    for rec in records:
        seq = rec.get("origin_sequence", "") or ""
        if not seq:
            continue

        s_count, s_fraction = serine_stats(seq)
        motif_metrics = greedy_sxn_runs(seq, max_n=max_n, min_n=min_n)

        taxonomy_lineage = (
            rec.get("taxonomy_from_order")
            or rec.get("taxonomy_from_araneae")
            or rec.get("taxonomy_full")
            or []
        )
        suborder = taxonomy_lineage[1] if len(taxonomy_lineage) > 1 else "Unknown suborder"
        superfamily = taxonomy_lineage[2] if len(taxonomy_lineage) > 2 else "Unknown superfamily"
        family = taxonomy_lineage[3] if len(taxonomy_lineage) > 3 else "Unknown family"

        result = {
            "accession": rec.get("accession", ""),
            "title": rec.get("title", ""),
            "organism": rec.get("organism_name", "Unknown organism"),
            "taxonomy_from_order": rec.get("taxonomy_from_order", []),
            "taxonomy_from_araneae": rec.get("taxonomy_from_araneae", []),
            "taxonomy_full": rec.get("taxonomy_full", []),
            "suborder": suborder,
            "superfamily": superfamily,
            "family": family,
            "partial_full": rec.get("partial_full", "") or "unknown",
            "type": rec.get("type", "") or "unknown",
            "length": len(seq),
            "serine_count": s_count,
            "serine_fraction": s_fraction,
            **motif_metrics,
        }
        analyzed.append(result)
        type_counts[result["type"]] += 1
        pf_key = result["partial_full"] or "unknown"
        type_pf_counts[result["type"]][pf_key] += 1
        species_pf_counts[result["organism"]][pf_key] += 1

        info = species_info.setdefault(
            result["organism"],
            {"taxonomy": result["taxonomy_from_araneae"] or result["taxonomy_full"], "types": set()},
        )
        info["types"].add(result["type"])

    for v in species_info.values():
        v["types"] = sorted(v["types"])

    return {
        "analyzed_records": analyzed,
        "type_counts": dict(type_counts),
        "type_partial_full_counts": {t: dict(pf) for t, pf in type_pf_counts.items()},
        "species_partial_full_counts": {s: dict(pf) for s, pf in species_pf_counts.items()},
        "species_info": species_info,
        "max_n": max_n,
        "min_n": min_n,
    }


def write_json(root: Path, data: Dict[str, Any], out_name: str) -> Path:
    out_path = root / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return out_path


def write_md(root: Path, summary: Dict[str, Any], json_path: Path, out_name: str) -> Path:
    """Render the markdown summary similar to the fibroin SXn reports."""
    lines: List[str] = []
    lines.append("# Serine + [SX]_n Motif Analysis (Taxonomic Grouping)")
    lines.append("")
    lines.append(
        "This analysis identifies serine composition and continuous runs of the **SX** unit "
        "($\\mathbf{[SX]_n}$) in sequences. Results are grouped by taxonomic hierarchy "
        "(Suborder, Superfamily, and Family)."
    )
    lines.append(f"- JSON data: `{json_path.name}`")
    lines.append(f"- n range: {summary.get('min_n', 0)}..{summary.get('max_n', 0)}")
    lines.append("---")

    all_records = summary.get("analyzed_records", [])
    if not all_records:
        out_path = root / out_name
        out_path.write_text("\n".join(lines + ["", "_No records found._"]), encoding="utf-8")
        return out_path

    all_n_values = set()
    all_x_residues = set()
    for record in all_records:
        all_n_values.update(record.get("motif_runs", {}).keys())
        for x_counts in record.get("x_residue_counts", {}).values():
            all_x_residues.update(x_counts.keys())

    sorted_n_values = sorted(all_n_values, reverse=True)
    sorted_x_residues = sorted(all_x_residues)

    lines.append("## Detailed Run Counts by Taxonomic Group (Overall)")
    records_sorted = sorted(
        all_records,
        key=lambda r: (
            r.get("suborder", "Unknown suborder"),
            r.get("superfamily", "Unknown superfamily"),
            r.get("family", "Unknown family"),
            r.get("organism", ""),
            r.get("type", ""),
            r.get("partial_full", ""),
        ),
    )

    current_suborder = None
    current_superfamily = None
    current_family = None

    header_cols = [
        "Organism",
        "Type",
        "Partial/Full",
        "Accession ID",
        "Length",
        "Serine Count",
        "Serine Fraction",
        "Total Runs",
        "Residues in Motifs",
        "Residues Fraction",
    ]
    n_cols = [f"n={n} Count" for n in sorted_n_values]
    header_line = "| " + " | ".join(header_cols + n_cols) + " |"
    align_line = "| " + " :--- |" * len(header_cols) + " :---: |" * len(n_cols)

    for rec in records_sorted:
        suborder = rec.get("suborder", "Unknown suborder")
        superfamily = rec.get("superfamily", "Unknown superfamily")
        family = rec.get("family", "Unknown family")

        if suborder != current_suborder:
            current_suborder = suborder
            lines.append(f"\n# Suborder: {suborder}\n")
            lines.append("===")
            current_superfamily = None

        if superfamily != current_superfamily:
            current_superfamily = superfamily
            lines.append(f"\n## Superfamily: {superfamily}\n")
            lines.append("---")
            current_family = None

        if family != current_family:
            current_family = family
            lines.append(f"\n### Family: {family}")
            lines.append(header_line)
            lines.append(align_line)

        row = [
            rec.get("organism", "Unknown"),
            rec.get("type", "Unknown"),
            rec.get("partial_full", "unknown"),
            f"`{rec.get('accession', '')}`",
            str(rec.get("length", 0)),
            str(rec.get("serine_count", 0)),
            f"**{rec.get('serine_fraction', 0):.2f}%**",
            f"**{rec.get('total_runs', 0)}**",
            str(rec.get("total_motif_residues", 0)),
            f"**{rec.get('fraction_motif_residues', 0):.2f}%**",
        ]
        for n in sorted_n_values:
            row.append(str(rec.get("motif_runs", {}).get(n, 0)))
        lines.append("| " + " | ".join(row) + " |")

    lines.append("\n## X-Residue Composition Analysis by Run Length ($\\mathbf{[SX]_n}$)")
    lines.append(
        "This table shows the total count for each amino acid (X) found in the "
        "$\\mathbf{X}$ position across all runs of a specific length $\\mathbf{n}$."
    )

    current_suborder = None
    current_superfamily = None
    current_family = None

    x_header_base = ["Organism/Type/ID", "n (Run Length)", "Total X Residues"]
    x_cols = [f"X={x}" for x in sorted_x_residues]
    x_header_line = "| " + " | ".join(x_header_base + x_cols) + " |"
    x_align_line = "| :--- | :---: | :---: |" + " :---: |" * len(x_cols)

    for rec in records_sorted:
        suborder = rec.get("suborder", "Unknown suborder")
        superfamily = rec.get("superfamily", "Unknown superfamily")
        family = rec.get("family", "Unknown family")

        if suborder != current_suborder:
            current_suborder = suborder
            lines.append(f"\n# Suborder: {suborder}\n")
            lines.append("===")
            current_superfamily = None

        if superfamily != current_superfamily:
            current_superfamily = superfamily
            lines.append(f"\n## Superfamily: {superfamily}\n")
            lines.append("---")
            current_family = None

        if family != current_family:
            current_family = family
            lines.append(f"\n### Family: {family} (X-Residue Detail)")

        lines.append(x_header_line)
        lines.append(x_align_line)

        x_counts_by_n = rec.get("x_residue_counts", {})
        for n in sorted_n_values:
            x_counts = x_counts_by_n.get(n)
            if not x_counts:
                continue
            total_x = sum(x_counts.values())
            row = [
                f"{rec.get('organism', 'Unknown')} ({rec.get('type', '')}) - `{rec.get('accession', '')}`",
                str(n),
                f"**{total_x}**",
            ]
            for aa in sorted_x_residues:
                row.append(str(x_counts.get(aa, 0)))
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    lines.append("---")
    lines.append("**Key Metrics:**")
    lines.append(r"* **$n$ (Run Count):** Number of times the base $\mathbf{SX}$ unit is repeated (e.g., $SXSXSX$ has $n=3$).")
    lines.append(r"* **Total Runs:** Sum of all continuous $\mathbf{[SX]_n}$ segments found in the sequence.")
    lines.append(r"* **Residues in Motifs:** Total amino acids covered by all $\mathbf{[SX]_n}$ runs.")
    lines.append(r"* **Residues Fraction:** Percentage of the protein length covered by $\mathbf{[SX]_n}$ runs.")
    lines.append(r"* **X-Residue Composition:** Counts of amino acids substituting for $\mathbf{X}$ in $\mathbf{SX}$.")

    out_path = root / out_name
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def analyze_and_save(
    records: List[Dict[str, Any]],
    root: Path,
    json_out: str = "serine_sxn_analysis.json",
    md_out: str = "serine_sxn_analysis.md",
    max_n: int = 50,
    min_n: int = 2,
) -> Tuple[Path, Path]:
    summary = analyze_records(records, max_n=max_n, min_n=min_n)
    json_path = write_json(root, summary, json_out)
    md_path = write_md(root, summary, json_path, md_out)
    return json_path, md_path


def load_records_from_root(root: Path, skip_tree: bool = True) -> List[Dict[str, Any]]:
    """
    Recursively load all JSON files under `root` that look like record files
    (must contain origin_sequence). Skips phylo_tree.json by default.
    """
    records: List[Dict[str, Any]] = []
    for json_path in root.rglob("*.json"):
        if skip_tree and json_path.name.lower() == "phylo_tree.json":
            continue
        try:
            data = json.loads(json_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        candidates = data if isinstance(data, list) else [data]
        for rec in candidates:
            if not isinstance(rec, dict) or not rec.get("origin_sequence"):
                continue

            # Infer partial/type from path if missing: .../<partial|full>/<type>/<file.json>
            parts = json_path.parent.parts
            if len(parts) >= 2:
                partial_candidate = parts[-2]
                type_candidate = parts[-1]
                if partial_candidate in ("partial", "full") and not rec.get("partial_full"):
                    rec["partial_full"] = partial_candidate
                if partial_candidate in ("partial", "full") and not rec.get("type"):
                    rec["type"] = type_candidate

            records.append(rec)
    return records


def dedupe_records_by_accession(records: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Return records with duplicate accession IDs removed (keeps first seen)."""
    seen: set[str] = set()
    deduped: List[Dict[str, Any]] = []
    for rec in records:
        accession = str(rec.get("accession", "") or "").strip()
        if accession:
            if accession in seen:
                continue
            seen.add(accession)
        deduped.append(rec)
    return deduped


def _record_length_for_filters(rec: Dict[str, Any]) -> int:
    seq_len = rec.get("sequence_length")
    if isinstance(seq_len, (int, float)) and seq_len > 0:
        return int(seq_len)
    if isinstance(seq_len, str) and seq_len.strip().isdigit():
        return int(seq_len.strip())
    seq = rec.get("origin_sequence", "") or ""
    return len(seq)


def filter_records(
    records: Sequence[Dict[str, Any]],
    taxonomy_terms: Optional[Iterable[str]] = None,
    protein_types: Optional[Iterable[str]] = None,
    partial_full: Optional[str] = None,
    length_range: Optional[Tuple[int, int]] = None,
    length_threshold: Optional[int] = None,
    length_mode: str = "ge",
    longest_factor: Optional[float] = None,
    longest_factor_scope: str = "species",
) -> List[Dict[str, Any]]:
    filtered = list(records)

    if taxonomy_terms:
        terms_lower = {str(t).strip().lower() for t in taxonomy_terms if str(t).strip()}
        tmp: List[Dict[str, Any]] = []
        for rec in filtered:
            lineage = (
                rec.get("taxonomy_from_order")
                or rec.get("taxonomy_from_araneae")
                or rec.get("taxonomy_full")
                or []
            )
            lineage_lower = {str(x).lower() for x in lineage}
            if lineage_lower.intersection(terms_lower):
                tmp.append(rec)
        filtered = tmp

    if protein_types:
        type_terms = [str(t).strip().lower() for t in protein_types if str(t).strip()]
        tmp = []
        for rec in filtered:
            rec_type = str(rec.get("type", "")).strip().lower()
            if any(term in rec_type for term in type_terms):
                tmp.append(rec)
        filtered = tmp

    if partial_full is not None:
        pf_target = partial_full.lower()
        filtered = [rec for rec in filtered if str(rec.get("partial_full", "")).lower() == pf_target]

    if length_range is not None:
        min_len, max_len = length_range
        tmp = []
        for rec in filtered:
            L = _record_length_for_filters(rec)
            if (min_len is None or L >= min_len) and (max_len is None or L <= max_len):
                tmp.append(rec)
        filtered = tmp

    if length_threshold is not None:
        mode = length_mode.lower()
        tmp = []
        for rec in filtered:
            L = _record_length_for_filters(rec)
            if mode == "ge" and L >= length_threshold:
                tmp.append(rec)
            elif mode == "le" and L <= length_threshold:
                tmp.append(rec)
        filtered = tmp

    if longest_factor is not None and longest_factor > 1 and filtered:
        scope = (longest_factor_scope or "species").lower()
        if scope not in {"species", "global"}:
            raise ValueError("longest_factor_scope must be 'species' or 'global'.")

        def _org_name(rec: Dict[str, Any]) -> str:
            return rec.get("organism_name") or rec.get("organism") or "Unknown"

        if scope == "global":
            max_len = max(_record_length_for_filters(rec) for rec in filtered)
            cutoff = max_len - (max_len / float(longest_factor))
            filtered = [rec for rec in filtered if _record_length_for_filters(rec) >= cutoff]
        else:
            max_by_org: Dict[str, int] = {}
            for rec in filtered:
                org = _org_name(rec)
                L = _record_length_for_filters(rec)
                if L > max_by_org.get(org, 0):
                    max_by_org[org] = L
            tmp = []
            for rec in filtered:
                org = _org_name(rec)
                L = _record_length_for_filters(rec)
                max_len = max_by_org.get(org, 0)
                cutoff = max_len - (max_len / float(longest_factor))
                if L >= cutoff:
                    tmp.append(rec)
            filtered = tmp

    return filtered


def _analysis_output_dir(
    root: Path,
    taxonomy_terms: Optional[Iterable[str]] = None,
    partial_full: Optional[str] = None,
    protein_types: Optional[Iterable[str]] = None,
    length_range: Optional[Tuple[int, int]] = None,
    length_threshold: Optional[int] = None,
    length_mode: str = "ge",
    longest_factor: Optional[float] = None,
) -> Path:
    base = root / "data_analysis" / "sxn_analysis"

    if taxonomy_terms:
        tax_part = "__".join(_slugify(str(t)) for t in taxonomy_terms)
    else:
        tax_part = "all_taxonomy"

    pf_part = _slugify(partial_full or "all_partial_full")

    if protein_types:
        type_part = "__".join(_slugify(str(t)) for t in protein_types)
    else:
        type_part = "all_types"

    length_parts: List[str] = []
    if length_range is not None:
        length_parts.append(f"{length_range[0]}-{length_range[1]}")
    if length_threshold is not None:
        length_parts.append(f"{length_mode}_{length_threshold}")
    if longest_factor is not None:
        length_parts.append(f"longest_factor_{longest_factor}")
    if not length_parts:
        length_parts = ["all_lengths"]

    out_dir = base / tax_part / pf_part / type_part
    for part in length_parts:
        out_dir = out_dir / part
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def analyze_from_root(
    root: Path,
    json_out: str = "serine_sxn_analysis.json",
    md_out: str = "serine_sxn_analysis.md",
    max_n: int = 50,
    min_n: int = 2,
    taxonomy_terms: Optional[Iterable[str]] = None,
    protein_types: Optional[Iterable[str]] = None,
    partial_full: Optional[str] = None,
    length_range: Optional[Tuple[int, int]] = None,
    length_threshold: Optional[int] = None,
    length_mode: str = "ge",
    longest_factor: Optional[float] = None,
    longest_factor_scope: str = "species",
) -> Tuple[Path, Path]:
    if not root.exists():
        raise FileNotFoundError(f"Root directory not found: {root}")
    records = load_records_from_root(root)
    records = dedupe_records_by_accession(records)
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
    output_root = _analysis_output_dir(
        root,
        taxonomy_terms=taxonomy_terms,
        partial_full=partial_full,
        protein_types=protein_types,
        length_range=length_range,
        length_threshold=length_threshold,
        length_mode=length_mode,
        longest_factor=longest_factor,
    )
    return analyze_and_save(
        filtered_records,
        output_root,
        json_out=json_out,
        md_out=md_out,
        max_n=max_n,
        min_n=min_n,
    )
