"""Amino-acid group motif labeling utilities."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


DEFAULT_GROUPS: Dict[str, Sequence[str]] = {
    "π": ["A", "V", "L", "I", "M", "P"],  # nonpolar
    "α": ["F", "W", "Y"],  # aromatic
    "η": ["S", "T"],  # hydroxyl
    "θ": ["C", "N", "Q", "G"],  # polar
    "ψ": ["D", "E"],  # negative charged
    "φ": ["K", "R", "H"],  # positive charged
}


GROUP_LEGEND = {
    "π": "nonpolar",
    "α": "aromatic",
    "η": "hydroxyl",
    "θ": "polar",
    "ψ": "negative charged",
    "φ": "positive charged",
}


RANK_ORDER = [
    "class",
    "subclass",
    "infraclass",
    "cohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species",
    "subspecies",
]
def build_aa_to_groups(groups: Mapping[str, Iterable[str]]) -> Dict[str, List[str]]:
    """Build a mapping of amino-acid letter -> list of group names."""
    aa_to_groups: Dict[str, List[str]] = {}
    for group_name, letters in groups.items():
        for letter in letters:
            up = str(letter).upper()
            aa_to_groups.setdefault(up, []).append(group_name)
    return aa_to_groups


def dedupe_records_by_accession(records: Sequence[Mapping[str, object]]) -> List[Mapping[str, object]]:
    """Return records with duplicate accession IDs removed (keeps first seen)."""
    seen: set[str] = set()
    deduped: List[Mapping[str, object]] = []
    for rec in records:
        accession = str(rec.get("accession", "") or "").strip()
        if accession:
            if accession in seen:
                continue
            seen.add(accession)
        deduped.append(rec)
    return deduped


def _record_length(rec: Mapping[str, object]) -> int:
    seq_len = rec.get("sequence_length")
    if isinstance(seq_len, (int, float)) and seq_len > 0:
        return int(seq_len)
    if isinstance(seq_len, str) and seq_len.strip().isdigit():
        return int(seq_len.strip())
    seq = rec.get("origin_sequence", "") or ""
    return len(seq)


def select_longest_records(
    records: Sequence[Mapping[str, object]],
    scope: str = "species",
) -> List[Mapping[str, object]]:
    """Select the longest record(s) by species or globally."""
    if not records:
        return []
    scope_norm = (scope or "species").lower()
    if scope_norm == "global":
        max_len = max(_record_length(rec) for rec in records)
        return [rec for rec in records if _record_length(rec) == max_len]
    if scope_norm != "species":
        raise ValueError("scope must be 'species' or 'global'.")

    by_species: Dict[str, Mapping[str, object]] = {}
    for rec in records:
        name = rec.get("organism_name") or rec.get("organism") or "Unknown"
        if name not in by_species or _record_length(rec) > _record_length(by_species[name]):
            by_species[name] = rec
    return list(by_species.values())


def characterize_sequence(
    seq: str,
    groups: Mapping[str, Iterable[str]] | None = None,
    multi_group_joiner: str = "|",
    unknown_label: str = "unknown",
) -> List[str]:
    """
    Return per-residue group labels for a sequence.

    If a residue belongs to multiple groups (e.g., Y), the labels are joined
    with `multi_group_joiner`.
    """
    group_map = build_aa_to_groups(groups or DEFAULT_GROUPS)
    labels: List[str] = []
    for ch in seq:
        if not ch.isalpha():
            continue
        up = ch.upper()
        hit = group_map.get(up)
        if not hit:
            labels.append(unknown_label)
        elif len(hit) == 1:
            labels.append(hit[0])
        else:
            labels.append(multi_group_joiner.join(hit))
    return labels


def motif_string(
    seq: str,
    groups: Mapping[str, Iterable[str]] | None = None,
    separator: str = ",",
    **kwargs,
) -> str:
    """Return a compact comma-separated motif label string."""
    return separator.join(characterize_sequence(seq, groups=groups, **kwargs))


def _slugify(text: str) -> str:
    """Make a filesystem-friendly slug."""
    return "".join(ch if ch.isalnum() else "_" for ch in text).strip("_") or "none"


def _load_phylo_tree(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _build_accession_lineages(
    tree: Any, accessions_filter: Optional[set[str]] = None
) -> Dict[str, List[Tuple[str, str]]]:
    accessions: Dict[str, List[Tuple[str, str]]] = {}

    def walk(node: Any, lineage: List[Tuple[str, str]]) -> None:
        if isinstance(node, dict):
            for name, child in node.items():
                if name == "__rank__":
                    continue
                if isinstance(child, dict):
                    rank = child.get("__rank__")
                    next_lineage = lineage
                    if rank:
                        next_lineage = lineage + [(str(rank).lower(), str(name))]
                    walk(child, next_lineage)
                elif isinstance(child, list):
                    for item in child:
                        if isinstance(item, dict):
                            acc = str(item.get("accession", "") or "").strip()
                            if acc and (accessions_filter is None or acc in accessions_filter):
                                accessions[acc] = list(lineage)
        elif isinstance(node, list):
            for item in node:
                walk(item, lineage)

    walk(tree, [])
    return accessions


def load_phylo_tree_lineages(
    root: Path,
    phylo_tree_path: Optional[Path] = None,
    accessions_filter: Optional[set[str]] = None,
) -> Dict[str, List[Tuple[str, str]]]:
    path = phylo_tree_path or (root / "phylo_tree.json")
    tree = _load_phylo_tree(path)
    if not tree:
        return {}
    return _build_accession_lineages(tree, accessions_filter=accessions_filter)


def resolve_analysis_output_dir(root: Path, analysis_subdir: str = "aa_motifs") -> Path:
    """Return the output directory under <root>/data_analysis/<analysis_subdir>."""
    out_dir = root / "data_analysis" / analysis_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def resolve_aa_motifs_output_dir(
    root: Path,
    taxonomy_terms: Iterable[str] | None = None,
    partial_full: str | None = None,
    protein_types: Iterable[str] | None = None,
    length_range: tuple[int, int] | None = None,
    length_threshold: int | None = None,
    length_mode: str = "ge",
    longest_factor: float | None = None,
    min_motif_length: int | None = None,
    max_motif_length: int | None = None,
) -> Path:
    """
    Build output path under:
      <root>/aa_motifs/aa_composition/<taxonomy>/<partial_full>/<type>/<length_filters>
    """
    base = root / "data_analysis" / "aa_motifs"

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
    if min_motif_length is not None and max_motif_length is not None:
        out_dir = out_dir / f"{min_motif_length}_to_{max_motif_length}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def motif_frequency_from_counts(
    counts: Mapping[str, int],
    lengths: Mapping[str, int],
    sequence_length: int,
) -> Dict[str, float]:
    """Compute relative frequency percentages for motifs using window length."""
    if sequence_length <= 0:
        return {motif: 0.0 for motif in counts}
    return {
        motif: (count * lengths.get(motif, 0) / sequence_length * 100.0)
        for motif, count in counts.items()
    }


def unique_motif_frequency_from_counts(
    counts: Mapping[str, int],
    motif_length: int,
    sequence_length: int,
) -> Dict[str, float]:
    """Compute relative frequency percentages for fixed-length motifs."""
    if sequence_length <= 0:
        return {motif: 0.0 for motif in counts}
    return {
        motif: (count * motif_length / sequence_length * 100.0)
        for motif, count in counts.items()
    }


def count_group_motifs(
    labels: Sequence[str],
    min_length: int | None,
    max_length: int | None,
    unknown_label: str = "unknown",
    skip_unknown: bool = True,
    labels_str: str | None = None,
    label_offsets: Sequence[int] | None = None,
) -> Dict[str, int]:
    if min_length is None or max_length is None:
        return {}
    if min_length < 1 or max_length < min_length:
        raise ValueError("Invalid motif length range.")

    counts: Dict[str, int] = {}
    L = len(labels)
    unknown_prefix = _build_unknown_prefix(labels, unknown_label) if skip_unknown else None
    if labels_str is None or label_offsets is None:
        labels_str = "".join(labels)
        label_offsets = _build_label_offsets(labels)
    for n in range(min_length, max_length + 1):
        if n > L:
            break
        for i in range(L - n + 1):
            if skip_unknown and unknown_prefix and (unknown_prefix[i + n] - unknown_prefix[i]) > 0:
                continue
            motif = labels_str[label_offsets[i] : label_offsets[i + n]]
            counts[motif] = counts.get(motif, 0) + 1
    return counts


def count_group_motif_lengths(
    labels: Sequence[str],
    min_length: int | None,
    max_length: int | None,
    unknown_label: str = "unknown",
    skip_unknown: bool = True,
    labels_str: str | None = None,
    label_offsets: Sequence[int] | None = None,
) -> Dict[str, int]:
    if min_length is None or max_length is None:
        return {}
    if min_length < 1 or max_length < min_length:
        raise ValueError("Invalid motif length range.")
    lengths: Dict[str, int] = {}
    L = len(labels)
    unknown_prefix = _build_unknown_prefix(labels, unknown_label) if skip_unknown else None
    if labels_str is None or label_offsets is None:
        labels_str = "".join(labels)
        label_offsets = _build_label_offsets(labels)
    for n in range(min_length, max_length + 1):
        if n > L:
            break
        for i in range(L - n + 1):
            if skip_unknown and unknown_prefix and (unknown_prefix[i + n] - unknown_prefix[i]) > 0:
                continue
            motif = labels_str[label_offsets[i] : label_offsets[i + n]]
            if motif not in lengths:
                lengths[motif] = n
    return lengths


def count_group_motifs_with_lengths(
    labels: Sequence[str],
    min_length: int | None,
    max_length: int | None,
    unknown_label: str = "unknown",
    skip_unknown: bool = True,
    labels_str: str | None = None,
    label_offsets: Sequence[int] | None = None,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    if min_length is None or max_length is None:
        return {}, {}
    if min_length < 1 or max_length < min_length:
        raise ValueError("Invalid motif length range.")

    counts: Dict[str, int] = {}
    lengths: Dict[str, int] = {}
    L = len(labels)
    unknown_prefix = _build_unknown_prefix(labels, unknown_label) if skip_unknown else None
    if labels_str is None or label_offsets is None:
        labels_str = "".join(labels)
        label_offsets = _build_label_offsets(labels)
    for n in range(min_length, max_length + 1):
        if n > L:
            break
        for i in range(L - n + 1):
            if skip_unknown and unknown_prefix and (unknown_prefix[i + n] - unknown_prefix[i]) > 0:
                continue
            motif = labels_str[label_offsets[i] : label_offsets[i + n]]
            counts[motif] = counts.get(motif, 0) + 1
            if motif not in lengths:
                lengths[motif] = n
    return counts, lengths


def collect_unique_motifs(
    labels: Sequence[str],
    motif_length: int | None,
    unknown_label: str = "unknown",
    skip_unknown: bool = True,
    labels_str: str | None = None,
    label_offsets: Sequence[int] | None = None,
) -> Dict[str, int]:
    if motif_length is None or motif_length < 1:
        return {}
    counts: Dict[str, int] = {}
    L = len(labels)
    if motif_length > L:
        return {}
    unknown_prefix = _build_unknown_prefix(labels, unknown_label) if skip_unknown else None
    if labels_str is None or label_offsets is None:
        labels_str = "".join(labels)
        label_offsets = _build_label_offsets(labels)
    for i in range(L - motif_length + 1):
        if skip_unknown and unknown_prefix and (unknown_prefix[i + motif_length] - unknown_prefix[i]) > 0:
            continue
        key = labels_str[label_offsets[i] : label_offsets[i + motif_length]]
        counts[key] = counts.get(key, 0) + 1
    return counts


def _build_unknown_prefix(labels: Sequence[str], unknown_label: str) -> List[int]:
    prefix = [0] * (len(labels) + 1)
    for idx, label in enumerate(labels):
        prefix[idx + 1] = prefix[idx] + (1 if label == unknown_label else 0)
    return prefix


def _build_label_offsets(labels: Sequence[str]) -> List[int]:
    offsets = [0] * (len(labels) + 1)
    total = 0
    for idx, label in enumerate(labels):
        total += len(label)
        offsets[idx + 1] = total
    return offsets




def is_repeated_motif(motif: str) -> bool:
    """Return True if motif is a pure repetition of a shorter motif."""
    L = len(motif)
    for k in range(1, L // 2 + 1):
        if L % k != 0:
            continue
        unit = motif[:k]
        if unit * (L // k) == motif:
            return True
    return False


def _is_combined_motif(motif: str, motifs: set[str]) -> bool:
    for i in range(1, len(motif)):
        if motif[:i] in motifs and motif[i:] in motifs:
            return True
    return False


def _prune_combined_counts(counts: Dict[str, int]) -> Dict[str, int]:
    motifs = set(counts.keys())
    combined = {
        motif
        for motif in motifs
        if is_repeated_motif(motif) or _is_combined_motif(motif, motifs)
    }
    return {motif: count for motif, count in counts.items() if motif not in combined}


def _prune_combined_unique(counts: Dict[str, int]) -> Dict[str, int]:
    motif_set = set(counts.keys())
    combined = {
        motif
        for motif in motif_set
        if is_repeated_motif(motif) or _is_combined_motif(motif, motif_set)
    }
    return {motif: count for motif, count in counts.items() if motif not in combined}


def _taxonomy_chain_from_lineage(
    lineage: List[Tuple[str, str]],
    taxonomy_terms: Optional[Iterable[str]] = None,
) -> List[Tuple[str, str]]:
    rank_to_name: Dict[str, str] = {}
    for rank, name in lineage:
        rank_to_name[str(rank).lower()] = str(name)

    chain: List[Tuple[str, str]] = []
    for rank in RANK_ORDER:
        if rank in rank_to_name:
            chain.append((rank, rank_to_name[rank]))

    if not taxonomy_terms:
        return chain

    terms_lower = {str(t).strip().lower() for t in taxonomy_terms if str(t).strip()}
    cutoff_idx = None
    for idx, (_, name) in enumerate(chain):
        if name.lower() in terms_lower:
            cutoff_idx = idx
    if cutoff_idx is None:
        return chain
    return chain[cutoff_idx:]


def _shared_from_records(
    records: Sequence[Dict[str, Any]],
    counts_key: str,
    freq_key: str,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    sets = [set((rec.get(counts_key) or {}).keys()) for rec in records]
    if not sets:
        return {}, {}
    shared = set.intersection(*sets)
    avg_freqs: Dict[str, float] = {}
    for motif in shared:
        freqs = [rec.get(freq_key, {}).get(motif, 0.0) for rec in records]
        avg_freqs[motif] = sum(freqs) / len(freqs) if freqs else 0.0
    return {motif: 0.0 for motif in shared}, avg_freqs


def compute_hierarchical_shared_motifs(
    records: Sequence[Dict[str, Any]],
    taxonomy_terms: Optional[Iterable[str]],
    phylo_lineages: Dict[str, List[Tuple[str, str]]],
) -> Dict[str, Any]:
    nodes: Dict[Tuple[str, str], Dict[str, Any]] = {}
    rank_index = {rank: idx for idx, rank in enumerate(RANK_ORDER)}

    def _node(key: Tuple[str, str]) -> Dict[str, Any]:
        if key not in nodes:
            nodes[key] = {"records": [], "children": set(), "parent": None}
        return nodes[key]

    for rec in records:
        acc = str(rec.get("accession", "") or "").strip()
        lineage = phylo_lineages.get(acc, [])
        if not lineage:
            continue
        chain = _taxonomy_chain_from_lineage(lineage, taxonomy_terms=taxonomy_terms)
        if not chain:
            continue
        lowest = chain[-1]
        _node(lowest)["records"].append(rec)
        for idx in range(len(chain) - 1):
            parent = chain[idx]
            child = chain[idx + 1]
            _node(parent)["children"].add(child)
            _node(child)["parent"] = parent

    # Bottom-up aggregation
    shared: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for rank in reversed(RANK_ORDER):
        for key, info in list(nodes.items()):
            if key[0] != rank:
                continue
            children = info["children"]
            if children:
                child_entries = [shared.get(child) for child in children if child in shared]
                if not child_entries:
                    continue
                motif_sets = [set(entry["shared_motifs"].keys()) for entry in child_entries]
                unique_sets = [set(entry["shared_unique_motifs"].keys()) for entry in child_entries]
                motifs = set.intersection(*motif_sets) if motif_sets else set()
                uniques = set.intersection(*unique_sets) if unique_sets else set()

                motif_avgs: Dict[str, float] = {}
                for motif in motifs:
                    freqs = [entry["shared_motifs"].get(motif, 0.0) for entry in child_entries]
                    motif_avgs[motif] = sum(freqs) / len(freqs) if freqs else 0.0

                unique_avgs: Dict[str, float] = {}
                for motif in uniques:
                    freqs = [entry["shared_unique_motifs"].get(motif, 0.0) for entry in child_entries]
                    unique_avgs[motif] = sum(freqs) / len(freqs) if freqs else 0.0

                shared[key] = {
                    "rank": key[0],
                    "name": key[1],
                    "num_records": sum(entry.get("num_records", 0) for entry in child_entries),
                    "shared_motifs": motif_avgs,
                    "shared_unique_motifs": unique_avgs,
                }
            else:
                recs = info["records"]
                if len(recs) < 2:
                    continue
                _, motif_avgs = _shared_from_records(recs, "motif_counts", "motif_frequencies")
                _, unique_avgs = _shared_from_records(recs, "unique_motifs", "unique_motif_frequencies")
                shared[key] = {
                    "rank": key[0],
                    "name": key[1],
                    "num_records": len(recs),
                    "shared_motifs": motif_avgs,
                    "shared_unique_motifs": unique_avgs,
                }

    # Sort for output
    output: List[Dict[str, Any]] = []
    for key, entry in shared.items():
        output.append(entry)
    output.sort(key=lambda e: rank_index.get(e["rank"], -1), reverse=True)
    return {"nodes": output}


def analyze_group_motifs(
    records: Sequence[Dict[str, Any]],
    groups: Mapping[str, Iterable[str]] | None = None,
    multi_group_joiner: str = "|",
    unknown_label: str = "unknown",
    include_motif_string: bool = False,
    include_motif_counts: bool = False,
    min_motif_length: int | None = None,
    max_motif_length: int | None = None,
    skip_unknown_motifs: bool = True,
    include_unique_motifs: bool = False,
    unique_motif_length: int | None = None,
    skip_combined_motifs: bool = False,
    include_reflected_sequence: bool = False,
    include_similar_motifs: bool = False,
    include_taxonomy_shared_motifs: bool = False,
    taxonomy_terms: Optional[Iterable[str]] = None,
    phylo_tree_path: Optional[Path] = None,
    root: Optional[Path] = None,
) -> Dict[str, Any]:
    group_def = groups or DEFAULT_GROUPS
    group_map = build_aa_to_groups(group_def)
    group_names = list(group_def.keys())

    analyzed: List[Dict[str, Any]] = []
    global_group_counts = {name: 0 for name in group_names}
    global_label_counts: Dict[str, int] = {}
    total_length = 0

    for rec in records:
        seq = rec.get("origin_sequence", "") or ""
        if not seq:
            continue

        label_counts: Dict[str, int] = {}
        group_counts = {name: 0 for name in group_names}

        motif_labels: List[str] = []
        group_labels: List[str] = []
        seq_len = 0
        for ch in seq:
            if not ch.isalpha():
                continue
            seq_len += 1
            up = ch.upper()
            hit = group_map.get(up)
            if not hit:
                label = unknown_label
            elif len(hit) == 1:
                label = hit[0]
                group_counts[label] += 1
            else:
                label = multi_group_joiner.join(hit)
                for name in hit:
                    group_counts[name] += 1
            label_counts[label] = label_counts.get(label, 0) + 1
            group_labels.append(label)
            if include_motif_string:
                motif_labels.append(label)

        total_length += seq_len
        for name in group_names:
            global_group_counts[name] += group_counts[name]
        for label, count in label_counts.items():
            global_label_counts[label] = global_label_counts.get(label, 0) + count

        group_fractions = {
            name: (group_counts[name] / seq_len * 100.0) if seq_len else 0.0
            for name in group_names
        }

        entry: Dict[str, Any] = {
            "accession": rec.get("accession", ""),
            "title": rec.get("title", ""),
            "organism": rec.get("organism_name", "Unknown organism"),
            "taxonomy_from_order": rec.get("taxonomy_from_order", []),
            "taxonomy_from_araneae": rec.get("taxonomy_from_araneae", []),
            "taxonomy_full": rec.get("taxonomy_full", []),
            "partial_full": rec.get("partial_full", "") or "unknown",
            "type": rec.get("type", "") or "unknown",
            "length": seq_len,
            "group_counts": group_counts,
            "group_fractions": group_fractions,
            "label_counts": label_counts,
        }
        if include_motif_string:
            entry["motif_string"] = ",".join(motif_labels)
        if include_motif_counts:
            labels_str = "".join(group_labels)
            label_offsets = _build_label_offsets(group_labels)
            motif_counts, motif_lengths = count_group_motifs_with_lengths(
                group_labels,
                min_length=min_motif_length,
                max_length=max_motif_length,
                unknown_label=unknown_label,
                skip_unknown=skip_unknown_motifs,
                labels_str=labels_str,
                label_offsets=label_offsets,
            )
            entry["motif_counts"] = motif_counts
            entry["motif_lengths"] = motif_lengths
            if skip_combined_motifs:
                entry["motif_counts"] = _prune_combined_counts(entry["motif_counts"])
                entry["motif_lengths"] = {
                    motif: length
                    for motif, length in entry["motif_lengths"].items()
                    if motif in entry["motif_counts"]
                }
            entry["motif_frequencies"] = motif_frequency_from_counts(
                entry["motif_counts"], entry["motif_lengths"], seq_len
            )
        if include_unique_motifs:
            labels_str = labels_str if include_motif_counts else "".join(group_labels)
            label_offsets = label_offsets if include_motif_counts else _build_label_offsets(group_labels)
            entry["unique_motifs"] = collect_unique_motifs(
                group_labels,
                motif_length=unique_motif_length,
                unknown_label=unknown_label,
                skip_unknown=skip_unknown_motifs,
                labels_str=labels_str,
                label_offsets=label_offsets,
            )
            if skip_combined_motifs:
                entry["unique_motifs"] = _prune_combined_unique(entry["unique_motifs"])
            entry["unique_motif_frequencies"] = unique_motif_frequency_from_counts(
                entry["unique_motifs"], unique_motif_length or 0, seq_len
            )
        if include_reflected_sequence:
            entry["reflected_sequence"] = "".join(group_labels)
        analyzed.append(entry)

    similar_motifs: Dict[str, Any] = {}
    if include_similar_motifs:
        by_species: Dict[str, List[Dict[str, Any]]] = {}
        for rec in analyzed:
            species = rec.get("organism", "Unknown") or "Unknown"
            by_species.setdefault(species, []).append(rec)

        for species, recs in by_species.items():
            if len(recs) < 2:
                continue
            motif_sets = [
                set((rec.get("motif_counts") or {}).keys())
                for rec in recs
                if rec.get("motif_counts")
            ]
            unique_sets = [
                set((rec.get("unique_motifs") or {}).keys())
                for rec in recs
                if rec.get("unique_motifs")
            ]
            shared_motifs = set.intersection(*motif_sets) if motif_sets else set()
            shared_unique = set.intersection(*unique_sets) if unique_sets else set()

            motif_avg_freq = {}
            for motif in shared_motifs:
                freqs = [rec.get("motif_frequencies", {}).get(motif, 0.0) for rec in recs]
                motif_avg_freq[motif] = sum(freqs) / len(freqs) if freqs else 0.0

            unique_avg_freq = {}
            for motif in shared_unique:
                freqs = [rec.get("unique_motif_frequencies", {}).get(motif, 0.0) for rec in recs]
                unique_avg_freq[motif] = sum(freqs) / len(freqs) if freqs else 0.0

            similar_motifs[species] = {
                "shared_motifs": motif_avg_freq,
                "shared_unique_motifs": unique_avg_freq,
                "num_records": len(recs),
            }

    taxonomy_shared: Dict[str, Any] = {}
    if include_taxonomy_shared_motifs:
        resolved_root = root or Path(".")
        accession_filter = {
            str(rec.get("accession", "") or "").strip()
            for rec in analyzed
            if str(rec.get("accession", "") or "").strip()
        }
        phylo_lineages = load_phylo_tree_lineages(
            resolved_root,
            phylo_tree_path=phylo_tree_path,
            accessions_filter=accession_filter or None,
        )
        taxonomy_shared = compute_hierarchical_shared_motifs(
            analyzed,
            taxonomy_terms=taxonomy_terms,
            phylo_lineages=phylo_lineages,
        )

    global_group_fractions = {
        name: (global_group_counts[name] / total_length * 100.0) if total_length else 0.0
        for name in group_names
    }

    return {
        "groups": {k: list(v) for k, v in group_def.items()},
        "group_legend": GROUP_LEGEND,
        "group_order": group_names,
        "multi_group_joiner": multi_group_joiner,
        "unknown_label": unknown_label,
        "include_motif_string": include_motif_string,
        "include_motif_counts": include_motif_counts,
        "min_motif_length": min_motif_length,
        "max_motif_length": max_motif_length,
        "skip_unknown_motifs": skip_unknown_motifs,
        "include_unique_motifs": include_unique_motifs,
        "unique_motif_length": unique_motif_length,
        "skip_combined_motifs": skip_combined_motifs,
        "include_reflected_sequence": include_reflected_sequence,
        "include_similar_motifs": include_similar_motifs,
        "include_taxonomy_shared_motifs": include_taxonomy_shared_motifs,
        "analyzed_records": analyzed,
        "similar_motifs": similar_motifs,
        "taxonomy_shared_motifs": taxonomy_shared,
        "global_group_counts": global_group_counts,
        "global_group_fractions": global_group_fractions,
        "global_label_counts": global_label_counts,
        "num_records": len(analyzed),
        "total_length": total_length,
    }


def write_json(root: Path, data: Dict[str, Any], out_name: str) -> Path:
    out_path = root / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return out_path


def write_md(root: Path, summary: Dict[str, Any], out_name: str) -> Path:
    group_order = summary.get("group_order", [])
    include_motif_counts = summary.get("include_motif_counts", False)
    include_unique_motifs = summary.get("include_unique_motifs", False)
    skip_combined = summary.get("skip_combined_motifs", False)
    lines: List[str] = []
    lines.append("# Amino Acid Group Motif Analysis")
    lines.append("")
    lines.append("This report summarizes residue group composition per sequence.")
    lines.append(f"- Groups: {', '.join(group_order)}")
    lines.append(
        "- Group legend: "
        + ", ".join(f"{symbol} = {meaning}" for symbol, meaning in GROUP_LEGEND.items())
    )
    if summary.get("include_motif_string"):
        lines.append("- Motif strings: included in JSON output")
    lines.append("---")

    lines.append("## Global Summary")
    lines.append("")
    lines.append(f"- Number of sequences analyzed: **{summary.get('num_records', 0)}**")
    lines.append(f"- Total residues across all sequences: **{summary.get('total_length', 0)}**")
    lines.append("")
    lines.append("| Group | Total Count | Global Percentage |")
    lines.append("| :--- | ---: | ---: |")
    for name in group_order:
        count = summary["global_group_counts"].get(name, 0)
        frac = summary["global_group_fractions"].get(name, 0.0)
        lines.append(f"| {name} | {count} | {frac:.2f}% |")
    lines.append("---")

    records = summary.get("analyzed_records", [])
    if not records:
        out_path = root / out_name
        lines.append("")
        lines.append("_No records remained after filtering._")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return out_path

    header_cols = ["Organism", "Type", "Partial/Full", "Accession ID", "Length"]
    group_cols = [f"{name} %" for name in group_order]
    lines.append("## Per-Record Summary")
    lines.append("")
    lines.append("| " + " | ".join(header_cols + group_cols) + " |")
    lines.append("| " + " :--- |" * len(header_cols) + " ---: |" * len(group_cols))

    for rec in records:
        row = [
            rec.get("organism", "Unknown"),
            rec.get("type", "unknown"),
            rec.get("partial_full", "unknown"),
            f"`{rec.get('accession', '')}`",
            str(rec.get("length", 0)),
        ]
        fracs = rec.get("group_fractions", {})
        for name in group_order:
            row.append(f"{fracs.get(name, 0.0):.2f}%")
        lines.append("| " + " | ".join(row) + " |")

    if include_motif_counts:
        min_len = summary.get("min_motif_length")
        max_len = summary.get("max_motif_length")
        lines.append("---")
        lines.append("## Motif Counts")
        lines.append("")
        lines.append(f"- Motif length range: **{min_len}** to **{max_len}**")
        if summary.get("skip_unknown_motifs"):
            lines.append("- Motifs containing unknown labels are skipped.")
        if skip_combined:
            lines.append("- Motifs that are pure repeats or concatenations of smaller motifs are skipped.")

        for rec in records:
            motif_counts = rec.get("motif_counts") or {}
            motif_freqs = rec.get("motif_frequencies") or {}
            if not motif_counts:
                continue
            organism = rec.get("organism", "Unknown")
            accession = rec.get("accession", "")
            lines.append("")
            lines.append(f"### {organism} `{accession}`")
            lines.append("| Motif | Count | Relative % |")
            lines.append("| :--- | ---: | ---: |")
            for motif, count in sorted(
                motif_counts.items(), key=lambda kv: (-kv[1], kv[0])
            ):
                rel = motif_freqs.get(motif, 0.0)
                lines.append(f"| {motif} | {count} | {rel:.4f}% |")

    if include_unique_motifs:
        motif_len = summary.get("unique_motif_length")
        lines.append("---")
        lines.append("## Unique Motifs")
        lines.append("")
        lines.append(f"- Unique motif length: **{motif_len}**")
        if summary.get("skip_unknown_motifs"):
            lines.append("- Motifs containing unknown labels are skipped.")
        if skip_combined:
            lines.append("- Motifs that are pure repeats or concatenations of smaller motifs are skipped.")

        for rec in records:
            unique_counts = rec.get("unique_motifs") or {}
            unique_freqs = rec.get("unique_motif_frequencies") or {}
            if not unique_counts:
                continue
            organism = rec.get("organism", "Unknown")
            accession = rec.get("accession", "")
            lines.append("")
            lines.append(f"### {organism} `{accession}`")
            lines.append(f"- Unique motifs found: **{len(unique_counts)}**")
            lines.append("| Motif | Count | Relative % |")
            lines.append("| :--- | ---: | ---: |")
            for motif, count in sorted(unique_counts.items(), key=lambda kv: (-kv[1], kv[0])):
                rel = unique_freqs.get(motif, 0.0)
                lines.append(f"| {motif} | {count} | {rel:.4f}% |")

    if summary.get("include_similar_motifs"):
        similar = summary.get("similar_motifs", {})
        if similar:
            lines.append("---")
            lines.append("## Similar Motifs By Species")
            for species, data in similar.items():
                lines.append("")
                lines.append(f"### {species}")
                lines.append(f"- Number of records: **{data.get('num_records', 0)}**")

                shared_motifs = data.get("shared_motifs") or {}
                if shared_motifs:
                    lines.append("")
                    lines.append("#### Shared Motif Counts (avg relative %)")
                    lines.append("| Motif | Avg Relative % |")
                    lines.append("| :--- | ---: |")
                    for motif, freq in sorted(shared_motifs.items(), key=lambda kv: (-kv[1], kv[0])):
                        lines.append(f"| {motif} | {freq:.4f}% |")

                shared_unique = data.get("shared_unique_motifs") or {}
                if shared_unique:
                    lines.append("")
                    lines.append("#### Shared Unique Motifs (avg relative %)")
                    lines.append("| Motif | Avg Relative % |")
                    lines.append("| :--- | ---: |")
                    for motif, freq in sorted(shared_unique.items(), key=lambda kv: (-kv[1], kv[0])):
                        lines.append(f"| {motif} | {freq:.4f}% |")

    if summary.get("include_taxonomy_shared_motifs"):
        nodes = (summary.get("taxonomy_shared_motifs") or {}).get("nodes", [])
        if nodes:
            lines.append("---")
            lines.append("## Hierarchical Shared Motifs (Taxonomy)")
            current_rank = None
            for node in nodes:
                rank = node.get("rank", "unknown")
                name = node.get("name", "Unknown")
                if rank != current_rank:
                    current_rank = rank
                    lines.append("")
                    lines.append(f"### {rank}")
                lines.append("")
                lines.append(f"#### {name}")
                lines.append(f"- Number of records: **{node.get('num_records', 0)}**")

                shared_motifs = node.get("shared_motifs") or {}
                if shared_motifs:
                    lines.append("##### Motif")
                    lines.append("| Motif | Avg Relative % |")
                    lines.append("| :--- | ---: |")
                    for motif, freq in sorted(shared_motifs.items(), key=lambda kv: (-kv[1], kv[0])):
                        lines.append(f"| {motif} | {freq:.4f}% |")

                shared_unique = node.get("shared_unique_motifs") or {}
                if shared_unique:
                    lines.append("##### Unique Motif")
                    lines.append("| Unique Motif | Avg Relative % |")
                    lines.append("| :--- | ---: |")
                    for motif, freq in sorted(shared_unique.items(), key=lambda kv: (-kv[1], kv[0])):
                        lines.append(f"| {motif} | {freq:.4f}% |")
                if not shared_motifs and not shared_unique:
                    lines.append("- No shared motifs found.")

    out_path = root / out_name
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def write_reflected_json(root: Path, summary: Dict[str, Any], out_name: str) -> Path:
    records = summary.get("analyzed_records", []) or []
    reflected = []
    for rec in records:
        reflected.append(
            {
                "accession": rec.get("accession", ""),
                "organism": rec.get("organism", ""),
                "type": rec.get("type", ""),
                "partial_full": rec.get("partial_full", ""),
                "length": rec.get("length", 0),
                "reflected_sequence": rec.get("reflected_sequence", ""),
            }
        )
    payload = {
        "groups": summary.get("groups", {}),
        "group_legend": summary.get("group_legend", {}),
        "group_order": summary.get("group_order", []),
        "num_records": summary.get("num_records", 0),
        "total_length": summary.get("total_length", 0),
        "records": reflected,
    }
    out_path = root / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return out_path


def write_reflected_md(root: Path, summary: Dict[str, Any], out_name: str) -> Path:
    group_order = summary.get("group_order", [])
    lines: List[str] = []
    lines.append("# Amino Acid Group Motif Analysis")
    lines.append("")
    lines.append("This report summarizes residue group composition per sequence.")
    lines.append(f"- Groups: {', '.join(group_order)}")
    lines.append(
        "- Group legend: "
        + ", ".join(f"{symbol} = {meaning}" for symbol, meaning in GROUP_LEGEND.items())
    )
    lines.append("---")

    lines.append("## Global Summary")
    lines.append("")
    lines.append(f"- Number of sequences analyzed: **{summary.get('num_records', 0)}**")
    lines.append(f"- Total residues across all sequences: **{summary.get('total_length', 0)}**")
    lines.append("")
    lines.append("| Group | Total Count | Global Percentage |")
    lines.append("| :--- | ---: | ---: |")
    for name in group_order:
        count = summary["global_group_counts"].get(name, 0)
        frac = summary["global_group_fractions"].get(name, 0.0)
        lines.append(f"| {name} | {count} | {frac:.2f}% |")
    lines.append("---")

    records = summary.get("analyzed_records", [])
    if not records:
        out_path = root / out_name
        lines.append("")
        lines.append("_No records remained after filtering._")
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return out_path

    header_cols = ["Organism", "Type", "Partial/Full", "Accession ID", "Length"]
    group_cols = [f"{name} %" for name in group_order]
    lines.append("## Per-Record Summary")
    lines.append("")
    lines.append("| " + " | ".join(header_cols + group_cols) + " |")
    lines.append("| " + " :--- |" * len(header_cols) + " ---: |" * len(group_cols))

    for rec in records:
        row = [
            rec.get("organism", "Unknown"),
            rec.get("type", "unknown"),
            rec.get("partial_full", "unknown"),
            f"`{rec.get('accession', '')}`",
            str(rec.get("length", 0)),
        ]
        fracs = rec.get("group_fractions", {})
        for name in group_order:
            row.append(f"{fracs.get(name, 0.0):.2f}%")
        lines.append("| " + " | ".join(row) + " |")

    lines.append("---")
    lines.append("## Reflected Sequences")
    for rec in records:
        reflected = rec.get("reflected_sequence")
        if not reflected:
            continue
        organism = rec.get("organism", "Unknown")
        accession = rec.get("accession", "")
        lines.append("")
        lines.append(f"### {organism} `{accession}`")
        lines.append("```text")
        lines.append(reflected)
        lines.append("```")

    out_path = root / out_name
    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


if __name__ == "__main__":
    example = "AVLWWSTYDK"
    print(example)
    print(motif_string(example))
