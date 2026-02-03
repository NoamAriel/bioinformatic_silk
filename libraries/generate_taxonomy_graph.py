import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


MARK_ORDER = [
    ("Full (Heavy)", "F-H"),
    ("Full (Light)", "F-L"),
    ("Full (Other)", "F-O"),
    ("Partial (Heavy)", "P-H"),
    ("Partial (Light)", "P-L"),
    ("Partial (Other)", "P-O"),
    ("Full (MaSp1)", "F-MaSp1"),
    ("Full (MaSp2)", "F-MaSp2"),
    ("Full (MaSp3)", "F-MaSp3"),
    ("Full (Major_ampullate_unspecified)", "F-MaSp_Unsp"),
    ("Full (MiSp)", "F-MiSp"),
    ("Full (AcSp)", "F-AcSp"),
    ("Full (TuSp)", "F-TuSp"),
    ("Full (Flag)", "F-Flag"),
    ("Full (PySp)", "F-PySp"),
    ("Full (AgSp)", "F-AgSp"),
    ("Full (CrSp)", "F-CrSp"),
    ("Full (Ecribellate_glue)", "F-EcGl"),
    ("Full (CySp)", "F-CySp"),
    ("Full (ASp_unknown)", "F-ASp_Unkn"),
    ("Full (unknown_type)", "F-Unkn"),
    ("Full (fibroin 1)", "F-fib1"),
    ("Full (fibroin 2)", "F-fib2"),
    ("Full (fibroin 3)", "F-fib3"),
    ("Full (fibroin 4)", "F-fib4"),
    ("Partial (MaSp1)", "P-MaSp1"),
    ("Partial (MaSp2)", "P-MaSp2"),
    ("Partial (MaSp3)", "P-MaSp3"),
    ("Partial (Major_ampullate_unspecified)", "P-MaSp_Unsp"),
    ("Partial (MiSp)", "P-MiSp"),
    ("Partial (AcSp)", "P-AcSp"),
    ("Partial (TuSp)", "P-TuSp"),
    ("Partial (Flag)", "P-Flag"),
    ("Partial (PySp)", "P-PySp"),
    ("Partial (AgSp)", "P-AgSp"),
    ("Partial (CrSp)", "P-CrSp"),
    ("Partial (Ecribellate_glue)", "P-EcGl"),
    ("Partial (CySp)", "P-CySp"),
    ("Partial (ASp_unknown)", "P-ASp_Unkn"),
    ("Partial (unknown_type)", "P-Unkn"),
    ("Partial (fibroin 1)", "P-fib1"),
    ("Partial (fibroin 2)", "P-fib2"),
    ("Partial (fibroin 3)", "P-fib3"),
    ("Partial (fibroin 4)", "P-fib4"),
]

RANK_ORDER = [
    "class",
    "Subclass",
    "Infraclass",
    "Cohort",
    "Superorder",
    "Order",
    "Suborder",
    "Infraorder",
    "parvorder",
    "Superfamily",
    "Family",
    "Subfamily",
    "Tribe",
    "Subtribe",
    "Genus",
    "Subgenus",
    "Species",
    "Subspecies",
]

RANK_ORDER_MAP = {rank.lower(): idx for idx, rank in enumerate(RANK_ORDER)}
RANK_Y_SPACING = 1.35
NODE_X_SPACING = 1.6


def _parse_int(value: str) -> int:
    cleaned = value.replace("**", "").strip()
    return int(cleaned) if cleaned.isdigit() else 0


def _normalize_species_key(name: str) -> str:
    return " ".join(name.replace("_", " ").strip().lower().split())


def parse_species_index(md_path: Path) -> dict[str, list[str]]:
    text = md_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    species_marks: dict[str, list[str]] = {}
    in_table = False
    header_map: dict[str, int] = {}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|") and "Species Name" in stripped:
            in_table = True
            header_cells = [c.strip() for c in stripped.strip("|").split("|")]
            header_map = {name: idx for idx, name in enumerate(header_cells)}
            continue
        if not in_table:
            continue
        if stripped.startswith("---"):
            break
        if stripped.startswith("|") and "---" in stripped:
            continue
        if not stripped.startswith("|"):
            continue

        cells = [c.strip() for c in stripped.strip("|").split("|")]
        if len(cells) < 2:
            continue

        species = cells[0]
        normalized_header_map = {
            "".join(ch for ch in key.lower() if ch.isalnum()): idx
            for key, idx in header_map.items()
        }
        counts = {}
        for label, _abbr in MARK_ORDER:
            normalized_label = "".join(ch for ch in label.lower() if ch.isalnum())
            idx = normalized_header_map.get(normalized_label)
            if idx is None or idx >= len(cells):
                counts[label] = 0
            else:
                counts[label] = _parse_int(cells[idx])

        marks = [abbr for key, abbr in MARK_ORDER if counts.get(key, 0) > 0]
        species_marks[_normalize_species_key(species)] = marks

    return species_marks


def parse_taxonomy_tree(tree: dict) -> dict:
    def build_node(name: str, payload: object) -> dict:
        rank = None
        children: list[dict] = []
        if isinstance(payload, dict):
            rank = payload.get("__rank__")
            if isinstance(rank, str):
                rank = rank.strip().lower()
            for child_name, child_data in payload.items():
                if child_name in ("__rank__", "full", "partial"):
                    continue
                if str(child_name).startswith("__"):
                    continue
                children.append(build_node(child_name, child_data))
        return {"name": name, "rank": rank, "children": children, "key": ""}

    roots = [
        build_node(name, payload)
        for name, payload in tree.items()
        if not str(name).startswith("__")
    ]
    if not roots:
        raise ValueError("No taxonomy nodes found in tree.")

    root = roots[0] if len(roots) == 1 else {"name": "ROOT", "rank": None, "children": roots, "key": ""}
    assign_keys(root)
    return root


def assign_keys(node: dict, parent_key: str = "") -> None:
    key = node["name"] if not parent_key else f"{parent_key}/{node['name']}"
    node["key"] = key
    for child in node["children"]:
        assign_keys(child, key)


def gather_rank_labels(root: dict) -> dict[int, str]:
    present = set()
    for node in iter_nodes(root):
        rank = node.get("rank")
        if isinstance(rank, str) and rank in RANK_ORDER_MAP:
            present.add(RANK_ORDER_MAP[rank])
    return {idx: RANK_ORDER[idx] for idx in sorted(present)}


def layout_tree(root: dict) -> tuple[dict[str, tuple[float, float]], int, int, int]:
    positions: dict[str, tuple[float, float]] = {}
    leaf_index = 0
    max_depth = 0
    min_rank_index: int | None = None
    present_rank_indices: set[int] = set()

    def walk(node: dict, depth: int) -> float:
        nonlocal leaf_index, max_depth
        max_depth = max(max_depth, depth)
        children = node["children"]
        if not children:
            x = float(leaf_index) * NODE_X_SPACING
            leaf_index += 1
            positions[node["key"]] = (x, -float(depth) * RANK_Y_SPACING)
            return x

        child_xs = [walk(child, depth + 1) for child in children]
        x = sum(child_xs) / len(child_xs)
        positions[node["key"]] = (x, -float(depth) * RANK_Y_SPACING)
        return x

    walk(root, 0)

    max_rank_index = -1
    for node in iter_nodes(root):
        rank = node.get("rank")
        if isinstance(rank, str) and rank in RANK_ORDER_MAP:
            x, _ = positions[node["key"]]
            rank_idx = RANK_ORDER_MAP[rank]
            present_rank_indices.add(rank_idx)
            positions[node["key"]] = (x, -float(rank_idx) * RANK_Y_SPACING)
            max_rank_index = max(max_rank_index, RANK_ORDER_MAP[rank])
            if min_rank_index is None:
                min_rank_index = rank_idx
            else:
                min_rank_index = min(min_rank_index, rank_idx)

    rank_offset = min_rank_index or 0
    if present_rank_indices:
        sorted_present = sorted(present_rank_indices)
        compact_positions = {idx: i for i, idx in enumerate(sorted_present)}
        for node in iter_nodes(root):
            rank = node.get("rank")
            if isinstance(rank, str) and rank in RANK_ORDER_MAP:
                rank_idx = RANK_ORDER_MAP[rank]
                compact_idx = compact_positions[rank_idx]
                x, _ = positions[node["key"]]
                positions[node["key"]] = (x, -float(compact_idx) * RANK_Y_SPACING)
        rank_offset = 0
        max_rank_index = len(sorted_present) - 1
    elif rank_offset:
        y_offset = float(rank_offset) * RANK_Y_SPACING
        for key, (x, y) in positions.items():
            positions[key] = (x, y + y_offset)

    depth_count = max(max_depth + 1 - rank_offset, max_rank_index - rank_offset + 1)
    return positions, leaf_index, depth_count, rank_offset


def iter_nodes(root: dict) -> list[dict]:
    nodes: list[dict] = []

    def walk(node: dict) -> None:
        nodes.append(node)
        for child in node["children"]:
            walk(child)

    walk(root)
    return nodes


def draw_edges(ax, node: dict, positions: dict[str, tuple[float, float]]) -> None:
    parent_key = node["key"]
    x1, y1 = positions[parent_key]
    for child in node["children"]:
        child_key = child["key"]
        x2, y2 = positions[child_key]
        ax.plot([x1, x2], [y1, y2], color="#4a4a4a", linewidth=1.0)
        draw_edges(ax, child, positions)


def draw_labels(
    ax,
    node: dict,
    positions: dict[str, tuple[float, float]],
    species_marks: dict[str, list[str]],
    species_label_offset: float,
    species_fontsize: int,
    name_rank_map: dict[str, str],
    label_seen: set[tuple[str, str]],
) -> None:
    key = node["key"]
    name = node["name"]
    x, y = positions[key]
    is_leaf = not node["children"]
    label = name
    name_key = _normalize_species_key(name)
    is_species = is_leaf or node.get("rank") == "species" or name_key in species_marks
    if is_species:
        marks = species_marks.get(name_key, [])
        if marks:
            label = f"{name}\n[{', '.join(marks)}]"
    rank_key = node.get("rank") or name_rank_map.get(name_key)
    if rank_key:
        label_group = f"rank:{rank_key}"
    else:
        label_group = f"node:{key}"
    label_key = (label_group, name_key)
    if label_key not in label_seen:
        label_seen.add(label_key)
        if is_species:
            ax.text(
                x,
                y - species_label_offset,
                label,
                ha="right",
                va="center",
                fontsize=species_fontsize,
                rotation=45,
                color="#2f2f2f",
            )
        else:
            ax.text(
                x,
                y,
                label,
                ha="center",
                va="center",
                fontsize=species_fontsize,
                color="#2f2f2f",
            )
    for child in node["children"]:
        draw_labels(
            ax,
            child,
            positions,
            species_marks,
            species_label_offset,
            species_fontsize,
            name_rank_map,
            label_seen,
        )


def render_tree(
    root: dict,
    species_marks: dict[str, list[str]],
    outpath: Path,
    fmt: str = "png",
) -> Path:
    positions, leaf_count, depth_count, rank_offset = layout_tree(root)
    depth_labels = gather_rank_labels(root)
    name_rank_map: dict[str, str] = {}
    for node in iter_nodes(root):
        rank = node.get("rank")
        name = node.get("name", "")
        if not rank or not name:
            continue
        key = name.strip().lower()
        if key in name_rank_map and name_rank_map[key] != rank:
            name_rank_map[key] = ""
        elif key not in name_rank_map:
            name_rank_map[key] = rank
    name_rank_map = {k: v for k, v in name_rank_map.items() if v}
    width = max(8.0, leaf_count * 1.4 * NODE_X_SPACING)
    height = max(4.5, depth_count * RANK_Y_SPACING)
    species_label_offset = 0.6 * RANK_Y_SPACING
    species_fontsize = 9

    fig, ax = plt.subplots(figsize=(width, height))
    draw_edges(ax, root, positions)
    draw_labels(
        ax,
        root,
        positions,
        species_marks,
        species_label_offset,
        species_fontsize,
        name_rank_map,
        set(),
    )

    min_x = min(x for x, _ in positions.values())
    max_x = max(x for x, _ in positions.values())
    min_y = min(y for _, y in positions.values())
    max_y = max(y for _, y in positions.values())
    rank_x = min_x - 1.2

    if depth_labels:
        sorted_present = sorted(depth_labels.keys())
        compact_labels = {i: depth_labels[idx] for i, idx in enumerate(sorted_present)}
    else:
        compact_labels = {}

    for depth in range(depth_count):
        label = compact_labels.get(depth)
        if not label:
            continue
        ax.text(
            rank_x,
            -float(depth) * RANK_Y_SPACING,
            label,
            ha="right",
            va="center",
            fontsize=species_fontsize,
            color="#2f2f2f",
        )

    species_names = {
        _normalize_species_key(node.get("name", ""))
        for node in iter_nodes(root)
        if node.get("rank") == "species"
    }
    if species_names:
        used_marks = {
            mark
            for name, marks in species_marks.items()
            if name in species_names
            for mark in marks
        }
        marked_species = sum(1 for name in species_names if species_marks.get(name))
        print(
            "Taxonomy graph marks:",
            f"species={len(species_names)}",
            f"marked_species={marked_species}",
            f"used_marks={len(used_marks)}",
        )
        legend_handles = []
        for label, abbr in MARK_ORDER:
            if abbr not in used_marks:
                continue
            legend_handles.append(Patch(facecolor="white", edgecolor="#4a4a4a", label=f"{abbr} = {label}"))
        if legend_handles:
            legend = ax.legend(
                handles=legend_handles,
                loc="upper right",
                bbox_to_anchor=(1.0, 1.0),
                borderaxespad=0.0,
                frameon=False,
                fontsize=species_fontsize,
            )
            for text in legend.get_texts():
                text.set_color("#2f2f2f")

    ax.set_xlim(rank_x - 0.8, max_x + 0.8)
    ax.set_ylim(min_y - 0.8 - species_label_offset, max_y + 0.8)
    ax.set_axis_off()

    output_file = outpath.with_suffix(f".{fmt}")
    fig.savefig(output_file, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a taxonomy hierarchy graph from phylo_tree.json with protein-type marks."
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=Path("caddisfly/phylo_tree.json"),
        help="Path to phylo_tree.json",
    )
    parser.add_argument(
        "--index",
        type=Path,
        default=Path("caddisfly/Species_Index.md"),
        help="Path to Species_Index.md",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("caddisfly/taxonomy_tree_with_protein_marks"),
        help="Output path without extension",
    )
    parser.add_argument(
        "--format",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output format",
    )

    args = parser.parse_args()
    tree = json.loads(args.tree.read_text(encoding="utf-8"))
    species_marks = parse_species_index(args.index)
    root = parse_taxonomy_tree(tree)
    output_file = render_tree(
        root=root,
        species_marks=species_marks,
        outpath=args.out,
        fmt=args.format,
    )
    print(f"Wrote: {output_file}")


if __name__ == "__main__":
    main()
