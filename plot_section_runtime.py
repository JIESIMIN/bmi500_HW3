#!/usr/bin/env python
"""
Aggregate the instrumented-profiling results and plot per-section runtime
against dataset size as a stacked area chart.

Run: python plot_section_runtime.py
"""

import glob
import os
import re
import tarfile

import numpy as np
import matplotlib
matplotlib.use("Agg")         
import matplotlib.pyplot as plt


SURFACE = "#fcfcfb"
INK     = "#0b0b0b"
INK_2   = "#52514e"
GRID    = "#e3e2de"
SLOTS   = ["#56B4E9", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#800020", "#4a3aa7"]
OTHER   = "#9a9992"     

MAX_SERIES = 7          # + "Other" = 8 bands, the categorical cap

# Pipeline sections
PIPELINE = [
    "startup_argparse", "read_10x_mtx", "filter_cells_genes",
    "normalize_log1p", "highly_variable_genes", "scale", "pca",
    "neighbors", "louvain", "umap", "write_h5ad", "rank_genes_groups",
]

# rows of the instrumented table, e.g. "umap    8.323 s  ( 36.3%)"
ROW_RE = re.compile(r"^(\w+)\s+([0-9.]+)\s+s")
SUFFIX = "_instrumented_time.txt"


def n_cells(dataset):
    """Cell count from the MatrixMarket header (line 3: genes cells nnz)."""
    mtx = f"data/{dataset}/filtered_gene_bc_matrices/matrix.mtx"
    if os.path.exists(mtx):
        with open(mtx) as fh:
            return int(fh.readlines()[2].split()[1])
    # not extracted -- stream the header straight out of the tarball
    with tarfile.open(f"data/{dataset}.tgz") as tf:
        fh = tf.extractfile(f"{dataset}/filtered_gene_bc_matrices/matrix.mtx")
        for i, line in enumerate(fh):
            if i == 2:
                return int(line.split()[1])
    raise SystemExit(f"cannot determine cell count for {dataset}")


def load():
    """{dataset: {'n_cells': int, 'sections': {name: seconds}}}"""
    out = {}
    for path in sorted(glob.glob(f"profiles/*{SUFFIX}")):
        ds = os.path.basename(path)[:-len(SUFFIX)]
        sections = {}
        with open(path) as fh:
            for line in fh:
                m = ROW_RE.match(line)
                if m and m.group(1) != "TOTAL":
                    sections[m.group(1)] = float(m.group(2))
        if sections:
            out[ds] = {"n_cells": n_cells(ds), "sections": sections}
    return out


def fold(data):
    """Keep the MAX_SERIES largest sections; bucket the rest into 'Other'."""
    biggest = max(data, key=lambda d: data[d]["n_cells"])
    ranked = sorted(data[biggest]["sections"],
                    key=lambda s: data[biggest]["sections"][s], reverse=True)
    keep = set(ranked[:MAX_SERIES])
    series = [s for s in PIPELINE if s in keep]          # pipeline order
    folded = {}
    for ds, d in data.items():
        vals = {s: d["sections"].get(s, 0.0) for s in series}
        vals["Other"] = sum(v for k, v in d["sections"].items() if k not in keep)
        folded[ds] = vals
    return series + ["Other"], folded


def main():
    data = load()
    if not data:
        raise SystemExit(f"No profiles/*{SUFFIX} found.")

    order = sorted(data, key=lambda d: data[d]["n_cells"])

    series, folded = fold(data)
    colors = [OTHER if s == "Other" else SLOTS[i] for i, s in enumerate(series)]

    cells = np.array([data[d]["n_cells"] for d in order], dtype=float)
    mat = np.array([[folded[d][s] for d in order] for s in series])

    fig, ax = plt.subplots(figsize=(9.5, 6.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    ax.stackplot(cells, mat, colors=colors, labels=series,
                 edgecolor=SURFACE, linewidth=2)   # 2px surface gap between bands

    ax.set_xticks(cells)
    ax.set_xticklabels([f"{d}\n{int(c):,} cells" for d, c in zip(order, cells)])
    ax.set_xlim(cells.min(), cells.max())
    ax.set_ylim(0, None)
    ax.set_xlabel("Dataset")
    ax.set_ylabel("Cumulative runtime (seconds)")

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9.5, length=0)
    ax.grid(True, axis="y", color=GRID, linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)

    ax.set_title("Per-section runtime vs dataset size",fontsize=14)

    handles, names = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], names[::-1], loc="upper left", frameon=False,
              fontsize=9.5, labelcolor=INK_2, bbox_to_anchor=(1.01, 1.0))

    fig.tight_layout()
    out = "profiles/section_runtime_vs_dataset.png"
    os.makedirs("profiles", exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.3)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
