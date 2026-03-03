#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import csv

import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from matplotlib.ticker import FuncFormatter

def _tick_no_trailing_zero(x, pos):
    if abs(x - int(x)) < 1e-9:
        return f"{int(x)}"
    return f"{x:g}"

valid_nts = set("ACGT")
nucleotide_order = ["A", "T", "G", "C"]
default_point_color = "#76a89f"

def parse_genotype_file(path: Path) -> dict[str, str]:
    d: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            rsid = parts[0].strip()
            gt = parts[3].strip().upper()
            if rsid:
                d[rsid] = gt
    return d

def load_variant_list_any(path: Path) -> list[str]:
    vids: list[str] = []

    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            r = csv.DictReader(f)
            if r.fieldnames is None:
                raise SystemExit(f"[error] CSV has no header: {path}")

            candidates = ["rsid", "RSID", "Rsid", "snp", "SNP", "marker", "Marker", "id", "ID"]
            rsid_col = None
            for c in candidates:
                if c in r.fieldnames:
                    rsid_col = c
                    break
            if rsid_col is None:
                rsid_col = r.fieldnames[0]

            for row in r:
                v = (row.get(rsid_col) or "").strip()
                if not v:
                    continue
                if v.lower() in ("rsid", "variant_id", "id", "snp"):
                    continue
                vids.append(v)
    else:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if not parts:
                    continue
                v = parts[0].strip()
                if not v:
                    continue
                if v.lower() in ("variant_id", "rsid", "id", "snp"):
                    continue
                vids.append(v)

    seen = set()
    out: list[str] = []
    for v in vids:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out

def is_missing_genotype(gt: str) -> bool:
    if not gt or gt == "--":
        return True
    return not any(ch in valid_nts for ch in gt)

def genotype_to_counts(gt: str) -> dict[str, int]:
    counts = {nt: 0 for nt in nucleotide_order}
    if is_missing_genotype(gt):
        return counts
    for ch in gt:
        if ch in valid_nts:
            counts[ch] += 1
    return counts

def get_prefix(dirname: str) -> str:
    base = Path(dirname).name
    parts = base.split("-")
    if len(parts) >= 2:
        return f"{parts[0]}-{parts[1]}"
    return base

def _is_33_dataset(dataset_name: str) -> bool:
    n = dataset_name.strip()
    return n.startswith("33-") or n.startswith("33_") or n == "33" or n.startswith("33")

def _infer_group_from_stem(stem: str) -> str | None:
    parts = stem.split("_")
    if len(parts) < 3:
        return None

    raw = parts[2].strip().lower()

    if raw == "white":
        return "EUR"
    if raw == "asian":
        return "EAS"
    if raw == "black":
        return "AFR"
    if raw == "blackwhite":
        return "AFR"
    if raw in ("hispanic", "americanindian"):
        return "AMR"

    return None

def _style_axes_left_bottom_only(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.tick_params(axis="both", which="both", direction="out", length=4)

def run_pca_and_plot(matrix_tsv: Path, out_scatter_png: Path, out_var_png: Path, dataset_name: str):
    persons: list[str] = []
    rows: list[list[float]] = []

    with matrix_tsv.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        if not header:
            raise SystemExit(f"[error] Empty matrix file: {matrix_tsv}")

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            persons.append(parts[0])
            rows.append([float(x) for x in parts[1:]])

    X = np.asarray(rows, dtype=float)
    if X.shape[0] < 2:
        raise SystemExit("[error] Need at least 2 individuals for PCA.")

    # Standardize features for PCA
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

    # Fit PCA; keep several PCs for variance plot
    ncomp = min(10, Xs.shape[0], Xs.shape[1])
    pca = PCA(n_components=ncomp, random_state=0)
    pcs = pca.fit_transform(Xs)

    pc1 = pcs[:, 0]
    pc2 = pcs[:, 1]

    is33 = _is_33_dataset(dataset_name)

    pal = {
        "EUR": "#7c3eb3",
        "EAS": "#2871a8",
        "AFR": "#76a879",
        "AMR": "#e3aa00",
    }

    s = 48
    edge_lw = 0.9
    alpha = 0.95

    fig_w = 5 + (1.0 if is33 else 0.0)
    fig_h = 5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.xaxis.set_major_formatter(FuncFormatter(_tick_no_trailing_zero))
    ax.yaxis.set_major_formatter(FuncFormatter(_tick_no_trailing_zero))

    if is33:
        groups: list[str] = []
        colors: list[str] = []
        for name in persons:
            g = _infer_group_from_stem(name) or "EUR"
            groups.append(g)
            colors.append(pal[g])

        ax.scatter(
            pc1, pc2,
            s=s,
            c=colors,
            edgecolors="white",
            linewidth=edge_lw,
            alpha=alpha,
            zorder=2,
        )

        from matplotlib.lines import Line2D
        handles = [
            Line2D(
                [0], [0],
                marker="o",
                linestyle="None",
                markerfacecolor=pal[k],
                markeredgecolor="white",
                markeredgewidth=edge_lw,
                markersize=7,
                label=k,
            )
            for k in ["EUR", "EAS", "AFR", "AMR"]
        ]
        leg = ax.legend(
            title='Self-reported\nglobal population',
            handles=handles,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )
        leg.get_title().set_ha("right")
        leg._legend_box.align = "right"
    else:
        ax.scatter(
            pc1, pc2, s=s, c=default_point_color,
            edgecolors="white", linewidth=edge_lw, alpha=alpha, zorder=2,
        )

    _style_axes_left_bottom_only(ax)
    ax.set_xlabel("PC1", fontsize=12)
    ax.set_ylabel("PC2", fontsize=12)

    plt.tight_layout()
    plt.savefig(out_scatter_png, dpi=300)
    plt.close(fig)

    var = pca.explained_variance_ratio_
    fig2, ax2 = plt.subplots(figsize=(5, 3))
    xs = np.arange(1, len(var) + 1)
    ax2.bar(xs, var)

    _style_axes_left_bottom_only(ax2)
    ax2.set_xlabel("Principal component")
    ax2.set_ylabel("Variance explained")
    ax2.xaxis.set_major_formatter(FuncFormatter(_tick_no_trailing_zero))
    ax2.yaxis.set_major_formatter(FuncFormatter(_tick_no_trailing_zero))
    ax2.set_xticks(xs)

    plt.tight_layout()
    plt.savefig(out_var_png, dpi=300)
    plt.close(fig2)

def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} <genotype_directory>")

    genotype_dir = Path(sys.argv[1])
    prefix = get_prefix(genotype_dir.name)

    out_matrix = Path(f"{prefix}_genotype_matrix.tsv")
    out_hist = Path(f"{prefix}_missing_hist.png")
    out_zero = Path(f"{prefix}_zero_coverage_snps.tsv")
    out_pca = Path(f"{prefix}_pca_pc1_pc2.png")
    out_pca_var = Path(f"{prefix}_pca_variance_explained.png")

    snp_file = Path("Seldin_AIMS_SNPs.csv")
    files = sorted([p for p in genotype_dir.iterdir() if p.is_file() and p.suffix == ".txt"])
    individuals = [(p.stem, parse_genotype_file(p)) for p in files]
    n_individuals = len(individuals)

    # load seldin snps
    variant_ids = load_variant_list_any(snp_file)
    if not variant_ids:
        raise SystemExit(f"[error] SNP list empty: {snp_file}")
    n_snps = len(variant_ids)

    # missingness/no coverage snps
    missing_per_variant: list[int] = []
    zero_coverage: list[str] = []

    for vid in variant_ids:
        miss = sum(1 for _, m in individuals if is_missing_genotype(m.get(vid, "")))
        missing_per_variant.append(miss)
        if miss == n_individuals:
            zero_coverage.append(vid)

    # Write gtm
    header = ["Person"] + [f"{vid}{nt}" for vid in variant_ids for nt in nucleotide_order]
    with out_matrix.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for iid, m in individuals:
            row = [iid]
            for vid in variant_ids:
                counts = genotype_to_counts(m.get(vid, ""))
                for nt in nucleotide_order:
                    row.append(str(counts[nt]))
            w.writerow(row)

    # report zero-coverage SNPs
    with out_zero.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rsid"])
        for vid in zero_coverage:
            w.writerow([vid])

    # missingness histogram
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.hist(missing_per_variant, bins=30)
    _style_axes_left_bottom_only(ax)
    ax.set_xlabel("Number of individuals missing genotype (per SNP)")
    ax.set_ylabel("Number of SNPs")
    plt.tight_layout()
    plt.savefig(out_hist, dpi=300)
    plt.close(fig)

    run_pca_and_plot(out_matrix, out_pca, out_pca_var, dataset_name=genotype_dir.name)

    n_rows = n_individuals
    n_cols = 1 + 4 * n_snps

    xs = sorted(missing_per_variant)

    def quantile(p: float):
        if not xs:
            return None
        i = int(round(p * (len(xs) - 1)))
        return xs[i]

    thr = max(1, int(round(0.10 * n_individuals)))
    n_good = sum(1 for x in missing_per_variant if x <= thr)

    print(f"[dir] {genotype_dir}  prefix={prefix}")
    print(f"[individuals] {n_individuals}")
    print(f"[snps] {n_snps} (from {snp_file})")
    print(f"[matrix] {out_matrix}  size={n_rows} x {n_cols}")
    print(f"[histogram] {out_hist}")
    print(f"[zero_coverage] {len(zero_coverage)} SNPs -> {out_zero}")
    print(f"[pca] {out_pca}")
    print(f"[pca_var] {out_pca_var}")
    print(f"[missing] min={xs[0]} median={quantile(0.5)} p90={quantile(0.9)} max={xs[-1]}")
    print(f"[coverage] <= 10% missing (<= {thr} people): {n_good}/{n_snps} = {n_good/n_snps:.3f}")

    if n_good / n_snps >= 0.5:
        print("Majority of SNPs have good coverage.")
    else:
        print("Majority of SNPs do not have good coverage.")

if __name__ == "__main__":
    main()