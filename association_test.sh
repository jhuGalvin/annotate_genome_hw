#!/usr/bin/env bash

PREFIX="relatedness_qc/HapMap_3_r3_1_relatedness_filtered"
COVAR_FILE="kg_harmonization/mds_covariates.txt"
COVAR_NAMES="C1,C2,C3,C4"
OUTDIR="gwas_assoc"

mkdir -p "${OUTDIR}/plots"

echo "Running basic association test..."
plink \
  --bfile "${PREFIX}" \
  --assoc \
  --allow-no-sex \
  --out "${OUTDIR}/assoc_basic"

echo "Running logistic regression with MDS covariates..."
plink \
  --bfile "${PREFIX}" \
  --covar "${COVAR_FILE}" \
  --covar-name ${COVAR_NAMES} \
  --logistic \
  --allow-no-sex \
  --out "${OUTDIR}/assoc_mds"

echo "Filtering ADD model only..."
awk 'NR==1 || $5=="ADD"' "${OUTDIR}/assoc_mds.assoc.logistic" > "${OUTDIR}/assoc_mds.ADD.txt"

echo "Running multiple testing correction on basic association results..."
plink \
  --bfile "${PREFIX}" \
  --assoc \
  --adjust \
  --allow-no-sex \
  --out "${OUTDIR}/assoc_adjusted"

echo "Creating SNP subset for permutation testing..."
awk 'BEGIN{srand(2026)} {print rand(), $2}' "${PREFIX}.bim" | sort -k1,1g | head -n 5000 | cut -d' ' -f2 > "${OUTDIR}/perm_snps.txt"

echo "Running permutation-based association test on SNP subset..."
plink \
  --bfile "${PREFIX}" \
  --extract "${OUTDIR}/perm_snps.txt" \
  --assoc \
  --mperm 1000 \
  --allow-no-sex \
  --out "${OUTDIR}/assoc_perm"

echo "Generating Manhattan and QQ plots..."
python3 << 'EOF'
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

outdir = Path("gwas_assoc/plots")
outdir.mkdir(parents=True, exist_ok=True)

from matplotlib import font_manager as fm
import matplotlib.pyplot as plt

REG_PATH = "/Users/jtglvn/helvetica_neue_extracted/HelveticaNeue-Regular.ttf"
ITA_PATH = "/Users/jtglvn/helvetica_neue_extracted/HelveticaNeue-Italic.ttf"
BOLD_PATH = "/Users/jtglvn/helvetica_neue_extracted/HelveticaNeue-Bold.ttf"

for fp in [REG_PATH, ITA_PATH, BOLD_PATH]:
    fm.fontManager.addfont(fp)

plt.rcParams.update({
    "font.family": "Helvetica Neue",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
})

def style_axes(ax, remove_bottom=False, x_tick_width=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if remove_bottom:
        ax.spines["bottom"].set_visible(False)

    ax.spines["left"].set_linewidth(1)
    if not remove_bottom:
        ax.spines["bottom"].set_linewidth(1)

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        width=1,
        length=4,
        labelsize=16,
        pad=5
    )

    if x_tick_width is not None:
        ax.tick_params(axis="x", which="both", width=x_tick_width)


def read_assoc(path):
    df = pd.read_csv(path, sep=r"\s+")
    df = df[df["P"].notna()].copy()
    df = df[df["P"] > 0].copy()
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["BP"] = pd.to_numeric(df["BP"], errors="coerce")
    df = df.dropna(subset=["CHR", "BP"])

    # drop PLINK chromosome 25 (XY)
    df = df[df["CHR"] != 25].copy()

    df = df.sort_values(["CHR", "BP"]).copy()
    return df


def manhattan(df, fname):
    df = df.copy()
    df["minus_log10_p"] = -np.log10(df["P"])

    chroms = sorted(df["CHR"].unique())
    x_positions = []
    tick_positions = []
    tick_labels = []

    offset = 0
    gap = 0

    for chrom in chroms:
        sub = df[df["CHR"] == chrom].copy()
        xpos = sub["BP"].values + offset
        x_positions.extend(xpos)

        tick_positions.append((xpos.min() + xpos.max()) / 2)
        if int(chrom) == 23:
            tick_labels.append("X")
        else:
            tick_labels.append(str(int(chrom)))

        offset = xpos.max() + gap

    df["x"] = x_positions

    chrom_colors = ["#9ba9cc", "#445580"]

    fig, ax = plt.subplots(figsize=(12, 5))

    for i, chrom in enumerate(chroms):
        sub = df[df["CHR"] == chrom]
        ax.scatter(
            sub["x"],
            sub["minus_log10_p"],
            c=chrom_colors[i % 2],
            s=8,
            linewidths=0,
            alpha=1,
            rasterized=True,
        )

    ax.axhline(-np.log10(5e-8), color="0.7", linestyle=(0, (2, 2)), linewidth=1.5)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Chromosome", fontsize=20, labelpad=14)
    ax.set_ylabel(r"$-log_{10}(\it{P})$", fontsize=20, labelpad=14)
    ax.set_ylabel(r"$-\log_{10}(\mathit{P})$", fontsize=20, labelpad=14)

    ax.set_ylim(0, 8)
    ax.set_yticks(np.arange(0, 9, 2))

    style_axes(ax, remove_bottom=True, x_tick_width=0)
    
    ax.tick_params(
        axis="x", which="both", direction="out", width=0, length=4, labelsize=12, pad=5
    )
    plt.tight_layout()
    #ax.set_xlim(df["x"].min() - 20, df["x"].max() + 20)
    plt.savefig(outdir / fname, dpi=300, bbox_inches="tight")
    plt.close()


def qqplot(pvals, fname):
    pvals = np.asarray(pvals)
    pvals = pvals[np.isfinite(pvals)]
    pvals = pvals[(pvals > 0) & (pvals <= 1)]
    pvals = np.sort(pvals)

    n = len(pvals)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(pvals)

    lim = max(exp.max(), obs.max())

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(
        exp, obs, s=8, c="0.2", linewidths=0, alpha=1, rasterized=True,
    )

    ax.set_yticks(np.arange(0, lim, 2))
    ax.set_xticks(np.arange(0, lim, 2))

    ax.plot([0, lim], [0, lim], color="0.7", linestyle='--', linewidth=1.2)

    ax.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=20, labelpad=14)
    ax.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=20, labelpad=14)

    style_axes(ax, remove_bottom=False)

    ax.spines["left"].set_position(("outward", 5))
    ax.spines["bottom"].set_position(("outward", 5))
    ax.spines["bottom"].set_bounds(0, lim)
    ax.spines["left"].set_bounds(0, lim)

    plt.tight_layout()
    plt.savefig(outdir / fname, dpi=300, bbox_inches="tight")
    plt.close()

basic = read_assoc("gwas_assoc/assoc_basic.assoc")
manhattan(basic, "manhattan_basic.png")
qqplot(basic["P"], "qq_basic.png")

log = read_assoc("gwas_assoc/assoc_mds.ADD.txt")
manhattan(log, "manhattan_mds.png")
qqplot(log["P"], "qq_mds.png")
EOF

cat > "${OUTDIR}/summary.txt" << 'EOF'

1a. Basic association test:
plink --bfile relatedness_qc/HapMap_3_r3_1_relatedness_filtered --assoc

1b. Logistic regression adjusted for MDS covariates:
plink --bfile relatedness_qc/HapMap_3_r3_1_relatedness_filtered \
      --covar kg_harmonization/mds_covariates.txt \
      --covar-name C1,C2,C3,C4 \
      --logistic

1c. Post-processing:
- Retained only the ADD rows from the logistic regression output.
- Excluded missing or invalid p-values before plotting.
- Interpreted additive SNP effects for downstream comparison.

1d. Why adjusted and unadjusted results may differ:
- The covariate-adjusted model accounts for ancestry/population structure using MDS axes.
- This can reduce confounding and spurious associations caused by stratification.
- Therefore p-values and top SNPs may differ from the unadjusted association test.

2a. Multiple testing correction:
plink --bfile relatedness_qc/HapMap_3_r3_1_relatedness_filtered --assoc --adjust

2b. Permutation testing:
- A subset of 5000 SNPs was sampled from the .bim file.
- Permutation-based association was then run on that subset:
  plink --bfile relatedness_qc/HapMap_3_r3_1_relatedness_filtered \
        --extract gwas_assoc/perm_snps.txt \
        --assoc --mperm 1000
EOF

echo "Done. Results are in ${OUTDIR}/"