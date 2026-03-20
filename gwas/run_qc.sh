#!/usr/bin/env bash

set -euo pipefail

PREFIX="HapMap_3_r3_1"

########################################
# Part 1: Missingness
########################################
OUTDIR="missingness_qc"
OUTPREFIX="${OUTDIR}/${PREFIX}_missing"

mkdir -p "${OUTDIR}"

echo "Computing missingness with PLINK..."
plink --bfile "${PREFIX}" --missing --out "${OUTPREFIX}"

echo "Plotting histograms..."
python <<'PY'
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FuncFormatter

prefix = "HapMap_3_r3_1"
outdir = Path("missingness_qc")
outprefix = outdir / f"{prefix}_missing"

BAR_COLOR = "#70b8a1"

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.bottom"] = False

def pct_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:g}"

imiss = pd.read_csv(f"{outprefix}.imiss", sep=r"\s+")
lmiss = pd.read_csv(f"{outprefix}.lmiss", sep=r"\s+")

imiss_pct = imiss["F_MISS"] * 100.0
lmiss_pct = lmiss["F_MISS"] * 100.0

# Individual missingness
fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(imiss_pct, bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1)
ax.set_xlabel("SNPs missing in dataset (%)", fontsize=20, labelpad=14)
ax.set_ylabel("Number of individuals", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1.0)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_major_formatter(FuncFormatter(pct_tick))
ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

plt.tight_layout()
plt.savefig(outdir / "individual_missingness_hist.png", dpi=300, bbox_inches="tight")
plt.close()

# SNP missingness
fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(lmiss_pct, bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1.0)

ax.set_xlabel("Individuals missing in dataset (%)", fontsize=20, labelpad=14)
ax.set_ylabel("Number of SNPs (×10⁶)", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1.0)
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_major_formatter(FuncFormatter(pct_tick))

ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

yticks = [0, 2e5, 4e5, 6e5, 8e5, 1e6]
ax.set_yticks(yticks)
ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])

plt.tight_layout()
plt.savefig(outdir / "snp_missingness_hist.png", dpi=300, bbox_inches="tight")
plt.close()
PY

########################################
# Part 2: Sex check
########################################
SEXDIR="sexcheck_qc"
mkdir -p "${SEXDIR}"

echo "Running sex check..."

plink --bfile "${PREFIX}" --chr X --make-bed --out "${SEXDIR}/${PREFIX}_chrX"
plink --bfile "${SEXDIR}/${PREFIX}_chrX" --check-sex --out "${SEXDIR}/${PREFIX}_sexcheck"

awk 'NR>1 && $5 == "PROBLEM" {print $1, $2}' \
  "${SEXDIR}/${PREFIX}_sexcheck.sexcheck" \
  > "${SEXDIR}/flagged_samples.txt"

echo "Plotting F histograms..."

python <<'PY'
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

prefix = "HapMap_3_r3_1"
outdir = Path("sexcheck_qc")

BAR_COLOR = "#70b8a1"

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.bottom"] = False

def clean_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:g}"

df = pd.read_csv(outdir / f"{prefix}_sexcheck.sexcheck", sep=r"\s+")

def plot(df_sub, fname):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.hist(
        df_sub["F"].dropna(), bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1
    )
    ax.set_xlabel("Inbreeding coefficient", fontsize=20, labelpad=14)
    ax.set_ylabel("Number of samples", fontsize=20, labelpad=14)

    ax.spines["left"].set_linewidth(1.0)
    ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
    ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
    ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

    plt.tight_layout()
    plt.savefig(outdir / fname, dpi=300, bbox_inches="tight")
    plt.close()

# all samples
plot(df, "F_all.png")

# reported males (PEDSEX == 1)
plot(df[df["PEDSEX"] == 1], "F_males.png")

# reported females (PEDSEX == 2)
plot(df[df["PEDSEX"] == 2], "F_females.png")
PY

########################################
# Part 2c: Conditional imputation
########################################
N_PROBLEMS=$(awk 'NR>1 && $5 == "PROBLEM" {count++} END {print count+0}' \
  "${SEXDIR}/${PREFIX}_sexcheck.sexcheck")

echo "Number of sex discrepancies: ${N_PROBLEMS}"

echo "Flagged sex discrepancy sample(s):"
awk 'NR==1 || $5 == "PROBLEM"' \
  "${SEXDIR}/${PREFIX}_sexcheck.sexcheck"

if [ "${N_PROBLEMS}" -gt 0 ]; then
  echo "Imputing sex..."

  plink \
    --bfile "${SEXDIR}/${PREFIX}_chrX" \
    --impute-sex \
    --make-bed \
    --out "${SEXDIR}/${PREFIX}_chrX_imputed"

  awk '{print $1, $2, $5}' \
    "${SEXDIR}/${PREFIX}_chrX_imputed.fam" \
    > "${SEXDIR}/imputed_sex_update.txt"

  plink \
    --bfile "${PREFIX}" \
    --update-sex "${SEXDIR}/imputed_sex_update.txt" \
    --make-bed \
    --out "${SEXDIR}/${PREFIX}_updated_sex"

else
  echo "No discrepancies found. Skipping imputation."
fi

########################################
# Part 3: Minor Allele Frequency (MAF) Filtering
########################################
MAFDIR="maf_qc"
mkdir -p "${MAFDIR}"

if [ -f "${SEXDIR}/${PREFIX}_updated_sex.bed" ]; then
  INPUT_PREFIX="${SEXDIR}/${PREFIX}_updated_sex"
else
  INPUT_PREFIX="${PREFIX}"
fi

MAF_THRESHOLD="0.01"
PRE_MAF_PREFIX="${MAFDIR}/${PREFIX}_pre_maf"
POST_MAF_PREFIX="${MAFDIR}/${PREFIX}_post_maf"
FILTERED_PREFIX="${MAFDIR}/${PREFIX}_maf_filtered"

echo "Computing allele frequencies before MAF filtering..."
plink \
  --bfile "${INPUT_PREFIX}" \
  --freq \
  --out "${PRE_MAF_PREFIX}"

echo "Filtering SNPs with MAF < ${MAF_THRESHOLD}..."
plink \
  --bfile "${INPUT_PREFIX}" \
  --maf "${MAF_THRESHOLD}" \
  --make-bed \
  --out "${FILTERED_PREFIX}"

echo "Computing allele frequencies after MAF filtering..."
plink \
  --bfile "${FILTERED_PREFIX}" \
  --freq \
  --out "${POST_MAF_PREFIX}"

echo "Plotting MAF histograms..."
python <<'PY'
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

prefix = "HapMap_3_r3_1"
outdir = Path("maf_qc")

pre = pd.read_csv(outdir / f"{prefix}_pre_maf.frq", sep=r"\s+")
post = pd.read_csv(outdir / f"{prefix}_post_maf.frq", sep=r"\s+")

BAR_COLOR = "#70b8a1"

def thousands_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{int(x/1000)}"

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.bottom"] = False

def clean_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:g}"

pre_maf = pre["MAF"].dropna()
post_maf = post["MAF"].dropna()

# Before filtering
fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(pre_maf, bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1)
ax.set_xlabel("Minor allele frequency", fontsize=20, labelpad=14)
ax.set_ylabel("Number of SNPs (×10³)", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1)
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
ax.yaxis.set_major_formatter(FuncFormatter(thousands_tick))
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

plt.tight_layout()
plt.savefig(outdir / "maf_hist_before_filtering.png", dpi=300, bbox_inches="tight")
plt.close()

# After filtering
fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(post_maf, bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1)
ax.set_xlabel("Minor allele frequency", fontsize=20, labelpad=14)
ax.set_ylabel("Number of SNPs (×10³)", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1)
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
ax.yaxis.set_major_formatter(FuncFormatter(thousands_tick))
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

plt.tight_layout()
plt.savefig(outdir / "maf_hist_after_filtering.png", dpi=300, bbox_inches="tight")
plt.close()
PY

echo "MAF filtering complete."
echo "Threshold used: ${MAF_THRESHOLD}"

########################################
# Part 4: Hardy-Weinberg Equilibrium (HWE) Filtering
########################################
HWEDIR="hwe_qc"
mkdir -p "${HWEDIR}"

HWE_INPUT_PREFIX="${MAFDIR}/${PREFIX}_maf_filtered"
HWE_THRESHOLD="1e-6"

PRE_HWE_PREFIX="${HWEDIR}/${PREFIX}_pre_hwe"
FILTERED_HWE_PREFIX="${HWEDIR}/${PREFIX}_hwe_filtered"
POST_HWE_PREFIX="${HWEDIR}/${PREFIX}_post_hwe"

echo "Computing HWE p-values before filtering..."
plink \
  --bfile "${HWE_INPUT_PREFIX}" \
  --hardy \
  --out "${PRE_HWE_PREFIX}"

echo "Filtering SNPs failing HWE (p < ${HWE_THRESHOLD})..."
plink \
  --bfile "${HWE_INPUT_PREFIX}" \
  --hwe "${HWE_THRESHOLD}" \
  --make-bed \
  --out "${FILTERED_HWE_PREFIX}"

echo "Computing HWE p-values after filtering..."
plink \
  --bfile "${FILTERED_HWE_PREFIX}" \
  --hardy \
  --out "${POST_HWE_PREFIX}"

echo "Plotting HWE p-value histogram after filtering..."
python <<'PY'
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

prefix = "HapMap_3_r3_1"
outdir = Path("hwe_qc")

BAR_COLOR = "#70b8a1"

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.bottom"] = False

def clean_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:g}"

hwe = pd.read_csv(outdir / f"{prefix}_post_hwe.hwe", sep=r"\s+")

# Keep only the SNP test rows
hwe = hwe[hwe["TEST"] == "ALL"].copy()

pvals = pd.to_numeric(hwe["P"], errors="coerce").dropna()

fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(pvals, bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1)

ax.set_xlabel(r"HWE $P$ value", fontsize=20, labelpad=14)
ax.set_ylabel("Number of SNPs (×10³)", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1.0)
ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)

def thousands_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{int(x/1000)}"

ax.yaxis.set_major_formatter(FuncFormatter(thousands_tick))
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

plt.tight_layout()
plt.savefig(outdir / "hwe_pval_hist_after_filtering.png", dpi=300, bbox_inches="tight")
plt.close()
PY

echo "HWE filtering complete."
echo "HWE threshold used: ${HWE_THRESHOLD}"

########################################
# Part 5: Heterozygosity Rate Filtering
########################################
HETDIR="heterozygosity_qc"
mkdir -p "${HETDIR}"

HET_INPUT_PREFIX="${HWEDIR}/${PREFIX}_hwe_filtered"

PRUNE_PREFIX="${HETDIR}/${PREFIX}_pruned"
HET_PREFIX="${HETDIR}/${PREFIX}_het"
OUTLIER_FILE="${HETDIR}/heterozygosity_outliers.txt"
FILTERED_PREFIX="${HETDIR}/${PREFIX}_het_filtered"
POST_HET_PREFIX="${HETDIR}/${PREFIX}_post_het"

echo "LD-pruning SNPs for heterozygosity estimation..."
plink \
  --bfile "${HET_INPUT_PREFIX}" \
  --indep-pairwise 50 5 0.2 \
  --out "${PRUNE_PREFIX}"

echo "Computing heterozygosity rates..."
plink \
  --bfile "${HET_INPUT_PREFIX}" \
  --extract "${PRUNE_PREFIX}.prune.in" \
  --het \
  --out "${HET_PREFIX}"

echo "Identifying heterozygosity outliers..."
python <<'PY'
import pandas as pd
from pathlib import Path

prefix = "HapMap_3_r3_1"
outdir = Path("heterozygosity_qc")

het = pd.read_csv(outdir / f"{prefix}_het.het", sep=r"\s+")

# Heterozygosity rate
het["HET_RATE"] = (het["N(NM)"] - het["O(HOM)"]) / het["N(NM)"]

mean_het = het["HET_RATE"].mean()
sd_het = het["HET_RATE"].std()

lower = mean_het - 3 * sd_het
upper = mean_het + 3 * sd_het

outliers = het[(het["HET_RATE"] < lower) | (het["HET_RATE"] > upper)].copy()
outliers[["FID", "IID"]].to_csv(
    outdir / "heterozygosity_outliers.txt",
    sep="\t",
    index=False,
    header=False,
)

summary = pd.DataFrame({
    "mean_het": [mean_het],
    "sd_het": [sd_het],
    "lower_bound": [lower],
    "upper_bound": [upper],
    "n_outliers": [len(outliers)],
})
summary.to_csv(outdir / "heterozygosity_summary.tsv", sep="\t", index=False)
PY

N_HET_OUTLIERS=$(wc -l < "${OUTLIER_FILE}" || echo 0)
echo "Number of heterozygosity outliers: ${N_HET_OUTLIERS}"

echo "Removed individual(s) (FID IID):"
if [ "${N_HET_OUTLIERS}" -gt 0 ]; then
  cat "${OUTLIER_FILE}"

  echo "Removing heterozygosity outliers..."
  plink \
    --bfile "${HET_INPUT_PREFIX}" \
    --remove "${OUTLIER_FILE}" \
    --make-bed \
    --out "${FILTERED_PREFIX}"
else
  echo "None"
  plink \
    --bfile "${HET_INPUT_PREFIX}" \
    --make-bed \
    --out "${FILTERED_PREFIX}"
fi

echo "Recomputing heterozygosity after filtering..."
plink \
  --bfile "${FILTERED_PREFIX}" \
  --extract "${PRUNE_PREFIX}.prune.in" \
  --het \
  --out "${POST_HET_PREFIX}"

echo "Plotting heterozygosity histogram after filtering..."
python <<'PY'
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

prefix = "HapMap_3_r3_1"
outdir = Path("heterozygosity_qc")

BAR_COLOR = "#70b8a1"

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.bottom"] = False

def clean_tick(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:g}"

het = pd.read_csv(outdir / f"{prefix}_post_het.het", sep=r"\s+")
het["HET_RATE"] = (het["N(NM)"] - het["O(HOM)"]) / het["N(NM)"]

fig, ax = plt.subplots(figsize=(6, 5))
ax.hist(het["HET_RATE"].dropna(), bins=30, color=BAR_COLOR, edgecolor="white", linewidth=1)

ax.set_xlabel("Heterozygosity rate", fontsize=20, labelpad=14)
ax.set_ylabel("Number of individuals", fontsize=20, labelpad=14)

ax.spines["left"].set_linewidth(1.0)
ax.xaxis.set_major_formatter(FuncFormatter(clean_tick))
ax.tick_params(axis="x", length=4, width=0, labelsize=16, pad=5)
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

from matplotlib.ticker import MultipleLocator
ax.yaxis.set_major_locator(MultipleLocator(4))
ax.tick_params(axis="y", length=4, width=1, labelsize=16, pad=5)

plt.tight_layout()
plt.savefig(outdir / "heterozygosity_hist_after_filtering.png", dpi=300, bbox_inches="tight")
plt.close()
PY

echo "Heterozygosity filtering complete."


########################################
# Part 6: Relatedness Filtering
########################################
RELDIR="relatedness_qc"
mkdir -p "${RELDIR}"

REL_INPUT_PREFIX="${HETDIR}/${PREFIX}_het_filtered"

REL_PRUNE_PREFIX="${RELDIR}/${PREFIX}_pruned"
REL_GENOME_PREFIX="${RELDIR}/${PREFIX}_relatedness"
REL_PAIRS_FILE="${RELDIR}/related_pairs.txt"
REL_REMOVE_FILE="${RELDIR}/related_samples_to_remove.txt"
REL_FILTERED_PREFIX="${RELDIR}/${PREFIX}_relatedness_filtered"

PIHAT_THRESHOLD="0.185"

echo "LD-pruning SNPs for relatedness estimation..."
plink \
  --bfile "${REL_INPUT_PREFIX}" \
  --indep-pairwise 50 5 0.2 \
  --out "${REL_PRUNE_PREFIX}"

echo "Computing pairwise relatedness..."
plink \
  --bfile "${REL_INPUT_PREFIX}" \
  --extract "${REL_PRUNE_PREFIX}.prune.in" \
  --genome \
  --out "${REL_GENOME_PREFIX}"

echo "Extracting related pairs with PI_HAT > ${PIHAT_THRESHOLD}..."
awk -v thr="${PIHAT_THRESHOLD}" '
  NR==1 || $10 > thr
' "${REL_GENOME_PREFIX}.genome" > "${REL_PAIRS_FILE}"

echo "Choosing one sample from each related pair to remove..."
awk -v thr="${PIHAT_THRESHOLD}" '
  NR > 1 && $10 > thr {print $1, $2}
' "${REL_GENOME_PREFIX}.genome" | sort -u > "${REL_REMOVE_FILE}"

N_REMOVE=$(wc -l < "${REL_REMOVE_FILE}" || echo 0)
echo "Number of samples marked for removal: ${N_REMOVE}"

echo "Removed related sample(s) (FID IID):"
if [ "${N_REMOVE}" -gt 0 ]; then
  cat "${REL_REMOVE_FILE}"

  echo "Removing related individuals..."
  plink \
    --bfile "${REL_INPUT_PREFIX}" \
    --remove "${REL_REMOVE_FILE}" \
    --make-bed \
    --out "${REL_FILTERED_PREFIX}"
else
  echo "None"
  plink \
    --bfile "${REL_INPUT_PREFIX}" \
    --make-bed \
    --out "${REL_FILTERED_PREFIX}"
fi

echo "Computing final sample and SNP counts..."
python <<'PY'
from pathlib import Path

prefix = "HapMap_3_r3_1"
rel_filtered = Path("relatedness_qc") / f"{prefix}_relatedness_filtered"

fam_file = rel_filtered.with_suffix(".fam")
bim_file = rel_filtered.with_suffix(".bim")

n_samples = sum(1 for _ in open(fam_file))
n_snps = sum(1 for _ in open(bim_file))

with open(Path("relatedness_qc") / "final_counts.txt", "w") as f:
    f.write(f"Final number of samples: {n_samples}\n")
    f.write(f"Final number of SNPs: {n_snps}\n")

print(f"Final number of samples: {n_samples}")
print(f"Final number of SNPs: {n_snps}")
PY

echo "Relatedness filtering complete."