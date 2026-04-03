#!/usr/bin/env bash
set -euo pipefail

BFILE="../gwas_assignment/HapMap_3_r3_1"
SUMSTATS="discovery_gwas_sim_sumstats_cleaned.tsv"
WEIGHTS="prs_weights_harmonized.tsv"
PHENO="pheno.tsv"
COVAR="../gwas_assignment/kg_harmonization/mds_covariates.txt"

OUTDIR="prs_ct_results"
mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/scores"
mkdir -p "${OUTDIR}/weights"
mkdir -p "${OUTDIR}/snplists"
mkdir -p "${OUTDIR}/plots"

echo "=== Creating phenotype file ==="

{
  printf "FID\tIID\tPHENO\n"
  awk '$6 != -9 {print $1 "\t" $2 "\t" $6}' "${BFILE}.fam"
} > "${PHENO}"

echo "Phenotype file created:"
head "${PHENO}"

for f in "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam" "${SUMSTATS}" "${WEIGHTS}" "${PHENO}" "${COVAR}"; do
    if [[ ! -f "$f" ]]; then
        echo "Missing required file: $f" >&2
        exit 1
    fi
done

# -----------------------------
# 4.1 LD clumping
# -----------------------------
echo "=== Step 4.1: LD clumping ==="

plink \
  --bfile "${BFILE}" \
  --clump "${SUMSTATS}" \
  --clump-snp-field SNP \
  --clump-field P \
  --clump-p1 1 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out "${OUTDIR}/clump"

awk 'NR>1 && $3 != "" {print $3}' "${OUTDIR}/clump.clumped" > "${OUTDIR}/clumped.snps"

echo "Number of clumped SNPs:"
wc -l "${OUTDIR}/clumped.snps"

# -----------------------------
# Filter harmonized weights to clumped SNPs
# -----------------------------
awk 'NR==FNR {keep[$1]=1; next} FNR==1 || ($1 in keep)' \
  "${OUTDIR}/clumped.snps" "${WEIGHTS}" > "${OUTDIR}/prs_weights_harmonized_clumped.tsv"

echo "Rows in clumped harmonized weight file:"
wc -l "${OUTDIR}/prs_weights_harmonized_clumped.tsv"

# -----------------------------
# 4.2 C+T p-value thresholds
# -----------------------------
echo "=== Step 4.2: C+T scoring across p-value thresholds ==="

THRESHOLDS=(
  "5e-8" "1e-6" "1e-4" "1e-3" "1e-2" "0.05" "0.1" "0.5" "1"
)

for PT in "${THRESHOLDS[@]}"; do
    SAFE_PT=$(echo "${PT}" | tr '.' '_' )
    WFILE="${OUTDIR}/weights/prs_weights.PT_${SAFE_PT}.tsv"
    SNPFILE="${OUTDIR}/snplists/prs_${SAFE_PT}.snps"
    OUTPREFIX="${OUTDIR}/scores/prs.PT_${SAFE_PT}"

    awk -v pt="${PT}" 'BEGIN{FS=OFS="\t"} NR==1 || $4 <= pt' \
      "${OUTDIR}/prs_weights_harmonized_clumped.tsv" > "${WFILE}"

    awk 'NR>1 {print $1}' "${WFILE}" > "${SNPFILE}"

    NROWS=$(($(wc -l < "${WFILE}") - 1))
    NSNPS=$(wc -l < "${SNPFILE}")
    echo "PT=${PT} SNPs=${NROWS} list_entries=${NSNPS}"

    if [[ "${NROWS}" -le 0 ]]; then
        echo "Skipping PT=${PT}, no SNPs retained"
        continue
    fi

    plink \
      --bfile "${BFILE}" \
      --score "${WFILE}" 1 2 3 header \
      --out "${OUTPREFIX}"
done

# -----------------------------
# 4.3 Logistic regression in Python
# -----------------------------
echo "=== Step 4.3: Logistic regression with MDS covariates ==="

python - <<'PYCODE'
import os
import glob
import math
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

outdir = "prs_ct_results"
score_dir = os.path.join(outdir, "scores")
snplist_dir = os.path.join(outdir, "snplists")
plot_dir = os.path.join(outdir, "plots")

pheno_file = "pheno.tsv"
covar_file = "../gwas_assignment/kg_harmonization/mds_covariates.txt"

# -----------------------------
# Helpers
# -----------------------------
def parse_pt_label(pt_str):
    s = str(pt_str).strip()
    try:
        return float(s)
    except ValueError:
        return np.nan

def format_pt_label(x):
    if x >= 1e-2:
        return f"{x:g}"
    return f"{x:.0e}".replace("e-0", "e-").replace("e+0", "e+")

# -----------------------------
# Load phenotype and covariates
# -----------------------------
pheno = pd.read_csv(pheno_file, sep="\t")
covar = pd.read_csv(covar_file, sep="\t")

pheno["FID"] = pheno["FID"].astype(str)
pheno["IID"] = pheno["IID"].astype(str)

covar["FID"] = covar["FID"].astype(str)
covar["IID"] = covar["IID"].astype(str)

required_pheno = ["FID", "IID", "PHENO"]
required_covar = ["FID", "IID", "C1", "C2", "C3", "C4"]

missing_pheno = [c for c in required_pheno if c not in pheno.columns]
missing_covar = [c for c in required_covar if c not in covar.columns]

if missing_pheno:
    raise ValueError(f"Missing phenotype columns: {missing_pheno}")
if missing_covar:
    raise ValueError(f"Missing covariate columns: {missing_covar}")

# PLINK convention: 1 = control, 2 = case
pheno = pheno[pheno["PHENO"].isin([1, 2])].copy()
pheno["CASE"] = (pheno["PHENO"] == 2).astype(int)

results = []

profile_files = sorted(glob.glob(os.path.join(score_dir, "prs.PT_*.profile")))
if len(profile_files) == 0:
    raise ValueError("No .profile files found in prs_ct_results/scores")

for pf in profile_files:
    prs = pd.read_csv(pf, sep=r"\s+")
    prs["FID"] = prs["FID"].astype(str)
    prs["IID"] = prs["IID"].astype(str)

    score_col = None
    for candidate in ["SCORE", "SCORESUM", "SCORE1_AVG", "SCORE1_SUM"]:
        if candidate in prs.columns:
            score_col = candidate
            break

    if score_col is None:
        raise ValueError(
            f"Could not find PRS score column in {pf}. "
            f"Columns were: {prs.columns.tolist()}"
        )

    df = prs.merge(
        pheno[["FID", "IID", "PHENO", "CASE"]], on=["FID", "IID"], how="inner"
    )
    df = df.merge(
        covar[["FID", "IID", "C1", "C2", "C3", "C4"]], on=["FID", "IID"], how="inner"
    )

    prs_std = df[score_col].std(ddof=0)
    if prs_std == 0 or pd.isna(prs_std):
        print(f"Skipping {pf}: PRS variance is zero")
        continue

    df["PRS_Z"] = (df[score_col] - df[score_col].mean()) / prs_std

    X = df[["PRS_Z", "C1", "C2", "C3", "C4"]].copy()
    X = sm.add_constant(X)
    y = df["CASE"]

    try:
        model = sm.Logit(y, X).fit(disp=False)
    except Exception as e:
        print(f"Skipping {pf}: logistic regression failed with error: {e}")
        continue

    beta = model.params["PRS_Z"]
    se = model.bse["PRS_Z"]
    p = model.pvalues["PRS_Z"]
    or_val = np.exp(beta)
    lci = np.exp(beta - 1.96 * se)
    uci = np.exp(beta + 1.96 * se)

    pseudo_r2 = getattr(model, "prsquared", np.nan)

    n = len(df)
    n_case = int(df["CASE"].sum())
    n_control = int((df["CASE"] == 0).sum())

    pt = os.path.basename(pf).replace("prs.PT_", "").replace(".profile", "")
    pt = pt.replace("_", ".")

    safe_pt = os.path.basename(pf).replace("prs.PT_", "").replace(".profile", "")
    snpfile = os.path.join(snplist_dir, f"prs_{safe_pt}.snps")
    n_snps = 0
    if os.path.exists(snpfile):
        with open(snpfile) as fh:
            n_snps = sum(1 for _ in fh)

    results.append({
        "PT": pt,
        "PT_num": parse_pt_label(pt),
        "N_SNPs": n_snps,
        "N": n,
        "N_case": n_case,
        "N_control": n_control,
        "Beta_PRS_Z": beta,
        "SE_PRS_Z": se,
        "OR_per_SD": or_val,
        "OR_L95": lci,
        "OR_U95": uci,
        "P_PRS_Z": p,
        "minus_log10_P": -np.log10(p) if p > 0 else np.inf,
        "Pseudo_R2": pseudo_r2,
        "Score_col": score_col
    })

results_df = pd.DataFrame(results)

if len(results_df) == 0:
    raise ValueError("No successful logistic regression results were produced.")

# numeric sort for plotting / reporting
results_df = results_df.sort_values("PT_num", ascending=True).reset_index(drop=True)

results_file = os.path.join(outdir, "prs_logistic_results.tsv")
results_df.to_csv(results_file, sep="\t", index=False)

# best threshold by smallest p-value
best_idx = results_df["P_PRS_Z"].astype(float).idxmin()
best = results_df.loc[best_idx]

with open(os.path.join(outdir, "best_prs_result.txt"), "w") as f:
    f.write("Best PRS threshold result\n")
    f.write(f"PT\t{best['PT']}\n")
    f.write(f"N_SNPs\t{best['N_SNPs']}\n")
    f.write(f"N\t{best['N']}\n")
    f.write(f"N_case\t{best['N_case']}\n")
    f.write(f"N_control\t{best['N_control']}\n")
    f.write(f"Beta_PRS_Z\t{best['Beta_PRS_Z']}\n")
    f.write(f"SE_PRS_Z\t{best['SE_PRS_Z']}\n")
    f.write(f"OR_per_SD\t{best['OR_per_SD']}\n")
    f.write(f"OR_L95\t{best['OR_L95']}\n")
    f.write(f"OR_U95\t{best['OR_U95']}\n")
    f.write(f"P_PRS_Z\t{best['P_PRS_Z']}\n")
    f.write(f"minus_log10_P\t{best['minus_log10_P']}\n")
    f.write(f"Pseudo_R2\t{best['Pseudo_R2']}\n")

# short interpretation text
sig_text = "statistically significant" if best["P_PRS_Z"] < 0.05 else "not statistically significant"
direction_text = "higher" if best["Beta_PRS_Z"] > 0 else "lower"

with open(os.path.join(outdir, "prs_interpretation.txt"), "w") as f:
    f.write("PRS assignment interpretation\n")
    f.write("============================\n\n")
    f.write(f"Best-performing p-value threshold: {best['PT']}\n")
    f.write(f"Number of SNPs retained at this threshold: {int(best['N_SNPs'])}\n")
    f.write(f"Sample size used in regression: {int(best['N'])} "
            f"({int(best['N_case'])} cases, {int(best['N_control'])} controls)\n")
    f.write(f"PRS beta (per SD): {best['Beta_PRS_Z']:.6g}\n")
    f.write(f"Odds ratio per 1 SD increase in PRS: {best['OR_per_SD']:.6g} "
            f"(95% CI {best['OR_L95']:.6g} to {best['OR_U95']:.6g})\n")
    f.write(f"P-value for PRS: {best['P_PRS_Z']:.6g}\n")
    f.write(f"McFadden pseudo-R2: {best['Pseudo_R2']:.6g}\n\n")
    f.write(
        f"Interpretation: at the best threshold, a 1 SD {direction_text} PRS is associated "
        f"with an odds ratio of {best['OR_per_SD']:.4f} for case status. "
        f"This result is {sig_text} at alpha = 0.05.\n"
    )

# -----------------------------
# Plot: dot + line
# x-axis = p-value threshold
# y-axis = -log10(p)
# -----------------------------
plot_df = results_df.dropna(subset=["PT_num", "minus_log10_P"]).copy()
plot_df = plot_df.sort_values("PT_num")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(plot_df["PT_num"], plot_df["minus_log10_P"], marker="o", linewidth=1.5)
ax.set_xscale("log")
ax.set_xlabel("P-value threshold")
ax.set_ylabel("-log10(PRS p-value)")
ax.set_title("PRS performance across C+T thresholds")

# horizontal line at p = 0.05
sig_line = -math.log10(0.05)
ax.axhline(sig_line, linestyle="--", linewidth=1)

# mark best threshold
ax.scatter([best["PT_num"]], [best["minus_log10_P"]], s=60, zorder=3)
ax.annotate(
    f"best PT={best['PT']}",
    xy=(best["PT_num"], best["minus_log10_P"]),
    xytext=(8, 8),
    textcoords="offset points"
)

# use clean threshold tick labels
xticks = plot_df["PT_num"].tolist()
xlabels = [format_pt_label(x) for x in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=45, ha="right")

fig.tight_layout()
plot_path = os.path.join(plot_dir, "prs_threshold_performance.png")
fig.savefig(plot_path, dpi=300)
plt.close(fig)

print("\n=== PRS logistic regression results ===")
print(results_df.sort_values('P_PRS_Z', ascending=True).to_string(index=False))
print(f"\nWrote: {results_file}")
print(f"Wrote: {os.path.join(outdir, 'best_prs_result.txt')}")
print(f"Wrote: {os.path.join(outdir, 'prs_interpretation.txt')}")
print(f"Wrote: {plot_path}")
PYCODE

echo "=== Done ==="
echo "Main outputs:"
echo "  ${OUTDIR}/clump.clumped"
echo "  ${OUTDIR}/clumped.snps"
echo "  ${OUTDIR}/prs_weights_harmonized_clumped.tsv"
echo "  ${OUTDIR}/snplists/prs_*.snps"
echo "  ${OUTDIR}/weights/prs_weights.PT_*.tsv"
echo "  ${OUTDIR}/scores/*.profile"
echo "  ${OUTDIR}/prs_logistic_results.tsv"
echo "  ${OUTDIR}/best_prs_result.txt"
echo "  ${OUTDIR}/prs_interpretation.txt"
echo "  ${OUTDIR}/plots/prs_threshold_performance.png"