#!/usr/bin/env bash
set -euo pipefail

VCF="ALL.2of4intersection.20100804.genotypes.vcf.gz"
PANEL="20100804.ALL.panel.txt"
HAPMAP_PREFIX="relatedness_qc/HapMap_3_r3_1_relatedness_filtered"

OUTDIR="kg_harmonization"
mkdir -p "${OUTDIR}"

have_plink_prefix () {
  local pfx="$1"
  [[ -f "${pfx}.bed" && -f "${pfx}.bim" && -f "${pfx}.fam" ]]
}

count_variants () {
  local pfx="$1"
  wc -l < "${pfx}.bim"
}

count_samples () {
  local pfx="$1"
  wc -l < "${pfx}.fam"
}

########################################
# 1) Preprocess input genotype data
########################################

if have_plink_prefix "${OUTDIR}/1KG_raw"; then
  echo "==> ${OUTDIR}/1KG_raw exists; skipping VCF -> PLINK conversion"
else
  echo "==> Converting 1KG VCF to PLINK format with filled-in variant IDs"
  plink \
    --vcf "${VCF}" \
    --set-missing-var-ids @:#:\$1:\$2 \
    --make-bed \
    --out "${OUTDIR}/1KG_raw"
fi

RAW_VARIANTS=$(count_variants "${OUTDIR}/1KG_raw")
RAW_SAMPLES=$(count_samples "${OUTDIR}/1KG_raw")

########################################
# 2) QC and SNP Harmonization of 1000 Genomes Data
########################################

if have_plink_prefix "${OUTDIR}/1KG_qc"; then
  echo "==> ${OUTDIR}/1KG_qc exists; skipping 1KG QC"
else
  echo "==> Running QC on 1KG dataset"
  plink \
    --bfile "${OUTDIR}/1KG_raw" \
    --geno 0.05 \
    --maf 0.01 \
    --make-bed \
    --out "${OUTDIR}/1KG_qc"
fi

QC_VARIANTS=$(count_variants "${OUTDIR}/1KG_qc")
QC_SAMPLES=$(count_samples "${OUTDIR}/1KG_qc")

if [[ -f "${OUTDIR}/shared.snps" ]]; then
  echo "==> ${OUTDIR}/shared.snps exists; skipping shared SNP discovery"
else
  echo "==> Finding shared SNPs between QCed 1KG and HapMap"
  cut -f2 "${OUTDIR}/1KG_qc.bim" | sort > "${OUTDIR}/1KG_qc.snps"
  cut -f2 "${HAPMAP_PREFIX}.bim" | sort > "${OUTDIR}/HapMap_qc.snps"
  comm -12 "${OUTDIR}/1KG_qc.snps" "${OUTDIR}/HapMap_qc.snps" > "${OUTDIR}/shared.snps"
fi

SHARED_SNPS=$(wc -l < "${OUTDIR}/shared.snps")

if have_plink_prefix "${OUTDIR}/1KG_qc_shared"; then
  echo "==> ${OUTDIR}/1KG_qc_shared exists; skipping 1KG shared extraction"
else
  echo "==> Restricting QCed 1KG to shared SNP set"
  plink \
    --bfile "${OUTDIR}/1KG_qc" \
    --extract "${OUTDIR}/shared.snps" \
    --make-bed \
    --out "${OUTDIR}/1KG_qc_shared"
fi

if have_plink_prefix "${OUTDIR}/HapMap_shared"; then
  echo "==> ${OUTDIR}/HapMap_shared exists; skipping HapMap shared extraction"
else
  echo "==> Restricting HapMap to shared SNP set"
  plink \
    --bfile "${HAPMAP_PREFIX}" \
    --extract "${OUTDIR}/shared.snps" \
    --make-bed \
    --out "${OUTDIR}/HapMap_shared"
fi

KG_SHARED_VARIANTS=$(count_variants "${OUTDIR}/1KG_qc_shared")
KG_SHARED_SAMPLES=$(count_samples "${OUTDIR}/1KG_qc_shared")
HAPMAP_SHARED_VARIANTS=$(count_variants "${OUTDIR}/HapMap_shared")
HAPMAP_SHARED_SAMPLES=$(count_samples "${OUTDIR}/HapMap_shared")

########################################
# 3) Harmonize genomic build
########################################

if [[ -f "${OUTDIR}/hapmap_chr_update.txt" && -f "${OUTDIR}/hapmap_pos_update.txt" ]]; then
  echo "==> HapMap update files exist; skipping coordinate file generation"
else
  echo "==> Building coordinate update files from HapMap shared BIM"
  awk 'BEGIN{OFS="\t"} {print $2, $1}' "${OUTDIR}/HapMap_shared.bim" > "${OUTDIR}/hapmap_chr_update.txt"
  awk 'BEGIN{OFS="\t"} {print $2, $4}' "${OUTDIR}/HapMap_shared.bim" > "${OUTDIR}/hapmap_pos_update.txt"
fi

if have_plink_prefix "${OUTDIR}/1KG_qc_shared_build_updated"; then
  echo "==> ${OUTDIR}/1KG_qc_shared_build_updated exists; skipping build update"
else
  echo "==> Updating 1KG shared SNP coordinates to HapMap build"
  plink \
    --bfile "${OUTDIR}/1KG_qc_shared" \
    --update-chr "${OUTDIR}/hapmap_chr_update.txt" 1 2 \
    --update-map "${OUTDIR}/hapmap_pos_update.txt" 1 2 \
    --make-bed \
    --out "${OUTDIR}/1KG_qc_shared_build_updated"
fi

UPDATED_VARIANTS=$(count_variants "${OUTDIR}/1KG_qc_shared_build_updated")
UPDATED_SAMPLES=$(count_samples "${OUTDIR}/1KG_qc_shared_build_updated")

########################################
# 4) Harmonize prior to merging
########################################

if [[ -f "${OUTDIR}/strand_mismatch_snps.txt" ]]; then
  echo "==> ${OUTDIR}/strand_mismatch_snps.txt exists; skipping first mismatch discovery"
else
  echo "==> First merge attempt to identify allele/strand mismatches"
  set +e
  plink \
    --bfile "${OUTDIR}/HapMap_shared" \
    --bmerge "${OUTDIR}/1KG_qc_shared_build_updated.bed" \
             "${OUTDIR}/1KG_qc_shared_build_updated.bim" \
             "${OUTDIR}/1KG_qc_shared_build_updated.fam" \
    --make-bed \
    --out "${OUTDIR}/merge_attempt_1"
  MERGE1_STATUS=$?
  set -e

  if [[ -f "${OUTDIR}/merge_attempt_1-merge.missnp" ]]; then
    cp "${OUTDIR}/merge_attempt_1-merge.missnp" "${OUTDIR}/strand_mismatch_snps.txt"
  else
    : > "${OUTDIR}/strand_mismatch_snps.txt"
  fi
fi

STRAND_MISMATCH_COUNT=$(wc -l < "${OUTDIR}/strand_mismatch_snps.txt")

if have_plink_prefix "${OUTDIR}/1KG_qc_shared_build_updated_flipped"; then
  echo "==> ${OUTDIR}/1KG_qc_shared_build_updated_flipped exists; skipping strand flip"
else
  echo "==> Flipping strand-mismatched SNPs in 1KG"
  if [[ "${STRAND_MISMATCH_COUNT}" -gt 0 ]]; then
    plink \
      --bfile "${OUTDIR}/1KG_qc_shared_build_updated" \
      --flip "${OUTDIR}/strand_mismatch_snps.txt" \
      --make-bed \
      --out "${OUTDIR}/1KG_qc_shared_build_updated_flipped"
  else
    plink \
      --bfile "${OUTDIR}/1KG_qc_shared_build_updated" \
      --make-bed \
      --out "${OUTDIR}/1KG_qc_shared_build_updated_flipped"
  fi
fi

if [[ -f "${OUTDIR}/unreconciled_snps.txt" ]]; then
  echo "==> ${OUTDIR}/unreconciled_snps.txt exists; skipping second mismatch discovery"
else
  echo "==> Second merge attempt after strand correction"
  set +e
  plink \
    --bfile "${OUTDIR}/HapMap_shared" \
    --bmerge "${OUTDIR}/1KG_qc_shared_build_updated_flipped.bed" \
             "${OUTDIR}/1KG_qc_shared_build_updated_flipped.bim" \
             "${OUTDIR}/1KG_qc_shared_build_updated_flipped.fam" \
    --make-bed \
    --out "${OUTDIR}/merge_attempt_2"
  MERGE2_STATUS=$?
  set -e

  if [[ -f "${OUTDIR}/merge_attempt_2-merge.missnp" ]]; then
    cp "${OUTDIR}/merge_attempt_2-merge.missnp" "${OUTDIR}/unreconciled_snps.txt"
  else
    : > "${OUTDIR}/unreconciled_snps.txt"
  fi
fi

UNRECONCILED_COUNT=$(wc -l < "${OUTDIR}/unreconciled_snps.txt")

if have_plink_prefix "${OUTDIR}/HapMap_shared_clean" && have_plink_prefix "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean"; then
  echo "==> Clean harmonized datasets exist; skipping problematic SNP removal"
else
  echo "==> Removing unreconciled SNPs from both datasets"
  if [[ "${UNRECONCILED_COUNT}" -gt 0 ]]; then
    plink \
      --bfile "${OUTDIR}/HapMap_shared" \
      --exclude "${OUTDIR}/unreconciled_snps.txt" \
      --make-bed \
      --out "${OUTDIR}/HapMap_shared_clean"

    plink \
      --bfile "${OUTDIR}/1KG_qc_shared_build_updated_flipped" \
      --exclude "${OUTDIR}/unreconciled_snps.txt" \
      --make-bed \
      --out "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean"
  else
    plink \
      --bfile "${OUTDIR}/HapMap_shared" \
      --make-bed \
      --out "${OUTDIR}/HapMap_shared_clean"

    plink \
      --bfile "${OUTDIR}/1KG_qc_shared_build_updated_flipped" \
      --make-bed \
      --out "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean"
  fi
fi

HAPMAP_CLEAN_VARIANTS=$(count_variants "${OUTDIR}/HapMap_shared_clean")
HAPMAP_CLEAN_SAMPLES=$(count_samples "${OUTDIR}/HapMap_shared_clean")
KG_CLEAN_VARIANTS=$(count_variants "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean")
KG_CLEAN_SAMPLES=$(count_samples "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean")

########################################
# 5) Merge HapMap and 1000 Genomes
########################################

if have_plink_prefix "${OUTDIR}/HapMap_1KG_merged"; then
  echo "==> ${OUTDIR}/HapMap_1KG_merged exists; skipping final merge"
  MERGE_SUCCESS="Yes"
else
  echo "==> Merging harmonized HapMap and 1KG datasets"
  set +e
  plink \
    --bfile "${OUTDIR}/HapMap_shared_clean" \
    --bmerge "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.bed" \
             "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.bim" \
             "${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.fam" \
    --make-bed \
    --out "${OUTDIR}/HapMap_1KG_merged"
  MERGE_FINAL_STATUS=$?
  set -e

  if [[ ${MERGE_FINAL_STATUS} -eq 0 ]] && have_plink_prefix "${OUTDIR}/HapMap_1KG_merged"; then
    MERGE_SUCCESS="Yes"
  else
    MERGE_SUCCESS="No"
  fi
fi

if have_plink_prefix "${OUTDIR}/HapMap_1KG_merged"; then
  MERGED_VARIANTS=$(count_variants "${OUTDIR}/HapMap_1KG_merged")
  MERGED_SAMPLES=$(count_samples "${OUTDIR}/HapMap_1KG_merged")
else
  MERGED_VARIANTS=0
  MERGED_SAMPLES=0
fi

########################################
# 6) Perform MDS and visualize population structure
########################################

MDS_VARIANTS=0
MDS_SAMPLES=0

if [[ "${MERGE_SUCCESS}" == "Yes" ]]; then
  if [[ -f "${OUTDIR}/HapMap_1KG_merged.genome" ]]; then
    echo "==> ${OUTDIR}/HapMap_1KG_merged.genome exists; skipping --genome"
  else
    echo "==> Computing pairwise genetic distances"
    plink \
      --bfile "${OUTDIR}/HapMap_1KG_merged" \
      --genome \
      --out "${OUTDIR}/HapMap_1KG_merged"
  fi

  if [[ -f "${OUTDIR}/HapMap_1KG_merged_mds.mds" ]]; then
    echo "==> ${OUTDIR}/HapMap_1KG_merged_mds.mds exists; skipping MDS"
  else
    echo "==> Running MDS"
    plink \
      --bfile "${OUTDIR}/HapMap_1KG_merged" \
      --read-genome "${OUTDIR}/HapMap_1KG_merged.genome" \
      --cluster \
      --mds-plot 4 \
      --out "${OUTDIR}/HapMap_1KG_merged_mds"
  fi

  if [[ -f "${OUTDIR}/HapMap_1KG_merged_mds.mds" ]]; then
    MDS_SAMPLES=$(tail -n +2 "${OUTDIR}/HapMap_1KG_merged_mds.mds" | wc -l)
    MDS_VARIANTS=${MERGED_VARIANTS}
  fi

  echo "==> Preparing population labels"
  awk 'BEGIN{OFS="\t"} {print $1, $1, $2}' "${PANEL}" > "${OUTDIR}/1KG_population_labels.txt"
  awk 'BEGIN{OFS="\t"} {print $1, $2, "STUDY"}' "${HAPMAP_PREFIX}.fam" > "${OUTDIR}/study_labels.txt"
  cat "${OUTDIR}/1KG_population_labels.txt" "${OUTDIR}/study_labels.txt" > "${OUTDIR}/all_labels.txt"

  echo "==> Joining labels to MDS output"
  python3 - <<'PY'
import pandas as pd

outdir = "kg_harmonization"

mds = pd.read_csv(
    f"{outdir}/HapMap_1KG_merged_mds.mds",
    sep=r"\s+",
    engine="python",
)

labels = pd.read_csv(
    f"{outdir}/all_labels.txt",
    sep=r"\s+",
    engine="python",
    header=None,
    names=["FID", "IID", "POP"],
)

plot_df = mds.merge(labels, on=["FID", "IID"], how="left")
plot_df["POP"] = plot_df["POP"].fillna("UNKNOWN")
plot_df.to_csv(
    f"{outdir}/HapMap_1KG_merged_mds_annotated.tsv",
    sep="\t",
    index=False,
)
PY

  echo "==> Plotting MDS results"
  python3 - <<'PY'
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

outdir = "kg_harmonization"
df = pd.read_csv(f"{outdir}/HapMap_1KG_merged_mds_annotated.tsv", sep="\t")

study = df[df["POP"] == "STUDY"].copy()
ref = df[df["POP"] != "STUDY"].copy()

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial"]

TEAL = "#5BA4A4"

pop_colors = {
    "CEU": "#4E79A7",
    "TSI": "#A0CBE8",
    "GBR": "#F28E2B",
    "FIN": "#FFBE7D",
    "CHB": "#59A14F",
    "CHS": "#8CD17D",
    "JPT": "#B6992D",
    "YRI": "#E15759",
    "LWK": "#FF9D9A",
    "ASW": "#B07AA1",
    "MXL": "#9C755F",
    "PUR": "#ebc14d",
    "UNKNOWN": "0.7",
}

pop_colors = {
    "CEU": "#4E79A7",
    "TSI": "#A0CBE8",
    "GBR": "#2691bf",
    "FIN": "#88a6bf",
    "CHB": "#E15759",
    "CHS": "#FF9D9A",
    "JPT": "#B07AA1",
    "YRI": "#069c3b",
    "LWK": "#168c6d",
    "ASW": "#88bf50",
    "MXL": "#FFBE7D",
    "PUR": "#ebc14d",
    "UNKNOWN": "0.7",
}

def clean_fmt(x, pos):
    if abs(x) < 1e-12:
        return "0"
    return f"{x:.1f}"

fig, ax = plt.subplots(figsize=(6, 6))

pop_order = [
    p for p in [
        "CEU", "TSI", "GBR", "FIN", "CHB", "CHS",
        "JPT", "YRI", "LWK", "ASW", "MXL", "PUR", "UNKNOWN"
    ]
    if p in set(ref["POP"])
]

for pop in pop_order:
    sub = ref[ref["POP"] == pop]
    if len(sub) == 0:
        continue
    ax.scatter(
        sub["C1"],
        sub["C2"],
        s=16,
        color=pop_colors.get(pop, "0.7"),
        edgecolor="white",
        linewidth=0.4,
        alpha=1,
        label=pop,
    )

ax.scatter(
    study["C1"],
    study["C2"],
    s=16,
    color=TEAL,
    edgecolor="white",
    linewidth=0.4,
    alpha=1,
    label="Study",
    zorder=10,
)

ax.set_xlabel("MDS1", fontsize=20, labelpad=14)
ax.set_ylabel("MDS2", fontsize=20, labelpad=14)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.0)
ax.spines["bottom"].set_linewidth(1.0)

ax.set_aspect("equal", adjustable="box")
ax.tick_params(axis="both", length=4, width=1, labelsize=16, pad=5)

ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_major_formatter(FuncFormatter(clean_fmt))
ax.yaxis.set_major_formatter(FuncFormatter(clean_fmt))

leg = ax.legend(
    frameon=False,
    fontsize=11,
    loc="lower right",
    handletextpad=0.4,
    borderpad=1,
    markerscale=1.5,
    ncol=2,
)

for h in leg.legend_handles:
    try:
        h.set_edgecolor("white")
        h.set_linewidth(0.8)
    except Exception:
        pass

plt.tight_layout()
plt.savefig(f"{outdir}/HapMap_1KG_MDS.png", dpi=300, bbox_inches="tight")
plt.close()
PY

  ########################################
  # 7) Generate MDS-based covariates
  ########################################

  echo "==> Generating MDS covariate file"
  python3 - <<'PY'
import pandas as pd

outdir = "kg_harmonization"

mds = pd.read_csv(
    f"{outdir}/HapMap_1KG_merged_mds.mds",
    sep=r"\s+",
    engine="python",
)

cov = mds[["FID", "IID", "C1", "C2", "C3", "C4"]].copy()
cov.to_csv(f"{outdir}/mds_covariates.txt", sep="\t", index=False)
PY
fi

########################################
# Final report
########################################

echo
echo "========================================"
echo "Assignment report"
echo "========================================"
echo
echo "1) Preprocess input genotype data"
echo "a) VCF -> PLINK command:"
echo "   plink --vcf ${VCF} --set-missing-var-ids @:#:\$1:\$2 --make-bed --out ${OUTDIR}/1KG_raw"
echo
echo "b) Missing variant ID assignment command:"
echo "   --set-missing-var-ids @:#:\$1:\$2"
echo
echo "c) Resulting PLINK dataset:"
echo "   Variants: ${RAW_VARIANTS}"
echo "   Samples:  ${RAW_SAMPLES}"
echo
echo "2) QC and SNP Harmonization of 1000 Genomes Data"
echo "a) QC command:"
echo "   plink --bfile ${OUTDIR}/1KG_raw --geno 0.05 --maf 0.01 --make-bed --out ${OUTDIR}/1KG_qc"
echo
echo "b) Remaining after QC:"
echo "   Variants: ${QC_VARIANTS}"
echo "   Samples:  ${QC_SAMPLES}"
echo
echo "c) Shared SNP set obtained by intersecting SNP IDs from:"
echo "   ${OUTDIR}/1KG_qc.bim"
echo "   ${HAPMAP_PREFIX}.bim"
echo "   Shared SNPs: ${SHARED_SNPS}"
echo
echo "d-i) Filtered 1KG dataset after restricting to shared SNPs:"
echo "   Variants: ${KG_SHARED_VARIANTS}"
echo "   Samples:  ${KG_SHARED_SAMPLES}"
echo
echo "d-ii) HapMap dataset after extracting the same SNP set:"
echo "   Variants: ${HAPMAP_SHARED_VARIANTS}"
echo "   Samples:  ${HAPMAP_SHARED_SAMPLES}"
echo
echo "3) Harmonize Genomic Build Between HapMap and 1000 Genomes"
echo "a) Coordinate information source:"
echo "   Chromosome and base-pair positions were taken from ${OUTDIR}/HapMap_shared.bim"
echo "   for the shared SNPs, then used to update the 1KG shared dataset."
echo
echo "b) Commands used to update SNP positions:"
echo "   plink --bfile ${OUTDIR}/1KG_qc_shared \\"
echo "         --update-chr ${OUTDIR}/hapmap_chr_update.txt 1 2 \\"
echo "         --update-map ${OUTDIR}/hapmap_pos_update.txt 1 2 \\"
echo "         --make-bed --out ${OUTDIR}/1KG_qc_shared_build_updated"
echo
echo "c) Updated 1KG dataset:"
echo "   Variants: ${UPDATED_VARIANTS}"
echo "   Samples:  ${UPDATED_SAMPLES}"
echo
echo "4) Harmonize HapMap and 1000 Genomes Data Prior to Merging"
echo "a) Reference allele coding was matched by attempting a PLINK merge and using"
echo "   the resulting .missnp file to identify SNPs with allele/strand mismatches."
echo
echo "b) Commands used to identify and correct strand mismatches:"
echo "   plink --bfile ${OUTDIR}/HapMap_shared --bmerge ${OUTDIR}/1KG_qc_shared_build_updated.bed ${OUTDIR}/1KG_qc_shared_build_updated.bim ${OUTDIR}/1KG_qc_shared_build_updated.fam --make-bed --out ${OUTDIR}/merge_attempt_1"
echo "   plink --bfile ${OUTDIR}/1KG_qc_shared_build_updated --flip ${OUTDIR}/strand_mismatch_snps.txt --make-bed --out ${OUTDIR}/1KG_qc_shared_build_updated_flipped"
echo "   Strand-mismatched SNPs identified: ${STRAND_MISMATCH_COUNT}"
echo
echo "c) SNPs that could not be reconciled after strand correction were identified"
echo "   from ${OUTDIR}/merge_attempt_2-merge.missnp after the second merge attempt."
echo "   Unreconciled SNPs: ${UNRECONCILED_COUNT}"
echo
echo "d) These problematic SNPs were removed with:"
echo "   plink --bfile ${OUTDIR}/HapMap_shared --exclude ${OUTDIR}/unreconciled_snps.txt --make-bed --out ${OUTDIR}/HapMap_shared_clean"
echo "   plink --bfile ${OUTDIR}/1KG_qc_shared_build_updated_flipped --exclude ${OUTDIR}/unreconciled_snps.txt --make-bed --out ${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean"
echo
echo "e) Counts after harmonization and filtering:"
echo "   HapMap clean variants: ${HAPMAP_CLEAN_VARIANTS}"
echo "   HapMap clean samples:  ${HAPMAP_CLEAN_SAMPLES}"
echo "   1KG clean variants:    ${KG_CLEAN_VARIANTS}"
echo "   1KG clean samples:     ${KG_CLEAN_SAMPLES}"
echo
echo "5) Merge HapMap and 1000 Genomes Data"
echo "a) Merge command:"
echo "   plink --bfile ${OUTDIR}/HapMap_shared_clean --bmerge ${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.bed ${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.bim ${OUTDIR}/1KG_qc_shared_build_updated_flipped_clean.fam --make-bed --out ${OUTDIR}/HapMap_1KG_merged"
echo
echo "b) Merged dataset:"
echo "   Variants: ${MERGED_VARIANTS}"
echo "   Samples:  ${MERGED_SAMPLES}"
echo
echo "c) Was the merge successful?"
echo "   ${MERGE_SUCCESS}"
echo

if [[ "${MERGE_SUCCESS}" == "Yes" ]]; then
  echo "6) Perform MDS Analysis and Visualize Population Structure"
  echo "a) Commands used:"
  echo "   plink --bfile ${OUTDIR}/HapMap_1KG_merged --genome --out ${OUTDIR}/HapMap_1KG_merged"
  echo "   plink --bfile ${OUTDIR}/HapMap_1KG_merged --read-genome ${OUTDIR}/HapMap_1KG_merged.genome --cluster --mds-plot 4 --out ${OUTDIR}/HapMap_1KG_merged_mds"
  echo
  echo "b) Included in MDS:"
  echo "   Variants: ${MDS_VARIANTS}"
  echo "   Samples:  ${MDS_SAMPLES}"
  echo
  echo "c) MDS results were annotated by joining the .mds output to 1KG population labels"
  echo "   from ${PANEL}; HapMap samples were labeled as STUDY."
  echo
  echo "d) Final MDS plot:"
  echo "   ${OUTDIR}/HapMap_1KG_MDS.png"
  echo
  echo "e) Interpretation:"
  echo "   Inspect where the STUDY samples cluster relative to the 1KG reference populations."
  echo
  echo "7) Generate MDS-Based Covariates for Downstream Analysis"
  echo "a) Covariate file was created by extracting FID, IID, and the first four MDS components"
  echo "   from ${OUTDIR}/HapMap_1KG_merged_mds.mds."
  echo
  echo "b) Included components:"
  echo "   C1, C2, C3, C4"
  echo
  echo "c) These covariates are included to account for population structure in downstream analyses."
  echo
  echo "Output covariate file:"
  echo "   ${OUTDIR}/mds_covariates.txt"
else
  echo "MDS and covariate generation were skipped because the merge was not successful."
fi

echo
echo "Done."