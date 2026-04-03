import pandas as pd
import numpy as np

# ---------------------------
# Input files
# ---------------------------
SUMSTATS_FILE = "discovery_gwas_sim_sumstats_noisy.tsv"
BIM_FILE = "../gwas_assignment/HapMap_3_r3_1.bim"

CLEANED_FILE = "discovery_gwas_sim_sumstats_cleaned.tsv"
HARMONIZED_FILE = "prs_weights_harmonized.tsv"
QCD_FILE = "prs_weights_harmonized_qcd.tsv"

# ---------------------------
# Helper
# ---------------------------
def is_ambiguous(a1, a2):
    pair = (str(a1).upper(), str(a2).upper())
    return pair in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# ---------------------------
# 1) Discovery summary-stat QC
# ---------------------------
df = pd.read_csv(SUMSTATS_FILE, sep="\t")

required_cols = ["SNP", "A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP"]
missing_cols = [c for c in required_cols if c not in df.columns]
if len(missing_cols) > 0:
    raise ValueError(f"Missing required columns: {missing_cols}")

print("=== STEP 1: DISCOVERY SUMMARY-STAT QC ===")
print(f"Initial rows: {len(df)}")

# standardize allele case
df["A1"] = df["A1"].astype(str).str.upper()
df["A2"] = df["A2"].astype(str).str.upper()

# Missingness check
before = len(df)
df = df.dropna(subset=["SNP", "A1", "A2", "BETA", "SE", "P", "N"])
after = len(df)
print(f"Removed {before - after} rows due to missing SNP/A1/A2/BETA/SE/P/N")

# Range checks
before = len(df)
df = df[(df["P"] >= 0) & (df["P"] <= 1)]
print(f"Removed {before - len(df)} rows with invalid P")

before = len(df)
df = df[df["SE"] > 0]
print(f"Removed {before - len(df)} rows with SE <= 0")

beta_threshold = 10
extreme = df[np.abs(df["BETA"]) > beta_threshold]
print(f"Extreme |BETA| rows flagged: {len(extreme)}")

# Duplicate SNP IDs: keep smallest P
dup_count = int(df["SNP"].duplicated().sum())
print(f"Duplicate SNP IDs before deduplication: {dup_count}")

df = df.sort_values("P").drop_duplicates(subset="SNP", keep="first")
print(f"Rows after deduplication: {len(df)}")

df.to_csv(CLEANED_FILE, sep="\t", index=False)
print(f"Wrote cleaned summary stats: {CLEANED_FILE}")

# ---------------------------
# 2) Harmonization with target BIM
# ---------------------------
print("\n=== STEP 2: HARMONIZATION WITH TARGET BIM ===")

bim = pd.read_csv(
    BIM_FILE,
    sep=r"\s+",
    header=None,
    names=["CHR_bim", "SNP", "CM", "BP_bim", "A1_bim", "A2_bim"]
)

bim["A1_bim"] = bim["A1_bim"].astype(str).str.upper()
bim["A2_bim"] = bim["A2_bim"].astype(str).str.upper()

# SNP overlap
merged = df.merge(bim, on="SNP", how="inner")
overlap_count = len(merged)
pct_retained = 100 * overlap_count / len(df) if len(df) > 0 else 0

print(f"Overlap SNP count with target BIM: {overlap_count}")
print(f"Percent retained after overlap: {pct_retained:.2f}%")

# Ambiguous SNP removal
ambig_mask = merged.apply(lambda r: is_ambiguous(r["A1"], r["A2"]), axis=1)
n_ambig = int(ambig_mask.sum())
merged = merged.loc[~ambig_mask].copy()
print(f"Removed ambiguous A/T and C/G SNPs: {n_ambig}")

# Allele matching and flipping
matched_rows = []
flipped_count = 0
dropped_mismatch_count = 0

for _, row in merged.iterrows():
    ss_a1 = row["A1"]
    ss_a2 = row["A2"]
    bim_a1 = row["A1_bim"]
    bim_a2 = row["A2_bim"]

    # same orientation
    if ss_a1 == bim_a1 and ss_a2 == bim_a2:
        row["flip_status"] = "not_flipped"
        matched_rows.append(row)

    # swapped orientation
    elif ss_a1 == bim_a2 and ss_a2 == bim_a1:
        row["BETA"] = -row["BETA"]
        row["A1"] = bim_a1
        row["A2"] = bim_a2
        row["flip_status"] = "flipped"
        flipped_count += 1
        matched_rows.append(row)

    # mismatch
    else:
        dropped_mismatch_count += 1

harmonized = pd.DataFrame(matched_rows)

print(f"Flipped SNPs: {flipped_count}")
print(f"Dropped overlapping SNPs due to allele mismatch: {dropped_mismatch_count}")
print(f"Final harmonized SNP count: {len(harmonized)}")

# Output PRS weights file
prs_weights = harmonized[["SNP", "A1", "BETA", "P"]].copy()
prs_weights.to_csv(HARMONIZED_FILE, sep="\t", index=False)
print(f"Wrote harmonized PRS weights: {HARMONIZED_FILE}")

# ---------------------------
# 3) PRS-ready output sanity checks
# ---------------------------
print("\n=== STEP 3: PRS-READY OUTPUT SANITY CHECKS ===")

# a) Final file integrity
n_unique = prs_weights["SNP"].nunique()
n_rows = len(prs_weights)
n_missing_beta = int(prs_weights["BETA"].isna().sum())

print(f"Final file rows: {n_rows}")
print(f"Unique SNP IDs: {n_unique}")
print(f"Missing BETA values: {n_missing_beta}")

if n_unique != n_rows:
    print("WARNING: final file does not have unique SNP IDs")
else:
    print("PASS: final file has unique SNP IDs")

if n_missing_beta != 0:
    print("WARNING: final file has missing BETA values")
else:
    print("PASS: final file has no missing BETA values")

# b) Random ~10 SNP manual spot-check table
spot_n = min(10, len(harmonized))
spotcheck = harmonized.sample(n=spot_n, random_state=42).copy()

spotcheck_out = spotcheck[
    ["SNP", "A1", "A2", "A1_bim", "A2_bim", "BETA", "P", "flip_status"]
].copy()

print("\nRandom spot-check SNPs:")
print(spotcheck_out.to_string(index=False))

# Save qcd file
# This includes the final harmonized PRS weights plus flip status for traceability.
qcd = harmonized[["SNP", "A1", "BETA", "P", "flip_status"]].copy()
qcd.to_csv(QCD_FILE, sep="\t", index=False)
print(f"\nWrote QCd harmonized output: {QCD_FILE}")

# Optional: save the manual spot-check table too
spotcheck_out.to_csv("prs_weights_harmonized_spotcheck.tsv", sep="\t", index=False)
print("Wrote manual spot-check table: prs_weights_harmonized_spotcheck.tsv")