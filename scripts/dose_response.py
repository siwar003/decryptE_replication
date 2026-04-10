"""
Dose-response pipeline (Python port of the R workflow).
4-parameter log-logistic model (LL.4 / 4PL) implemented directly
using scipy.optimize.curve_fit — no drc-equivalent package used.
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ═══════════════════════════════════════════════════════════════
# ── 1. USER SETTINGS  ← only edit this block ─────────────────
# ═══════════════════════════════════════════════════════════════
DRUG_NAME     = "Vincristine"          # ← change to any drug in the mapping sheet

INPUT_FILE    = "data/proteinGroups_fdr0.01_example.txt"
MAPPING_FILE  = "data/Mapping_Sheet.xlsx"
OUTPUT_DIR    = "results"

DOSES_NM      = [1, 10, 100, 1000, 10000]   # nM — fixed for all drugs
DMSO_PER_PLATE = 6                           # 6 DMSO controls per plate

# Quality / threshold settings
FC_UP_THRESHOLD   = 1.5
FC_DOWN_THRESHOLD = 0.7
MIN_R2            = 0.8
MIN_VALID_PTS     = 4
MIN_EFFECT        = 0.15

# ═══════════════════════════════════════════════════════════════
# ── 2. 4-PARAMETER LOG-LOGISTIC MODEL (LL.4 / 4PL) ──────────
# ═══════════════════════════════════════════════════════════════
#
#   f(x) = bottom + (top - bottom) / (1 + (x / EC50)^slope)
#
# Parameters:
#   slope  – Hill slope (steepness); negative = decreasing curve
#   bottom – lower asymptote (response at x → ∞ when slope > 0)
#   top    – upper asymptote (response at x → 0)
#   EC50   – dose giving half-maximal response

def ll4(x, slope, bottom, top, ec50):
    """Numerically stable 4PL."""
    x = np.asarray(x, dtype=float)
    
    # Prevent division issues
    ec50 = np.maximum(ec50, 1e-12)
    
    # Work in log space to avoid overflow
    log_ratio = np.log(x) - np.log(ec50)
    
    # Clip before exponentiation
    log_term = slope * log_ratio
    log_term = np.clip(log_term, -700, 700)  # avoids exp overflow
    
    term = np.exp(log_term)
    
    return bottom + (top - bottom) / (1.0 + term)


def fit_ll4(doses, responses, dose_min_guard, dose_max_guard):
    """
    Fit the LL.4 model to observed dose-response data.

    Returns (params_dict, fitted_values) on success, or (None, fail_reason).
    """
    # Initial guesses
    top0    = responses[np.argmin(doses)]   # response at lowest dose ≈ top
    bottom0 = responses[np.argmax(doses)]   # response at highest dose ≈ bottom
    ec50_0  = np.sqrt(doses.min() * doses.max())  # geometric mean
    slope0  = 1.0

    # Try multiple initial slope signs
    best_fit = None
    best_cost = np.inf

    for s0 in [1.0, -1.0, 0.5, -0.5, 2.0, -2.0]:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, _ = curve_fit(
                    ll4, doses, responses,
                    p0=[s0, bottom0, top0, ec50_0],
                    bounds=([-50, -10, -10, doses.min() * 0.001],
                            [ 50,  10,  10, doses.max() * 1000]),
                    maxfev=10000,
                    method="trf",
                )
            residuals = responses - ll4(doses, *popt)
            cost = np.sum(residuals ** 2)
            if cost < best_cost:
                best_cost = cost
                best_fit = popt
        except (RuntimeError, ValueError):
            continue

    if best_fit is None:
        return None, "curve_fit error"

    slope_val, bottom_val, top_val, ec50_val = best_fit

    # Validate EC50 range
    if not np.isfinite(ec50_val) or ec50_val <= 0:
        return None, f"EC50 not finite ({ec50_val:.2e} nM)"
    if ec50_val < dose_min_guard or ec50_val > dose_max_guard:
        return None, f"EC50 outside range ({ec50_val:.2e} nM)"

    # Check minimum effect
    if abs(top_val - bottom_val) < MIN_EFFECT:
        return None, f"|top-bottom| = {abs(top_val - bottom_val):.3f} < {MIN_EFFECT:.2f} (flat)"

    fitted_vals = ll4(doses, *best_fit)
    params = {
        "slope":  slope_val,
        "bottom": bottom_val,
        "top":    top_val,
        "EC50_nM": ec50_val,
    }
    return params, fitted_vals


# ═══════════════════════════════════════════════════════════════
# ── 3. RESOLVE DRUG METADATA FROM MAPPING SHEET ──────────────
# ═══════════════════════════════════════════════════════════════
print(f"Reading mapping sheet: {MAPPING_FILE}")
wb = openpyxl.load_workbook(MAPPING_FILE, read_only=True)
ws = wb.active

# Read all rows
mapping_rows = []
header = None
for i, row in enumerate(ws.iter_rows(values_only=True)):
    if i == 0:
        header = [str(c).strip().lower() if c else "" for c in row]
        continue
    mapping_rows.append(row)
wb.close()

# Build a DataFrame-like structure from the relevant columns
# Columns: sample (idx 0), drug (idx 2), ID (idx 3)
col_sample = 0
col_drug   = header.index("drug")
col_id     = header.index("id")

# ── Find drug ID from the 'drug' column ──────────────────────
drug_id = None
for row in mapping_rows:
    d = row[col_drug]
    if d and str(d).strip().lower() == DRUG_NAME.strip().lower():
        drug_id = str(int(row[col_id]) if isinstance(row[col_id], (int, float)) else row[col_id])
        break

if drug_id is None:
    available = sorted({str(r[col_drug]).strip() for r in mapping_rows if r[col_drug]})
    raise SystemExit(
        f"Drug '{DRUG_NAME}' not found in mapping sheet.\n"
        f"Available drugs:\n  " + "\n  ".join(available)
    )

# ── Find plate number by matching drug name in the 'sample' column ──
# Sample format: Plate<N>_<DrugName>_<dose> nM
# We search the full sample column (not the drug-matched row) because
# the sample and drug columns are independent sequences.
escaped_drug = re.escape(DRUG_NAME)
sample_pattern = re.compile(rf'(?i)^plate(\d+)_{escaped_drug}_\d')

matched_samples = []
plate_nums_found = []
for row in mapping_rows:
    s = str(row[col_sample]) if row[col_sample] else ""
    m = sample_pattern.search(s)
    if m:
        matched_samples.append(s)
        plate_nums_found.append(int(m.group(1)))

if not plate_nums_found:
    raise SystemExit(
        f"Could not find any sample rows matching drug '{DRUG_NAME}' "
        f"in the sample column.\n"
        f"  Searched for pattern: Plate<N>_{DRUG_NAME}_<dose>"
    )

print(f"  Matched sample rows (first 5):")
for s in matched_samples[:5]:
    print(f"    {s}")

unique_plates = sorted(set(plate_nums_found))
if len(unique_plates) > 1:
    print(f"  WARNING: Drug '{DRUG_NAME}' appears on multiple plates: "
          f"{unique_plates} -- using first: {unique_plates[0]}")
plate_num = unique_plates[0]

DRUG_ID = drug_id
PLATE_NUM = plate_num

print(f"  Drug: {DRUG_NAME}  |  ID: {DRUG_ID}")
print(f"  Plate: {PLATE_NUM}")

# ── Compute DMSO column names for this plate ────────────────
dmso_start = (PLATE_NUM - 1) * DMSO_PER_PLATE + 1
dmso_end   = PLATE_NUM * DMSO_PER_PLATE
DMSO_COLS  = [f"LFQ Intensity DMSO{i}" for i in range(dmso_start, dmso_end + 1)]
print(f"  DMSO columns: LFQ Intensity DMSO{dmso_start} - DMSO{dmso_end}")

# ── Build treatment column names ─────────────────────────────
TREATMENT_COLS = [f"LFQ Intensity {DRUG_ID} {d}" for d in DOSES_NM]
print("  Treatment columns:")
for tc in TREATMENT_COLS:
    print(f"    {tc}")

# ═══════════════════════════════════════════════════════════════
# ── 4. LOAD DATA ─────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
print(f"\nLoading: {INPUT_FILE}")
use_cols = ["Gene names"] + DMSO_COLS + TREATMENT_COLS
raw = pd.read_csv(INPUT_FILE, sep="\t", usecols=use_cols, low_memory=False)
print(f"  Rows: {len(raw)}  Cols: {raw.shape[1]}")

# ── 5. PREPARE GENE COLUMN ──────────────────────────────────
raw["gene"] = raw["Gene names"].astype(str).str.split(";").str[0].str.strip()

# Verify intensity columns exist
all_intensity_cols = DMSO_COLS + TREATMENT_COLS
missing_cols = [c for c in all_intensity_cols if c not in raw.columns]
if missing_cols:
    raise SystemExit(
        "Missing intensity columns in proteinGroups file:\n"
        + "\n".join(f"  - {c}" for c in missing_cols)
    )

# Convert 0 → NaN and coerce to numeric
for col in all_intensity_cols:
    raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw.loc[raw[col] == 0, col] = np.nan

dat = raw[raw["gene"].notna() & (raw["gene"] != "") & (raw["gene"] != "nan")].copy()
print(f"  Rows after removing empty gene names: {len(dat)}")

# Deduplicate: keep row with most valid intensities
if dat["gene"].duplicated().any():
    n_before = len(dat)
    dat["n_valid_total"] = dat[all_intensity_cols].notna().sum(axis=1)
    dat = dat.sort_values(["gene", "n_valid_total"], ascending=[True, False])
    dat = dat.drop_duplicates(subset="gene", keep="first")
    dat.drop(columns=["n_valid_total"], inplace=True)
    print(f"  Deduplicated gene names: {n_before} -> {len(dat)} rows")

dat = dat.reset_index(drop=True)

# ═══════════════════════════════════════════════════════════════
# ── 6. DMSO NORMALISATION ────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
dmso_mat  = dat[DMSO_COLS].values.astype(float)
treat_mat = dat[TREATMENT_COLS].values.astype(float)

dmso_mean    = np.nanmean(dmso_mat, axis=1)
n_dmso_valid = np.sum(~np.isnan(dmso_mat), axis=1)

keep = np.isfinite(dmso_mean) & (dmso_mean > 0)
dat          = dat.loc[keep].reset_index(drop=True)
dmso_mean    = dmso_mean[keep]
n_dmso_valid = n_dmso_valid[keep]
treat_mat    = treat_mat[keep]

# Relative intensity = treatment / mean(DMSO)
rel_mat = treat_mat / dmso_mean[:, np.newaxis]
n_treat_valid = np.sum(~np.isnan(rel_mat), axis=1)
log10_doses   = np.log10(DOSES_NM)
doses_arr     = np.array(DOSES_NM, dtype=float)

# ═══════════════════════════════════════════════════════════════
# ── 7. CURVE FITTING ─────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
dose_min_guard = min(DOSES_NM) * 0.1
dose_max_guard = max(DOSES_NM) * 10

n_prot = len(dat)
print(f"\nFitting LL.4 dose-response curves for {n_prot} genes...")

results = []

for i in range(n_prot):
    if (i + 1) % 500 == 0 or i == n_prot - 1:
        print(f"  {i + 1}/{n_prot}...", flush=True)

    rel_vals = rel_mat[i]
    valid_mask = np.isfinite(rel_vals)
    n_valid = int(valid_mask.sum())

    out = {
        "n_valid":            n_valid,
        "fit_success":        False,
        "regulated":          False,
        "fail_reason":        None,
        "EC50_nM":            np.nan,
        "pEC50":              np.nan,
        "slope":              np.nan,
        "top":                np.nan,
        "bottom":             np.nan,
        "AUC":                np.nan,
        "max_abs_effect":     np.nan,
        "fold_change_at_max": np.nan,
        "direction":          None,
        "r_squared":          np.nan,
        "mean_deviation":     np.nan,
        "linear_slope":       np.nan,
    }

    if n_valid < MIN_VALID_PTS:
        out["fail_reason"] = f"< {MIN_VALID_PTS} valid points"
        results.append(out)
        continue

    dose_fit = doses_arr[valid_mask]
    resp_fit = rel_vals[valid_mask]

    # Fit the 4PL model
    fit_result, fit_detail = fit_ll4(dose_fit, resp_fit, dose_min_guard, dose_max_guard)

    if fit_result is None:
        out["fail_reason"] = fit_detail  # fit_detail is the failure reason string
        results.append(out)
        continue

    slope_val  = fit_result["slope"]
    bottom_val = fit_result["bottom"]
    top_val    = fit_result["top"]
    ec50_val   = fit_result["EC50_nM"]
    fitted_vals = fit_detail  # on success, fit_detail is the fitted values array

    # pEC50: convert EC50 (in nM) to pEC50 = 9 - log10(EC50_nM)
    pec50_val = 9.0 - np.log10(ec50_val)

    # R²
    resids = resp_fit - fitted_vals
    ss_res = np.sum(resids ** 2)
    ss_tot = np.sum((resp_fit - np.mean(resp_fit)) ** 2)
    r_sq   = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    mean_dev = np.mean(np.abs(resids))

    out["fit_success"]    = True
    out["EC50_nM"]        = ec50_val
    out["pEC50"]          = pec50_val
    out["slope"]          = slope_val
    out["top"]            = top_val
    out["bottom"]         = bottom_val
    out["r_squared"]      = r_sq
    out["mean_deviation"] = mean_dev

    if np.isnan(r_sq) or r_sq < MIN_R2:
        out["fail_reason"] = f"R2 = {r_sq:.3f} < {MIN_R2:.2f}"
        results.append(out)
        continue

    # Fold change at max dose (predicted)
    fc_max = ll4(max(DOSES_NM), slope_val, bottom_val, top_val, ec50_val)
    is_up   = not np.isnan(fc_max) and fc_max > FC_UP_THRESHOLD
    is_down = not np.isnan(fc_max) and fc_max < FC_DOWN_THRESHOLD

    if not is_up and not is_down:
        out["fail_reason"] = (
            f"fold change at max dose = {fc_max:.3f} "
            f"(not > {FC_UP_THRESHOLD:.1f} or < {FC_DOWN_THRESHOLD:.1f})"
        )
        out["fold_change_at_max"] = fc_max
        results.append(out)
        continue

    # Regulated!
    out["regulated"]          = True
    out["fold_change_at_max"] = fc_max
    out["direction"]          = "up" if is_up else "down"
    out["max_abs_effect"]     = float(np.max(np.abs(resp_fit - 1.0)))

    # Linear slope over log10 doses
    log_d = np.log10(dose_fit)
    try:
        coeffs = np.polyfit(log_d, resp_fit, 1)
        out["linear_slope"] = coeffs[0]
    except Exception:
        pass

    # AUC via trapezoidal integration over fitted curve
    xs = np.linspace(log10_doses.min(), log10_doses.max(), 200)
    ys = ll4(10.0 ** xs, slope_val, bottom_val, top_val, ec50_val)
    out["AUC"] = float(trapezoid(ys, xs))

    results.append(out)

# ═══════════════════════════════════════════════════════════════
# ── 8. BUILD SUMMARY TABLE ───────────────────────────────────
# ═══════════════════════════════════════════════════════════════
results_df = pd.DataFrame(results)
summary_df = pd.DataFrame({
    "gene":          dat["gene"].values,
    "drug":          DRUG_NAME,
    "drug_id":       DRUG_ID,
    "plate":         PLATE_NUM,
    "n_dmso_valid":  n_dmso_valid.astype(int),
    "n_treat_valid": n_treat_valid.astype(int),
})
summary_df = pd.concat([summary_df, results_df], axis=1)

# Add relative intensity columns
for j, dose in enumerate(DOSES_NM):
    summary_df[f"rel_{dose}nM"] = rel_mat[:, j]

n_fit  = summary_df["fit_success"].sum()
n_reg  = summary_df["regulated"].sum()
n_up   = ((summary_df["regulated"]) & (summary_df["direction"] == "up")).sum()
n_down = ((summary_df["regulated"]) & (summary_df["direction"] == "down")).sum()

print(f"\n-- Fitting summary ------------------------------------------------")
print(f"  Drug                  : {DRUG_NAME}  (ID: {DRUG_ID}, Plate: {PLATE_NUM})")
print(f"  Total genes analysed  : {len(summary_df)}")
print(f"  Converged fits        : {n_fit}")
print(f"  Regulated             : {n_reg}  ({n_up} up / {n_down} down)")
print(f"  Fail/not-regulated    : {len(summary_df) - n_reg}")

fail_reasons = summary_df.loc[~summary_df["regulated"], "fail_reason"].dropna()
if len(fail_reasons) > 0:
    print("  Fail reasons:")
    for reason, count in fail_reasons.value_counts().items():
        print(f"    {reason:<55s} : {count}")

# ═══════════════════════════════════════════════════════════════
# ── 9. LONG FORMAT ───────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
long_records = []
genes = dat["gene"].values
for j, dose in enumerate(DOSES_NM):
    for i in range(len(genes)):
        long_records.append({
            "gene":               genes[i],
            "drug":               DRUG_NAME,
            "drug_id":            DRUG_ID,
            "plate":              PLATE_NUM,
            "dose_nM":            dose,
            "log10_dose":         np.log10(dose),
            "relative_intensity": rel_mat[i, j],
        })
long_df = pd.DataFrame(long_records)

# ═══════════════════════════════════════════════════════════════
# ── 10. WRITE OUTPUTS ────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
os.makedirs(OUTPUT_DIR, exist_ok=True)

safe_name = re.sub(r"[^A-Za-z0-9_]", "_", DRUG_NAME)
csv_path  = os.path.join(OUTPUT_DIR, f"summary_{safe_name}.csv")
xlsx_path = os.path.join(OUTPUT_DIR, f"summary_{safe_name}.xlsx")

print(f"\nWriting output files...")
summary_df.to_csv(csv_path, index=False)
print(f"  {csv_path}")

with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    long_df.to_excel(writer, sheet_name="Long_format", index=False)
print(f"  {xlsx_path}")

# ═══════════════════════════════════════════════════════════════
# ── 11. GENERATE PDF OF DOSE-RESPONSE CURVES ────────────────
# ═══════════════════════════════════════════════════════════════
# Include: regulated proteins + unregulated that only failed the FC cutoff
plot_mask = (
    summary_df["EC50_nM"].notna()
    & summary_df["pEC50"].notna()
    & summary_df["r_squared"].notna()
    & (summary_df["r_squared"] >= MIN_R2)
    & (
        summary_df["fail_reason"].isna()
        | summary_df["fail_reason"].str.startswith("fold change at max dose =", na=False)
    )
)
plot_df = summary_df[plot_mask].reset_index(drop=True)

# Build a lookup from gene → row index in dat
gene_to_idx = {g: i for i, g in enumerate(dat["gene"].values)}

print(f"\nGenerating PDF for {len(plot_df)} fitted proteins "
      f"(regulated + unregulated passing fit filters)...")

log10_plot_min = log10_doses.min() - 0.3
log10_plot_max = log10_doses.max() + 0.3
plot_doses_nM  = 10.0 ** np.linspace(log10_plot_min, log10_plot_max, 500)

pdf_path = os.path.join(OUTPUT_DIR, f"dose_response_curves_all_fitted_{safe_name}.pdf")

with PdfPages(pdf_path) as pdf:
    for idx in range(len(plot_df)):
        row = plot_df.iloc[idx]
        gene_id = row["gene"]
        data_idx = gene_to_idx.get(gene_id)
        if data_idx is None:
            continue

        rel_vals   = rel_mat[data_idx]
        valid_mask = np.isfinite(rel_vals)
        if valid_mask.sum() < MIN_VALID_PTS:
            continue

        dose_fit = doses_arr[valid_mask]
        resp_fit = rel_vals[valid_mask]

        # Re-fit for prediction (use stored parameters)
        s, b, t, e = row["slope"], row["bottom"], row["top"], row["EC50_nM"]
        pred_line = ll4(plot_doses_nM, s, b, t, e)

        if not np.all(np.isfinite(pred_line)):
            continue

        # Colours
        if row["direction"] == "up":
            curve_col, point_col = "#c0392b", "#e74c3c"
        elif row["direction"] == "down":
            curve_col, point_col = "#2980b9", "#3498db"
        else:
            curve_col, point_col = "#555555", "#777777"

        fig, ax = plt.subplots(figsize=(7, 5))

        # Fitted curve
        ax.plot(np.log10(plot_doses_nM), pred_line,
                color=curve_col, linewidth=1.0)
        # Observed points
        ax.scatter(np.log10(dose_fit), resp_fit,
                   color=point_col, s=30, zorder=5)
        # Reference lines
        ax.axhline(y=1, linestyle="--", color="grey", linewidth=0.4)
        ax.axvline(x=np.log10(row["EC50_nM"]), linestyle=":",
                   color=curve_col, linewidth=0.4, alpha=0.7)

        ax.set_xticks(log10_doses)
        ax.set_xticklabels(DOSES_NM)
        ax.set_xlim(log10_plot_min, log10_plot_max)
        ax.set_xlabel("Concentration (nM)")
        ax.set_ylabel("Relative intensity (treatment / mean DMSO)")

        # Y limits with padding
        all_y = np.concatenate([resp_fit, pred_line])
        y_min, y_max = np.nanmin(all_y), np.nanmax(all_y)
        y_pad = max((y_max - y_min) * 0.15, 0.05)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)

        # Title and annotations
        subtitle = (
            f"pEC50 = {row['pEC50']:.2f}  |  EC50 = {row['EC50_nM']:.1f} nM  |  "
            f"slope = {row['slope']:.2f}  |  top = {row['top']:.2f}  |  "
            f"bottom = {row['bottom']:.2f}  |  R2 = {row['r_squared']:.3f}"
        )
        direction_label = row["direction"] if row["direction"] else "unregulated"
        auc_label = f"{row['AUC']:.3f}" if not np.isnan(row["AUC"]) else "NA"
        fc_label  = f"{row['fold_change_at_max']:.2f}" if not np.isnan(row["fold_change_at_max"]) else "NA"

        caption = (
            f"Drug: {DRUG_NAME} (ID: {DRUG_ID}, Plate: {PLATE_NUM})  |  "
            f"direction: {direction_label}  |  FC at max dose: {fc_label}  |  "
            f"AUC: {auc_label}  |  n = {int(valid_mask.sum())} pts"
        )

        ax.set_title(gene_id, fontsize=12, fontweight="bold")
        ax.text(0.5, 1.01, subtitle, transform=ax.transAxes,
                fontsize=7.5, color="grey", ha="center", va="bottom")
        fig.text(0.5, 0.01, caption, fontsize=7, color="grey", ha="center")

        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

print(f"  {pdf_path}")

# ═══════════════════════════════════════════════════════════════
# ── 12. FINAL SUMMARY ────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════
print(f"\n========================================")
print(f"decryptE-like MaxQuant analysis complete")
print(f"========================================")
print(f"Drug    : {DRUG_NAME}  (ID: {DRUG_ID})")
print(f"Plate   : {PLATE_NUM}  (DMSO columns: DMSO{dmso_start}-DMSO{dmso_end})")
print(f"Output  : {OUTPUT_DIR}")
print(f"Files   :")
print(f"  summary_{safe_name}.csv")
print(f"  summary_{safe_name}.xlsx")
print(f"  dose_response_curves_all_fitted_{safe_name}.pdf")
print(f"")
print(f"Regulated proteins: {n_reg}  ({n_up} up / {n_down} down)")
print(f"FC thresholds: > {FC_UP_THRESHOLD} up  |  < {FC_DOWN_THRESHOLD} down")
print(f"========================================")
