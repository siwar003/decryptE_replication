# ============================================================
# decrypte_like_dose_response_maxquant.R
# Multi-drug version: user specifies a drug name, the script
# resolves its ID, plate, treatment columns, and DMSO columns
# automatically from the Mapping_Sheet.xlsx
# ============================================================

# ── 0. PACKAGE SETUP ──────────────────────────────────────────
required_pkgs <- c("drc", "writexl", "ggplot2",
                   "dplyr", "tidyr", "pracma", "data.table", "readxl")
missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install missing packages:\n",
       'install.packages(c(',
       paste0('"', missing_pkgs, '"', collapse = ", "), '))')
}
suppressPackageStartupMessages({
  library(drc);     library(dplyr);  library(tidyr)
  library(ggplot2); library(pracma); library(writexl)
  library(data.table); library(readxl)
})

# ══════════════════════════════════════════════════════════════
# ── 1. USER SETTINGS  ─────────────────────────────────────────
# ══════════════════════════════════════════════════════════════
DRUG_NAME     <- "Vincristine"    # ← change to any drug in the mapping sheet

INPUT_FILE    <- "data/proteinGroups_fdr0.01_example.txt"
MAPPING_FILE  <- "data/Mapping_Sheet.xlsx"
OUTPUT_DIR    <- "results"

DOSES_NM      <- c(1, 10, 100, 1000, 10000)   # nM — fixed for all drugs
DMSO_PER_PLATE <- 6                            # 6 DMSO controls per plate

# Quality / threshold settings
FC_UP_THRESHOLD   <- 1.5
FC_DOWN_THRESHOLD <- 0.7
MIN_R2            <- 0.8
MIN_VALID_PTS     <- 4
MIN_EFFECT        <- 0.15
# ══════════════════════════════════════════════════════════════

# ── 2. RESOLVE DRUG METADATA FROM MAPPING SHEET ───────────────
message("Reading mapping sheet: ", MAPPING_FILE)
mapping <- as.data.frame(read_excel(MAPPING_FILE))

# Normalise column names to lowercase for robust matching
colnames(mapping) <- tolower(trimws(colnames(mapping)))

required_map_cols <- c("sample", "drug", "id")
missing_map <- setdiff(required_map_cols, colnames(mapping))
if (length(missing_map) > 0) {
  stop("Mapping sheet is missing column(s): ",
       paste(missing_map, collapse = ", "),
       "\n  Found columns: ", paste(colnames(mapping), collapse = ", "))
}

# Filter to rows matching the requested drug (case-insensitive)
drug_rows <- mapping[tolower(trimws(mapping$drug)) == tolower(trimws(DRUG_NAME)), ]

if (nrow(drug_rows) == 0) {
  available <- sort(unique(trimws(mapping$drug)))
  stop("Drug '", DRUG_NAME, "' not found in mapping sheet.\n",
       "Available drugs:\n  ", paste(available, collapse = "\n  "))
}

# ── Extract drug ID ───────────────────────────────────────────
drug_ids <- unique(drug_rows$id)
if (length(drug_ids) > 1) {
  warning("Multiple IDs found for '", DRUG_NAME, "': ",
          paste(drug_ids, collapse = ", "), " — using first: ", drug_ids[1])
}
DRUG_ID <- as.character(drug_ids[1])
message(sprintf("  Drug: %s  |  ID: %s", DRUG_NAME, DRUG_ID))

# ── Extract plate number from the 'sample' column ─────────────
# Expected format: Plate#_drugName_dose nM  (e.g. "Plate1_Vincristine_100 nM")
# Grab the digit(s) right after "Plate" (case-insensitive)
plate_nums <- as.integer(
  regmatches(drug_rows$sample,
             regexpr("(?i)(?<=plate)\\d+", drug_rows$sample, perl = TRUE))
)
plate_nums <- plate_nums[!is.na(plate_nums)]

if (length(plate_nums) == 0) {
  stop("Could not parse plate number from 'sample' column for drug '", DRUG_NAME,
       "'.\n  Example sample values found:\n  ",
       paste(head(drug_rows$sample, 3), collapse = "\n  "),
       "\n  Expected format: Plate<number>_<drug>_<dose> nM")
}

unique_plates <- unique(plate_nums)
if (length(unique_plates) > 1) {
  warning("Drug '", DRUG_NAME, "' appears on multiple plates: ",
          paste(unique_plates, collapse = ", "),
          " — using first: ", unique_plates[1])
}
PLATE_NUM <- unique_plates[1]
message(sprintf("  Plate: %d", PLATE_NUM))

# ── Compute DMSO column indices for this plate ────────────────
# Plate 1 → DMSO 1–6, Plate 2 → DMSO 7–12, Plate N → DMSO (N-1)*6+1 : N*6
dmso_start <- (PLATE_NUM - 1) * DMSO_PER_PLATE + 1
dmso_end   <- PLATE_NUM * DMSO_PER_PLATE
DMSO_COLS  <- paste0("LFQ Intensity DMSO", dmso_start:dmso_end)
message(sprintf("  DMSO columns: LFQ Intensity DMSO%d – DMSO%d", dmso_start, dmso_end))

# ── Build treatment column names from drug ID + doses ─────────
# Format: "LFQ Intensity <ID> <dose_nM>"
TREATMENT_COLS <- paste0("LFQ Intensity ", DRUG_ID, " ", DOSES_NM)
message("  Treatment columns:")
for (tc in TREATMENT_COLS) message("    ", tc)

# ── 3. OUTPUT DIRECTORY ───────────────────────────────────────
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 4. LOAD DATA (memory-efficient) ──────────────────────────
message("\nLoading: ", INPUT_FILE)
raw <- data.table::fread(
  INPUT_FILE,
  select = c("Gene names", DMSO_COLS, TREATMENT_COLS)
)
message("  Rows: ", nrow(raw), "  Cols: ", ncol(raw))

# ── 5. PREPARE GENE COLUMN ────────────────────────────────────
raw$gene <- sapply(strsplit(raw[["Gene names"]], ";"), `[`, 1)
raw$gene <- trimws(raw$gene)

# Verify intensity columns exist
all_intensity_cols <- c(DMSO_COLS, TREATMENT_COLS)
missing_cols <- setdiff(all_intensity_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("Missing intensity columns in proteinGroups file:\n",
       paste("  •", missing_cols, collapse = "\n"),
       "\n\nCheck that the drug ID and column naming in the proteinGroups file",
       "\nmatch what is expected (LFQ Intensity <ID> <dose>).")
}

# Convert 0 → NA and coerce to numeric
for (col in all_intensity_cols) {
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
  raw[[col]][raw[[col]] == 0] <- NA_real_
}

dat <- raw[!is.na(raw$gene) & nchar(raw$gene) > 0, ]
message("  Rows after removing empty gene names: ", nrow(dat))

# Deduplicate: keep row with most valid intensities
if (anyDuplicated(dat$gene)) {
  n_before <- nrow(dat)
  dat$n_valid_total <- rowSums(!is.na(dat[, all_intensity_cols, with = FALSE]))
  dat <- dat[order(dat$gene, -dat$n_valid_total)]
  dat <- dat[!duplicated(dat$gene)]
  dat$n_valid_total <- NULL
  message("  Deduplicated gene names: ", n_before, " → ", nrow(dat), " rows")
}

# ── 6. DMSO NORMALISATION ─────────────────────────────────────
dmso_mat  <- as.matrix(dat[, ..DMSO_COLS])
treat_mat <- as.matrix(dat[, ..TREATMENT_COLS])

dmso_mean    <- rowMeans(dmso_mat, na.rm = TRUE)
n_dmso_valid <- rowSums(!is.na(dmso_mat))

keep      <- !is.na(dmso_mean) & is.finite(dmso_mean) & dmso_mean > 0
dat          <- dat[keep, ]
dmso_mean    <- dmso_mean[keep]
n_dmso_valid <- n_dmso_valid[keep]
treat_mat    <- treat_mat[keep, , drop = FALSE]

rel_mat <- treat_mat / dmso_mean
colnames(rel_mat) <- paste0("rel_", DOSES_NM, "nM")
n_treat_valid <- rowSums(!is.na(rel_mat))
log10_doses   <- log10(DOSES_NM)

# ── 7. CURVE FITTING ──────────────────────────────────────────
dose_min_guard <- min(DOSES_NM) * 0.1
dose_max_guard <- max(DOSES_NM) * 10

message("\nFitting LL.4 dose-response curves for ", nrow(dat), " genes...")

# Helper: linear slope over log10 doses
compute_linear_slope <- function(log_doses, responses) {
  tryCatch({
    fit <- lm(responses ~ log_doses)
    unname(coef(fit)[2])
  }, error = function(e) NA_real_)
}

# Helper: AUC via trapezoidal integration over the fitted curve
compute_auc_trapz <- function(fit_obj, log_dose_min, log_dose_max) {
  tryCatch({
    xs <- seq(log_dose_min, log_dose_max, length.out = 200)
    ys <- as.numeric(predict(fit_obj,
                             newdata = data.frame(dose_fit = 10^xs)))
    pracma::trapz(xs, ys)
  }, error = function(e) NA_real_)
}

n_prot  <- nrow(dat)
results <- vector("list", n_prot)
pb      <- txtProgressBar(min = 0, max = n_prot, style = 3)

for (i in seq_len(n_prot)) {
  rel_vals   <- rel_mat[i, ]
  valid_mask <- !is.na(rel_vals) & is.finite(rel_vals)
  n_valid    <- sum(valid_mask)
  out <- list(
    n_valid              = n_valid,
    fit_success          = FALSE,
    regulated            = FALSE,
    fail_reason          = NA_character_,
    EC50_nM              = NA_real_,
    pEC50                = NA_real_,
    slope                = NA_real_,
    top                  = NA_real_,
    bottom               = NA_real_,
    AUC                  = NA_real_,
    max_abs_effect       = NA_real_,
    fold_change_at_max   = NA_real_,
    direction            = NA_character_,
    r_squared            = NA_real_,
    mean_deviation       = NA_real_,
    linear_slope         = NA_real_
  )
  if (n_valid < MIN_VALID_PTS) {
    out$fail_reason <- paste0("< ", MIN_VALID_PTS, " valid points")
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  dose_fit <- DOSES_NM[valid_mask]
  resp_fit <- rel_vals[valid_mask]
  fit <- tryCatch(
    suppressWarnings(
      drm(resp_fit ~ dose_fit,
          fct = LL.4(names = c("slope", "bottom", "top", "EC50")),
          control = drmc(maxIt = 1000, errorm = FALSE, noMessage = TRUE))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    out$fail_reason <- "drc error"
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(cf) || length(cf) < 4 || anyNA(cf)) {
    out$fail_reason <- "NA coefficients"
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  slope_val  <- unname(cf[1]); bottom_val <- unname(cf[2])
  top_val    <- unname(cf[3]); ec50_val   <- unname(cf[4])
  if (!is.finite(ec50_val) || ec50_val <= 0 ||
      ec50_val < dose_min_guard || ec50_val > dose_max_guard) {
    out$fail_reason <- sprintf("EC50 outside range (%.2e nM)", ec50_val)
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  if (abs(top_val - bottom_val) < MIN_EFFECT) {
    out$fail_reason <- sprintf("|top-bottom| = %.3f < %.2f (flat)",
                               abs(top_val - bottom_val), MIN_EFFECT)
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  pec50_val   <- 9 - log10(ec50_val)
  fitted_vals <- tryCatch(fitted(fit),
                          error = function(e) rep(NA_real_, length(resp_fit)))
  resids  <- resp_fit - fitted_vals
  ss_res  <- sum(resids^2, na.rm = TRUE)
  ss_tot  <- sum((resp_fit - mean(resp_fit, na.rm = TRUE))^2, na.rm = TRUE)
  r_sq    <- if (is.finite(ss_tot) && ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  mean_dev <- mean(abs(resids), na.rm = TRUE)
  out$fit_success    <- TRUE
  out$EC50_nM        <- ec50_val; out$pEC50 <- pec50_val
  out$slope          <- slope_val; out$top  <- top_val; out$bottom <- bottom_val
  out$r_squared      <- r_sq;     out$mean_deviation <- mean_dev
  if (is.na(r_sq) || r_sq < MIN_R2) {
    out$fail_reason <- sprintf("R² = %.3f < %.2f", r_sq, MIN_R2)
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  fc_max <- tryCatch(
    as.numeric(predict(fit, newdata = data.frame(dose_fit = max(DOSES_NM)))),
    error = function(e) NA_real_
  )
  is_up   <- !is.na(fc_max) && fc_max > FC_UP_THRESHOLD
  is_down <- !is.na(fc_max) && fc_max < FC_DOWN_THRESHOLD
  if (!is_up && !is_down) {
    out$fail_reason <- sprintf(
      "fold change at max dose = %.3f (not > %.1f or < %.1f)",
      ifelse(is.na(fc_max), NA, fc_max), FC_UP_THRESHOLD, FC_DOWN_THRESHOLD)
    out$fold_change_at_max <- fc_max
    results[[i]] <- out; setTxtProgressBar(pb, i); next
  }
  out$regulated          <- TRUE
  out$fail_reason        <- NA_character_
  out$fold_change_at_max <- fc_max
  out$direction          <- if (is_up) "up" else "down"
  out$max_abs_effect     <- max(abs(resp_fit - 1), na.rm = TRUE)
  out$linear_slope       <- compute_linear_slope(log10_doses[valid_mask], resp_fit)
  out$AUC                <- compute_auc_trapz(fit, min(log10_doses), max(log10_doses))
  results[[i]] <- out
  setTxtProgressBar(pb, i)
}
close(pb)

# ── 8. BUILD SUMMARY TABLE ────────────────────────────────────
results_df <- dplyr::bind_rows(results)
summary_df <- data.frame(
  gene          = dat$gene,
  drug          = DRUG_NAME,
  drug_id       = DRUG_ID,
  plate         = PLATE_NUM,
  n_dmso_valid  = n_dmso_valid,
  n_treat_valid = n_treat_valid,
  stringsAsFactors = FALSE
)
summary_df <- cbind(summary_df, results_df)

n_fit  <- sum(summary_df$fit_success,  na.rm = TRUE)
n_reg  <- sum(summary_df$regulated,    na.rm = TRUE)
n_up   <- sum(summary_df$regulated & summary_df$direction == "up",   na.rm = TRUE)
n_down <- sum(summary_df$regulated & summary_df$direction == "down", na.rm = TRUE)

message("\n── Fitting summary ────────────────────────────────────")
message("  Drug                  : ", DRUG_NAME, "  (ID: ", DRUG_ID, ", Plate: ", PLATE_NUM, ")")
message("  Total genes analysed  : ", nrow(summary_df))
message("  Converged fits        : ", n_fit)
message("  Regulated             : ", n_reg,
        "  (", n_up, " up / ", n_down, " down)")
message("  Fail/not-regulated    : ", nrow(summary_df) - n_reg)
fail_tab <- sort(table(summary_df$fail_reason[!summary_df$regulated]), decreasing = TRUE)
if (length(fail_tab) > 0) {
  message("  Fail reasons:")
  for (nm in names(fail_tab))
    message(sprintf("    %-55s : %d", nm, fail_tab[nm]))
}

# ── 9. LONG FORMAT ────────────────────────────────────────────
long_df <- data.frame(
  gene               = rep(dat$gene, each = length(DOSES_NM)),
  drug               = DRUG_NAME,
  drug_id            = DRUG_ID,
  plate              = PLATE_NUM,
  dose_nM            = rep(DOSES_NM, times = nrow(dat)),
  log10_dose         = rep(log10_doses, times = nrow(dat)),
  relative_intensity = as.vector(t(rel_mat)),
  stringsAsFactors   = FALSE
)

# ── 10. WRITE OUTPUTS ─────────────────────────────────────────
# File names include the drug name so repeated runs don't overwrite each other
safe_name  <- gsub("[^A-Za-z0-9_]", "_", DRUG_NAME)   # filesystem-safe
csv_path   <- file.path(OUTPUT_DIR, paste0("summary_", safe_name, ".csv"))
xlsx_path  <- file.path(OUTPUT_DIR, paste0("summary_", safe_name, ".xlsx"))

message("\nWriting output files...")
write.csv(summary_df, csv_path, row.names = FALSE)
message("  ", csv_path)
writexl::write_xlsx(list("Summary" = summary_df, "Long_format" = long_df),
                    path = xlsx_path)
message("  ", xlsx_path)

# ── 11. PDF OF CURVES ───────────────────────────────
# Includes:
#   - regulated proteins (fail_reason = NA)
#   - unregulated proteins that only failed the FC cutoff
# Excludes:
#   - low-point, flat, bad-EC50, low-R², or fit-failure cases
 
plot_df <- summary_df[
  !is.na(summary_df$EC50_nM) &
    !is.na(summary_df$pEC50) &
    !is.na(summary_df$r_squared) &
    summary_df$r_squared >= MIN_R2 &
    (
      is.na(summary_df$fail_reason) |
      grepl("^fold change at max dose =", summary_df$fail_reason)
    ),
]
 
message("\nGenerating PDF for ", nrow(plot_df),
        " fitted proteins (regulated + unregulated passing fit filters)...")
 
log10_plot_min <- min(log10_doses) - 0.3
log10_plot_max <- max(log10_doses) + 0.3
plot_doses_nM  <- 10^seq(log10_plot_min, log10_plot_max, length.out = 500)
 
pdf_path_out <- file.path(
  OUTPUT_DIR,
  paste0("dose_response_curves_all_fitted_", safe_name, ".pdf")
)
pdf(pdf_path_out, width = 7, height = 5)
 
for (i in seq_len(nrow(plot_df))) {
  gene_id <- plot_df$gene[i]
  row_idx <- which(dat$gene == gene_id)[1]
  if (is.na(row_idx)) next
 
  rel_vals   <- rel_mat[row_idx, ]
  valid_mask <- !is.na(rel_vals) & is.finite(rel_vals)
  if (sum(valid_mask) < MIN_VALID_PTS) next
 
  dose_fit <- DOSES_NM[valid_mask]
  resp_fit <- rel_vals[valid_mask]
 
  fit <- tryCatch(
    suppressWarnings(
      drm(resp_fit ~ dose_fit,
          fct = LL.4(names = c("slope", "bottom", "top", "EC50")),
          control = drmc(maxIt = 1000, errorm = FALSE, noMessage = TRUE))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
 
  pred_line <- tryCatch(
    as.numeric(predict(fit, newdata = data.frame(dose_fit = plot_doses_nM))),
    error = function(e) NULL
  )
  if (is.null(pred_line) || anyNA(pred_line)) next
 
  obs_df  <- data.frame(log10_dose = log10(dose_fit), rel_int = resp_fit)
  line_df <- data.frame(log10_dose = log10(plot_doses_nM), rel_int = pred_line)
 
  curve_col <- if (identical(plot_df$direction[i], "up")) {
    "#c0392b"
  } else if (identical(plot_df$direction[i], "down")) {
    "#2980b9"
  } else {
    "#555555"
  }
 
  point_col <- if (identical(plot_df$direction[i], "up")) {
    "#e74c3c"
  } else if (identical(plot_df$direction[i], "down")) {
    "#3498db"
  } else {
    "#777777"
  }
 
  all_y  <- c(obs_df$rel_int, line_df$rel_int)
  y_rng  <- range(all_y, na.rm = TRUE)
  y_pad  <- max(diff(y_rng) * 0.15, 0.05)
  y_lims <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
 
  subtitle <- sprintf(
    "pEC50 = %.2f  |  EC50 = %.1f nM  |  slope = %.2f  |  top = %.2f  |  bottom = %.2f  |  R² = %.3f",
    plot_df$pEC50[i], plot_df$EC50_nM[i], plot_df$slope[i],
    plot_df$top[i], plot_df$bottom[i], plot_df$r_squared[i]
  )
 
  direction_label <- ifelse(is.na(plot_df$direction[i]), "unregulated", plot_df$direction[i])
 
  auc_label <- ifelse(
    is.na(plot_df$AUC[i]),
    "NA",
    sprintf("%.3f", plot_df$AUC[i])
  )
 
  fc_label <- ifelse(
    is.na(plot_df$fold_change_at_max[i]),
    "NA",
    sprintf("%.2f", plot_df$fold_change_at_max[i])
  )
 
  p <- ggplot() +
    geom_line(data = line_df, aes(x = log10_dose, y = rel_int),
              color = curve_col, linewidth = 1.0) +
    geom_point(data = obs_df, aes(x = log10_dose, y = rel_int),
               color = point_col, size = 3, shape = 16) +
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "grey55", linewidth = 0.4) +
    geom_vline(xintercept = log10(plot_df$EC50_nM[i]), linetype = "dotted",
               color = curve_col, linewidth = 0.4, alpha = 0.7) +
    scale_x_continuous(
      breaks = log10_doses,
      labels = DOSES_NM,
      limits = c(log10_plot_min, log10_plot_max),
      name = "Concentration (nM)"
    ) +
    scale_y_continuous(
      limits = y_lims,
      name = "Relative intensity (treatment / mean DMSO)"
    ) +
    labs(
      title    = gene_id,
      subtitle = subtitle,
      caption  = sprintf(
        "Drug: %s (ID: %s, Plate: %d)  |  direction: %s  |  FC at max dose: %s  |  AUC: %s  |  n = %d pts",
        DRUG_NAME, DRUG_ID, PLATE_NUM,
        direction_label, fc_label, auc_label, sum(valid_mask)
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 7.5, color = "grey30"),
      plot.caption     = element_text(size = 7, color = "grey50"),
      panel.grid.minor = element_blank()
    )
 
  print(p)
}
 
dev.off()
message("  ", pdf_path_out)
# ── 12. FINAL SUMMARY ─────────────────────────────────────────
message("\n========================================")
message("decryptE-like MaxQuant analysis complete")
message("========================================")
message("Drug    : ", DRUG_NAME, "  (ID: ", DRUG_ID, ")")
message("Plate   : ", PLATE_NUM,
        "  (DMSO columns: DMSO", dmso_start, "–DMSO", dmso_end, ")")
message("Output  : ", OUTPUT_DIR)
message("Files   :")
message("  summary_", safe_name, ".csv")
message("  summary_", safe_name, ".xlsx")
message("  dose_response_curves_", safe_name, ".pdf")
message("")
message("Regulated proteins: ", n_reg,
        "  (", n_up, " up / ", n_down, " down)")
message("FC thresholds: > ", FC_UP_THRESHOLD, " up  |  < ", FC_DOWN_THRESHOLD, " down")
message("========================================")
