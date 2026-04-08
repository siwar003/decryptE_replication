# decryptE Replication Using FragPipe and R

## Overview

This repository contains a reproducible workflow to replicate and extend the dose-response proteomics analysis described in the decryptE framework. The pipeline processes label-free quantification (LFQ) data generated from FragPipe and applies a plate-aware normalization strategy followed by sigmoidal curve fitting to identify drug-induced protein regulation.

The goal of this project is to develop a scalable and interpretable pipeline for analyzing large-scale drug screening datasets across multiple doses, time-points, and cell-cycle states.

---

## Key Features

* 4-parameter log-logistic (LL.4) dose-response modeling using the `drc` R package
* Extraction of key pharmacological metrics:
  * Fold change at highest dose
  * EC50 and pEC50
  * Hill slope
  * Top and bottom asymptotes
  * Area under the curve (AUC)
  * R² (goodness of fit)
* Filtering of biologically meaningful responses 

---

## Workflow

### 1. Data Input

* Input: MaxQuant's `protein.Groups_fdr0.01.xlsx` file containing LFQ intensities
** Each column corresponds to a treatment or DMSO control sample (144 treatments, 5 doses each, and 48 DMSO samples)


### 2. Column Naming Convention

To enable automated parsing, columns must follow this format:

**Treatment samples:**
LFQ Intensity drugID dose
Exmaple: 
LFQ Intensity 16 1 

**DMSO controls:**
LFQ Intensity DMSO#
Example: 
LFQ Intensity DMSO38

---

### 3. Data Preprocessing
* Map drug of interest to its respective DMSO control for normalization using the Mapping_Sheet
* Replace zero values with `NA`
* For rows with multiple gene names separated by semicolon, take the first gene name

---

### 4. Plate-Matched Normalization

Each treatment intensity is normalized to the mean DMSO intensity from the same plate:

relative intensity = treatment LFQ / mean(DMSO LFQ for that plate)

This preserves plate-specific variation and avoids cross-plate bias.

---

### 5. Dose-Response Curve Fitting

For each protein and drug:

* Fit a 4-parameter log-logistic model using the `drc` package:

  * slope
  * lower asymptote (bottom)
  * upper asymptote (top)
  * EC50

* Convert EC50 to pEC50:
  pEC50 = 9 − log10(EC50 in nM)

* Require at least 4 valid dose points per protein

---

### 6. Quality Filtering

Curves are retained if they satisfy:

* EC50 within plausible range relative to tested doses
* |top − bottom| ≥ 0.15 (non-flat response)
* R² ≥ 0.8 (goodness of fit)

---

### 7. Biological Filtering

Proteins are classified as regulated if:

* Upregulated: fold change at max dose > 1.5
* Downregulated: fold change at max dose < 0.7

---

### 8. Output

The pipeline generates:

* Summary table (CSV and Excel):

  * all extracted parameters per protein and drug
* Long-format dataset for visualization and downstream analysis
* Dose-response curve plots (PDF) for regulated proteins

---

## Repository Structure

scripts/ → analysis scripts (parsing, normalization, curve fitting)
metadata/ → naming conventions and plate templates
data/ → input data (not tracked in Git)
results/ → output files (not tracked in Git)
README.md → project documentation

---

## Requirements

R packages:

* drc
* readxl
* writexl
* dplyr
* tidyr
* ggplot2
* pracma

Install with:
install.packages(c("drc","readxl","writexl","dplyr","tidyr","ggplot2","pracma"))

---

## Notes

* EC50 values are computed in nM and converted to pEC50 in molar units
* The log-logistic model inherently accounts for logarithmic dose scaling
* Plotting is performed on a log10 concentration axis for interpretability

---

## Goals

* Extension to time- and cell-cycle-state-resolved datasets
* Incorporation of machine learning-based curve classification (e.g., decryptE RF model)
* Application to low-input or clinical proteomics datasets


