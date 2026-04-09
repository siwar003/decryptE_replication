# decryptE Replication Using FragPipe and R

# decryptE Replication (FragPipe + R)

## Overview  
Pipeline to replicate decryptE dose–response proteomics analysis using FragPipe LFQ data. Performs plate-aware normalization and 4-parameter log-logistic (LL.4) modeling to identify regulated proteins.

---

## Workflow  

### Input  
- MaxQuant `proteinGroups_fdr0.01.xlsx`  
- Column format:  
  - Treatment: `LFQ Intensity drugID dose`  
  - DMSO: `LFQ Intensity DMSO#`

### Preprocessing  
- Map treatments to plate-specific DMSOs (Mapping Sheet)  
- Replace 0 → NA  
- Keep first gene if multiple  

### Normalization  
- Relative intensity = treatment / mean(DMSO of same plate)

### Curve Fitting  
- LL.4 model (`drc`)  
- Extract: fold change (FC), EC50, pEC50, slope, top/bottom, AUC, R²  
- ≥4 valid dose points required  
- pEC50 = 9 − log10(EC50 in nM)

### Filtering  
- |top − bottom| ≥ 0.15  
- R² ≥ 0.8  
- Regulated: FC > 1.5 or < 0.7  

---

## Output  
- Summary table (Excel/CSV)  
- Long-format dataset  
- Dose-response plots (PDF, regulated proteins)

---

## Goals

* Extension to time- and cell-cycle-state-resolved datasets
* Incorporation of machine learning-based curve classification (e.g., decryptE RF model)
* Application to low-input or clinical proteomics datasets


