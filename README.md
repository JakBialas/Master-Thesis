# Master-Thesis
# Treatment Efficacy Score (TES) – Data Augmentation Analysis

This repository contains R code used in a master's thesis to analyze the effect of various data augmentation methods on the stability and clinical correlation of the **Treatment Efficacy Score (TES)** in the context of breast cancer trials.

---

## Project Structure

```
DaneTestowe/
├── Data for code testing/
│   └── Test_Data.csv          # Sample test data for running the scripts
└── Kod/
    ├── functions/             # Helper functions (sourced automatically)
    │   ├── CalculateTES.R         # TES and p-value computation
    │   ├── wKS.R                  # Weighted Kolmogorov-Smirnov test (primary method)
    │   ├── w_KS_stat.R / w_KS_stat_fast.R
    │   ├── DensDiff.R / DensRatio.R
    │   ├── RCB_diff.R / RCB_ratio.R / RCB_diff_fast.R / RCB_ratio_fast.R
    │   ├── augment_data_BootS.R       # Augmentation: Bootstrap
    │   ├── augment_data_BootSwLN.R    # Augmentation: Bootstrap + Log-Normal noise
    │   ├── augment_data_random.R      # Augmentation: random sampling
    │   ├── augment_data_smote.R       # Augmentation: SMOTE
    │   ├── generate_cox_data.R        # Data generation (Cox model)
    │   ├── generate_smoteCox_data.R   # Data generation (SMOTE + Cox)
    │   ├── generate_smote_data_2Ev.R
    │   ├── generate_smoteCox_data2Ev.R
    │   ├── calculate_rmst_single.R    # RMST computation for a single case
    │   ├── calculate_rmst_with_TES.R  # Joining TES results with RMST
    │   ├── w_eCDF.R                   # Weighted empirical CDF
    │   └── check_libs.R               # Library checking and installation
    │
    ├── BootS_TES.R            # Analysis: Bootstrap
    ├── BootS_TESv2.R
    ├── BootSwLN_TES.R         # Analysis: Bootstrap + Log-Normal
    ├── BootSwLN_TESv2.R
    ├── Cox_TES.R              # Analysis: Cox model
    ├── Hybrid_TES.R           # Analysis: Hybrid method
    ├── Random_TES.R           # Analysis: Random augmentation
    ├── Random_TESv2.R
    ├── Smote_TES.R            # Analysis: SMOTE
    ├── Smote_TESv2.R
    ├── Smote2E_TES.R          # Analysis: SMOTE two-event
    ├── Smote2E_TESv2.R
    ├── SmoteCox_TES2Evv2.R    # Analysis: SMOTE + Cox (two-event)
    └── SmoteCox_TESv2.R       # Analysis: SMOTE + Cox
```

---

## Method Overview

### Treatment Efficacy Score (TES)

TES is a statistic based on the **weighted Kolmogorov-Smirnov test (wKS)**, comparing the RCB (*Residual Cancer Burden*) score distributions between a control group and an experimental treatment group. A higher TES indicates greater treatment efficacy.

Steps for computing TES:
1. Compute the wKS statistic between the RCB distributions of the control and experimental groups.
2. Run a permutation test (default: `n_perm = 10,000`).
3. Adjust TES by subtracting the expected value from the permutation null distribution.
4. Compute a one-sided p-value.

Available methods for computing TES (controlled by the `method` parameter in `CalculateTES.R`):
- `wKS` – weighted KS test *(default and recommended)*
- `DensRatio` – density ratio
- `DensDiff` – density difference

### Data Augmentation Methods

Each main script tests a different augmentation method, scanning augmentation levels from 0% to 100% in steps of 1 percentage point:

| Method | Script | Description |
|---|---|---|
| Bootstrap | `BootS_TES.R` | Sampling with replacement from original data |
| Bootstrap + Log-Normal | `BootSwLN_TES.R` | Bootstrap with added log-normal noise |
| Random | `Random_TES.R` | Random sampling from a fitted distribution |
| SMOTE | `Smote_TES.R` | Synthetic augmentation using SMOTE |
| SMOTE + Cox | `SmoteCox_TESv2.R` | SMOTE with a Cox survival model |
| Hybrid | `Hybrid_TES.R` | Combination of multiple augmentation methods |
| Cox | `Cox_TES.R` | Data generation from a Cox proportional hazards model |

---

## Requirements

### R (≥ 4.0)

Required packages (installed automatically by `check_libs.R` if missing):

```r
c("ggplot2", "densratio", "parallel", "reshape", "pcg",
  "gridExtra", "dplyr", "writexl", "tidyr", "smotefamily", "survRM2")
```

---

## Running the Scripts

1. Set the working directory to the folder containing the script (`setwd(...)`).
2. Place the `Rscripts/` (or `functions/`) directory with helper functions in the same location.
3. Update the path to the input data file (`BC_Data.csv`) inside the script.
4. Source the desired script, e.g.:

```r
source("BootS_TES.R")
```

> **Note:** Scripts use parallel processing via `mclapply`. The default is 40 cores — adjust the `numCores` parameter to match your machine.

---

## Input Data Format

A CSV file with the following columns:

| Column | Description |
|---|---|
| `Subtype` | Breast cancer subtype (e.g. `HR+/HER2+`, `TNBC`) |
| `arm` | Trial arm: `Control` or treatment regimen name |
| `rcbindex` | Residual Cancer Burden (RCB) index |
| `efs.time` | Time to event (EFS) |
| `efs.ind` | Event indicator (0 = censored, 1 = event) |

The `Data for code testing/` directory contains sample test data (`Test_Data.csv`).

---

## Output

Each script produces:

- **PNG plots** (saved in a method-specific subdirectory, e.g. `plots_Bootstrap/`):
  - TES change as a function of augmentation level
  - TES variance as a function of augmentation level
  - Scatter plot: Delta RMST vs. TES for each augmentation level
  - Combined 2×2 grid: Pearson R, Spearman R, slope, and intercept vs. augmentation

- **An Excel file** (e.g. `Bootstrap.xlsx`) with a results table containing:
  - Mean TES and its variance for each (Subtype, Regimen, Percent) combination
  - Pearson and Spearman correlation coefficients between TES and Delta RMST
  - Linear regression parameters (slope and intercept)

---

## Author

Jakub Białas
