# Event-related potential signatures of auditory predictive processing in first-episode psychosis: a passive pattern-based study

**Muñoz-Caracuel et al. (2026)**

MATLAB and R code for EEG processing, ERP extraction, statistical analyses and figures. CNV, MMN and PINV are primary components; N1 and P2 following S1 are exploratory. Statistical models adjust for **age, sex and years of education**.

Participant data are not included. Authorized access to the EEG and phenotypic data is required to run the analyses.

## Files

| File | Purpose |
| --- | --- |
| `matlab/preprocessing_reference.m` | EEG preprocessing, artifact rejection and ERP averaging. |
| `matlab/binlist_S4.txt` | Event-bin definitions used during preprocessing. |
| `matlab/extract_erp_measures.m` | Extracts ERP measures and creates the input table for R. |
| `matlab/plot_figures.m` | Generates ERP waveform figures. |
| `matlab/erp_time_axis.m` | Shared time-axis checks; called automatically. |
| `run_analysis.R` | Runs the statistical analyses and generates Figure 2. |
| `R/functions.R` | Input checks, estimation, standardization and FDR correction. |
| `R/models.R` | Group comparisons, ERP–phenotype associations and medication models. |
| `R/figure2.R` | Figure 2 plotting code. |

The files in `R` are loaded by `run_analysis.R`; do not run them separately.

## Requirements

- **MATLAB R2025b**. Preprocessing requires EEGLAB, ERPLAB, ICLabel, DIPFIT and Parallel Computing Toolbox.
- **R 4.6.1**, with the following packages:

```r
install.packages(c("lmerTest", "emmeans", "ggplot2"))
```

Keep `Code` and the private `EEG` folder inside the same project directory. Replace the example paths below with your own paths. Keep data and generated results outside `Code`.

## Workflow

### 1. Preprocess EEG in MATLAB

Review and execute the relevant sections of `matlab/preprocessing_reference.m`, starting from EEGLAB `.set` datasets. Manual inspection is required; `runPreprocessing = false` prevents unattended execution. Use working copies, as preprocessing can overwrite intermediate files.

If averaged matrices are already available, skip preprocessing. To generate missing pre-S4-baselined PINV matrices for plotting, execute only the path configuration section and section 8.

### 2. Extract ERP measures

Prepare a phenotypic CSV with one row per participant and matching subject identifiers. Required covariates and outcome columns are defined in `R/functions.R`.

In MATLAB:

```matlab
addpath('/path/to/project/Code/matlab');
extract_erp_measures( ...
    '/path/to/project/EEG', ...
    '/path/to/project/private/phenotypes.csv', ...
    '/path/to/project/results/extracted');
```

This creates `analysis_input.csv`, extraction-window information and participant-level waveforms. The output directory must be new. If a compatible analysis input table is already available, proceed to step 3.

### 3. Run the statistical analyses

From a terminal:

```bash
cd /path/to/project/Code
Rscript run_analysis.R /path/to/project/results/extracted/analysis_input.csv /path/to/project/results/statistics
```

The output directory must be new. Results are exported as CSV files, together with Figure 2 (PNG, 300-ppi TIFF and PDF when available) and software-version information. The scripts export numerical results, not formatted Word or Excel tables.

FDR corrections are defined in the R scripts. For FCz associations in Table 2, use the dedicated `Table2_FCz_associations.csv` output, which applies the site-specific correction.

### 4. Plot ERP waveforms

Run `matlab/plot_figures.m` separately in MATLAB after preparing the matrices. It saves PNG/FIG waveform plots with mean ± SEM shading under `EEG/Figuras_PINV_principal_preS4BL`. This step is independent of the statistical analyses and does not generate Figure 2. Existing outputs are protected by default.

## Sharing

Upload only code and documentation. Do not upload participant-level data or generated individual-level outputs. Review the repository contents before publication; `.gitignore` is provided to help exclude data files.
