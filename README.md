# Event-related potential signatures of auditory predictive processing in first-episode psychosis: a passive pattern-based study

# Muñoz-Caracuel et al. (2026)

Code for EEG processing, ERP measure extraction, statistical analyses and figures for the study. Statistical models adjust for **age, sex and years of education**. CNV, MMN and PINV are primary components; N1 and P2 following S1 are exploratory.

Participant data are not distributed with this repository. Running the analyses requires authorized access to the EEG matrices and phenotypic data.

## Files and their roles

| File | Purpose | How to use it |
| --- | --- | --- |
| `matlab/preprocessing_reference.m` | EEG preprocessing, artifact rejection, participant-level bin averages and pre-S4-baselined PINV matrices. | Review and execute the relevant sections in MATLAB; manual stages require intervention. |
| `matlab/binlist_S4.txt` | Defines the six S4 event bins used by ERPLAB. | Read automatically during preprocessing; do not execute separately. |
| `matlab/extract_erp_measures.m` | Extracts ERP amplitudes and exploratory measures, joins phenotypic data and creates the input table for R. | Call as a MATLAB function with three arguments (step 2). |
| `matlab/plot_figures.m` | Generates ERP waveform plots with mean ± SEM shading. | Run separately after preparing the matrices (step 4). |
| `run_analysis.R` | Runs the statistical analyses, exports numerical results and generates Figure 2. | Execute from the `Code` directory (step 3). |
| `R/functions.R` | Defines input checks, outcome families, estimation, standardization and FDR functions. | Loaded automatically by `run_analysis.R`. |
| `R/models.R` | Fits group comparisons, ERP–phenotype associations, ERP-by-group interactions and medication models. | Loaded automatically by `run_analysis.R`. |
| `R/figure2.R` | Draws the three-panel Figure 2 from the fitted models. | Loaded automatically by `run_analysis.R`. |
| `.gitignore` | Excludes common data and output formats from Git tracking. | Keep in the repository; it does not replace checking files before publication. |

## Requirements and folder layout

- MATLAB R2025b. EEG preprocessing also requires EEGLAB, ERPLAB, ICLabel, DIPFIT channel coordinates and Parallel Computing Toolbox for the parallel sections. Add EEGLAB and its required plugins to the MATLAB path.
- R 4.6.1 with `lmerTest`, `emmeans` and `ggplot2`. Install the packages once from R:

```r
install.packages(c("lmerTest", "emmeans", "ggplot2"))
```

The MATLAB preprocessing and waveform scripts resolve paths from their location in `Code/matlab`. Keep this local layout, or adapt their path configuration before execution:

```text
project/
├── Code/                 # Public code only
│   ├── README.md
│   ├── run_analysis.R
│   ├── R/
│   └── matlab/
├── EEG/                  # Private EEG datasets and processed matrices
├── private/              # Private participant-level phenotypic CSV
└── results/              # Generated extraction and statistical outputs
```

Only the contents of `Code` should be uploaded to the public repository. Replace the example `/path/to/project` paths below with your own absolute paths. Use new output directories for extraction and statistical analyses.

## Workflow

**EEG preprocessing → processed matrices → ERP extraction → statistical analyses and Figure 2.** ERP waveform plotting is a separate branch from the processed matrices; it is not required to run the statistical models.

### 1. Prepare the processed EEG matrices in MATLAB

Open `matlab/preprocessing_reference.m`. Review the input directories, manual inspection stages and processing sections before execution. The script expects EEG datasets already converted to EEGLAB `.set` format and uses the dataset filename patterns specified in each section.

The `runPreprocessing = false` protection deliberately prevents unattended execution of the full pipeline. Enable it only after reviewing the processing stages and use private working copies: preprocessing writes intermediate datasets and may replace existing outputs. Do not rerun preprocessing merely to regenerate figures or statistical results.

The participant-level averages needed for extraction are:

```text
EEG/ERPtotal_CON_MMN_0p05Hz.mat
EEG/ERPtotal_PEP_MMN_0p05Hz.mat
EEG/ERPtotal_CON_CNV_PINV_0p05Hz.mat
EEG/ERPtotal_PEP_CNV_PINV_0p05Hz.mat
```

For waveform plotting, section 8 additionally saves:

```text
EEG/ERPtotal_CON_PINV_preS4BL_0p05Hz.mat
EEG/ERPtotal_PEP_PINV_preS4BL_0p05Hz.mat
```

If the four averaged matrices already exist, skip the earlier preprocessing stages. To create missing PINV plotting matrices, execute only `CONFIGURACION DE RUTAS` and section 8. This section applies the pre-S4 baseline to the saved long-epoch averages without repeating filtering, ICA or trial rejection, and refuses to overwrite existing PINV output files. Figures are not generated automatically by preprocessing.

### 2. Extract ERP measures and create the analysis input

Prepare a private CSV with one row per participant. `Subject` must match the identifiers generated by the extractor (`C01`, `P01`, etc.); if `GROUP` is supplied, use `CON` and `FEP`. The extractor currently parses participant numbers from matrix subject names containing `Sujeto<number>_`; adapt this parsing if your deidentified naming convention differs.

The phenotypic CSV must include `AGE`, `SEX`, `YEARS_EDUCATION`, `CPZ_EQ` and the outcome columns listed in `families` in `R/functions.R`. Use consistent two-level sex coding and numeric scores. PANSS, SFS and dose values are required for FEP participants; they may be missing for controls. The R input checks require complete ERP measures, age, education and pooled phenotypic outcomes.

Run the following in the MATLAB Command Window:

```matlab
addpath('/path/to/project/Code/matlab');
extract_erp_measures( ...
    '/path/to/project/EEG', ...
    '/path/to/project/private/phenotypes.csv', ...
    '/path/to/project/results/extracted');
```

The three arguments specify the matrix directory, phenotypic CSV and a **new** output directory. This file is a function, not a standalone script to run without arguments.

The extractor uses the four MMN/long-epoch matrices above and applies the PINV pre-S4 baseline internally. It does not require the two separate PINV plotting matrices. It creates:

- `analysis_input.csv`: one row per participant and frontocentral electrode, containing ERP measures, outcomes and covariates.
- `localizer_windows.csv`: group-collapsed FCz peak locations and amplitude windows for MMN, N1 and P2.
- `waveforms.mat`: participant-level waveforms supporting extraction; keep this file private.

The input contains Fz, FCz, Cz, FC3 and FC4 only. Do not add ROI rows: the R script calculates each participant's ROI measures by averaging the five electrode-level measures. If a complete, compatible `analysis_input.csv` is already available, proceed directly to step 3.

### 3. Run the statistical analyses and generate Figure 2

In a terminal, change to the `Code` directory and run:

```bash
cd /path/to/project/Code
Rscript run_analysis.R /path/to/project/results/extracted/analysis_input.csv /path/to/project/results/statistics
```

The output directory must not already exist. `run_analysis.R` loads the three files in `R` automatically; do not run them separately. It validates the input, fits the models and exports:

| Output | Contents |
| --- | --- |
| `Clinical.csv` | Adjusted group comparisons for clinical and cognitive outcomes. |
| `ERP_groups.csv` | Adjusted ERP group contrasts at each electrode and in the ROI. |
| `Group_electrode.csv` | Group-by-electrode interaction tests. |
| `Posthoc_electrode_contrasts.csv` | Holm-adjusted comparisons between electrode-specific group contrasts when an interaction survives FDR correction. |
| `ERP_phenotype.csv` | ERP–outcome associations. |
| `ERP_group_slopes.csv` | ERP-by-group interaction tests for pooled outcomes. |
| `Table2_FCz_associations.csv` | Primary FCz associations with their site-specific FDR correction. |
| `Medication.csv` | Dose–ERP associations and ERP–outcome models adjusted for dose, with corresponding estimates without dose. |
| `Figure_2_data.csv` | Estimates plotted in Figure 2. |
| `Figure_2.png`, `Figure_2.tiff` | Figure 2, including a 300-ppi TIFF. A PDF is also exported when Cairo is available. |
| `sessionInfo.txt` | R and package versions used for the run. |

The CSV files provide numerical results for the manuscript tables; the scripts do not create publication-formatted Word or Excel tables.

### 4. Generate the ERP waveform figures separately

Once the matrices listed in step 1 are available, open `matlab/plot_figures.m` and click **Run** in MATLAB. Run the saved file rather than evaluating the entire file as a selection, because it contains local functions.

The script loads the EEG matrices directly, checks their dimensions and PINV metadata, and saves PNG/FIG waveform plots and associated plotting matrices under:

```text
EEG/Figuras_PINV_principal_preS4BL/PNG/
EEG/Figuras_PINV_principal_preS4BL/FIG/
EEG/Plot_data_PINV_principal_preS4BL/
```

Shading represents the mean ± one between-participant SEM, not a confidence interval or jackknife estimate. Plotting does not rerun statistical models or create Figure 2. Existing outputs are protected by `overwriteExistingFigures = false`; use a new output location or explicitly review that setting before regenerating plots.

## Statistical conventions

All models include age, sex and years of education. Electrode-level ERP group comparisons use participant-specific random intercepts and fixed effects for group, electrode and their interaction. Pooled ERP–outcome models additionally adjust for group; symptom and functioning models use the FEP group only.

Benjamini–Hochberg FDR correction separates primary and exploratory families. Primary ERP group contrasts comprise 15 tests (three components × five electrodes); N1/P2 contrasts comprise 10. ROI and interaction families are corrected separately. For ERP–outcome associations, `q` uses the five-electrode family where applicable, whereas `q_within_site` uses the site-specific family. The FCz summary for Table 2 uses the latter. These values should not be interchanged. Confidence intervals are ordinary 95% intervals, not multiplicity-adjusted intervals.

Medication analyses are restricted to FEP participants with `CPZ_EQ > 0` and use `log(1 + CPZ_EQ)`. ERP–outcome models with and without dose use the same medicated sample and retain age, sex and education. No direct dose–phenotype models or P2−N1 analyses are included.

## Data protection and publication

Keep participant-level data and generated outputs outside the public code directory. Review all files before uploading; `.gitignore` does not remove files already tracked by Git. Before public release, select an author-approved license and add the manuscript citation and repository identifier. Verify exported estimates against the final manuscript tables before using the outputs in a submission.
