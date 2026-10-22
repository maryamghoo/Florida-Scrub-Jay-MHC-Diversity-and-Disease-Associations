# MHC associations with infection incidence and intensity

These two R scripts use `MHC_sample_alleles.csv` produced by `01_AmpliSAS_Postprocessing`. They share the same individual-level input but analyze different outcomes.

## Folder structure

```text
project/
├── MHC_Incidence_Models.R
├── MHC_Intensity_Models.R
├── README.md
├── r_input/
│   └── MHC_sample_alleles.csv
└── r_output/
    ├── incidence_models/
    └── intensity_models/
```

## Requirements

Install the R packages `dplyr`, `lme4`, `car`, `MuMIn`, `pROC`, and `lmtest`.

## Run

Set the working directory to `project/`, then run:

```bash
Rscript MHC_Incidence_Models.R
Rscript MHC_Intensity_Models.R
```

## Incidence analysis

Retains samples with binary `Incidence` (0 or 1), a recorded `Genetic_cluster`, and a cluster other than `Admixed`. The script calculates allele and, when available, supertype carrier counts, a numeric correlation matrix, and VIF screens. It defines 35 candidate binomial logistic mixed-effects models with a random intercept for `Year`; models requiring absent predictors or lacking variation are skipped, and fitting errors are caught.

Fitted models are ranked with `MuMIn::model.sel()`. **All fitted models with ΔAICc < 2 are selected automatically** for diagnostics and detailed summaries, rather than using a hardcoded list of model IDs. For each selected model, the script exports Pearson chi-square, residual degrees of freedom, the chi-square/df ratio, the chi-square upper-tail p-value, and a singular-fit indicator. It also produces an ROC plot and AUC CSV when ROC estimation succeeds. Fixed-effect estimates and odds ratios are exported for all fitted models.

## Intensity analysis

Starts from the same incidence/cluster filter, then keeps only infected samples (`Incidence == 1`) with finite `logIntensity`. It calculates allele and, when available, supertype carrier counts, a numeric correlation matrix, and VIF screens. It defines 90 candidate Gaussian linear mixed-effects models for `logIntensity`, with a random intercept for `Year`, fitted using maximum likelihood (`REML = FALSE`); models requiring absent predictors or lacking variation are skipped, and fitting errors are caught.

Fitted models are ranked with `MuMIn::model.sel()`. **All fitted models with ΔAICc < 2 are selected automatically** for diagnostics and detailed summaries. For each selected model, the script saves a residual Q–Q plot, a residual-versus-fitted plot, and a CSV containing the Breusch–Pagan test statistic, p-value, and singular-fit indicator. Fixed-effect estimates and exponentiated coefficients (`FoldChange`) are exported for all fitted models.

## Variables and model comparison

The scripts use `Incidence`, `Genetic_cluster`, and `Functional_Alleles` in place of the older names `Infection_Incidence`, `Genetic_Cluster`, and `Allele_count`. `Functional_Alleles` counts only functional alleles. `Genotype` categories occurring fewer than five times are grouped under `Rare` separately within each analysis dataset.

The input may lack `Supertype_count`, `SupertypeIIC`, or other `SupertypeII*` columns. Such variables are used only when present; candidate models requiring missing columns are omitted. No supertype values are generated.

Within each analysis, all successfully fitted candidate models use the same complete-case rows across the union of available candidate predictors, making their model-selection scores comparable on the same observations. The two outcomes have separate analysis datasets. Original model IDs are retained, including the duplicate intensity specifications (`mod31`–`mod33` repeat earlier candidate formulas). Model count therefore refers to *defined candidates*, not necessarily the number successfully fitted.

The file `model_selection_AIC.csv` retains its existing name, although the scripts select their best-supported models using the `delta` column of the `MuMIn::model.sel()` table (ΔAICc < 2). A model must have `delta < 2`, not `delta <= 2`, to be included in the diagnostic and summary exports.

## Outputs

Files are written to `r_output/incidence_models/` or `r_output/intensity_models/`, respectively:

| Output | Description |
| --- | --- |
| `allele_carrier_counts.csv` | Carrier totals for allele columns in the filtered analysis dataset. |
| `supertype_carrier_counts.csv` | Supertype carrier totals, only when supertype columns are available. |
| `correlation_matrix.csv` | Pairwise-complete correlations among eligible numeric variables. |
| `vif_results.csv` | VIF screening results for the full, reduced, and count/genotype predictor sets, when screening succeeds. |
| `model_selection_AIC.csv` | Ranked fitted-model comparison table, including `delta`. |
| `model_coefficients.csv` | Fixed-effect estimates for fitted models, with `OddsRatio` (incidence) or `FoldChange` (intensity). |
| `best_models_summary.txt` | Printed summaries of all fitted models with ΔAICc < 2. |

**Additional incidence outputs for each model with ΔAICc < 2:** `<model_id>_diagnostics.csv`, and, if ROC estimation succeeds, `<model_id>_ROC.png` and `<model_id>_AUC.csv`.

**Additional intensity outputs for each model with ΔAICc < 2:** `<model_id>_QQ.png`, `<model_id>_residuals.png`, and `<model_id>_diagnostics.csv`.

