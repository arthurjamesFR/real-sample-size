# The sample size you really need

Sample size for trials with a binary outcome when patients' baseline risk varies.

**Online calculator: https://arthurjamesfr.github.io/real-sample-size/**

Sample sizes are usually calculated by applying the expected treatment effect to the average risk of the population. If the treatment lowers each patient's odds of the outcome by the same proportion and patients' risks differ, it prevents fewer events than this calculation assumes, and a trial sized this way is underpowered. In 60,243 trauma patients from the French Traumabase registry, the usual calculation underestimated the sample size of a mortality trial nearly 4-fold.

## Contents

| File | Description |
| --- | --- |
| `index.html` | Online calculator. A single self-contained page; it contains only an aggregated distribution of predicted risk (share of patients in 200 risk bands), no individual data. Risks entered by users are processed in their browser and are not transmitted. |
| `R/analyse.R` | Complete analysis code of the article (data preparation, multiple imputation, simulations, tables and figures). |

## Calculator

The calculator gives three sample sizes per arm for a two-arm superiority trial with 1:1 allocation:

- the usual calculation, with the expected effect applied to the average risk;
- the sample size actually required for an unadjusted analysis, with the effect applied to each patient's risk;
- the sample size required for an analysis adjusted for the logit of predicted risk.

The risk distribution can be entered as its mean and SD, as individual predicted risks, or taken from the Traumabase. The method is described on the page and in the supplementary material of the article (ESM Text S4).

## Analysis code

Requirements: R ≥ 4.5 with the packages `readxl`, `tidyverse`, `mice` and `patchwork`.

The Traumabase data are not included. They are available upon reasonable request, after approval by the Traumabase Scientific Committee and the participating centers.

To run the code on an authorized extract:

1. set `DATA_PATH` to the extract and `OUTPUT_DIR` / `LOCAL_DIR` to the output folders (STEP 0);
2. set `MODE <- "TEST"` for a quick run, or `"FINAL"` for the full analysis (several hours);
3. run `source("R/analyse.R")`.

Only aggregated outputs are written to `OUTPUT_DIR`. The key linking center codes to centers is written to `LOCAL_DIR` and must not be shared.

## Citation

James A, et al. Applying the treatment effect to mean mortality underestimates the sample size of mortality trials: a multicenter trauma registry-based simulation study. Manuscript in preparation.

## License

Code released under the MIT License.

## Contact

Arthur James — arthur.james@aphp.fr
