# OLLMGo defective cure model — estimation and residual diagnostics

R code for maximum likelihood estimation and randomized residual diagnostics of the defective odd log-logistic modified Gompertz (OLLMGo) model and its modified Gompertz (MGo) submodel. The repository covers two applications from the accompanying manuscript.

## Repository contents

| Path | Description |
| --- | --- |
| `R/model_functions.R` | Shared likelihood, BFGS estimation, delta-method standard errors, cure-fraction prediction, AIC and BIC. |
| `R/residual_functions.R` | Shared randomized Cox–Snell and quantile residuals and two-panel QQ plots. |
| `R/application_escc.R` | Application 1: esophageal squamous cell carcinoma (ESCC), without covariates. |
| `R/application_covid.R` | Application 2: COVID-19 survival data, with and without covariate regression. |
| `data/escc.xlsx` | ESCC data: `time` (months) and `status`. |
| `data/covid19.xlsx` | COVID-19 data: `time` (days), `status`, `sex` and `age60`. |

The event coding is the same for both datasets: `status = 1` for an observed event and `status = 0` for right censoring. In the COVID-19 data, `sex = 0` denotes female and `sex = 1` male; `age60 = 0` denotes age below 60, and `age60 = 1` age 60 or older.

## Requirements

- R (base R functions are used for estimation)
- `readxl` for reading the Excel files
- `ggplot2` and `patchwork` for the two-panel residual QQ plots

Install the packages with `install.packages(c("readxl", "ggplot2", "patchwork"))`.

## Reproduce the results

Run these commands from the repository root:

```r
source("R/application_escc.R")
source("R/application_covid.R")
```

The scripts report maximum likelihood estimates, standard errors, AIC, BIC, and nested likelihood ratio tests of `H0: gamma = 1`. The COVID-19 script additionally reports fitted cure fractions for the four sex–age profiles and the sample-standardized mean cure fraction. Scripts do not refit models to individual sex or age groups. Each script also produces a two-panel randomized residual QQ plot in the `figures/` directory: `Figure5_ESCC.pdf` or `Figure9_residuals_COVID19.pdf`. When Cairo graphics are available, the corresponding EPS file is also saved.

### Expected information criteria

| Application | Model | Regression | AIC | BIC |
| --- | --- | --- | ---: | ---: |
| ESCC | MGo | Intercept only | 470.74 | 476.38 |
| ESCC | OLLMGo | Intercept only | 462.48 | 470.94 |
| COVID-19 | MGo | Intercept only | 400.93 | 406.80 |
| COVID-19 | OLLMGo | Intercept only | 371.35 | 380.15 |
| COVID-19 | MGo | Sex + age60 | 339.78 | 351.52 |
| COVID-19 | OLLMGo | Sex + age60 | 305.17 | 319.84 |

Expected LR statistics: ESCC approximately 10.2624; COVID-19 with covariates approximately 36.6166. Minor discrepancies may arise from optimizer tolerances.

## Model parameterization

The limiting survival probability is `eta_i`, with `logit(eta_i) = X_i xi`. The parameters `beta > 0` and `gamma > 0` are common across individuals. The defective modified Gompertz parameter is derived from the cure fraction as

`lambda_i = log[1 + exp{-(X_i xi)/gamma}]`, `alpha_i = beta * lambda_i`.

MGo is obtained at `gamma = 1`. The likelihood is evaluated under independent right censoring. Standard errors use the inverse observed-information matrix and the delta method. Numerical identifiability should be assessed for any new dataset.

## Data preparation and scope

The ESCC spreadsheet was transcribed from the analysis vectors supplied for this project (124 observations; 41 events). The standardized COVID-19 spreadsheet retains the 139 records from `data_china.xlsx` (33 events); its columns map from the original `time`, `cens`, `gender.c`, and `age.c`, respectively. The COVID-19 dataset was previously analyzed by Prataviera et al. (2024); see the manuscript for its source and citations.

This repository contains **core estimation code and randomized residual QQ plots**. It does not reproduce the survival/hazard figures or Monte Carlo simulations. Before publishing the datasets, verify that redistribution is permitted and include the corresponding original-data citations. Add an appropriate software license separately.

## Randomized residual diagnostics

Both applications use `source("R/residual_functions.R")` and the common
`randomized_residual_qq()` and `plot_residual_qq()` functions. For an observed
event, the probability-scale residual is `U = S(time)`; for an independently
right-censored observation, it is `U = V * S(time)`, with `V ~ Uniform(0, 1)`.
The randomized Cox–Snell and randomized quantile residuals are respectively
`-log(U)` and `qnorm(U, lower.tail = FALSE)`.

The QQ plots use the median of ordered residuals over **51 randomizations**
(with common uniform draws for MGo and OLLMGo). The plots are descriptive,
not formal goodness-of-fit tests; the reference distributions apply to
individual randomized residuals, not exactly to their rank-wise medians.
The ESCC plot uses the intercept-only fits; the COVID-19 plot uses the
**complete sex + age60 regression fits**, without stratified reestimation.
The random seeds are 23 (ESCC) and 579 (COVID-19). The figures reproduce the residual diagnostics described in the manuscript; minor graphical differences may arise from R/ggplot2 versions.
