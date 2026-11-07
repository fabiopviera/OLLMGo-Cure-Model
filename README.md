# OLLMGo defective cure model — R code and applications

This repository contains R code and data for maximum likelihood estimation and randomized residual diagnostics of the defective odd log-logistic modified Gompertz (OLLMGo) model and its modified Gompertz (MGo) submodel, as presented in the accompanying manuscript:

**Defective odd log-logistic modified Gompertz model with applications to medical data.**  
*Manuscript - full bibliographic details will be added when available.*

## Repository contents

| File | Description |
| --- | --- |
| `R/model_functions.R` | Shared model likelihoods, maximum likelihood estimation, standard errors, information criteria, and cure-fraction estimates. |
| `R/residual_functions.R` | Randomized Cox–Snell and quantile residuals, and QQ plots. |
| `R/application_escc.R` | Application 1: esophageal squamous cell carcinoma (ESCC). |
| `R/application_covid.R` | Application 2: COVID-19, including regression on sex and age group. |
| `data/escc.xlsx`, `data/covid19.xlsx` | Data used in the two applications. |

## Reproduce the analyses

Install the required R packages:

```r
install.packages(c("readxl", "ggplot2", "patchwork"))
```

From the repository root, run:

```r
source("R/application_escc.R")
source("R/application_covid.R")
```

The scripts report parameter estimates, standard errors, AIC, BIC, likelihood ratio tests, and randomized residual QQ plots (51 randomizations). Results are saved in `figures/` where applicable. The code does not include Monte Carlo simulations or survival/hazard plots.

The datasets use `time` and `status`, where `status = 1` denotes an observed event and `status = 0` right censoring. In the COVID-19 data, `sex` is coded 0 (female) or 1 (male), and `age60` is coded 0 (under 60) or 1 (60 or older).

## Applications and references

**Application 1 — ESCC (124 patients).** The original data were reported by Li et al. (2016) and subsequently analyzed using cure models by Martinez et al. (2022):

- Martínez, E. Z., Freitas, B. C. L., Achcar, J. A., Aragon, D. C., & Peres, M. V. O. (2022). *Exponentiated Weibull Models Applied to Medical Data in Presence of Right-censoring, Cure Fraction and Covariates.* Statistics, Optimization & Information Computing, 10(2), 548–571. https://doi.org/10.19139/soic-2310-5070-1266
- Li, J., et al. (2016). *Genetic Variants in the p14ARF/MDM2/TP53 Pathway Are Associated with the Prognosis of Esophageal Squamous Cell Carcinoma Patients Treated with Radical Resection.* PLOS ONE, 11(7), e0158613. https://doi.org/10.1371/journal.pone.0158613

**Application 2 — COVID-19 (139 patients).** These data were previously analyzed by Prataviera et al. (2024), who attributed the original data source to Dong et al. (2020):

- Prataviera, F., Hashimoto, E. M., Ortega, E. M. M., Cordeiro, G. M., Cancho, V. G., & Vila, R. (2024). *A new flexible regression model with application to recovery probability Covid-19 patients.* Journal of Applied Statistics, 51(5), 826–844. https://doi.org/10.1080/02664763.2022.2163229
- Dong, E., Du, H., & Gardner, L. (2020). *An interactive web-based dashboard to track COVID-19 in real time.* The Lancet Infectious Diseases, 20(5), 533–534. https://doi.org/10.1016/S1473-3099(20)30120-1

## Citation and reuse

When using the OLLMGo model or code, please cite the accompanying OLLMGo manuscript once its bibliographic details are available. When using the application data, also acknowledge the corresponding references above. Dataset redistribution permissions and a software license should be verified before public release.
