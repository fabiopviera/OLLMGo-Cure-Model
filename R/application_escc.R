# -----------------------------------------------------------------------------
# Application 1: Esophageal squamous cell carcinoma (ESCC)
# Fit the MGo and OLLMGo cure models without covariates.
# Run from the repository root: source("R/application_escc.R")
# Requires the readxl package to import the Excel dataset.
# -----------------------------------------------------------------------------

source("R/model_functions.R")
source("R/residual_functions.R")

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

data <- as.data.frame(readxl::read_excel("data/escc.xlsx", sheet = "ESCC"))
stopifnot(nrow(data) == 124L, sum(data$status == 1) == 41L,
          sum(data$status == 0) == 83L)

# Fit the MGo model (gamma fixed at 1).
mgo <- fit_cure_model(
  data, model = "MGo",
  start = list(c(qlogis(0.1), log(0.1)))
)

# Fit the OLLMGo model using multiple starting values for gamma.
ollmgo_starts <- lapply(c(0.7, 1.5, 3), function(gamma_start) {
  c(mgo$optim$par, log(gamma_start))
})
ollmgo <- fit_cure_model(data, model = "OLLMGo", start = ollmgo_starts)

# Parameter estimates and delta-method standard errors.
parameter_table <- rbind(
  data.frame(Model = "MGo", Parameter = names(mgo$estimates),
             Estimate = unname(mgo$estimates),
             SE = unname(mgo$standard_errors)),
  data.frame(Model = "OLLMGo", Parameter = names(ollmgo$estimates),
             Estimate = unname(ollmgo$estimates),
             SE = unname(ollmgo$standard_errors))
)

# Log-likelihoods and information criteria.
model_comparison <- data.frame(
  Model = c("MGo", "OLLMGo"),
  Parameters = c(2L, 3L),
  LogLik = c(mgo$logLik, ollmgo$logLik),
  AIC = c(mgo$AIC, ollmgo$AIC),
  BIC = c(mgo$BIC, ollmgo$BIC)
)

# Likelihood ratio test of H0: gamma = 1 (one degree of freedom).
LR <- 2 * (ollmgo$logLik - mgo$logLik)
if (LR < -1e-6) warning("Check model optimization: negative LR statistic.")
p_value <- pchisq(max(LR, 0), df = 1, lower.tail = FALSE)

cat("\nESCC: n =", nrow(data),
    "| events =", sum(data$status == 1),
    "| right-censored =", sum(data$status == 0), "\n")
cat("\nMaximum likelihood estimates (standard errors):\n")
print(transform(parameter_table,
                Estimate = round(Estimate, 4), SE = round(SE, 4)),
      row.names = FALSE)
cat("\nDerived alpha parameters:\n")
print(c(MGo = unname(mgo$alpha), OLLMGo = unname(ollmgo$alpha)))
cat("\nModel comparison:\n")
print(transform(model_comparison,
                LogLik = round(LogLik, 4),
                AIC = round(AIC, 4), BIC = round(BIC, 4)),
      row.names = FALSE)
cat("\nLR test (gamma = 1): statistic =", round(LR, 4),
    "| p-value =", signif(p_value, 5), "\n")

# -----------------------------------------------------------------------------
# Randomized residual QQ plots for the fitted intercept-only models.
# The same 51 randomizations are used for the two competing models.
# -----------------------------------------------------------------------------
qq_data <- randomized_residual_qq(
  data = data,
  fits = list(OLLMGo = ollmgo, MGo = mgo),
  B = 51L, seed = 23L
)
residual_figure <- plot_residual_qq(qq_data)
print(residual_figure)

# Save a vector graphic of the residual diagnostics.
##dir.create("figures", showWarnings = FALSE)
#ggplot2::ggsave("figures/Figure5_ESCC.pdf", residual_figure,
 #               width = 10, height = 4.2)
#if (capabilities("cairo")) {
#  ggplot2::ggsave("figures/Figure5_ESCC.eps", residual_figure,
 #                 width = 10, height = 4.2, device = grDevices::cairo_ps)
#}
