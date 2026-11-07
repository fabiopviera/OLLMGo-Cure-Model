# -----------------------------------------------------------------------------
# Application 2: COVID-19 survival data (China; n = 139)
# MGo and OLLMGo models with and without cure-fraction regression.
# Run from the repository root: source("R/application_covid.R")
# -----------------------------------------------------------------------------

source("R/model_functions.R")
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the 'readxl' package to import the Excel dataset.")
}

data <- as.data.frame(readxl::read_excel("data/covid19.xlsx", sheet = "COVID19"))
stopifnot(all(c("time", "status", "sex", "age60") %in% names(data)),
          nrow(data) == 139L, sum(data$status == 1L) == 33L,
          sum(data$status == 0L) == 106L,
          !anyNA(data), all(data$sex %in% 0:1), all(data$age60 %in% 0:1))

# sex: 0 = female, 1 = male; age60: 0 = under 60, 1 = 60 or older.
X <- model.matrix(~ sex + age60, data = data)
colnames(X) <- c("xi0", "xi_sex", "xi_age60")

# Intercept-only models.
mgo0 <- fit_cure_model(data, "MGo", start = list(
  c(qlogis(0.75), log(0.05)),
  c(qlogis(0.50), log(0.10))
))
ollmgo0 <- fit_cure_model(data, "OLLMGo", start = list(
  c(mgo0$optim$par, log(2)),
  c(qlogis(0.75), log(0.15), log(7))
))

# Models with sex and age-group effects on the cure fraction.
mgo <- fit_cure_model(data, "MGo", X = X, start = list(
  c(3, -1, -3, log(0.03)),
  c(2, -1, -3, log(0.06))
))
ollmgo <- fit_cure_model(data, "OLLMGo", X = X, start = list(
  c(mgo$optim$par, log(2)),
  c(3, -1, -4, log(0.15), log(7))
))

# Maximum likelihood estimates and standard errors.
fits <- list(
  "MGo: no covariates" = mgo0,
  "OLLMGo: no covariates" = ollmgo0,
  "MGo: sex + age60" = mgo,
  "OLLMGo: sex + age60" = ollmgo
)
for (name in names(fits)) {
  fit <- fits[[name]]
  cat("\n", name, "\n", sep = "")
  print(round(cbind(Estimate = fit$estimates,
                    SE = fit$standard_errors), 4))
}

# Information criteria and nested likelihood ratio tests (gamma = 1).
comparison <- data.frame(
  Model = c("MGo", "OLLMGo", "MGo", "OLLMGo"),
  Structure = c("Intercept", "Intercept", "Sex + age60", "Sex + age60"),
  logLik = vapply(list(mgo0, ollmgo0, mgo, ollmgo), `[[`, numeric(1), "logLik"),
  AIC = vapply(list(mgo0, ollmgo0, mgo, ollmgo), `[[`, numeric(1), "AIC"),
  BIC = vapply(list(mgo0, ollmgo0, mgo, ollmgo), `[[`, numeric(1), "BIC")
)
cat("\nModel comparison:\n")
print(comparison, row.names = FALSE, digits = 6)

LR_intercept <- 2 * (ollmgo0$logLik - mgo0$logLik)
LR_regression <- 2 * (ollmgo$logLik - mgo$logLik)
cat("\nLR test, intercept only: ", round(LR_intercept, 4),
    " (p = ", signif(pchisq(LR_intercept, 1, lower.tail = FALSE), 4), ")\n", sep = "")
cat("LR test, sex + age60: ", round(LR_regression, 4),
    " (p = ", signif(pchisq(LR_regression, 1, lower.tail = FALSE), 4), ")\n", sep = "")

# Cure fractions for the four covariate profiles (full regression models).
profiles <- data.frame(sex = c(0, 1, 0, 1), age60 = c(0, 0, 1, 1))
X_profiles <- model.matrix(~ sex + age60, data = profiles)
colnames(X_profiles) <- colnames(X)
eta_mgo <- predict_cure_fraction(mgo, X_profiles)
eta_ollmgo <- predict_cure_fraction(ollmgo, X_profiles)

cure_table <- data.frame(
  Group = c("Female <60", "Male <60", "Female >=60", "Male >=60"),
  n = vapply(seq_len(nrow(profiles)), function(i) {
    sum(data$sex == profiles$sex[i] & data$age60 == profiles$age60[i])
  }, integer(1)),
  MGo = eta_mgo$eta, SE_MGo = eta_mgo$SE,
  OLLMGo = eta_ollmgo$eta, SE_OLLMGo = eta_ollmgo$SE
)
cat("\nCure fractions by covariate profile:\n")
print(cure_table, row.names = FALSE, digits = 5)

# Sample-standardized cure fraction: mean of individual predictions.
cat("\nOverall estimated cure fractions:\n")
print(round(c(
  MGo = mean(predict_cure_fraction(mgo, X)$eta),
  OLLMGo = mean(predict_cure_fraction(ollmgo, X)$eta)
), 4))
