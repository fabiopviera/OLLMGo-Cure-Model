# -----------------------------------------------------------------------------
# Randomized residual diagnostics for the MGo and OLLMGo cure models
#
# RCS = -log(U) and RQ = Phi^{-1}(1 - U), where:
#   U = S(t)                 for an observed event (status = 1),
#   U = V * S(t), V ~ U(0,1) for independent right censoring (status = 0).
#
# For QQ plots, residuals are sorted within each randomization; the plotted
# sample quantile at each rank is the median across B randomizations.
# Both fitted models use the SAME uniform random numbers.
# These plots are descriptive diagnostics, not calibrated goodness-of-fit tests.
# -----------------------------------------------------------------------------

# Individual fitted survival probabilities at the observed times.
# Requires softplus() and log_expm1() from model_functions.R.
fitted_survival <- function(fit, data, X = NULL) {
  stopifnot(is.data.frame(data),
            all(c("time", "status") %in% names(data)),
            all(is.finite(data$time)), all(data$time > 0),
            all(data$status %in% 0:1), !anyNA(data$status))

  if (is.null(X)) {
    X <- matrix(1, nrow(data), 1L)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  stopifnot(nrow(X) == nrow(data), ncol(X) == fit$p,
            all(is.finite(X)), fit$model %in% c("MGo", "OLLMGo"))

  p <- fit$p
  par <- fit$optim$par
  linear_predictor <- drop(X %*% par[seq_len(p)])
  beta <- exp(par[p + 1L])
  gamma <- if (fit$model == "MGo") 1 else exp(par[p + 2L])
  lambda <- softplus(-linear_predictor / gamma)
  q <- lambda * (-expm1(-beta * data$time))

  log_S <- if (fit$model == "MGo") {
    -q
  } else {
    -softplus(gamma * log_expm1(q))
  }
  S <- exp(log_S)
  if (any(!is.finite(S)) || any(S <= 0) || any(S >= 1)) {
    stop("Invalid survival probabilities. Check model fit and input data.")
  }
  S
}

# Produce theoretical and sample quantiles for both randomized residuals.
# 'fits' must be a named list: list(OLLMGo = ollmgo, MGo = mgo).
# X = NULL for the intercept-only ESCC models; use the full design matrix
# for the COVID-19 models with sex and age-group covariates.
randomized_residual_qq <- function(data, fits, X = NULL,
                                   B = 51L, seed = 123L) {
  stopifnot(is.list(fits),
            setequal(names(fits), c("OLLMGo", "MGo")),
            length(B) == 1L, is.finite(B), B >= 1, B == as.integer(B),
            length(seed) == 1L, is.finite(seed),
            nrow(data) >= 2L)

  n <- nrow(data)
  set.seed(seed)
  V <- matrix(stats::runif(n * B), nrow = n, ncol = B)
  is_event <- data$status == 1L
  probabilities <- stats::ppoints(n)

  one_model <- function(fit) {
    S <- fitted_survival(fit, data, X)
    U <- sweep(V, 1L, S, `*`)
    U[is_event, ] <- S[is_event]

    # Avoid infinite normal or exponential quantiles from finite precision.
    U <- pmin(pmax(U, 1e-12), 1 - 1e-12)

    ordered_rcs <- vapply(seq_len(B),
                          function(b) sort(-log(U[, b])), numeric(n))
    ordered_rq <- vapply(seq_len(B),
                         function(b) sort(stats::qnorm(U[, b],
                                                       lower.tail = FALSE)),
                         numeric(n))
    list(
      RCS = data.frame(Theoretical = stats::qexp(probabilities),
                       Observed = apply(ordered_rcs, 1L, stats::median),
                       Model = fit$model),
      RQ = data.frame(Theoretical = stats::qnorm(probabilities),
                      Observed = apply(ordered_rq, 1L, stats::median),
                      Model = fit$model)
    )
  }

  results <- lapply(fits[c("OLLMGo", "MGo")], one_model)
  list(RCS = do.call(rbind, lapply(results, `[[`, "RCS")),
       RQ = do.call(rbind, lapply(results, `[[`, "RQ")),
       B = as.integer(B), seed = seed)
}

# Return a two-panel ggplot/patchwork figure (RCS on the left, RQ on the right).
plot_residual_qq <- function(qq_data) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Install 'ggplot2' and 'patchwork' to display residual QQ plots.")
  }

  plot_one <- function(df, title) {
    ggplot2::ggplot() +
      ggplot2::geom_abline(intercept = 0, slope = 1,
                           color = "gray40", linewidth = 0.7) +
      ggplot2::geom_point(
        data = subset(df, Model == "MGo"),
        ggplot2::aes(x = Theoretical, y = Observed, color = Model),
        size = 1.5, alpha = 0.85) +
      ggplot2::geom_point(
        data = subset(df, Model == "OLLMGo"),
        ggplot2::aes(x = Theoretical, y = Observed, color = Model),
        size = 1.5, alpha = 0.85) +
      ggplot2::scale_color_manual(
        values = c("OLLMGo" = "black", "MGo" = "gray60"),
        breaks = c("OLLMGo", "MGo")) +
      ggplot2::labs(title = title, x = "Theoretical quantiles",
                    y = "Sample quantiles") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 9),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.95),
        legend.justification = c(0, 1))
  }

  patchwork::wrap_plots(
    plot_one(qq_data$RCS, "Randomized Cox-Snell residuals"),
    plot_one(qq_data$RQ, "Randomized quantile residuals"),
    ncol = 2)
}
