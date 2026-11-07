# -----------------------------------------------------------------------------
# MGo and OLLMGo defective cure models: shared estimation functions
#
# Cure regression: logit(eta_i) = X_i %*% xi.
# beta > 0 and gamma > 0 are common across observations.
# lambda_i = log(1 + exp[-(X_i %*% xi)/gamma]); alpha_i = beta * lambda_i.
# MGo is nested in OLLMGo at gamma = 1.
# status = 1: observed event; status = 0: independent right censoring.
# -----------------------------------------------------------------------------

# Numerically stable log(1 + exp(x)).
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Numerically stable log(exp(x) - 1), defined for x > 0.
log_expm1 <- function(x) {
  out <- numeric(length(x))
  small <- x <= log(2)
  out[small] <- log(expm1(x[small]))
  out[!small] <- x[!small] + log1p(-exp(-x[!small]))
  out
}

# Negative log-likelihood for exact failures and independent right censoring.
# par = (xi, log_beta) for MGo; (xi, log_beta, log_gamma) for OLLMGo.
neg_log_likelihood <- function(par, time, status, X, model) {
  p <- ncol(X)
  lp <- drop(X %*% par[seq_len(p)])
  beta <- exp(par[p + 1L])
  gamma <- if (model == "MGo") 1 else exp(par[p + 2L])

  if (!all(is.finite(c(beta, gamma, lp))) || beta <= 0 || gamma <= 0) {
    return(1e100)
  }

  lambda <- softplus(-lp / gamma)
  if (any(!is.finite(lambda)) || any(lambda <= 0)) return(1e100)
  q <- lambda * (-expm1(-beta * time))
  if (any(!is.finite(q)) || any(q <= 0)) return(1e100)

  if (model == "MGo") {
    log_S <- -q
    log_h <- log(beta) + log(lambda) - beta * time
  } else {
    log_z <- log_expm1(q)
    log_S <- -softplus(gamma * log_z)
    log_h <- log(gamma) + log(beta) + log(lambda) - beta * time +
      q + (gamma - 1) * log_z + log_S
  }

  value <- -sum(log_S + status * log_h)
  if (is.finite(value)) value else 1e100
}

# Maximum likelihood estimation by BFGS with multiple starting values.
# X = NULL fits an intercept-only cure fraction; otherwise supply a design
# matrix with one row per observation (including the intercept column).
fit_cure_model <- function(data, model = c("MGo", "OLLMGo"), start, X = NULL) {
  model <- match.arg(model)
  stopifnot(is.data.frame(data), all(c("time", "status") %in% names(data)),
            all(is.finite(data$time)), all(data$time > 0),
            !anyNA(data$status), all(data$status %in% 0:1))

  if (is.null(X)) {
    X <- matrix(1, nrow(data), 1L, dimnames = list(NULL, "xi0"))
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (is.null(colnames(X))) colnames(X) <- paste0("xi", seq_len(ncol(X)) - 1L)
  stopifnot(nrow(X) == nrow(data), ncol(X) >= 1L,
            all(is.finite(X)))

  p <- ncol(X)
  k <- p + if (model == "MGo") 1L else 2L
  stopifnot(nrow(data) > k, is.list(start), length(start) >= 1L,
            all(vapply(start, length, integer(1)) == k),
            all(vapply(start, function(v) all(is.finite(v)), logical(1))))

  fits <- lapply(start, function(initial) {
    optim(par = initial, fn = neg_log_likelihood,
          time = data$time, status = data$status, X = X, model = model,
          method = "BFGS", hessian = TRUE,
          control = list(maxit = 5000, reltol = 1e-10))
  })

  valid <- which(vapply(fits, function(f) {
    f$convergence == 0L && is.finite(f$value)
  }, logical(1)))
  if (length(valid) == 0L) stop("No optimization run converged.")
  best <- valid[which.min(vapply(fits[valid], `[[`, numeric(1), "value"))]
  fit <- fits[[best]]

  xi <- fit$par[seq_len(p)]
  beta <- exp(fit$par[p + 1L])
  gamma <- if (model == "MGo") 1 else exp(fit$par[p + 2L])
  lp <- drop(X %*% xi)
  eta <- plogis(lp)
  alpha <- beta * softplus(-lp / gamma)
  intercept_only <- p == 1L && all(X[, 1L] == 1)

  # Return eta rather than xi for an intercept-only model; regression
  # coefficients remain on the logit scale when covariates are present.
  if (intercept_only) {
    estimates <- c(eta = eta[1L], beta = beta,
                   if (model == "OLLMGo") c(gamma = gamma))
    jacobian_diag <- c(eta[1L] * (1 - eta[1L]), beta,
                       if (model == "OLLMGo") gamma)
  } else {
    estimates <- c(setNames(xi, colnames(X)), beta = beta,
                   if (model == "OLLMGo") c(gamma = gamma))
    jacobian_diag <- c(rep(1, p), beta,
                       if (model == "OLLMGo") gamma)
  }

  # Inverse observed information; delta method on the reported scale.
  H <- (fit$hessian + t(fit$hessian)) / 2
  eig <- if (all(is.finite(H))) {
    eigen(H, symmetric = TRUE, only.values = TRUE)$values
  } else NA_real_

  if (any(!is.finite(eig)) || min(eig) <= 0) {
    warning("Observed information is not positive definite; SEs are NA.")
    covariance <- matrix(NA_real_, k, k)
    standard_errors <- rep(NA_real_, k)
  } else {
    covariance <- solve(H)
    J <- diag(jacobian_diag, nrow = k)
    standard_errors <- sqrt(diag(J %*% covariance %*% J))
  }
  names(standard_errors) <- names(estimates)

  log_likelihood <- -fit$value
  n <- nrow(data)
  list(
    model = model,
    estimates = estimates,
    standard_errors = standard_errors,
    alpha = alpha,
    p = p,
    covariance = covariance,
    logLik = log_likelihood,
    AIC = -2 * log_likelihood + 2 * k,
    BIC = -2 * log_likelihood + k * log(n),
    optim = fit
  )
}

# Cure fraction and its delta-method standard error for new covariate profiles.
# X_new must contain the same design-matrix columns as the fitted model.
predict_cure_fraction <- function(fit, X_new) {
  X_new <- as.matrix(X_new)
  stopifnot(ncol(X_new) == fit$p, all(is.finite(X_new)))
  xi <- fit$optim$par[seq_len(fit$p)]
  eta <- plogis(drop(X_new %*% xi))
  gradient <- sweep(X_new, 1L, eta * (1 - eta), `*`)
  Vxi <- fit$covariance[seq_len(fit$p), seq_len(fit$p), drop = FALSE]
  se <- sqrt(rowSums((gradient %*% Vxi) * gradient))
  data.frame(eta = as.numeric(eta), SE = as.numeric(se))
}
