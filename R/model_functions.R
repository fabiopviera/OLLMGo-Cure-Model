# -----------------------------------------------------------------------------
# MGo and OLLMGo cure models: likelihood and maximum likelihood estimation
#
# Parameterization:
#   eta   = logistic(xi), with 0 < eta < 1
#   beta  = exp(log_beta), with beta > 0
#   gamma = exp(log_gamma), with gamma > 0 (OLLMGo only)
#   lambda = log[1 + exp(-xi / gamma)]
#   alpha = beta * lambda (derived parameter)
#
# MGo is the special case gamma = 1.
# Event indicator: status = 1 for event, 0 for right censoring.
# -----------------------------------------------------------------------------

# Stable evaluation of log(1 + exp(x)).
softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Stable evaluation of log(exp(x) - 1), for x > 0.
log_expm1 <- function(x) {
  out <- numeric(length(x))
  small <- x <= log(2)
  out[small] <- log(expm1(x[small]))
  out[!small] <- x[!small] + log1p(-exp(-x[!small]))
  out
}

# Negative log-likelihood under independent right censoring.
# par = (xi, log_beta) for MGo, or (xi, log_beta, log_gamma) for OLLMGo.
neg_log_likelihood <- function(par, time, status) {
  beta <- exp(par[2])
  gamma <- if (length(par) == 2L) 1 else exp(par[3])
  lambda <- softplus(-par[1] / gamma)

  if (!all(is.finite(c(beta, gamma, lambda))) ||
      beta <= 0 || gamma <= 0 || lambda <= 0) {
    return(1e100)
  }

  q <- lambda * (-expm1(-beta * time))
  if (any(!is.finite(q)) || any(q <= 0)) return(1e100)

  log_z <- log_expm1(q)
  log_S <- -softplus(gamma * log_z)
  log_h <- log(gamma) + log(lambda) + par[2] - beta * time +
    q + (gamma - 1) * log_z + log_S

  value <- -sum(log_S + status * log_h)
  if (is.finite(value)) value else 1e100
}

# Fit either model by BFGS and transform standard errors using the delta method.
# start: a list of starting vectors on the unconstrained parameter scale.
fit_cure_model <- function(data, model = c("MGo", "OLLMGo"), start) {
  model <- match.arg(model)
  k <- if (model == "MGo") 2L else 3L

  stopifnot(is.data.frame(data), all(c("time", "status") %in% names(data)),
            nrow(data) > k, all(is.finite(data$time)),
            all(data$time > 0), all(data$status %in% c(0, 1)))
  stopifnot(is.list(start), length(start) > 0L,
            all(vapply(start, length, integer(1)) == k))

  fits <- lapply(start, function(initial) {
    optim(par = initial, fn = neg_log_likelihood,
          time = data$time, status = data$status,
          method = "BFGS", hessian = TRUE,
          control = list(maxit = 5000, reltol = 1e-10))
  })

  valid <- which(vapply(fits, function(f) {
    f$convergence == 0L && is.finite(f$value)
  }, logical(1)))
  if (length(valid) == 0L) stop("No optimization run converged.")

  fit <- fits[[valid[which.min(vapply(fits[valid], `[[`, numeric(1), "value"))]]]
  xi <- fit$par[1]
  beta <- exp(fit$par[2])
  gamma <- if (k == 2L) 1 else exp(fit$par[3])
  eta <- plogis(xi)
  alpha <- beta * softplus(-xi / gamma)

  estimates <- if (k == 2L) {
    c(eta = eta, beta = beta)
  } else {
    c(eta = eta, beta = beta, gamma = gamma)
  }

  # The observed information is the Hessian of the negative log-likelihood.
  H <- (fit$hessian + t(fit$hessian)) / 2
  eigenvalues <- if (all(is.finite(H))) {
    eigen(H, symmetric = TRUE, only.values = TRUE)$values
  } else {
    NA_real_
  }

  if (any(!is.finite(eigenvalues)) || min(eigenvalues) <= 0) {
    warning("Observed information is not positive definite; SEs set to NA.")
    se <- setNames(rep(NA_real_, k), names(estimates))
  } else {
    J <- diag(c(eta * (1 - eta), beta, if (k == 3L) gamma))
    V <- J %*% solve(H) %*% J
    se <- setNames(sqrt(diag(V)), names(estimates))
  }

  log_likelihood <- -fit$value
  n <- nrow(data)
  list(
    model = model,
    estimates = estimates,
    standard_errors = se,
    alpha = alpha,
    logLik = log_likelihood,
    AIC = -2 * log_likelihood + 2 * k,
    BIC = -2 * log_likelihood + k * log(n),
    optim = fit
  )
}
