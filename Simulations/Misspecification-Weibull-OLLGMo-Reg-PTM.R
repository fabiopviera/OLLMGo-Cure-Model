# ==============================================================================
# SIMULAÇÃO DE MONTE CARLO - OLL-MGo COM REGRESSÃO NA FRAÇÃO DE CURA
# DGP: Weibull de Tempo de Promoção | eta_i = exp(-exp(X %*% gammas))
# Focado em: xi0, xi1 e xi2
# ==============================================================================

rm(list = ls())
set.seed(123)

# ------------------------------------------------------------------------------
# 1. Gerador de Dados Estável (Weibull de Tempo de Promoção)
# ------------------------------------------------------------------------------
simudata_weibull_promotion_stable <- function(n, gamma, cp, shape = 1.5, scale_weibull = 1.0) {
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  X  <- cbind(1, x1, x2)
  
  lin_pred <- pmin(pmax(drop(X %*% gamma[1:3]), -5), 5)
  theta    <- exp(lin_pred)
  eta      <- exp(-theta)
  
  M <- rpois(n, lambda = theta)
  
  SurvT <- rep(Inf, n)
  susceptible <- which(M > 0)
  
  for (i in susceptible) {
    latent_times <- rweibull(n = M[i], shape = shape, scale = scale_weibull)
    SurvT[i] <- min(latent_times)
  }
  
  CensT  <- rexp(n, rate = cp)
  times  <- pmin(SurvT, CensT)
  status <- as.numeric(SurvT <= CensT)
  
  data.frame(
    times  = times,
    status = status,
    x1     = x1,
    x2     = x2,
    y      = as.numeric(M == 0)
  )
}

# ------------------------------------------------------------------------------
# 2. Log-Verossimilhança Estável do Modelo OLL-MGo com Regressão
# ------------------------------------------------------------------------------
ollmgo_reg_loglik_stable <- function(param, data) {
  gammas <- param[1:3]
  beta   <- param[4]
  g      <- param[5]
  
  if (beta <= 1e-4 || g <= 1e-4) return(Inf)
  
  X <- cbind(1, data$x1, data$x2)
  lin_pred <- pmin(pmax(drop(X %*% gammas), -5), 5)
  theta    <- exp(lin_pred)
  eta      <- exp(-theta)
  eta      <- pmin(pmax(eta, 1e-12), 1 - 1e-12)
  
  a_i <- beta * log(1 + ((1 - eta) / eta)^(1 / g))
  t   <- data$times
  dj  <- data$status
  
  z  <- 1 - exp(-beta * t)
  S0 <- exp(-(a_i / beta) * z)
  G  <- pmax(pmin(1 - S0, 1 - 1e-12), 1e-12)
  S0 <- pmax(pmin(S0, 1 - 1e-12), 1e-12)
  
  g_pdf <- a_i * exp(-beta * t) * S0
  g_pdf <- pmax(g_pdf, 1e-12)
  
  G_g  <- G^g
  S0_g <- S0^g
  Den  <- G_g + S0_g
  
  S_t <- S0_g / Den
  f_t <- (g * g_pdf * (G^(g - 1)) * (S0^(g - 1))) / (Den^2)
  
  f_t <- pmax(f_t, 1e-15)
  S_t <- pmax(S_t, 1e-15)
  
  loglik <- sum(dj * log(f_t) + (1 - dj) * log(S_t))
  if (!is.finite(loglik)) return(Inf)
  return(-loglik)
}

# ------------------------------------------------------------------------------
# 3. Harness de Monte Carlo Completo (Apenas xi0, xi1 e xi2)
# ------------------------------------------------------------------------------
run_monte_carlo_diagnostics <- function(ns, param_true, cp, mc.samples) {
  p <- length(param_true)
  
  Table <- matrix(NA_real_, nrow = length(ns), ncol = 22)
  colnames(Table) <- c(
    "n", "%cens", "%cure", 
    "xi0.mc", "xi1.mc", "xi2.mc", 
    "RB.xi0", "RB.xi1", "RB.xi2", 
    "MCse.xi0", "MCse.xi1", "MCse.xi2", 
    "Ase.xi0", "Ase.xi1", "Ase.xi2", 
    "RMSE.xi0", "RMSE.xi1", "RMSE.xi2", 
    "CP.xi0", "CP.xi1", "CP.xi2", 
    "%Conv"
  )
  Table[, 1] <- ns
  
  for (r in seq_along(ns)) {
    n <- ns[r]
    
    parhats      <- matrix(NA_real_, mc.samples, p)
    stdhats      <- matrix(NA_real_, mc.samples, p)
    CoverageProb <- matrix(NA_real_, mc.samples, p)
    
    censor_sample <- numeric(mc.samples)
    cure_sample   <- numeric(mc.samples)
    conv_sample   <- numeric(mc.samples)
    
    starting_grids <- list(
      param_true,
      param_true * c(0.9, 1.1, 0.95, 1.05, 1.0),
      param_true * c(1.1, 0.9, 1.05, 0.95, 1.1)
    )
    
    for (j in seq_len(mc.samples)) {
      valid_fit <- FALSE
      data_sim <- simudata_weibull_promotion_stable(n = n, gamma = param_true[1:3], cp = cp)
      
      fitted_solutions <- list()
      
      for (g_idx in seq_along(starting_grids)) {
        init_p <- starting_grids[[g_idx]]
        
        fit_try <- tryCatch(
          optim(
            par     = init_p,
            fn      = ollmgo_reg_loglik_stable,
            data    = data_sim,
            method  = "L-BFGS-B",
            lower   = c(-5, -5, -5, 1e-3, 1e-3),
            upper   = c(5, 5, 5, 10, 10),
            hessian = TRUE
          ),
          error = function(e) NULL
        )
        
        if (!is.null(fit_try) && fit_try$convergence == 0 && all(is.finite(fit_try$par)) && !is.null(fit_try$hessian)) {
          fitted_solutions[[length(fitted_solutions) + 1]] <- fit_try
        }
      }
      
      if (length(fitted_solutions) > 0) {
        vals <- sapply(fitted_solutions, function(x) x$value)
        best_fit <- fitted_solutions[[which.min(vals)]]
        
        H <- (best_fit$hessian + t(best_fit$hessian)) / 2
        mystd_try <- tryCatch({
          evals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
          if (any(evals <= 1e-6)) stop("Não positiva definida")
          sqrt(diag(solve(H)))
        }, error = function(e) NULL)
        
        if (!is.null(mystd_try) && length(mystd_try) == p && all(is.finite(mystd_try)) && all(mystd_try > 0)) {
          parhats[j, ] <- best_fit$par
          stdhats[j, ] <- mystd_try
          
          Lower.CI <- best_fit$par - qnorm(0.975) * mystd_try
          Upper.CI <- best_fit$par + qnorm(0.975) * mystd_try
          CoverageProb[j, ] <- as.numeric(Lower.CI <= param_true & param_true <= Upper.CI)
          
          conv_sample[j]   <- 1
          censor_sample[j] <- 100 * mean(1 - data_sim$status)
          cure_sample[j]   <- 100 * mean(data_sim$y)
          valid_fit        <- TRUE
        }
      }
      
      if (!valid_fit) {
        conv_sample[j]   <- 0
      }
    }
    
    # Extrai apenas as métricas referentes aos 3 primeiros parâmetros (xi0, xi1, xi2)
    col_mean      <- colMeans(parhats, na.rm = TRUE)[1:3]
    col_mc_se     <- apply(parhats, 2, sd, na.rm = TRUE)[1:3]
    col_ase       <- colMeans(stdhats, na.rm = TRUE)[1:3]
    relative_bias <- (col_mean - param_true[1:3]) / abs(param_true[1:3])
    rmse          <- sqrt(colMeans((sweep(parhats[, 1:3, drop = FALSE], 2, param_true[1:3], "-"))^2, na.rm = TRUE))
    col_cp        <- colMeans(CoverageProb, na.rm = TRUE)[1:3]
    
    Table[r, 2]     <- round(mean(censor_sample, na.rm = TRUE), 3)
    Table[r, 3]     <- round(mean(cure_sample, na.rm = TRUE), 3)
    Table[r, 4:6]   <- round(col_mean, 3)
    Table[r, 7:9]   <- round(relative_bias, 3)
    Table[r, 10:12] <- round(col_mc_se, 3)
    Table[r, 13:15] <- round(col_ase, 3)
    Table[r, 16:18] <- round(rmse, 3)
    Table[r, 19:21] <- round(col_cp, 3)
    Table[r, 22]    <- round(100 * mean(conv_sample, na.rm = TRUE), 2)
    
    cat(sprintf("Concluído bloco n = %d\n", n))
  }
  
  return(Table)
}

# ------------------------------------------------------------------------------
# 4. Execução da Simulação
# ------------------------------------------------------------------------------
ssizes <- c(100, 300, 500, 800, 1000)
param_true <- c(-1.5, -1.5, 0.3, 2, 1.49)
cp1 <- 0.71

Results_Diagnostics <- run_monte_carlo_diagnostics(
  ns = ssizes,
  param_true = param_true,
  cp = cp1,
  mc.samples = 1000)

print(Results_Diagnostics)
