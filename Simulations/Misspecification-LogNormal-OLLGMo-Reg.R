# ==============================================================================
# EXPERIMENTO DE MÁ ESPECIFICAÇÃO: DGP Log-Normal -> Ajuste OLLMGo
# Relato Focado Exclusivamente nos Coeficientes de Regressão (xi0, xi1, xi2)
# ==============================================================================

rm(list = ls())
set.seed(123)

# ------------------------------------------------------------------------------
# 1. Funções do Modelo OLLMGo
# ------------------------------------------------------------------------------

## Density Function
den.func <- function(data, param) {
  X <- cbind(1, data$x1, data$x2)
  gammas <- param[1:3]
  b <- param[4]
  g <- param[5]
  
  eta <- plogis(drop(X %*% gammas))
  a_i <- b * log(1 + ((1 - eta) / eta)**(1 / g))
  
  t <- data$times
  z <- 1 - exp(-b * t)
  A <- exp(-(a_i / b) * z)
  B <- exp(-(g * a_i / b) * z)
  
  Num <- g * a_i * exp(-b * t) * B * (1 - A)**(g - 1)
  Den <- ((1 - A)**g + B)**2
  
  density <- Num / Den
  return(density)
}

## Survival Function
surv.func <- function(data, param) {
  X <- cbind(1, data$x1, data$x2)
  gammas <- param[1:3]
  b <- param[4]
  g <- param[5]
  
  eta <- plogis(drop(X %*% gammas))
  a_i <- b * log(1 + ((1 - eta) / eta)**(1 / g))
  
  t <- data$times
  z <- 1 - exp(-b * t)
  
  Num <- exp(-(g * a_i / b) * z)
  Den <- (1 - exp(-(a_i / b) * z))**g + exp(-(g * a_i / b) * z)
  
  surv <- Num / Den
  return(surv)
}

## Log-likelihood Function
LogLik <- function(data, param) {
  dj <- data$status
  p_cure <- plogis(drop(cbind(1, data$x1, data$x2) %*% param[1:3]))
  
  if ((param[4] > 0) & (param[5] > 0)) {
    f_val <- pmax(den.func(data, param), 1e-15)
    s_val <- pmax(surv.func(data, param), 1e-15)
    p_cure <- pmax(pmin(p_cure, 1 - 1e-12), 1e-12)
    
    term1 <- dj * (log(1 - p_cure) + log(f_val))
    term2 <- (1 - dj) * log(p_cure + (1 - p_cure) * s_val)
    
    val <- sum(term1 + term2)
    if (!is.finite(val)) return(Inf)
    return((-1) * val)
  } else {
    return(Inf)
  }
}

# ------------------------------------------------------------------------------
# 2. Gerador de Dados: Modelo de Mistura Padrão com Log-Normal e Regressão (DGP)
# ------------------------------------------------------------------------------
simudata_lognormal_reg <- function(n, gamma, cp, meanlog = 0.0, sdlog = 1.0) {
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  X  <- cbind(1, x1, x2)
  
  lin_pred <- pmin(pmax(drop(X %*% gamma[1:3]), -5), 5)
  p_cure   <- plogis(lin_pred)
  
  cured <- rbinom(n, 1, 1-p_cure)
  
  t_event <- rep(Inf, n)
  susceptible <- which(cured == 0)
  
  if (length(susceptible) > 0) {
    t_event[susceptible] <- rlnorm(
      n = length(susceptible),
      meanlog = meanlog,
      sdlog = sdlog
    )
  }
  
  cens_times <- rexp(n, rate = cp)
  times  <- pmin(t_event, cens_times)
  status <- as.numeric(t_event <= cens_times)
  
  data.frame(
    times  = times,
    status = status,
    x1     = x1,
    x2     = x2,
    y      = cured
  )
}

# ------------------------------------------------------------------------------
# 3. Harness de Monte Carlo Focado em (xi0, xi1, xi2)
# ------------------------------------------------------------------------------
run_monte_carlo_regression_only <- function(ns, param_true, cp, meanlog_val, sdlog_val, mc.samples) {
  p <- length(param_true)
  
  # Tabela restrita e limpa contendo apenas as métricas para xi0, xi1 e xi2
  Table <- matrix(NA_real_, nrow = length(ns), ncol = 23)
  colnames(Table) <- c(
    "n", "%cens", "%cure", 
    "xi0.mc", "xi1.mc", "xi2.mc", 
    "RB.xi0", "RB.xi1", "RB.xi2", 
    "MCse.xi0", "MCse.xi1", "MCse.xi2", 
    "Ase.xi0", "Ase.xi1", "Ase.xi2", 
    "RMSE.xi0", "RMSE.xi1", "RMSE.xi2", 
    "CP.xi0", "CP.xi1", "CP.xi2", 
    "%Conv", "error.blg"
  )
  Table[, 1] <- ns
  
  starting_grids <- list(
    param_true,
    param_true * c(1.0, 1.0, 1.0, 0.8, 1.2),
    param_true * c(1.0, 1.0, 1.0, 1.2, 0.8),
    param_true + c(0.1, -0.1, 0.1, 0.0, 0.0)
  )
  
  for (r in seq_along(ns)) {
    n <- ns[r]
    
    parhats      <- matrix(NA_real_, mc.samples, p)
    stdhats      <- matrix(NA_real_, mc.samples, p)
    CoverageProb <- matrix(NA_real_, mc.samples, p)
    
    censor_sample <- numeric(mc.samples)
    cure_sample   <- numeric(mc.samples)
    conv_sample   <- numeric(mc.samples)
    err.alg       <- 0
    
    for (j in seq_len(mc.samples)) {
      valid_fit <- FALSE
      
      data_sim <- simudata_lognormal_reg(
        n = n, 
        gamma = param_true[1:3], 
        cp = cp, 
        meanlog = meanlog_val, 
        sdlog = sdlog_val
      )
      
      for (attempt in seq_along(starting_grids)) {
        init_par <- starting_grids[[attempt]]
        
        fit_try <- tryCatch(
          optim(
            par     = init_par,
            fn      = LogLik,
            data    = data_sim,
            method  = "L-BFGS-B",
            lower   = c(-10, -10, -10, 1e-4, 1e-4),
            upper   = c(10, 10, 10, Inf, Inf),
            hessian = TRUE
          ),
          error = function(e) NULL
        )
        
        if (!is.null(fit_try) && fit_try$convergence == 0 && all(is.finite(fit_try$par)) && !is.null(fit_try$hessian)) {
          
          mystd_try <- tryCatch({
            H <- (fit_try$hessian + t(fit_try$hessian)) / 2
            evals <- eigen(H, symmetric = TRUE)$values
            if (any(evals <= 1e-8)) stop("Hessiana não positiva definida")
            sqrt(diag(solve(H)))
          }, error = function(e) NULL)
          
          if (!is.null(mystd_try) && length(mystd_try) == p && all(is.finite(mystd_try)) && all(mystd_try > 0)) {
            parhats[j, ]      <- fit_try$par
            stdhats[j, ]      <- mystd_try
            
            lower_ci <- fit_try$par - qnorm(0.975) * mystd_try
            upper_ci <- fit_try$par + qnorm(0.975) * mystd_try
            CoverageProb[j, ] <- as.numeric(lower_ci <= param_true & param_true <= upper_ci)
            
            conv_sample[j]   <- 1
            censor_sample[j] <- 100 * mean(1 - data_sim$status)
            cure_sample[j]   <- 100 * mean(data_sim$y)
            valid_fit        <- TRUE
            break
          }
        }
      }
      
      if (!valid_fit) {
        conv_sample[j] <- 0
        err.alg        <- err.alg + 1
      }
      
      cat(sprintf("Log-Normal -> OLLMGo (Regressão) | n = %d | Réplica %d/%d\r", n, j, mc.samples))
    }
    cat("\n")
    
    col_mean  <- colMeans(parhats, na.rm = TRUE)
    bias_ord  <- (col_mean - param_true) / abs(param_true)         
    col_mc_se <- apply(parhats, 2, sd, na.rm = TRUE) 
    col_ase   <- colMeans(stdhats, na.rm = TRUE)     
    col_rmse  <- sqrt(colMeans((sweep(parhats, 2, param_true, "-"))^2, na.rm = TRUE))
    col_cp    <- colMeans(CoverageProb, na.rm = TRUE)
    
    # Preenchimento restrito às colunas de interesse (n, %cens, %cure e índices 1:3 para xi0, xi1, xi2)
    Table[r, 2]     <- round(mean(censor_sample, na.rm = TRUE), 3)
    Table[r, 3]     <- round(mean(cure_sample, na.rm = TRUE), 3)
    Table[r, 4:6]   <- round(col_mean[1:3], 3)
    Table[r, 7:9]   <- round(bias_ord[1:3], 3)
    Table[r, 10:12] <- round(col_mc_se[1:3], 3)
    Table[r, 13:15] <- round(col_ase[1:3], 3)
    Table[r, 16:18] <- round(col_rmse[1:3], 3)
    Table[r, 19:21] <- round(col_cp[1:3], 3)
    Table[r, 22]    <- round(100 * mean(conv_sample, na.rm = TRUE), 2)
    Table[r, 23]    <- err.alg
  }
  
  return(Table)
}

# ------------------------------------------------------------------------------
# 4. Execução da Simulação
# ------------------------------------------------------------------------------
param_true  <- c(-1.5, -1.5, 0.3, 1.0, 1.2) # [xi0, xi1, xi2, b, g]
meanlog_val <- 0.0
sdlog_val   <- 1.0
cp_gen      <- 0.71
ssizes      <- c(100, 300, 500, 800, 1000)

Results_Regression_Only <- run_monte_carlo_regression_only(
  ns          = ssizes,
  param_true  = param_true,
  cp          = cp_gen,
  meanlog_val = meanlog_val,
  sdlog_val   = sdlog_val,
  mc.samples  = 1000 # Ajuste para 1000 na rodada definitiva
)

print("--- RESULTADOS FINAIS: COEFICIENTES DE REGRESSÃO (xi0, xi1, xi2) ---")
print(Results_Regression_Only)
