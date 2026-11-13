# ==============================================================================
# SIMULAÇÃO DE MONTE CARLO - OLL-MGo (Foco exclusivo nos Coeficientes de Regressão)
# ==============================================================================

rm(list = ls())
set.seed(123)
library(parallel)

# ------------------------------------------------------------------------------
# 1. Gerador de Dados Otimizado
# ------------------------------------------------------------------------------
simudata_smcm_weibull_reg <- function(n, gamma, cp, shape = 1.5, scale_weibull = 1.0) {
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  X  <- cbind(1, x1, x2)
  
  lin_pred <- pmin(pmax(drop(X %*% gamma[1:3]), -10), 10)
  eta      <- 1 / (1 + exp(-lin_pred))
  
  cured <- rbinom(n, 1, 1-eta)
  t_event <- rep(Inf, n)
  susceptible <- which(cured == 0)
  
  if (length(susceptible) > 0) {
    t_event[susceptible] <- rweibull(
      n = length(susceptible),
      shape = shape,
      scale = scale_weibull
    )
  }
  
  CensT  <- rexp(n, rate = cp)
  times  <- pmin(t_event, CensT)
  status <- as.numeric(t_event <= CensT)
  
  list(
    times  = times,
    status = status,
    X      = X,
    y      = cured
  )
}

# ------------------------------------------------------------------------------
# 2. Função de Log-Verossimilhança Direta na Escala Logarítmica
# ------------------------------------------------------------------------------
LogLik_fast <- function(param, data) {
  b <- param[4]
  g <- param[5]
  
  if (b <= 1e-5 || g <= 1e-5) return(1e10)
  
  X      <- data$X
  times  <- data$times
  status <- data$status
  
  eta <- 1 / (1 + exp(-drop(X %*% param[1:3])))
  eta <- pmin(pmax(eta, 1e-7), 1 - 1e-7)
  
  a_i <- b * log1p(((1 - eta) / eta)^(1 / g))
  
  exp_bt <- exp(-b * times)
  z      <- 1 - exp_bt
  
  ai_b_z <- (a_i / b) * z
  A      <- exp(-ai_b_z)
  B      <- exp(-g * ai_b_z)
  
  one_minus_A <- pmax(1 - A, 1e-12)
  Den         <- pmax(one_minus_A^g + B, 1e-12)
  
  log_f <- log(g) + log(a_i) - b * times - g * ai_b_z + (g - 1) * log(one_minus_A) - 2 * log(Den)
  log_s <- -g * ai_b_z - log(Den)
  
  val <- sum(status * log_f + (1 - status) * log_s)
  if (!is.finite(val)) return(1e10)
  
  return(-val)
}

# ------------------------------------------------------------------------------
# 3. Harness Paralelo de Monte Carlo (Apenas Coeficientes de Regressão)
# ------------------------------------------------------------------------------
run_monte_carlo_parallel_reg <- function(ns, param_true, cp, mc.samples, num_cores = NULL) {
  if (is.null(num_cores)) {
    num_cores <- max(1, detectCores() - 1)
  }
  
  # Apenas os 3 coeficientes de regressão (gama0, gama1, gama2)
  eval_params <- c("xi0", "xi1", "xi2")
  true_eval   <- param_true[1:3]
  p_eval      <- length(eval_params)
  
  # Construção dinâmica das colunas da tabela resumida
  col_names <- c("n", "%cens", "%cure")
  for (m in c("mc", "RB", "MCse", "Ase", "RMSE", "CP")) {
    col_names <- c(col_names, paste(m, eval_params, sep = "."))
  }
  col_names <- c(col_names, "%Conv", "%PosDef", "Total.Errors")
  
  Table <- matrix(NA_real_, nrow = length(ns), ncol = length(col_names))
  colnames(Table) <- col_names
  Table[, 1] <- ns
  
  cat(sprintf("Iniciando Cluster Paralelo com %d núcleos...\n", num_cores))
  cl <- makeCluster(num_cores)
  
  clusterExport(cl, c("simudata_smcm_weibull_reg", "LogLik_fast"))
  
  starting_grids <- list(
    param_true,
    param_true * c(0.9, 1.1, 0.95, 1.05, 1.0),
    param_true * c(1.1, 0.9, 1.05, 0.95, 1.1)
  )
  
  for (r in seq_along(ns)) {
    n <- ns[r]
    clusterExport(cl, c("n", "param_true", "cp", "starting_grids"), envir = environment())
    
    results_list <- parLapply(cl, seq_len(mc.samples), function(j) {
      set.seed(12345 + j + n * 10)
      data_sim <- simudata_smcm_weibull_reg(n = n, gamma = param_true[1:3], cp = cp)
      
      fitted_solutions <- list()
      
      for (g_idx in seq_along(starting_grids)) {
        init_p <- starting_grids[[g_idx]]
        
        fit_try <- tryCatch(
          optim(
            par     = init_p,
            fn      = LogLik_fast,
            data    = data_sim,
            method  = "L-BFGS-B",
            lower   = c(-10, -10, -10, 1e-4, 1e-4),
            upper   = c(10, 10, 10, 20, 20),
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
          if (any(evals <= 1e-6)) stop("Hessiana Não Positiva Definida")
          sqrt(diag(solve(H)))
        }, error = function(e) NULL)
        
        if (!is.null(mystd_try) && length(mystd_try) == length(param_true) && all(is.finite(mystd_try)) && all(mystd_try > 0)) {
          # Extrai apenas os 3 primeiros parâmetros (coeficientes de regressão)
          par_eval_j <- best_fit$par[1:3]
          std_eval_j <- mystd_try[1:3]
          
          Lower.CI <- par_eval_j - 1.959964 * std_eval_j
          Upper.CI <- par_eval_j + 1.959964 * std_eval_j
          cp_vec   <- as.numeric(Lower.CI <= true_eval & true_eval <= Upper.CI)
          
          return(list(
            par  = par_eval_j,
            std  = std_eval_j,
            cp   = cp_vec,
            cens = 100 * mean(1 - data_sim$status),
            cure = 100 * mean(data_sim$y),
            valid = TRUE
          ))
        }
      }
      
      return(list(valid = FALSE))
    })
    
    parhats       <- matrix(NA_real_, mc.samples, p_eval)
    stdhats       <- matrix(NA_real_, mc.samples, p_eval)
    CoverageProb  <- matrix(NA_real_, mc.samples, p_eval)
    censor_sample <- numeric(mc.samples)
    cure_sample   <- numeric(mc.samples)
    conv_count    <- 0
    
    for (j in seq_len(mc.samples)) {
      res <- results_list[[j]]
      if (res$valid) {
        parhats[j, ]    <- res$par
        stdhats[j, ]    <- res$std
        CoverageProb[j, ] <- res$cp
        censor_sample[j]  <- res$cens
        cure_sample[j]    <- res$cure
        conv_count        <- conv_count + 1
      }
    }
    
    col_mean      <- colMeans(parhats, na.rm = TRUE)
    col_mc_se     <- apply(parhats, 2, sd, na.rm = TRUE)
    col_ase       <- colMeans(stdhats, na.rm = TRUE)
    relative_bias <- (col_mean - true_eval) / abs(true_eval)
    rmse          <- sqrt(colMeans((sweep(parhats, 2, true_eval, "-"))^2, na.rm = TRUE))
    col_cp        <- colMeans(CoverageProb, na.rm = TRUE)
    
    Table[r, "n"]     <- n
    Table[r, "%cens"] <- round(mean(censor_sample, na.rm = TRUE), 3)
    Table[r, "%cure"] <- round(mean(cure_sample, na.rm = TRUE), 3)
    
    Table[r, paste0("mc.", eval_params)]   <- round(col_mean, 3)
    Table[r, paste0("RB.", eval_params)]   <- round(relative_bias, 3)
    Table[r, paste0("MCse.", eval_params)] <- round(col_mc_se, 3)
    Table[r, paste0("Ase.", eval_params)]  <- round(col_ase, 3)
    Table[r, paste0("RMSE.", eval_params)] <- round(rmse, 3)
    Table[r, paste0("CP.", eval_params)]   <- round(col_cp, 3)
    
    Table[r, "%Conv"]       <- round(100 * (conv_count / mc.samples), 2)
    Table[r, "%PosDef"]     <- round(100 * (conv_count / mc.samples), 2)
    Table[r, "Total.Errors"] <- mc.samples - conv_count
    
    cat(sprintf("Concluído bloco n = %d em paralelo!\n", n))
  }
  
  stopCluster(cl)
  return(as.data.frame(Table))
}

# ------------------------------------------------------------------------------
# 4. Execução
# ------------------------------------------------------------------------------
ssizes     <- c(100, 300, 500, 800, 1000)
param_true <- c(-1.5, -1.5, 0.3, 2, 1.49)
cp1        <- 0.71

Results_Diagnostics_Reg <- run_monte_carlo_parallel_reg(
  ns = ssizes,
  param_true = param_true,
  cp = cp1,
  mc.samples = 1000)

print(Results_Diagnostics_Reg)

