rm(list = ls())

## Density Function (Modelo Ajustado: OLL-MGo)
den.func <- function(data, param) {
  X <- cbind(1, data$x1, data$x2)
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  
  eta <- plogis(drop(X %*% xis))
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

## Survival Function (Modelo Ajustado: OLL-MGo)
surv.func <- function(data, param) {
  X <- cbind(1, data$x1, data$x2)
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  
  eta <- plogis(drop(X %*% xis))
  a_i <- b * log(1 + ((1 - eta) / eta)**(1 / g))
  
  t <- data$times
  z <- 1 - exp(-b * t)
  
  Num <- exp(-(g * a_i / b) * z)
  Den <- (1 - exp(-(a_i / b) * z))**g + exp(-(g * a_i / b) * z)
  
  surv <- Num / Den
  return(surv)
}

## Hazard Function
haz.func <- function(data, param) {
  haz <- den.func(data, param) / surv.func(data, param)
  return(haz)
}

## Função Geradora de Dados (DGP)
simudata <- function(n, param, cp) {
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  
  eta <- plogis(drop(X %*% xis))
  a_i <- b * log(1 + ((1 - eta) / eta)**(1 / g))
  
  cured <- rbinom(n, 1, eta)
  SurvT <- rep(Inf, n)
  
  ind <- which(cured == 0)
  
  if (length(ind) > 0) {
    U.susc <- eta[ind] + runif(length(ind)) * (1 - eta[ind])
    C.star <- ((1 - U.susc) / U.susc)^(1 / g)
    SurvT[ind] <- -(1 / b) * log(1 - (b / a_i[ind]) * log(1 + C.star))
  }
  
  CensT <- rexp(n, rate = cp)
  times <- pmin(SurvT, CensT)
  status <- as.numeric(SurvT <= CensT)
  
  data.frame(times = times, status = status, x1 = x1, x2 = x2, y = cured)
}

## Log-likelihood Function
LogLik <- function(data, param) {
  dj <- data$status
  if ((param[4] > 0) & (param[5] > 0)) {
    d_vals <- pmax(den.func(data, param), 1e-15)
    s_vals <- pmax(surv.func(data, param), 1e-15)
    return((-1) * sum(dj * log(d_vals) + (1 - dj) * log(s_vals)))
  } else {
    return((-1) * Inf)
  }
}

## Convergence Check
CheckCov <- function(fit, LogLik = NULL, data = NULL, tol_eigen = 1e-6) {
  optimizer_ok <- 0
  CheckHessian_ok <- 0
  min_eigenvalue <- NA_real_
  
  if (is.null(fit) || inherits(fit, "try-error") || is.null(fit$convergence)) {
    return(list(converged = FALSE, optimizer_ok = optimizer_ok, CheckHessian_ok = CheckHessian_ok, min_eigenvalue = min_eigenvalue))
  }
  
  optimizer_ok <- as.numeric(length(fit$convergence) == 1L && is.finite(fit$convergence) && fit$convergence == 0)
  
  if (!is.null(fit$hessian) && is.matrix(fit$hessian) && nrow(fit$hessian) == ncol(fit$hessian) && all(is.finite(fit$hessian))) {
    H <- (fit$hessian + t(fit$hessian)) / 2
    eig <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
    
    if (!is.null(eig) && all(is.finite(eig))) {
      min_eigenvalue <- min(eig)
      CheckHessian_ok <- as.numeric(min_eigenvalue > tol_eigen)
    }
  }
  
  converged <- (optimizer_ok == 1 && CheckHessian_ok == 1)
  return(list(converged = converged, optimizer_ok = optimizer_ok, CheckHessian_ok = CheckHessian_ok, min_eigenvalue = min_eigenvalue))
}

## Monte Carlo Simulation Function (Apenas xi0, xi1, xi2 na Tabela)
MCOLLGD <- function(param, ns, cp, mc.samples, max_attempts = 500) {
  p <- length(param)
  parameter_names <- c("xi0", "xi1", "xi2", "b", "g")
  
  if (p != 5) {
    stop("Este código pressupõe exatamente cinco parâmetros.")
  }
  
  # Matriz formatada com 23 colunas contendo apenas xi0, xi1 e xi2
  Table <- matrix(NA_real_, nrow = length(ns), ncol = 23)
  Table[, 1] <- ns
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
  
  ParhatList <- vector("list", length(ns))
  BiashatList <- vector("list", length(ns))
  CoverageList <- vector("list", length(ns))
  ConvList <- vector("list", length(ns))
  
  for (r in seq_along(ns)) {
    n <- ns[r]
    
    parhats <- matrix(NA_real_, mc.samples, p)
    stdhats <- matrix(NA_real_, mc.samples, p)
    RMSEhats <- matrix(NA_real_, mc.samples, p)
    Biashats <- matrix(NA_real_, mc.samples, p)
    
    Lower.CI <- matrix(NA_real_, mc.samples, p)
    Upper.CI <- matrix(NA_real_, mc.samples, p)
    CoverageProb <- matrix(NA_real_, mc.samples, p)
    
    ConvSample <- numeric(mc.samples)
    censor_sample <- numeric(mc.samples)
    cure_sample <- numeric(mc.samples)
    err.alg <- 0
    
    for (j in seq_len(mc.samples)) {
      valid_fit <- FALSE
      myfit <- NULL
      mystd <- rep(NA_real_, p)
      data <- NULL
      
      for (attempt in seq_len(max_attempts)) {
        data <- simudata(n, param, cp)
        
        fit_try <- tryCatch(
          optim(
            par = param,
            fn = LogLik,
            data = data,
            hessian = TRUE,
            method = ifelse(attempt == 1L, "BFGS", "L-BFGS-B"),
            lower = c(rep(-Inf, 3), 1e-5, 1e-5)
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit_try) || is.null(fit_try$par) || is.null(fit_try$hessian) || 
            length(fit_try$par) != p || !all(is.finite(fit_try$par)) || any(!is.finite(fit_try$hessian))) {
          err.alg <- err.alg + 1
          next
        }
        
        H <- (fit_try$hessian + t(fit_try$hessian)) / 2
        mystd_try <- tryCatch({
          V <- solve(H)
          variances <- diag(V)
          if (any(!is.finite(variances)) || any(variances <= 0)) stop("Variância inválida")
          sqrt(variances)
        }, error = function(e) NULL)
        
        if (is.null(mystd_try) || length(mystd_try) != p || any(!is.finite(mystd_try))) {
          err.alg <- err.alg + 1
          next
        }
        
        myfit <- fit_try
        mystd <- mystd_try
        valid_fit <- TRUE
        break
      }
      
      if (!valid_fit) {
        ConvSample[j] <- 0
        next
      }
      
      check <- CheckCov(fit = myfit, LogLik = LogLik, data = data)
      ConvSample[j] <- as.numeric(isTRUE(check$converged))
      
      parhats[j, ] <- myfit$par
      stdhats[j, ] <- mystd
      
      Lower.CI[j, ] <- myfit$par - qnorm(0.975) * mystd
      Upper.CI[j, ] <- myfit$par + qnorm(0.975) * mystd
      
      CoverageProb[j, ] <- as.numeric(Lower.CI[j, ] <= param & param <= Upper.CI[j, ])
      Biashats[j, ] <- (myfit$par - param) / abs(param)
      RMSEhats[j, ] <- (myfit$par - param)^2
      
      censor_sample[j] <- 100 * mean(1 - data$status)
      cure_sample[j] <- 100 * mean(data$y)
      
      cat("Bloco", r, "| Iteração", j, "\r")
    }
    
    colnames(parhats) <- parameter_names
    colnames(stdhats) <- parameter_names
    colnames(Biashats) <- parameter_names
    colnames(RMSEhats) <- parameter_names
    colnames(CoverageProb) <- parameter_names
    
    # Salva apenas a submatriz dos parâmetros xi (colunas 1 a 3) nas listas de saída
    ParhatList[[r]]   <- parhats[, 1:3]
    BiashatList[[r]]  <- Biashats[, 1:3]
    CoverageList[[r]] <- CoverageProb[, 1:3]
    ConvList[[r]]     <- ConvSample
    
    col_mean <- function(x) {
      vapply(seq_len(ncol(x)), function(k) {
        z <- x[, k]
        if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
      }, numeric(1))
    }
    
    col_sd <- function(x) {
      vapply(seq_len(ncol(x)), function(k) {
        z <- x[, k]
        z <- z[is.finite(z)]
        if (length(z) < 2L) NA_real_ else sd(z)
      }, numeric(1))
    }
    
    Table[r, 2] <- if (all(is.na(censor_sample))) NA_real_ else round(mean(censor_sample, na.rm = TRUE), 3)
    Table[r, 3] <- if (all(is.na(cure_sample))) NA_real_ else round(mean(cure_sample, na.rm = TRUE), 3)
    
    # Preenchimento exclusivo dos resultados referentes a xi0, xi1 e xi2 (índices 1 a 3)
    Table[r, 4:6]   <- round(col_mean(parhats[, 1:3]), 3)
    Table[r, 7:9]   <- round(col_mean(Biashats[, 1:3]), 3)
    Table[r, 10:12] <- round(col_sd(parhats[, 1:3]), 3)
    Table[r, 13:15] <- round(col_mean(stdhats[, 1:3]), 3)
    Table[r, 16:18] <- round(sqrt(col_mean(RMSEhats[, 1:3])), 3)
    Table[r, 19:21] <- round(col_mean(CoverageProb[, 1:3]), 3)
    Table[r, 22]    <- round(mean(ConvSample), 3)
    Table[r, 23]    <- err.alg
    cat("\n")
  }
  
  return(list(
    Table = Table,
    CoverageProb = CoverageList,
    BiashatList = BiashatList,
    ParhatList = ParhatList,
    ConvSample = ConvList
  ))
}

## Execução de Exemplo
set.seed(123)
ssizes <- c(100, 300, 500, 800, 1000)
param_sc1 <- c(-1.5, -1.5, 0.3, 0.5, 1.0) # xis (3), b (1), g (1)
cp_val <- 0.71

Res_MGo_DGP <- MCOLLGD(param_sc1, ssizes, cp_val, mc.samples = 1000)
print(Res_MGo_DGP$Table)

# Salvamento automático
dir_path <- "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\Misspecification"
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
save(Res_MGo_DGP, file = file.path(dir_path, "Res-MGo-Reg-OLLMGo.RData"))