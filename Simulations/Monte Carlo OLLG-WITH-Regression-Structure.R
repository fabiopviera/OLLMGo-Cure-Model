rm(list = ls())

## Density Function
den.func <- function(data, param) {
  
  X <- cbind(1, data$x1, data$x2)
  ##
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  ##
  ##--------------------------------##
  ##        Cure fraction           ##
  ##--------------------------------##
  ##  
  eta <- plogis(drop(X %*% xis))
  ##
  ##--------------------------------##
  ##        Individual parameter a_i ##
  ##--------------------------------##
  ##
  a_i <- b*log(1+((1-eta)/eta)**(1/g))
  ##
  t <- data$times
  ##
  z <- 1 - exp(-b*t)
  A <- exp(-(a_i/b)*z)
  B <- exp(-(g*a_i/b)*z)
  ##
  Num <- g*a_i*exp(-b*t)*B*(1-A)**(g-1)
  Den <- ((1-A)**g+B)**2
  ##
  density <- Num/Den
  return(density)
}


## Survival Function
surv.func <- function(data, param) {
  
  X <- cbind(1, data$x1, data$x2)
  ## 
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  ##
  ##--------------------------------##
  ##        Cure fraction           ##
  ##--------------------------------##
  ##  
  eta <- plogis(drop(X %*% xis))
  ##
  ##--------------------------------##
  ##        Individual parameter a_i ##
  ##--------------------------------##
  ##
  a_i <- b*log(1+((1-eta)/eta)**(1/g))
  ##
  t <- data$times
  z <- 1-exp(-b*t)
  ##
  Num <- exp(-(g*a_i/b)*z)
  Den <- (1-exp(-(a_i/b)*z))**g + exp(-(g*a_i/b)*z)
  ##
  surv <- Num/Den
  return(surv)
}

## Hazard Function
haz.func <- function(data, param){
  haz <- den.func(data, param)/surv.func(data, param)
  return(haz)
}

## Simulation Function
simudata <- function(n, param, cp) {
  ##
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  ##
  xis <- param[1:3]
  b <- param[4]
  g <- param[5]
  ##
  ##--------------------------------##
  ##        Cure fraction           ##
  ##--------------------------------##
  ##  
  eta <- plogis(drop(X %*% xis))
  ##
  ##--------------------------------##
  ##        Individual parameter a_i ##
  ##--------------------------------##
  ##
  a_i <- b*log(1+((1-eta)/eta)**(1/g))
  
  ##
  ##--------------------------------##
  ##        Cure indicator          ##
  ##--------------------------------##
  cured <- rbinom(n, 1, eta)
  ##
  ##--------------------------------##
  ##        Survival times          ##
  ##--------------------------------##
  SurvT <- rep(Inf, n)
  
  ##
  ## Individuals susceptible to the event
  ##
  ind <- which(cured == 0)
  
  if (length(ind) > 0) {
    U.susc <- eta[ind] + runif(length(ind))*(1-eta[ind])
    ##
    ## Inverse transformation
    ##
    C.star <- ((1 - U.susc) / U.susc)^(1 / g)
    ##
    ## OLLMGo survival time
    ##
    SurvT[ind]<- -(1/b)*log(1-(b/a_i[ind])*log(1+C.star))
  }
  
  ##
  ##--------------------------------##
  ##        Censoring times         ##
  ##--------------------------------##
  CensT <- rexp(n, rate = cp)
  ##
  ##--------------------------------##
  ##        Observed times          ##
  ##--------------------------------##
  times <- pmin(SurvT, CensT)
  ##
  ##--------------------------------##
  ##        Event indicator         ##
  ##--------------------------------##
  status <- as.numeric(SurvT <= CensT)
  ##
  data.frame(times = times, status = status, x1 = x1, x2 = x2, y = cured)
}

## Log-likelihood Function
LogLik <- function(data, param){
  dj <- data$status
  
  if((param[4] > 0) & (param[5] > 0)){
    return((-1) * sum(dj * log(den.func(data, param)) + (1 - dj) * log(surv.func(data, param))))
  } else {
    return((-1) * Inf)
  }
}

##-------------------##
## Convergence Check ##
##-------------------##
CheckCov <- function(fit,
                     LogLik = NULL,
                     data = NULL,
                     tol_eigen = 1e-6) {
  
  optimizer_ok <- 0
  CheckHessian_ok <- 0
  min_eigenvalue <- NA_real_
  
  if (is.null(fit) ||
      inherits(fit, "try-error") ||
      is.null(fit$convergence)) {
    
    return(list(
      converged = FALSE,
      optimizer_ok = optimizer_ok,
      CheckHessian_ok = CheckHessian_ok,
      min_eigenvalue = min_eigenvalue
    ))
  }
  
  optimizer_ok <- as.numeric(
    length(fit$convergence) == 1L &&
      is.finite(fit$convergence) &&
      fit$convergence == 0
  )
  
  if (!is.null(fit$hessian) &&
      is.matrix(fit$hessian) &&
      nrow(fit$hessian) == ncol(fit$hessian) &&
      all(is.finite(fit$hessian))) {
    
    H <- (fit$hessian + t(fit$hessian)) / 2
    
    eig <- tryCatch(
      eigen(
        H,
        symmetric = TRUE,
        only.values = TRUE
      )$values,
      error = function(e) NULL
    )
    
    if (!is.null(eig) && all(is.finite(eig))) {
      min_eigenvalue <- min(eig)
      CheckHessian_ok <- as.numeric(
        min_eigenvalue > tol_eigen
      )
    }
  }
  
  converged <- (
    optimizer_ok == 1 &&
      CheckHessian_ok == 1
  )
  
  return(list(
    converged = converged,
    optimizer_ok = optimizer_ok,
    CheckHessian_ok = CheckHessian_ok,
    min_eigenvalue = min_eigenvalue
  ))
}

##-----------------------##
## MONTE CARLO SIMULATION##
##-----------------------##
MCOLLGD <- function(param, ns, cp, mc.samples = 1000, max_attempts = 100, seed = NULL) {
  
  # Inicializa a semente de forma isolada se fornecida
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  p <- length(param)
  parameter_names <- c("xi0", "xi1", "xi2", "b", "g")
  
  if (p != 5) {
    stop("Este código pressupõe exatamente cinco parâmetros.")
  }
  
  Table <- matrix(NA_real_, nrow = length(ns), ncol = 35)
  Table[, 1] <- ns
  colnames(Table) <- c(
    "n", "%cens", "%cure",
    "xi0.mc", "xi1.mc", "xi2.mc", "b.mc", "g.mc",
    "RB.xi0", "RB.xi1", "RB.xi2", "RB.b", "RB.g",
    "MCse.xi0", "MCse.xi1", "MCse.xi2", "MCse.b", "MCse.g",
    "Ase.xi0", "Ase.xi1", "Ase.xi2", "Ase.b", "Ase.g",
    "RMSE.xi0", "RMSE.xi1", "RMSE.xi2", "RMSE.b", "RMSE.g",
    "CP.xi0", "CP.xi1", "CP.xi2", "CP.b", "CP.g",
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
            method = if (attempt == 1L) "BFGS" else "L-BFGS-B"
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit_try) ||
            is.null(fit_try$par) ||
            is.null(fit_try$hessian) ||
            length(fit_try$par) != p ||
            !all(is.finite(fit_try$par)) ||
            length(fit_try$value) != 1L ||
            !is.finite(fit_try$value) ||
            any(!is.finite(fit_try$hessian))) {
          
          err.alg <- err.alg + 1
          next
        }
        
        H <- (fit_try$hessian + t(fit_try$hessian)) / 2
        
        mystd_try <- tryCatch({
          V <- solve(H)
          variances <- diag(V)
          
          if (any(!is.finite(variances)) || any(variances <= 0)) {
            stop("Variâncias inválidas.")
          }
          
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
    
    ParhatList[[r]] <- parhats
    BiashatList[[r]] <- Biashats
    CoverageList[[r]] <- CoverageProb
    ConvList[[r]] <- ConvSample
    
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
    Table[r, 4:8] <- round(col_mean(parhats), 3)
    Table[r, 9:13] <- round(col_mean(Biashats), 3)
    Table[r, 14:18] <- round(col_sd(parhats), 3)
    Table[r, 19:23] <- round(col_mean(stdhats), 3)
    Table[r, 24:28] <- round(sqrt(col_mean(RMSEhats)), 3)
    Table[r, 29:33] <- round(col_mean(CoverageProb), 3)
    Table[r, 34] <- round(mean(ConvSample), 3)
    Table[r, 35] <- err.alg
    
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

##--------------------##
## The Main Scenarios ##
##--------------------##
ssizes <- c(100, 300)

##----------------------------------##
## Scenario 1: %Cen=30 and %Cure=10 ##
##----------------------------------##
param1 <- c(-1.5, -1.5, 0.3, 2, 1.49); cp1 <- 0.71
Res1 <- MCOLLGD(param1, ssizes, cp1, seed = 101)
Res1$Table
save(Res1, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc1.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc1.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
num_rows <- length(ssizes)
mat <- matrix(0, nrow = num_rows, ncol = 6)
idx <- 1
for (i in 1:num_rows) {
  for (j in 1:5) {
    mat[i, j] <- idx
    idx <- idx + 1
  }
  mat[i, 6] <- idx
  idx <- idx + 1
}
layout(mat, widths = c(rep(1, 5), 0.25))

labels_par <- c(expression(hat(xi)[0]), expression(hat(xi)[1]), expression(hat(xi)[2]), expression(hat(beta)), expression(hat(gamma)))

for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res1$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param1[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##----------------------------------##
## Scenario 2: %Cen=30 and %Cure=20 ##
##----------------------------------##
param2 <- c(-1.5, -0.3, 0.9, 1.23, 0.71); cp2 <- 0.46
Res2 <- MCOLLGD(param2, ssizes, cp2, seed = 102)
Res2$Table
save(Res2, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc2.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc2.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res2$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param2[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 3: %Cen=50 and %Cure=10 ##
##-----------------------------------##
param3 <- c(-1.5, -1.5, -0.3, 0.97, 1.74); cp3 <- 0.97
Res3 <- MCOLLGD(param3, ssizes, cp3, seed = 103)
Res3$Table
save(Res3, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc3.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc3.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res3$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param3[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 4: %Cen=50 and %Cure=20 ##
##-----------------------------------##
param4 <- c(-0.9, -1.5, 0.3, 0.97, 1.23); cp4 <- 0.97
Res4 <- MCOLLGD(param4, ssizes, cp4, seed = 104)
Res4$Table
save(Res4, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc4.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc4.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res4$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param4[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 5: %Cen=50 and %Cure=30 ##
##-----------------------------------##
param5 <- c(-0.9, -0.3, 0.9, 0.71, 0.97); cp5 <- 0.46
Res5 <- MCOLLGD(param5, ssizes, cp5, seed = 105)
Res5$Table
save(Res5, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc5.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc5.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res5$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param5[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 6: %Cen=50 and %Cure=40 ##
##-----------------------------------##
param6 <- c(-0.9, 0.9, 0.3, 1.23, 0.71); cp6 <- 0.46
Res6 <- MCOLLGD(param6, ssizes, cp6, seed = 106)
Res6$Table
save(Res6, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc6.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc6.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res6$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param6[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 7: %Cen=60 and %Cure=30 ##
##-----------------------------------##
param7 <- c(-0.3, -1.5, 0.3, 1.74, 1.49); cp7 <- 1.49
Res7 <- MCOLLGD(param7, ssizes, cp7, seed = 107)
Res7$Table
save(Res7, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc7.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc7.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res7$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param7[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()

##-----------------------------------##
##  Scenario 8: %Cen=60 and %Cure=40 ##
##-----------------------------------##
param8 <- c(-0.3, -0.3, -0.3, 1.74, 1.74); cp8 <- 0.71
Res8 <- MCOLLGD(param8, ssizes, cp8, seed = 108)
Res8$Table
save(Res8, file = "C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Res-sc8.RData")

postscript("C:\\Users\\UnB\\OneDrive - unb.br\\Area de Trabalho\\MATERIAL_ProfUNB\\DEFECTIVE MODELS\\CÓDIGOS\\RCODE-REVISED-PAPER\\SIMURES_Alpha\\Hist_sc8.eps", width = 18, height = 12, paper = "special", horizontal = FALSE)
layout(mat, widths = c(rep(1, 5), 0.25))
for (g in 1:num_rows) {
  for (j in 1:5) {
    par(mai = c(0.9, 0.8, 0.3, 0.1)) 
    Parhat <- Res8$ParhatList[[g]][, j]
    hist(Parhat, main = "", prob = TRUE, cex.axis = 2, cex.main = 2, cex.lab = 2, font.lab = 2, col = "lightblue2", xlab = labels_par[j], border = "white", mgp = c(4.5, 1, 0)) 
    x <- seq(min(Parhat), max(Parhat), length.out = 200)
    y <- dnorm(x, mean(Parhat), sd(Parhat))
    lines(x, y, lwd = 2)
    abline(v = param8[j], lty = 2, lwd = 2)
  }
  par(mai = c(0.8, 0.1, 0.3, 0.2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, 0.1, 1, 0.9, col = "gray85", border = "gray50", lwd = 1.5)
  text(0.5, 0.5, paste("n =", ssizes[g]), cex = 2, font = 2, col = "black", srt = 90)
}
dev.off()