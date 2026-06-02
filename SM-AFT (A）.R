# Load required packages
library("survival")
library("flexsurv")
library("numDeriv")

# A.1 LR_int function
LR_int = function(y1, len1, l1){ 
  if(y1 > 0 & y1 <= l1){
    a = c(.Machine$double.eps, l1)
  } else {
    k = as.integer((y1 - l1)/len1) + 1
    a = c(l1 + ((k-1)*len1), l1 + (k*len1))
  }
  return(a)
}

# Variance computation using Hessian (robust)
compute_var_se <- function(negLogLik, params, param_names, 
                           eps = 1e-6,              
                           method.args = list(d = 0.001, r = 4)) {  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    install.packages("numDeriv")
  }
  library(numDeriv)
  hess <- try(hessian(negLogLik, params, method.args = method.args), silent = TRUE)
  if (inherits(hess, "try-error")) {
    method.args$d <- 0.01   
    hess <- try(hessian(negLogLik, params, method.args = method.args), silent = TRUE)
  }
  if (inherits(hess, "try-error")) {
    se <- rep(999, length(params))
    names(se) <- param_names
    return(list(se = se, hessian = NULL))
  }
  hess_reg <- hess + diag(eps, nrow(hess))
  inv_hess <- try(MASS::ginv(hess_reg), silent = TRUE)
  if (inherits(inv_hess, "try-error")) {
    hess_reg2 <- hess + diag(eps * 10, nrow(hess))
    inv_hess <- try(MASS::ginv(hess_reg2), silent = TRUE)
    if (inherits(inv_hess, "try-error")) {
      se <- rep(999, length(params))
      names(se) <- param_names
      return(list(se = se, hessian = hess))
    }
  }
  
  var_diag <- pmax(diag(inv_hess), 1e-8)
  se <- sqrt(var_diag)
  names(se) <- param_names
  
  return(list(se = se, hessian = hess))
}

###### Data generation #####################################
N <- 1000
gamma1 <- c(1,2)
gamma2 <- c(1.5,1.5)
xx <- cbind(rnorm(N, 0, 1), rbinom(N, 1, 0.5))
zz <- cbind(rnorm(N, 0, 2), rbinom(N, 1, 0.5))
a1 <- 2
beta1 <- c(1, 1.5)
b1 <- 1
a2 <- 1.5
beta2 <- c(1,1)
b2 <- 2
cens <- 10

data.int <- function(N, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, cens, zz, xx,
                     baseline = "Weibull",
                     # baseline = "lognormal",
                     # baseline = "loglogistic",
                     link = "logit"){
  
  LR_int = function(y1, len1, l1){ 
    if(y1 > 0 & y1 <= l1){
      a = c(.Machine$double.eps, l1)
    } else {
      k = as.integer((y1 - l1)/len1) + 1
      a = c(l1 + ((k-1)*len1), l1 + (k*len1))
    }
    return(a)
  } 
  
  if (link == "logit") {
    temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
    temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
    temp0 <- 1 - temp1 - temp2
  } 
  
  u <- sapply(1:N, function(i) {
    sample(c(0,1,2), size = 1, prob = c(temp0[i], temp1[i], temp2[i]))
  })
  
  xbeta1 <- xx %*% beta1
  xbeta2 <- xx %*% beta2
  
  if (baseline == "Weibull") {
    t.true1 <- rweibull(N, shape = 1/b1, scale = exp(xbeta1 + a1))
    t.true2 <- rweibull(N, shape = 1/b2, scale = exp(xbeta2 + a2))
  }
  # For lognormal or loglogistic, uncomment if needed
  # if (baseline == "lognormal") { ... }
  # if (baseline == "loglogistic") { ... }
  
  C <- runif(N, 0, cens)
  L <- rep(NA, N)  
  R <- rep(NA, N) 
  Status <- rep(NA, N) 
  
  for (i in 1:N) {
    if (u[i] == 0) { 
      L[i] <- C[i]
      R[i] <- Inf
      Status[i] <- 0
    } else if (u[i] == 1) { 
      if (t.true1[i] >= C[i] + 40) {  
        L[i] <- C[i] + 40
        R[i] <- Inf
        Status[i] <- 0
      } else {  
        len <- runif(1, 0.2, 0.7)
        l <- runif(1, 0, 20)
        ans <- LR_int(t.true1[i], len, l)
        L[i] <- ans[1]
        R[i] <- ans[2]
        Status[i] <- 1
      }
    } else if (u[i] == 2) {  
      if (t.true2[i] >= C[i] + 40) {  
        L[i] <- C[i] + 40
        R[i] <- Inf
        Status[i] <- 0
      } else {  
        len <- runif(1, 0.2, 0.7)
        l <- runif(1, 0, 20)
        ans <- LR_int(t.true2[i], len, l)
        L[i] <- ans[1]
        R[i] <- ans[2]
        Status[i] <- 1
      }
    }
  }
  cure <- sum(u == 0)
  uncure1 <- sum(u == 1)
  uncure2 <- sum(u == 2)
  sum(Status)
  censorrate <- 1 - sum(Status) / N
  mydata <- data.frame(L = L, R = R, Status = Status, t.true1 = t.true1, t.true2 = t.true2,
                       u = u, xx = xx, zz = zz, cure = cure, uncure1 = uncure1, uncure2 = uncure2, censorrate = censorrate)
  return(list(mydata = mydata, xx = xx, zz = zz, cure = cure, uncure1 = uncure1, uncure2 = uncure2, censorrate = censorrate))
}

data <- data.int(N, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, cens, zz, xx,
                 baseline = "Weibull",
                 # baseline = "lognormal",
                 # baseline = "loglogistic",
                 link = "logit")

# Find censoring time to achieve ~20% censoring rate
c22 <- seq(1, 1e+04, by = 100)
c11 <- 0
j <- 0
xx <- cbind(rnorm(N, 0, 1), rbinom(N, 1, 0.5))
zz <- cbind(rnorm(N, 0, 1), rbinom(N, 1, 0.5))
rate <- numeric()
for(i in 1:length(c22)){
  result <- data.int(N, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, c22[i], zz, xx,
                     baseline = "Weibull",
                     # baseline = "lognormal",
                     # baseline = "loglogistic",
                     link = "logit")
  rate[i] <- result$censorrate
  if(rate[i] > 0.175 & rate[i] < 0.225){
    c11 <- c11 + c22[i]
    j <- j + 1
  }
}
cens <- c11 / j
data <- data.int(N, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, cens, zz, xx,
                 baseline = "Weibull",
                 # baseline = "lognormal",
                 # baseline = "loglogistic",
                 link = "logit")

### EM algorithm #############################
em_ic <- function(data) {
  mydata <- data$mydata
  xx <- data$xx
  zz <- data$zz
  L = mydata$L
  R = mydata$R
  t.true1 = mydata$t.true1
  t.true2 = mydata$t.true2
  Status = mydata$Status
  u = mydata$u
  gamma1hat = c(0.1,0.1)
  gamma2hat = c(0.5,0.5)
  a1hat = 0.1
  beta1hat = c(0.1,0.5)
  b1hat = 0.1
  a2hat = 0.5
  beta2hat = c(0.1,0.1)
  b2hat = 0.1
  convergence = 1000; iter = 1; emmax = 100; eps = 1e-3
  
  while (convergence >= eps & iter <= emmax) {
    # E-step
    em_temp1 <- exp(cbind(zz) %*% gamma1hat) / (exp(cbind(zz) %*% gamma1hat) + exp(cbind(zz) %*% gamma2hat) + 1)
    em_temp2 <- exp(cbind(zz) %*% gamma2hat) / (exp(cbind(zz) %*% gamma1hat) + exp(cbind(zz) %*% gamma2hat) + 1)
    em_temp0 <- 1 - em_temp1 - em_temp2
    
    xbeta1hat <- xx %*% beta1hat
    xbeta2hat <- xx %*% beta2hat
    sur1 <- 1 - pweibull(L, shape = 1/b1hat, scale = exp(xbeta1hat + a1hat))
    sur2 <- 1 - pweibull(L, shape = 1/b2hat, scale = exp(xbeta2hat + a2hat))
    # For lognormal or loglogistic, uncomment as needed
    # sur1 <- 1 - plnorm(L, meanlog = xbeta1hat + a1hat, sdlog = b1hat)
    # sur2 <- 1 - plnorm(L, meanlog = xbeta2hat + a2hat, sdlog = b2hat)
    # sur1 <- 1 - pllogis(L, shape = 1/b1hat, scale = exp(xbeta1hat + a1hat))
    # sur2 <- 1 - pllogis(L, shape = 1/b2hat, scale = exp(xbeta2hat + a2hat))
    sur0 <- em_temp0 + sur1 * em_temp1 + sur2 * em_temp2
    
    p10 <- Status * as.numeric(u == 0) + (1 - Status) * (em_temp0 / sur0)
    p11 <- Status * as.numeric(u == 1) + (1 - Status) * (em_temp1 * sur1 / sur0)
    p12 <- Status * as.numeric(u == 2) + (1 - Status) * (em_temp2 * sur2 / sur0)
    
    # Negative log-likelihood for gamma (cure probabilities)
    lc1 <- function(gamma1, gamma2, a1, b1, beta1) {
      temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
      temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
      temp0 <- 1 - temp1 - temp2
      loglikelihood <- sum(as.numeric(u == 1) * log(temp1) + 
                             as.numeric(u == 2) * log(temp2) + 
                             (1 - Status) * as.numeric(u == 0) * log(pmax(temp0, 1e-10)), na.rm = TRUE)
      return(-loglikelihood)
    }
    
    # Negative log-likelihood for event type 1 (default)
    lc2 <- function(a1, b1, beta1) {
      xbeta1 <- xx %*% beta1
      L <- ifelse(L == 0, 1e-3, L)
      R <- ifelse(R == 0, 1e-3, R)
      S1L <- 1 - pweibull(L, shape = 1/b1, scale = exp(xbeta1 + a1))
      S1R <- 1 - pweibull(R, shape = 1/b1, scale = exp(xbeta1 + a1))
      # For lognormal or loglogistic, uncomment as needed
      # S1L <- 1 - plnorm(L, meanlog = xbeta1 + a1, sdlog = b1)
      # S1R <- 1 - plnorm(R, meanlog = xbeta1 + a1, sdlog = b1)
      # S1L <- 1 - pllogis(L, shape = 1/b1, scale = exp(xbeta1 + a1))
      # S1R <- 1 - pllogis(R, shape = 1/b1, scale = exp(xbeta1 + a1))
      
      S1L <- ifelse(is.infinite(S1L), 1e+10, S1L)
      S1R <- ifelse(is.infinite(S1R), 1e+10, S1R)  
      dan <- S1L - S1R
      dan[dan <= 0] <- 1e-10  
      log_dan <- log(dan)
      loglikelihood <- sum(Status * as.numeric(u == 1) * log_dan + 
                             (1 - Status) * as.numeric(u == 1) * log(pmax(S1L, 1e-10))) 
      return(-loglikelihood)
    }
    
    # Negative log-likelihood for event type 2 (prepayment)
    lc3 <- function(a2, b2, beta2) {
      xbeta2 <- xx %*% beta2
      L <- ifelse(L == 0, 1e-3, L)
      R <- ifelse(R == 0, 1e-3, R)
      S2L <- 1 - pweibull(L, shape = 1/b2, scale = exp(xbeta2 + a2))
      S2R <- 1 - pweibull(R, shape = 1/b2, scale = exp(xbeta2 + a2))
      # For lognormal or loglogistic, uncomment as needed
      # S2L <- 1 - plnorm(L, meanlog = xbeta2 + a2, sdlog = b2)
      # S2R <- 1 - plnorm(R, meanlog = xbeta2 + a2, sdlog = b2)
      # S2L <- 1 - pllogis(L, shape = 1/b2, scale = exp(xbeta2 + a2))
      # S2R <- 1 - pllogis(R, shape = 1/b2, scale = exp(xbeta2 + a2))
      
      S2L <- ifelse(is.infinite(S2L), 1e+10, S2L)
      S2R <- ifelse(is.infinite(S2R), 1e+10, S2R) 
      ddf <- S2L - S2R
      ddf[ddf <= 0] <- 1e-10  
      log_ddf <- log(ddf)
      loglikelihood <- sum(Status * as.numeric(u == 2) * log_ddf + 
                             (1 - Status) * as.numeric(u == 2) * log(pmax(S2L, 1e-10))) 
      return(-loglikelihood)
    }
    
    # M-step: group-wise optimization
    # Group 1: gamma1, gamma2
    gamma_result <- optim(par = c(gamma1hat, gamma2hat),
                          fn = function(params) {
                            p <- length(gamma1hat)
                            gamma1 <- params[1:p]
                            gamma2 <- params[(p+1):(2*p)]
                            lc1(gamma1, gamma2, a1 = a1hat, b1 = b1hat, beta1 = beta1hat)
                          },
                          method = "Nelder-Mead",
                          control = list(reltol = 1e-4, maxit = 500))
    
    p_gamma <- length(gamma1hat)
    update_gamma1 <- gamma_result$par[1:p_gamma]
    update_gamma2 <- gamma_result$par[(p_gamma+1):(2*p_gamma)]
    
    # Group 2: a1, beta1, b1
    params1_result <- optim(par = c(a1hat, beta1hat, b1hat),
                            fn = function(params) {
                              p_beta1 <- length(beta1hat)
                              a1 <- params[1]
                              beta1 <- params[2:(1+p_beta1)]
                              b1 <- params[2+p_beta1]
                              lc2(a1 = a1, b1 = b1, beta1 = beta1)
                            },
                            method = "L-BFGS-B",
                            lower = c(-Inf, rep(-Inf, length(beta1hat)), 0.1),
                            upper = c(Inf, rep(Inf, length(beta1hat)), 10),
                            control = list(factr = 1e4, pgtol = 1e-4, maxit = 500))
    
    update_a1 <- params1_result$par[1]
    update_beta1 <- params1_result$par[2:(1+length(beta1hat))]
    update_b1 <- params1_result$par[2+length(beta1hat)]
    
    # Group 3: a2, beta2, b2
    params2_result <- optim(par = c(a2hat, beta2hat, b2hat),
                            fn = function(params) {
                              p_beta2 <- length(beta2hat)
                              a2 <- params[1]
                              beta2 <- params[2:(1+p_beta2)]
                              b2 <- params[2+p_beta2]
                              lc3(a2 = a2, b2 = b2, beta2 = beta2)
                            },
                            method = "L-BFGS-B",
                            lower = c(-Inf, rep(-Inf, length(beta2hat)), -Inf),
                            upper = c(Inf, rep(Inf, length(beta2hat)), Inf),
                            control = list(factr = 1e4, pgtol = 1e-4, maxit = 500))
    
    update_a2 <- params2_result$par[1]
    update_beta2 <- params2_result$par[2:(1+length(beta2hat))]
    update_b2 <- params2_result$par[2+length(beta2hat)]
    
    convergence <- sum(abs(update_gamma1 - gamma1hat),
                       abs(update_gamma2 - gamma2hat),
                       abs(update_a1 - a1hat),
                       abs(update_beta1 - beta1hat),
                       abs(update_b1 - b1hat),
                       abs(update_a2 - a2hat),
                       abs(update_beta2 - beta2hat),
                       abs(update_b2 - b2hat))
    
    gamma1hat <- update_gamma1
    gamma2hat <- update_gamma2
    a1hat <- update_a1
    beta1hat <- update_beta1
    b1hat <- update_b1
    a2hat <- update_a2
    beta2hat <- update_beta2
    b2hat <- update_b2
    
    iter <- iter + 1
  }
  
  variance_results <- calculate_variance_em(gamma1hat, gamma2hat, a1hat, beta1hat, b1hat, 
                                            a2hat, beta2hat, b2hat, data)
  list(
    gamma1hat = gamma1hat,
    gamma2hat = gamma2hat,
    a1hat = a1hat,
    beta1hat = beta1hat,
    b1hat = b1hat,
    a2hat = a2hat,
    beta2hat = beta2hat,
    b2hat = b2hat,
    variance = variance_results,
    convergence = convergence,
    iterations = iter
  )
}

calculate_variance_em <- function(gamma1hat, gamma2hat, a1hat, beta1hat, b1hat, 
                                  a2hat, beta2hat, b2hat, data) {
  
  mydata <- data$mydata
  xx <- data$xx
  zz <- data$zz
  L = mydata$L
  R = mydata$R
  Status = mydata$Status
  u = mydata$u
  
  # Complete data log-likelihood for Hessian
  complete_loglik <- function(all_params) {
    p_gamma <- length(gamma1hat)
    p_beta1 <- length(beta1hat)
    p_beta2 <- length(beta2hat)
    
    # Unpack parameters
    gamma1 <- all_params[1:p_gamma]
    gamma2 <- all_params[(p_gamma+1):(2*p_gamma)]
    a1 <- all_params[2*p_gamma + 1]
    beta1 <- all_params[(2*p_gamma+2):(2*p_gamma+1+p_beta1)]
    b1 <- all_params[2*p_gamma+1+p_beta1+1]
    a2 <- all_params[2*p_gamma+1+p_beta1+2]
    beta2 <- all_params[(2*p_gamma+1+p_beta1+3):(2*p_gamma+1+p_beta1+2+p_beta2)]
    b2 <- all_params[2*p_gamma+1+p_beta1+2+p_beta2+1]
    
    ll <- 0
    
    # Multinomial part
    temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
    temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
    temp0 <- 1 - temp1 - temp2
    
    ll <- ll + sum(as.numeric(u == 1) * log(temp1) + 
                     as.numeric(u == 2) * log(temp2) + 
                     as.numeric(u == 0) * log(pmax(temp0, 1e-10)), na.rm = TRUE)
    
    # Event type 1 (default) part
    xbeta1 <- xx %*% beta1
    L_temp <- ifelse(L == 0, 1e-3, L)
    R_temp <- ifelse(R == 0, 1e-3, R)
    S1L <- 1 - pweibull(L, shape = 1/b1, scale = exp(xbeta1 + a1))
    S1R <- 1 - pweibull(R, shape = 1/b1, scale = exp(xbeta1 + a1))
    # For lognormal or loglogistic, uncomment as needed
    # S1L <- 1 - plnorm(L, meanlog = xbeta1 + a1, sdlog = b1)
    # S1R <- 1 - plnorm(R, meanlog = xbeta1 + a1, sdlog = b1)
    # S1L <- 1 - pllogis(L_temp, shape = 1/b1, scale = exp(xbeta1 + a1))
    # S1R <- 1 - pllogis(R_temp, shape = 1/b1, scale = exp(xbeta1 + a1))
    
    S1L <- ifelse(is.infinite(S1L), 1e+10, S1L)
    S1R <- ifelse(is.infinite(S1R), 1e+10, S1R)  
    dan <- S1L - S1R
    dan[dan <= 0] <- 1e-10
    
    ll <- ll + sum(Status * as.numeric(u == 1) * log(dan) + 
                     (1 - Status) * as.numeric(u == 1) * log(pmax(S1L, 1e-10)))
    
    # Event type 2 (prepayment) part
    xbeta2 <- xx %*% beta2
    S2L <- 1 - pweibull(L, shape = 1/b2, scale = exp(xbeta2 + a2))
    S2R <- 1 - pweibull(R, shape = 1/b2, scale = exp(xbeta2 + a2))
    # For lognormal or loglogistic, uncomment as needed
    # S2L <- 1 - plnorm(L, meanlog = xbeta2 + a2, sdlog = b2)
    # S2R <- 1 - plnorm(R, meanlog = xbeta2 + a2, sdlog = b2)
    # S2L <- 1 - pllogis(L_temp, shape = 1/b2, scale = exp(xbeta2 + a2))
    # S2R <- 1 - pllogis(R_temp, shape = 1/b2, scale = exp(xbeta2 + a2))
    
    S2L <- ifelse(is.infinite(S2L), 1e+10, S2L)
    S2R <- ifelse(is.infinite(S2R), 1e+10, S2R) 
    ddf <- S2L - S2R
    ddf[ddf <= 0] <- 1e-10
    
    ll <- ll + sum(Status * as.numeric(u == 2) * log(ddf) + 
                     (1 - Status) * as.numeric(u == 2) * log(pmax(S2L, 1e-10)))
    
    return(-ll)  
  }
  
  all_params <- c(gamma1hat, gamma2hat, a1hat, beta1hat, b1hat, a2hat, beta2hat, b2hat)
  tryCatch({
    hessian_matrix <- hessian(complete_loglik, all_params)
    fisher_info <- hessian_matrix
    cov_matrix <- solve(fisher_info)
    p_gamma <- length(gamma1hat)
    p_beta1 <- length(beta1hat)
    p_beta2 <- length(beta2hat)
    
    idx <- 1
    gamma1_se <- sqrt(diag(cov_matrix)[idx:(idx+p_gamma-1)]); idx <- idx + p_gamma
    gamma2_se <- sqrt(diag(cov_matrix)[idx:(idx+p_gamma-1)]); idx <- idx + p_gamma
    a1_se <- sqrt(cov_matrix[idx, idx]); idx <- idx + 1
    beta1_se <- sqrt(diag(cov_matrix)[idx:(idx+p_beta1-1)]); idx <- idx + p_beta1
    b1_se <- sqrt(cov_matrix[idx, idx]); idx <- idx + 1
    a2_se <- sqrt(cov_matrix[idx, idx]); idx <- idx + 1
    beta2_se <- sqrt(diag(cov_matrix)[idx:(idx+p_beta2-1)]); idx <- idx + p_beta2
    b2_se <- sqrt(cov_matrix[idx, idx])
    
    return(list(
      gamma1_se = gamma1_se,
      gamma2_se = gamma2_se,
      a1_se = a1_se,
      beta1_se = beta1_se,
      b1_se = b1_se,
      a2_se = a2_se,
      beta2_se = beta2_se,
      b2_se = b2_se,
      cov_matrix = cov_matrix
    ))
    
  }, error = function(e) {
    warning("Unable to compute variance matrix: ", e$message)
    return(list(
      gamma1_se = rep(NA, length(gamma1hat)),
      gamma2_se = rep(NA, length(gamma2hat)),
      a1_se = NA,
      beta1_se = rep(NA, length(beta1hat)),
      b1_se = NA,
      a2_se = NA,
      beta2_se = rep(NA, length(beta2hat)),
      b2_se = NA,
      cov_matrix = NULL
    ))
  })
}

### Simulation #############################
sim <- 200 
truevalue <- c(gamma1, gamma2, a1, beta1, b1, a2, beta2, b2)
tvl <- length(truevalue)

# Storage for estimates
simgamma1 <- simgamma2 <- matrix(0, nrow = sim, ncol = 2)
sima1 <- sima2 <- matrix(0, nrow = sim, ncol = 1)
simbeta1 <- simbeta2 <- matrix(0, nrow = sim, ncol = 2)
simb1 <- simb2 <- matrix(0, nrow = sim, ncol = 1)

# Storage for standard errors
se_gamma1 <- se_gamma2 <- matrix(0, nrow = sim, ncol = 2)
se_a1 <- se_a2 <- matrix(0, nrow = sim, ncol = 1)
se_beta1 <- se_beta2 <- matrix(0, nrow = sim, ncol = 2)
se_b1 <- se_b2 <- matrix(0, nrow = sim, ncol = 1)

cp_count <- rep(0, tvl)
start_time <- Sys.time()

for (w in 1:sim) {
  tryCatch({
    new_xx <- cbind(rnorm(N, 0, 1), rbinom(N, 1, 0.5))
    new_zz <- cbind(rnorm(N, 0, 1), rbinom(N, 1, 0.5))
    
    data <- data.int(N, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, cens,
                     new_zz, new_xx,
                     baseline = "Weibull",
                     # baseline = "lognormal",
                     # baseline = "loglogistic",
                     link = "logit")
    
    results <- em_ic(data)
    
    simgamma1[w,] <- results$gamma1hat
    simgamma2[w,] <- results$gamma2hat
    sima1[w] <- results$a1hat
    simbeta1[w,] <- results$beta1hat
    simb1[w] <- results$b1hat
    sima2[w] <- results$a2hat
    simbeta2[w,] <- results$beta2hat
    simb2[w] <- results$b2hat
    
    if (!is.null(results$variance)) {
      se_gamma1[w,] <- results$variance$gamma1_se
      se_gamma2[w,] <- results$variance$gamma2_se
      se_a1[w] <- results$variance$a1_se
      se_beta1[w,] <- results$variance$beta1_se
      se_b1[w] <- results$variance$b1_se
      se_a2[w] <- results$variance$a2_se
      se_beta2[w,] <- results$variance$beta2_se
      se_b2[w] <- results$variance$b2_se
      
      all_estimates <- c(results$gamma1hat, results$gamma2hat, results$a1hat, 
                         results$beta1hat, results$b1hat, results$a2hat, 
                         results$beta2hat, results$b2hat)
      all_se <- c(results$variance$gamma1_se, results$variance$gamma2_se, 
                  results$variance$a1_se, results$variance$beta1_se, 
                  results$variance$b1_se, results$variance$a2_se, 
                  results$variance$beta2_se, results$variance$b2_se)
      lower <- all_estimates - 1.96 * all_se
      upper <- all_estimates + 1.96 * all_se
      for (j in 1:tvl) {
        if (truevalue[j] >= lower[j] && truevalue[j] <= upper[j]) {
          cp_count[j] <- cp_count[j] + 1
        }
      }
    }
    
    print(paste("sim", w, "completed"))
    
  }, error = function(e) {
    print(paste("Error in simulation", w, ":", e$message))
  })
}

end_time <- Sys.time()
run_time <- end_time - start_time
print(paste("Total run time:", run_time))

### Summary statistics #############################
simdata <- cbind(simgamma1, simgamma2, sima1, simbeta1, simb1, sima2, simbeta2, simb2)
simmean <- apply(simdata, 2, mean, na.rm = TRUE)
sd_empirical <- apply(simdata, 2, sd, na.rm = TRUE)  # Empirical SE (ESE)
ase <- c(
  colMeans(se_gamma1, na.rm = TRUE),
  colMeans(se_gamma2, na.rm = TRUE),
  mean(se_a1, na.rm = TRUE),
  colMeans(se_beta1, na.rm = TRUE),
  mean(se_b1, na.rm = TRUE),
  mean(se_a2, na.rm = TRUE),
  colMeans(se_beta2, na.rm = TRUE),
  mean(se_b2, na.rm = TRUE)
)
bias <- simmean - truevalue
mse <- bias^2 + sd_empirical^2
cp_rate <- cp_count / sim
conf_level <- 0.95
z_value <- qnorm(1 - (1 - conf_level) / 2)
simlcil <- simmean - z_value * sd_empirical
simucil <- simmean + z_value * sd_empirical

results_summary <- data.frame(
  Parameter = c(paste0("gamma1_", 1:2), paste0("gamma2_", 1:2), 
                "a1", paste0("beta1_", 1:2), "b1", 
                "a2", paste0("beta2_", 1:2), "b2"),
  TrueValue = truevalue,
  Mean = round(simmean, 4),
  Bias = round(bias, 4),
  ESE = round(sd_empirical, 4),
  ASE = round(ase, 4),
  MSE = round(mse, 4),
  CP = round(cp_rate, 4),
  CI_Lower = round(simlcil, 4),
  CI_Upper = round(simucil, 4)
)

print("Simulation results summary:")
print(results_summary)

# Optional: write.csv(results_summary, "simulation_results.csv", row.names = FALSE)

cat("\nAdditional information:\n")
cat("Sample size: ", N, "\n")
cat("Number of simulations: ", sim, "\n")
cat("Total run time: ", run_time, "\n")
cat("Average time per simulation: ", as.numeric(run_time)/sim, "seconds\n")

ratio_ase_ese <- ase / sd_empirical
cat("\nASE/ESE ratio:\n")
print(round(ratio_ase_ese, 4))

convergence_check <- data.frame(
  Parameter = results_summary$Parameter,
  Bias_Relative = round(abs(bias / truevalue) * 100, 2),
  CP_Target = ifelse(cp_rate > 0.90 & cp_rate < 0.98, "Good", 
                     ifelse(cp_rate >= 0.85 & cp_rate <= 0.95, "Acceptable", "Poor"))
)

print("Convergence check:")
print(convergence_check)

# Plot survival curves
library(ggplot2)

weibull_survival <- function(t, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, xx, zz) {
  temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
  temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
  temp0 <- 1 - temp1 - temp2
  xbeta1 <- xx %*% beta1
  xbeta2 <- xx %*% beta2
  # For Weibull baseline (can change to lognormal or loglogistic as needed)
  S1 <- 1 - pllogis(t, shape = 1/b1, scale = exp(xbeta1 + a1))
  S2 <- 1 - pllogis(t, shape = 1/b2, scale = exp(xbeta2 + a2))
  surv <- temp0 + temp1 * S1 + temp2 * S2  
  return(surv)
}

xx <- cbind(0.1, 0.3)
zz <- cbind(0.1, 0.3)
time <- seq(0, 150, by = 0.1)
survival_results <- matrix(0, nrow = 200, ncol = length(time))

for (i in 1:200) {
  gamma1hat <- simdata[i, 1:2]
  gamma2hat <- simdata[i, 3:4]
  a1hat <- simdata[i, 5]
  beta1hat <- simdata[i, 6:7]
  b1hat <- simdata[i, 8]
  a2hat <- simdata[i, 9]
  beta2hat <- simdata[i, 10:11]
  b2hat <- simdata[i, 12]
  
  survival_results[i, ] <- sapply(time, function(t) weibull_survival(t, gamma1hat, gamma2hat, a1hat, beta1hat, b1hat,
                                                                     a2hat, beta2hat, b2hat, xx, zz))
}

true_surv <- sapply(time, function(t) weibull_survival(t, gamma1 = c(1,0.5), gamma2 = c(1,1), a1 = 2, a2 = 1.5, b1 = 1, b2 = 2, beta1 = c(1,1.5),
                                                       beta2 = c(1,1), xx = xx, zz = zz))
est_surv <- sapply(time, function(t) weibull_survival(t, gamma1 = simmean[1:2], gamma2 = simmean[3:4],
                                                      a1 = simmean[5], beta1 = simmean[6:7], b1 = simmean[8],
                                                      a2 = simmean[9], beta2 = simmean[10:11], b2 = simmean[12], xx = xx, zz = zz))
lcil_surv <- apply(survival_results, 2, quantile, probs = 0.025)
ucil_surv <- apply(survival_results, 2, quantile, probs = 0.975)

data <- data.frame(
  Time = rep(time, 4),
  Survival = c(true_surv, est_surv, lcil_surv, ucil_surv),
  Type = rep(c("True", "Estimate", "Lower CI", "Upper CI"), each = length(time))
)

data_clean <- na.omit(data)

ggplot(data_clean, aes(x = Time, y = Survival, color = Type, linetype = Type, linewidth = Type)) +
  geom_step() +
  scale_color_manual(values = c("True" = "blue", "Estimate" = "red", 
                                "Lower CI" = "darkgrey", "Upper CI" = "darkgrey")) +
  scale_linetype_manual(values = c("True" = "solid", "Estimate" = "dashed", 
                                   "Lower CI" = "dotted", "Upper CI" = "longdash")) +
  scale_linewidth_manual(values = c("True" = 0.7, "Estimate" = 1.2, 
                                    "Lower CI" = 0.8, "Upper CI" = 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +  
  labs(title = "",
       x = "Time",
       y = "Survival Probability") +
  theme_minimal() +
  guides(linewidth = "none") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = c(1, 1),           
    legend.justification = c(1, 1),      
    legend.background = element_blank(), 
    legend.key = element_blank(),      
    axis.title.x = element_text(face = "plain", size = 12, color = "black"),
    axis.title.y = element_text(face = "plain", size = 12, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 14),    
    legend.title = element_text(size = 16)   
  )

save_path <- "C:/Users/17001/Desktop/overall S"
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
ggsave(filename = file.path(save_path, "l.40.1000.png"),   
       plot = last_plot(), 
       width = 8, 
       height = 6, 
       dpi = 300)
