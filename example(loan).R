library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(themis)
library(recipes)

data <- read_csv("D:\\BaiduNetdiskDownload\\data_2017.csv")

# Keep only 36-month loans
data <- data %>%
  filter(grepl("36", term) | grepl("36 months", term))

# Remove currently active loans
data <- data %>%
  filter(!grepl("Current", loan_status, ignore.case = TRUE))

# Keep only grade E loans
data <- data %>%
  filter(grepl("E", sub_grade, ignore.case = TRUE))

head(data)

data_clean <- data 

# Remove columns with >75% missing values
missing_values_ratio <- sapply(data_clean, function(x) mean(is.na(x)))
columns_with_high_na <- names(missing_values_ratio[missing_values_ratio > 0.75])
data_clean <- data_clean %>% select(-all_of(columns_with_high_na))

# Recompute missing ratio after removal
missing_values_ratio_clean <- sapply(data_clean, function(x) mean(is.na(x)))

# Mean imputation for columns with 10-20% missing
columns_to_fill <- names(missing_values_ratio_clean[missing_values_ratio_clean > 0.10 & missing_values_ratio_clean <= 0.20])
for (col in columns_to_fill) {
  data_clean[[col]] <- ifelse(is.na(data_clean[[col]]), mean(data_clean[[col]], na.rm = TRUE), data_clean[[col]])
}

# Remove rows with missing values in columns that have <10% missing
columns_with_low_na <- names(missing_values_ratio_clean[missing_values_ratio_clean < 0.10])
data_clean <- data_clean %>% filter(!rowSums(is.na(.[, columns_with_low_na])) > 0)

# Standardize date formats
data_clean <- data_clean %>%
  mutate(last_pymnt_d = my(last_pymnt_d))

data_clean <- data_clean %>%
  mutate(issue_d = ymd(issue_d))

str(data_clean)
colSums(is.na(data_clean))

# Select variables for standardization
selected_vars <- data_clean %>%
  select(int_rate, loan_amnt, annual_inc, installment, dti)

# Convert character columns to numeric (remove %, commas, etc.)
selected_vars_numeric <- selected_vars %>%
  mutate(across(everything(), ~ {
    if(is.character(.)) {
      as.numeric(str_replace_all(., "[^0-9.-]", ""))
    } else {
      as.numeric(.)
    }
  }))

cat("Converted data types:\n")
print(sapply(selected_vars_numeric, class))

# Standardize variables (mean=0, sd=1)
standardized_vars <- selected_vars_numeric %>%
  mutate(across(everything(), ~ scale(.)[,1], .names = "std_{.col}"))

cat("Original variable summary:\n")
print(summary(selected_vars_numeric))

cat("\nStandardized variable summary:\n")
print(summary(standardized_vars))

# Add standardized variables to main dataset
data_final <- bind_cols(data_clean, standardized_vars)

cat("\nStructure of standardized variables:\n")
str(data_final %>% select(contains("std_")))

cat("\nMissing values in standardized variables:\n")
colSums(is.na(data_final %>% select(contains("std_"))))

# Build survival dataset (default vs prepayment vs censored)
survival_data <- data_final %>%
  mutate(
    issue_date      = as.Date(issue_d),
    last_pymnt_date = as.Date(last_pymnt_d),
    issue_year       = year(issue_date),
    issue_month      = month(issue_date),
    last_pymnt_year  = year(last_pymnt_date),
    last_pymnt_month = month(last_pymnt_date),
    actual_months    = (last_pymnt_year - issue_year) * 12 +
      (last_pymnt_month - issue_month),
    is_current_but_delinq = if_else(loan_status == "Current" & delinq_amnt > 0, 1, 0),
    # Event type: 0 = no event, 1 = default, 2 = prepayment
    u = case_when(
      loan_status %in% c("Default",
                         "Late (31-120 days)",
                         "Late (16-30 days)",
                         "Late (1-15 days)",
                         "In Grace Period",
                         "Charged Off")            ~ 1,
      is_current_but_delinq == 1                    ~ 1,
      loan_status == "Fully Paid" & actual_months < 36 ~ 2,
      TRUE ~ 0
    ),
    survival_time = case_when(
      u == 1 ~ actual_months,
      u == 2 ~ actual_months,
      loan_status == "Fully Paid" & actual_months >= 36 ~ 36,
      TRUE ~ pmin(actual_months, 36)
    ),
    survival_time = pmin(pmax(survival_time, 0.1), 36)
  )

# Generate final analysis set
survival_final <- survival_data %>%
  filter(!is.na(survival_time), !is.na(u))

# Event distribution after expanding defaults
event_summary <- survival_final %>%
  count(u) %>%
  mutate(
    event_type = case_when(
      u == 0 ~ "No event",
      u == 1 ~ "Default",
      u == 2 ~ "Prepayment"
    ),
    proportion = n / nrow(survival_final)
  ) %>%
  select(event_type, n, proportion)

cat("Event distribution after expanding defaults (Strategy A):\n")
print(event_summary)

# Monthly proportions of default and prepayment
monthly_proportions <- survival_final %>%
  filter(u %in% c(1, 2)) %>%
  mutate(month = round(survival_time)) %>%
  filter(month >= 1 & month <= 36) %>%
  group_by(month, u) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(event_type = ifelse(u == 1, "default", "prepay")) %>%
  mutate(proportion = count / nrow(survival_final)) %>%
  select(month, event_type, proportion)

all_months <- expand.grid(
  month = 1:36,
  event_type = c("default", "prepay")
)

monthly_proportions_complete <- all_months %>%
  left_join(monthly_proportions, by = c("month", "event_type")) %>%
  mutate(proportion = ifelse(is.na(proportion), 0, proportion))

# Build interval-censored data structure
library(survival)

z1 <- survival_final$std_loan_amnt
z2 <- survival_final$std_annual_inc
z3 <- survival_final$std_installment
z4 <- survival_final$std_int_rate
z5 <- survival_final$std_dti
zz <- cbind(z1, z2, z4)
xx <- cbind(z1, z2)

# Initial parameter values (reference)
gamma1 <- c(1, 1, 1)
gamma2 <- c(2, 2, 2)
a1 <- 2
beta1 <- c(1, 1.5)
b1 <- 1
a2 <- 1.5
beta2 <- c(1, 1)
b2 <- 2

create_interval_censored_data <- function(survival_final, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, zz, xx){
  
  N <- nrow(survival_final)
  survival_time <- survival_final$survival_time
  u <- survival_final$u
  
  L <- rep(NA, N)
  R <- rep(NA, N)
  Status <- rep(NA, N)
  
  for (i in 1:N) {
    if (u[i] == 0) {
      L[i] <- survival_time[i]
      R[i] <- Inf
      Status[i] <- 0
    } else if (u[i] == 1) {
      L[i] <- max(0, survival_time[i] - 1)
      R[i] <- survival_time[i] + 1
      Status[i] <- 1
    } else if (u[i] == 2) {
      L[i] <- max(0, survival_time[i] - 1)
      R[i] <- survival_time[i] + 1
      Status[i] <- 1
    }
  }
  
  cure <- sum(u == 0, na.rm = TRUE)
  uncure1 <- sum(u == 1, na.rm = TRUE)
  uncure2 <- sum(u == 2, na.rm = TRUE)
  censorrate <- cure / N
  
  mydata <- data.frame(
    L = L, R = R, Status = Status, u = u,
    survival_time = survival_time,
    z1 = z1, z2 = z2, z3 = z3, z4 = z4, z5 = z5
  )
  
  return(list(
    mydata = mydata, xx = xx, zz = zz,
    cure = cure, uncure1 = uncure1, uncure2 = uncure2,
    censorrate = censorrate
  ))
}

interval_censored_result <- create_interval_censored_data(
  survival_final, gamma1, gamma2, a1, beta1, b1, a2, beta2, b2, zz, xx
)

cat("Interval-censored data statistics:\n")
cat("Cured events:", interval_censored_result$cure, "\n")
cat("Default events:", interval_censored_result$uncure1, "\n")
cat("Prepayment events:", interval_censored_result$uncure2, "\n")
cat("Censoring rate:", round(interval_censored_result$censorrate * 100, 2), "%\n")

cat("\nFirst 10 rows of interval-censored data:\n")
print(head(interval_censored_result$mydata, 10))

# EM algorithm for mixture cure model (3 event types)
library("survival")
library("flexsurv")
library("numDeriv")

em_ic <- function(interval_censored_result) {
  mydata <- interval_censored_result$mydata
  xx <- interval_censored_result$xx
  zz <- interval_censored_result$zz
  L = mydata$L
  R = mydata$R
  Status = mydata$Status
  u = mydata$u
  
  # Initial values
  gamma1hat = c(1,1,1)
  gamma2hat = c(2,2,2)
  a1hat = 2
  beta1hat = c(1,1.5)
  b1hat = 1
  a2hat = 1.5
  beta2hat = c(1,1)
  b2hat = 2
  
  convergence = 1000; iter = 1; emmax = 100; eps = 1e-4
  
  while (convergence >= eps & iter <= emmax) {
    # E-step: compute posterior probabilities
    em_temp1 <- exp(cbind(zz) %*% gamma1hat) / (exp(cbind(zz) %*% gamma1hat) + exp(cbind(zz) %*% gamma2hat) + 1)
    em_temp2 <- exp(cbind(zz) %*% gamma2hat) / (exp(cbind(zz) %*% gamma1hat) + exp(cbind(zz) %*% gamma2hat) + 1)
    em_temp0 <- 1 - em_temp1 - em_temp2
    
    xbeta1hat <- xx %*% beta1hat
    xbeta2hat <- xx %*% beta2hat
    L <- ifelse(L <= 0, 1e-3, L)
    
    # Log-logistic baseline survival
    sur1 <- 1 - pllogis(L, shape = 1/b1hat, scale = exp(xbeta1hat + a1hat))
    sur2 <- 1 - pllogis(L, shape = 1/b2hat, scale = exp(xbeta2hat + a2hat))
    
    sur0 <- em_temp0 + sur1 * em_temp1 + sur2 * em_temp2
    
    p10 <- Status * as.numeric(u == 0) + (1 - Status) * (em_temp0 / sur0)
    p11 <- Status * as.numeric(u == 1) + (1 - Status) * (em_temp1 * sur1 / sur0)
    p12 <- Status * as.numeric(u == 2) + (1 - Status) * (em_temp2 * sur2 / sur0)
    
    # Negative log-likelihood for cure probabilities (gamma)
    lc1 <- function(gamma1, gamma2, a1, b1, beta1) {
      temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
      temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
      temp0 <- 1 - temp1 - temp2
      loglik <- sum(as.numeric(u == 1) * log(temp1) + 
                      as.numeric(u == 2) * log(temp2) + 
                      (1 - Status) * as.numeric(u == 0) * log(pmax(temp0, 1e-10)), na.rm = TRUE)
      return(-loglik)
    }
    
    # Negative log-likelihood for event type 1 (default)
    lc2 <- function(a1, b1, beta1) {
      xbeta1 <- xx %*% beta1
      L <- ifelse(L <= 0, 1e-3, L)
      R <- ifelse(R == 0, 1e-3, R)
      S1L <- 1 - pllogis(L, shape = 1/b1, scale = exp(xbeta1 + a1))
      S1R <- 1 - pllogis(R, shape = 1/b1, scale = exp(xbeta1 + a1))
      S1L <- ifelse(is.infinite(S1L), 1e+10, S1L)
      S1R <- ifelse(is.infinite(S1R), 1e+10, S1R)
      dan <- S1L - S1R
      dan[dan <= 0] <- 1e-10
      log_dan <- log(dan)
      loglik <- sum(Status * as.numeric(u == 1) * log_dan + 
                      (1 - Status) * as.numeric(u == 1) * log(pmax(S1L, 1e-10)))
      return(-loglik)
    }
    
    # Negative log-likelihood for event type 2 (prepayment)
    lc3 <- function(a2, b2, beta2) {
      xbeta2 <- xx %*% beta2
      L <- ifelse(L <= 0, 1e-3, L)
      R <- ifelse(R == 0, 1e-3, R)
      S2L <- 1 - pllogis(L, shape = 1/b2, scale = exp(xbeta2 + a2))
      S2R <- 1 - pllogis(R, shape = 1/b2, scale = exp(xbeta2 + a2))
      S2L <- ifelse(is.infinite(S2L), 1e+10, S2L)
      S2R <- ifelse(is.infinite(S2R), 1e+10, S2R)
      ddf <- S2L - S2R
      ddf[ddf <= 0] <- 1e-10
      log_ddf <- log(ddf)
      loglik <- sum(Status * as.numeric(u == 2) * log_ddf + 
                      (1 - Status) * as.numeric(u == 2) * log(pmax(S2L, 1e-10)))
      return(-loglik)
    }
    
    # M-step: group-wise optimization
    # Group 1: gamma1 and gamma2
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
                            lower = c(-Inf, rep(-Inf, length(beta2hat)), 0.1),
                            upper = c(Inf, rep(Inf, length(beta2hat)), 10),
                            control = list(factr = 1e4, pgtol = 1e-4, maxit = 500))
    
    update_a2 <- params2_result$par[1]
    update_beta2 <- params2_result$par[2:(1+length(beta2hat))]
    update_b2 <- params2_result$par[2+length(beta2hat)]
    
    # Convergence criterion
    convergence <- sum(abs(update_gamma1 - gamma1hat),
                       abs(update_gamma2 - gamma2hat),
                       abs(update_a1 - a1hat),
                       abs(update_beta1 - beta1hat),
                       abs(update_b1 - b1hat),
                       abs(update_a2 - a2hat),
                       abs(update_beta2 - beta2hat),
                       abs(update_b2 - b2hat))
    
    # Update parameters
    gamma1hat <- update_gamma1
    gamma2hat <- update_gamma2
    a1hat <- update_a1
    beta1hat <- update_beta1
    b1hat <- update_b1
    a2hat <- update_a2
    beta2hat <- update_beta2
    b2hat <- update_b2
    
    iter <- iter + 1
    print(iter)
    print(convergence)
  }
  
  # Compute final log-likelihood and AIC
  final_ll <- -lc1(gamma1hat, gamma2hat, a1hat, b1hat, beta1hat) -
    lc2(a1hat, b1hat, beta1hat) -
    lc3(a2hat, b2hat, beta2hat)
  
  k <- length(gamma1hat) + length(gamma2hat) + 1 + 1 + length(beta1hat) +
    1 + 1 + length(beta2hat)
  aic <- -2 * final_ll + 2 * k
  
  list(gamma1hat = gamma1hat, gamma2hat = gamma2hat,
       a1hat = a1hat, beta1hat = beta1hat, b1hat = b1hat,
       a2hat = a2hat, beta2hat = beta2hat, b2hat = b2hat,
       logLik = final_ll,
       params = k,
       AIC = aic)
}

results <- list(gamma1hat = gamma1hat, gamma2hat = gamma2hat,
                a1hat = a1hat, beta1hat = beta1hat, b1hat = b1hat,
                a2hat = a2hat, beta2hat = beta2hat, b2hat = b2hat,
                logLik = final_ll,
                params = k,
                AIC = aic)

# Bootstrap inference
B <- 200
n <- nrow(interval_censored_result$mydata)

param_names <- c(
  paste0("gamma1_", 1:3), 
  paste0("gamma2_", 1:3),
  "a1", "beta1_1", "beta1_2", "b1",
  "a2", "beta2_1", "beta2_2", "b2"
)

boot_results <- matrix(NA, nrow = B, ncol = length(param_names))
colnames(boot_results) <- param_names
start_time <- Sys.time()

for (b in 1:B) {
  set.seed(123 + b)
  ind <- sample(1:n, size = n, replace = TRUE)
  mydata <- interval_censored_result$mydata[ind, ]
  xx <- interval_censored_result$xx[ind, , drop = FALSE]
  zz <- interval_censored_result$zz[ind, , drop = FALSE]
  u <- mydata$u
  
  re.data <- list(
    mydata = mydata,
    Status = mydata$Status,
    L = mydata$L,
    R = mydata$R,
    xx = xx,
    zz = zz,
    u = u
  )
  tryCatch({
    fit <- em_ic(re.data)
    boot_results[b, ] <- c(
      fit$gamma1hat, fit$gamma2hat,
      fit$a1hat, fit$beta1hat, fit$b1hat,
      fit$a2hat, fit$beta2hat, fit$b2hat
    )
    cat(sprintf("Iteration %d completed\n", b))
  }, error = function(e) {
    cat(sprintf("Iteration %d failed: %s\n", b, e$message))
  })
}

boot_results <- boot_results[complete.cases(boot_results), ]
successful_boots <- nrow(boot_results)

bootstrap_means <- apply(boot_results, 2, mean)
bootstrap_se <- apply(boot_results, 2, sd)
bootstrap_ci <- apply(boot_results, 2, quantile, probs = c(0.025, 0.975))

final_results <- data.frame(
  Parameter = param_names,
  Estimate = bootstrap_means,
  SE = bootstrap_se,
  CI_lower = bootstrap_ci[1, ],
  CI_upper = bootstrap_ci[2, ]
)

print(final_results)

# Wald test using original fit and bootstrap SE
original_fit <- em_ic(interval_censored_result)
theta_hat <- c(
  original_fit$gamma1hat, original_fit$gamma2hat,
  original_fit$a1hat, original_fit$beta1hat, original_fit$b1hat,
  original_fit$a2hat, original_fit$beta2hat, original_fit$b2hat
)
stopifnot(names(theta_hat) == colnames(boot_results))

final_results$z_value <- theta_hat / final_results$SE
final_results$p_value <- 2 * pnorm(-abs(final_results$z_value))

final_results$Significance <- ifelse(final_results$p_value < 0.001, "***",
                                     ifelse(final_results$p_value < 0.01, "**",
                                            ifelse(final_results$p_value < 0.05, "*", "")))

knitr::kable(final_results, digits = 4,
             caption = "Parameter Estimates with Bootstrap Inference")

end_time <- Sys.time()
run_time <- end_time - start_time

# Plot survival curves for event type 1 (default)
library(ggplot2)

weibull_survival <- function(t, a1, b1, beta1, xx) {
  xbeta1 <- xx %*% beta1
  S1 <- 1 - pweibull(t, shape = 1/b1, scale = exp(xbeta1 + a1))
  return(S1)
}

xx <- cbind(1, 2)
time <- seq(0, 80, by = 0.1)
survival_results <- matrix(0, nrow = 199, ncol = length(time))

for (i in 1:199) {
  a1hat <- boot_results[i, 7]
  b1hat <- boot_results[i, 10]
  beta1hat <- boot_results[i, 8:9]
  survival_results[i, ] <- sapply(time, function(t) weibull_survival(t, a1hat, b1hat, beta1hat, xx))
}

lcil_surv <- apply(survival_results, 2, quantile, probs = 0.025)
ucil_surv <- apply(survival_results, 2, quantile, probs = 0.975)
est_surv <- sapply(time, function(t) weibull_survival(t, a1 = bootstrap_means[7],
                                                      b1 = bootstrap_means[10], beta1 = bootstrap_means[8:9], xx = xx))

data <- data.frame(
  Time = rep(time, 3),
  Survival = c(est_surv, lcil_surv, ucil_surv),
  Type = rep(c("Estimate", "Lower CI", "Upper CI"), each = length(time))
)
data_clean <- na.omit(data)

ggplot(data_clean, aes(x = Time, y = Survival, color = Type, linetype = Type, linewidth = Type)) +
  geom_step() +
  scale_color_manual(values = c("Estimate" = "red", "Lower CI" = "darkgrey", "Upper CI" = "darkgrey")) +
  scale_linetype_manual(values = c("Estimate" = "dashed", "Lower CI" = "dotted", "Upper CI" = "longdash")) +
  scale_linewidth_manual(values = c("Estimate" = 1.2, "Lower CI" = 0.8, "Upper CI" = 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(title = "", x = "Months", y = "Survival Probability") +
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
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

save_path <- "C:/Users/17001/Desktop/loan"
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

ggsave(filename = file.path(save_path, "w.CR.S1.png"), plot = last_plot(), width = 4, height = 5, dpi = 300)

# Plot survival curves for event type 2 (prepayment)
weibull_survival <- function(t, a2, b2, beta2, xx) {
  xbeta2 <- xx %*% beta2
  S2 <- 1 - pweibull(t, shape = 1/b2, scale = exp(xbeta2 + a2))
  return(S2)
}

xx <- cbind(1, 2)
time <- seq(0, 80, by = 0.1)
survival_results <- matrix(0, nrow = 199, ncol = length(time))

for (i in 1:199) {
  a2hat <- boot_results[i, 11]
  b2hat <- boot_results[i, 14]
  beta2hat <- boot_results[i, 12:13]
  survival_results[i, ] <- sapply(time, function(t) weibull_survival(t, a2hat, b2hat, beta2hat, xx))
}

lcil_surv <- apply(survival_results, 2, quantile, probs = 0.025)
ucil_surv <- apply(survival_results, 2, quantile, probs = 0.975)
est_surv <- sapply(time, function(t) weibull_survival(t, a2 = bootstrap_means[11],
                                                      b2 = bootstrap_means[14], beta2 = bootstrap_means[12:13], xx = xx))

data <- data.frame(
  Time = rep(time, 3),
  Survival = c(est_surv, lcil_surv, ucil_surv),
  Type = rep(c("Estimate", "Lower CI", "Upper CI"), each = length(time))
)
data_clean <- na.omit(data)

ggplot(data_clean, aes(x = Time, y = Survival, color = Type, linetype = Type, linewidth = Type)) +
  geom_step() +
  scale_color_manual(values = c("Estimate" = "red", "Lower CI" = "darkgrey", "Upper CI" = "darkgrey")) +
  scale_linetype_manual(values = c("Estimate" = "dashed", "Lower CI" = "dotted", "Upper CI" = "longdash")) +
  scale_linewidth_manual(values = c("Estimate" = 1.2, "Lower CI" = 0.8, "Upper CI" = 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(title = "", x = "Months", y = "Survival Probability") +
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
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(filename = file.path(save_path, "w.CR.S2.png"), plot = last_plot(), width = 4, height = 5, dpi = 300)

# Overall survival probability (cured + default + prepayment)
weibull_survival <- function(t, gamma1, gamma2, a1, b1, beta1, a2, b2, beta2, xx, zz) {
  temp1 <- exp(cbind(zz) %*% gamma1) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
  temp2 <- exp(cbind(zz) %*% gamma2) / (exp(cbind(zz) %*% gamma1) + exp(cbind(zz) %*% gamma2) + 1)
  temp0 <- 1 - temp1 - temp2
  xbeta1 <- xx %*% beta1
  xbeta2 <- xx %*% beta2
  S1 <- 1 - pweibull(t, shape = 1/b1, scale = exp(xbeta1 + a1))
  S2 <- 1 - pweibull(t, shape = 1/b2, scale = exp(xbeta2 + a2))
  surv <- temp0 + temp1 * S1 + temp2 * S2
  return(surv)
}

xx <- cbind(1, 2)
zz <- cbind(2, 2, 2)
time <- seq(0, 80, by = 0.1)
survival_results <- matrix(0, nrow = 199, ncol = length(time))

for (i in 1:199) {
  gamma1hat <- boot_results[i, 1:3]
  gamma2hat <- boot_results[i, 4:6]
  a1hat <- boot_results[i, 7]
  b1hat <- boot_results[i, 10]
  beta1hat <- boot_results[i, 8:9]
  a2hat <- boot_results[i, 11]
  b2hat <- boot_results[i, 14]
  beta2hat <- boot_results[i, 12:13]
  survival_results[i, ] <- sapply(time, function(t) weibull_survival(t, gamma1hat, gamma2hat, a1hat, b1hat, beta1hat,
                                                                     a2hat, b2hat, beta2hat, xx, zz))
}

lcil_surv <- apply(survival_results, 2, quantile, probs = 0.025)
ucil_surv <- apply(survival_results, 2, quantile, probs = 0.975)
est_surv <- sapply(time, function(t) weibull_survival(t, gamma1 = bootstrap_means[1:3], gamma2 = bootstrap_means[4:6],
                                                      a1 = bootstrap_means[7], b1 = bootstrap_means[10], beta1 = bootstrap_means[8:9],
                                                      a2 = bootstrap_means[11], b2 = bootstrap_means[14], beta2 = bootstrap_means[12:13],
                                                      xx = xx, zz = zz))

data <- data.frame(
  Time = rep(time, 3),
  Survival = c(est_surv, lcil_surv, ucil_surv),
  Type = rep(c("Estimate", "Lower CI", "Upper CI"), each = length(time))
)
data_clean <- na.omit(data)

ggplot(data_clean, aes(x = Time, y = Survival, color = Type, linetype = Type, linewidth = Type)) +
  geom_step() +
  scale_color_manual(values = c("Estimate" = "red", "Lower CI" = "darkgrey", "Upper CI" = "darkgrey")) +
  scale_linetype_manual(values = c("Estimate" = "dashed", "Lower CI" = "dotted", "Upper CI" = "longdash")) +
  scale_linewidth_manual(values = c("Estimate" = 1.2, "Lower CI" = 0.8, "Upper CI" = 0.8)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(title = "", x = "Months", y = "Survival Probability") +
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
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

ggsave(filename = file.path(save_path, "w.CR.S.png"), plot = last_plot(), width = 4, height = 5, dpi = 300)

