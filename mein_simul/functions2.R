library(stats)
library(dplyr)

simulate_data <- function(n1, n2, mean1, mean2, sd1, sd2, rho, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # 데이터 생성
  x1 <- rnorm(n1, mean = mean1, sd = sd1)
  x2 <- rnorm(n2, mean = mean2, sd = sd2)
  
  return(list(
    x1 = x1,
    x2 = x2,
    rho = rho,
    sdr = sd2 / sd1,
    mean1 = mean1,
    mean2 = mean2,
    sd1 = sd1,
    sd2 = sd2
  ))
}

pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

cohens_d <- function(mean1, mean2, sd1, sd2, n1, n2) {
  (mean1 - mean2) / pooled_sd(sd1, sd2, n1, n2)
}

calculate_sdr <- function(sd1, sd2) {
  return(sd2 / sd1)
}

calculate_rho <- function(var1, var2) {
  return(1 / var1 / (1 / var1 + 1 / var2))
}

get_true_model <- function(delta, rho, delta_threshold = 0.1, rho_threshold = 0.1) {
  effect <- abs(delta) > delta_threshold
  heterogeneity <- abs(rho - 0.5) > rho_threshold #rho가 변수이므로 임계값으로 설정
  
  case_when(
    !effect & !heterogeneity ~ 1,  # No effect, No heterogeneity
    !effect & heterogeneity ~ 2,   # No effect, Heterogeneity
    effect & !heterogeneity ~ 3,   # Effect, No heterogeneity
    effect & heterogeneity ~ 4     # Effect, Heterogeneity
  )
}

get_RMSE <- function(estimate, true) {
  return(sqrt(mean((estimate - true)^2)))
}

# Helper functions for simulation tracking and visualization

get_loop <- function(dir, file, max_loop) {
  tracker_file <- file.path(dir, paste0(file, ".txt"))
  
  if (file.exists(tracker_file)) {
    loop <- as.numeric(readLines(tracker_file))
  } else {
    loop <- 1
  }
  
  # Ensure loop doesn't exceed max_loop
  loop <- min(loop, max_loop)
  
  # Increment and save the loop counter
  next_loop <- min(loop + 1, max_loop)
  writeLines(as.character(next_loop), tracker_file)
  
  return(loop)
}

update_tracker <- function(home_dir, tracker) {
  # This function can remain empty
}

empty_plot <- function() {
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', 
       xlim = c(0,1), ylim = c(0,1), mar = c(0,0,0,0))
}

# Compile Stan model once
stan_model <- stan_model(model_code = "
data {
  int<lower=0> N1;
  int<lower=0> N2;
  vector[N1] y1;
  vector[N2] y2;
}
parameters {
  real mu;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real delta;
}
transformed parameters {
  real<lower=0> pooled_sigma;
  real alpha;
  pooled_sigma = sqrt((sigma1^2 * (N1 - 1) + sigma2^2 * (N2 - 1)) / (N1 + N2 - 2));
  alpha = delta * pooled_sigma;
}
model {
  mu ~ cauchy(0, 1);
  sigma1 ~ cauchy(0, 1);
  sigma2 ~ cauchy(0, 1);
  delta ~ cauchy(0, 1);
  
  y1 ~ normal(mu - alpha/2, sigma1);
  y2 ~ normal(mu + alpha/2, sigma2);
}
")

# Renamed and slightly modified Wetzels' MCMC t-test function
wetzels_ttest <- function(y1, y2, iter = 5000, chains = 2, warmup = 1000) {
  data_list <- list(N1 = length(y1), N2 = length(y2), y1 = y1, y2 = y2)
  fit <- sampling(stan_model, data = data_list, iter = iter, chains = chains, warmup = warmup)
  
  # Extract posterior samples
  posterior_samples <- extract(fit)$delta
  
  # Estimate posterior density
  posterior_density <- density(posterior_samples)
  
  # Estimate posterior density at δ = 0
  posterior_at_zero <- approx(posterior_density$x, posterior_density$y, xout = 0)$y
  
  # Prior density at δ = 0 (Cauchy(0, 1))
  prior_at_zero <- dcauchy(0, 0, 1)
  
  # Calculate Bayes Factor (BF10)
  BF10 <- prior_at_zero / posterior_at_zero
  
  return(list(
    fit = fit,
    BF10 = BF10,
    posterior_samples = posterior_samples,
    summary = summary(fit)
  ))
}
