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
