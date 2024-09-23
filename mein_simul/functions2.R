library(stats)

simulate_data <- function(n1, n2, delta, rho_alpha, rho_beta, seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  
  # rho를 베타 분포에서 샘플링
  rho <- rbeta(1, shape1 = rho_alpha, shape2 = rho_beta)
  
  # 전체 표준편차 설정
  grand_sd <- 1
  
  # rho를 이용해 그룹별 분산 계산
  var1 <- grand_sd^2 / rho
  var2 <- grand_sd^2 / (1 - rho)
  
  # 그룹별 표준편차 계산
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  
  # pooled standard deviation 계산
  pooled_sd <- sqrt((sd1^2 * (n1 - 1) + sd2^2 * (n2 - 1)) / (n1 + n2 - 2))
  
  # 그룹별 평균 계산
  grand_mean <- 0
  mean1 <- grand_mean + 0.5 * delta * pooled_sd
  mean2 <- grand_mean - 0.5 * delta * pooled_sd
  
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

pooled_sd <- function(sd1, sd2, n1, n2){
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

cohens_d <- function(mean1, mean2, sd1, sd2, n1, n2){
  (mean1 - mean2) / pooled_sd(sd1, sd2, n1, n2)
}

calculate_sdr <- function(sd1, sd2) {
  return(sd2 / sd1)
}

calculate_rho <- function(var1, var2) {
  return(1 / var1 / (1 / var1 + 1 / var2))
}

# 기존의 다른 헬퍼 함수들은 그대로 유지
get_loop <- function(dir, file){
  # 기존 코드 유지
}

get_true_model <- function(delta, rho){
  effect        <- delta != 0
  heterogeneity <- rho   != 0.5
  return(sapply(paste0(effect, heterogeneity), function(M) switch(
    M,
    "FALSEFALSE" = 1,
    "FALSETRUE"  = 2,
    "TRUEFALSE"  = 3,
    "TRUETRUE"   = 4
  )))
}

empty_plot <- function() {
  # 기존 코드 유지
}

get_RMSE <- function(estimate, true){
  return(sqrt(mean((estimate - true)^2)))
}

