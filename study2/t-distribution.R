# 중심 t분포 그래프 및 통계량 계산 함수 (Pooled variance 사용)
plot_t_distribution <- function(n1, n2, s1, s2, xlim = NULL, ylim_pdf = NULL, ylim_cdf = c(0, 1)) {
  
  # === 기본 계산 ===
  # 자유도 계산
  df1 <- n1 - 1
  df2 <- n2 - 1
  df_pooled <- df1 + df2
  
  # Pooled variance 계산
  pooled_var <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  pooled_sd <- sqrt(pooled_var)
  
  # 표준오차 계산
  se1 <- s1 / sqrt(n1)
  se2 <- s2 / sqrt(n2)
  se_pooled <- pooled_sd * sqrt(1/n1 + 1/n2)
  
  # === 모수 정보 출력 ===
  cat("=== 모수 정보 ===\n")
  cat(sprintf("n1 = %d , n2 = %d\n", n1, n2))
  cat(sprintf("s1 = %.1f , s2 = %.1f\n", s1, s2))
  cat(sprintf("자유도: df = %d (n1 + n2 - 2)\n", df_pooled))
  cat(sprintf("Pooled variance = %.4f\n", pooled_var))
  cat(sprintf("Pooled SD = %.4f\n", pooled_sd))
  cat(sprintf("s1/√n1 = %.4f\n", se1))
  cat(sprintf("s2/√n2 = %.4f\n", se2))
  cat(sprintf("SE(pooled) = %.4f\n", se_pooled))
  cat("\n")
  
  # === t분포 그래프 그리기 ===
  # x축 범위 설정 (사용자 지정 또는 기본값)
  if(is.null(xlim)) {
    x_range <- seq(-4 * se_pooled, 4 * se_pooled, length.out = 1000)
  } else {
    x_range <- seq(xlim[1], xlim[2], length.out = 1000)
  }
  
  # 표준화된 t값으로 변환
  t_values <- x_range / se_pooled
  
  # t분포 밀도와 누적확률 계산 (SE로 스케일링)
  y_density <- dt(t_values, df_pooled) / se_pooled
  y_cumulative <- pt(t_values, df_pooled)
  
  # 2x1 패널 설정
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  
  # 1. 밀도 함수 (PDF)
  plot(x_range, y_density, type = "l", col = "blue", lwd = 3,
       main = sprintf("t-Distribution Density Function\nt(%d), SE = %.4f", df_pooled, se_pooled),
       xlab = "Mean Difference (X̄₁ - X̄₂)", ylab = "Density", 
       ylim = if(is.null(ylim_pdf)) c(0, max(y_density) * 1.1) else ylim_pdf,
       xlim = if(is.null(xlim)) range(x_range) else xlim)
  
  # 평균선 추가
  abline(v = 0, col = "gray", lty = 2, lwd = 1)
  
  # 격자 추가
  grid(col = "lightgray", lty = 3)
  
  # 2. 누적분포함수 (CDF)
  plot(x_range, y_cumulative, type = "l", col = "red", lwd = 3,
       main = sprintf("t-Distribution Cumulative Distribution Function\nt(%d), SE = %.4f", df_pooled, se_pooled),
       xlab = "Mean Difference (X̄₁ - X̄₂)", ylab = "P(T ≤ t)", 
       ylim = ylim_cdf,
       xlim = if(is.null(xlim)) range(x_range) else xlim)
  
  # 50% 선 추가
  abline(h = 0.5, col = "gray", lty = 2, lwd = 1)
  abline(v = 0, col = "gray", lty = 2, lwd = 1)
  
  # 격자 추가
  grid(col = "lightgray", lty = 3)
  
  # 원래 레이아웃으로 복원
  par(mfrow = c(1, 1))
  
  # === 분위수 비교 (실제 평균차이 척도) ===
  cat("=== Key Quantiles for Mean Difference ===\n")
  cat(sprintf("Quantile   Mean Diff Scale   Standardized t\n"))
  cat("--------   --------------   --------------\n")
  
  quantiles <- c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975)
  q_labels <- c("2.5%", "5%", "10%", "90%", "95%", "97.5%")
  
  for(i in 1:length(quantiles)) {
    q_t_standardized <- qt(quantiles[i], df_pooled)
    q_mean_diff <- q_t_standardized * se_pooled
    
    cat(sprintf("%-8s %14.4f   %14.3f\n", q_labels[i], q_mean_diff, q_t_standardized))
  }
  
  # === 임계값 정보 ===
  cat("\n=== Two-tailed Test Critical Values (α = 0.05) ===\n")
  t_critical_std <- qt(0.975, df_pooled)
  mean_diff_critical <- t_critical_std * se_pooled
  cat(sprintf("Standardized t-value: ±%.3f\n", t_critical_std))
  cat(sprintf("Mean difference critical value: ±%.4f\n", mean_diff_critical))
  
  cat("\n=== One-tailed Test Critical Values (α = 0.05) ===\n")
  t_one_sided_std <- qt(0.95, df_pooled)
  mean_diff_one_sided <- t_one_sided_std * se_pooled
  cat(sprintf("Standardized t-value: %.3f\n", t_one_sided_std))
  cat(sprintf("Mean difference critical value: %.4f\n", mean_diff_one_sided))
  
  # 결과를 리스트로 반환
  results <- list(
    n1 = n1, n2 = n2, s1 = s1, s2 = s2,
    df_pooled = df_pooled,
    pooled_var = pooled_var, pooled_sd = pooled_sd,
    se1 = se1, se2 = se2, se_pooled = se_pooled,
    t_critical_two_sided = mean_diff_critical,
    t_critical_one_sided = mean_diff_one_sided,
    t_critical_std_two_sided = t_critical_std,
    t_critical_std_one_sided = t_one_sided_std
  )
  
  invisible(results)
}


#### noncentral t-distribution (Pooled variance 사용) ####
# 비중심 t분포 그래프 및 통계량 계산 함수 (Pooled variance 사용)
plot_noncentral_t <- function(n1, n2, s1, s2, delta = 0, xlim = NULL, ylim_pdf = NULL, ylim_cdf = c(0, 1), show_central = TRUE) {
  
  # === 기본 계산 ===
  # 자유도 계산
  df1 <- n1 - 1
  df2 <- n2 - 1
  df_pooled <- df1 + df2
  
  # Pooled variance 계산
  pooled_var <- ((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)
  pooled_sd <- sqrt(pooled_var)
  
  # 표준오차 계산
  se1 <- s1 / sqrt(n1)
  se2 <- s2 / sqrt(n2)
  se_pooled <- pooled_sd * sqrt(1/n1 + 1/n2)
  
  # 비중심 모수 계산 (noncentrality parameter)
  ncp <- delta / se_pooled
  
  # === 모수 정보 출력 ===
  cat("=== 모수 정보 ===\n")
  cat(sprintf("n1 = %d , n2 = %d\n", n1, n2))
  cat(sprintf("s1 = %.1f , s2 = %.1f\n", s1, s2))
  cat(sprintf("자유도: df = %d (n1 + n2 - 2)\n", df_pooled))
  cat(sprintf("Pooled variance = %.4f\n", pooled_var))
  cat(sprintf("Pooled SD = %.4f\n", pooled_sd))
  cat(sprintf("SE(pooled) = %.4f\n", se_pooled))
  cat(sprintf("Delta (δ) = %.2f\n", delta))
  cat(sprintf("Noncentrality parameter (δ) = %.4f\n", ncp))
  cat("\n")
  
  # === t분포 그래프 그리기 ===
  # x축 범위 설정 (사용자 지정 또는 기본값)
  if(is.null(xlim)) {
    x_range <- seq(delta - 4 * se_pooled, delta + 4 * se_pooled, length.out = 1000)
  } else {
    x_range <- seq(xlim[1], xlim[2], length.out = 1000)
  }
  
  # 표준화된 t값으로 변환
  t_values <- x_range / se_pooled
  
  # 비중심 t분포 밀도와 누적확률 계산 (SE로 스케일링)
  y_density_nc <- dt(t_values, df_pooled, ncp = ncp) / se_pooled
  y_cumulative_nc <- pt(t_values, df_pooled, ncp = ncp)
  
  # 중심 t분포도 계산 (비교용)
  y_density_c <- dt(t_values, df_pooled) / se_pooled
  y_cumulative_c <- pt(t_values, df_pooled)
  
  # 2x1 패널 설정
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # 1. 밀도 함수 (PDF)
  max_density <- max(c(y_density_nc, if(show_central) y_density_c else 0))
  
  plot(x_range, y_density_nc, type = "l", col = "red", lwd = 3,
       main = sprintf("Noncentral t-Distribution (δ = %.3f)\ndf = %d, Delta = %.2f", ncp, df_pooled, delta),
       xlab = "Mean Difference (X̄₁ - X̄₂)", ylab = "Density", 
       ylim = if(is.null(ylim_pdf)) c(0, max_density * 1.1) else ylim_pdf,
       xlim = if(is.null(xlim)) range(x_range) else xlim)
  
  # 중심 t분포 추가 (비교용)
  if(show_central) {
    lines(x_range, y_density_c, col = "blue", lwd = 2, lty = 2)
    legend("topright", 
           legend = c(sprintf("Noncentral t (δ=%.3f)", ncp), "Central t (δ=0)"),
           col = c("red", "blue"), lwd = c(3, 2), lty = c(1, 2),
           cex = 0.8, x.intersp = 0.5, y.intersp = 0.5, text.width=strwidth("Central t (δ=0)", cex=0.8))
  }
  
  # 0 기준선 추가
  abline(v = 0, col = "gray", lty = 2, lwd = 1)
  
  # 격자 추가
  grid(col = "lightgray", lty = 3)
  
  # 2. 누적분포함수 (CDF)
  plot(x_range, y_cumulative_nc, type = "l", col = "red", lwd = 3,
       main = sprintf("Noncentral t-Distribution CDF (δ = %.3f)\ndf = %d, Delta = %.2f", ncp, df_pooled, delta),
       xlab = "Mean Difference (X̄₁ - X̄₂)", ylab = "P(T ≤ t)", 
       ylim = ylim_cdf,
       xlim = if(is.null(xlim)) range(x_range) else xlim)
  
  # 중심 t분포 CDF 추가 (비교용)
  if(show_central) {
    lines(x_range, y_cumulative_c, col = "blue", lwd = 2, lty = 2)
    legend("bottomright", 
           legend = c(sprintf("Noncentral t (δ=%.3f)", ncp), "Central t (δ=0)"),
           col = c("red", "blue"), lwd = c(3, 2), lty = c(1, 2),
           cex = 0.8, x.intersp = 0.5, y.intersp = 0.5, text.width=strwidth("Central t (δ=0)", cex=0.8))
  }
  
  # 50% 선과 0 기준선 추가
  abline(h = 0.5, col = "gray", lty = 2, lwd = 1)
  abline(v = 0, col = "gray", lty = 2, lwd = 1)
  
  # 격자 추가
  grid(col = "lightgray", lty = 3)
  
  # 원래 레이아웃으로 복원
  par(mfrow = c(1, 1))
  
  # === 분위수 비교 (비중심 t분포) ===
  cat("=== Key Quantiles for Noncentral t-Distribution ===\n")
  cat(sprintf("Quantile   Mean Diff Scale   Standardized t\n"))
  cat("--------   --------------   --------------\n")
  
  quantiles <- c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975)
  q_labels <- c("2.5%", "5%", "10%", "90%", "95%", "97.5%")
  
  for(i in 1:length(quantiles)) {
    q_t_standardized <- qt(quantiles[i], df_pooled, ncp = ncp)
    q_mean_diff <- q_t_standardized * se_pooled
    
    cat(sprintf("%-8s %14.4f   %14.3f\n", q_labels[i], q_mean_diff, q_t_standardized))
  }
  
  # === 검정력 분석 ===
  cat("\n=== Power Analysis (Two-tailed test, α = 0.05) ===\n")
  
  # 중심 t분포의 임계값
  t_critical_std <- qt(0.975, df_pooled)
  mean_diff_critical_upper <- t_critical_std * se_pooled
  mean_diff_critical_lower <- -t_critical_std * se_pooled
  
  cat(sprintf("Critical values (central t): ±%.4f (±%.3f standardized)\n", 
              mean_diff_critical_upper, t_critical_std))
  
  # 비중심 t분포에서의 검정력 계산
  # P(T > t_critical | ncp) + P(T < -t_critical | ncp)
  power_upper <- 1 - pt(t_critical_std, df_pooled, ncp = ncp)
  power_lower <- pt(-t_critical_std, df_pooled, ncp = ncp)
  total_power <- power_upper + power_lower
  
  cat(sprintf("Power (reject H0 when delta = %.2f): %.4f (%.2f%%)\n", 
              delta, total_power, total_power * 100))
  cat(sprintf("  - Upper tail power: %.4f\n", power_upper))
  cat(sprintf("  - Lower tail power: %.4f\n", power_lower))
  
  # === 일측 검정 검정력 ===
  cat("\n=== Power Analysis (One-tailed test, α = 0.05) ===\n")
  t_one_sided_std <- qt(0.95, df_pooled)
  mean_diff_one_sided <- t_one_sided_std * se_pooled
  
  cat(sprintf("Critical value (one-sided): %.4f (%.3f standardized)\n", 
              mean_diff_one_sided, t_one_sided_std))
  
  # 우측 검정의 검정력 (H1: μ1 > μ2)
  power_right <- 1 - pt(t_one_sided_std, df_pooled, ncp = ncp)
  cat(sprintf("Power (right-tailed test): %.4f (%.2f%%)\n", 
              power_right, power_right * 100))
  
  # 좌측 검정의 검정력 (H1: μ1 < μ2)
  power_left <- pt(-t_one_sided_std, df_pooled, ncp = ncp)
  cat(sprintf("Power (left-tailed test): %.4f (%.2f%%)\n", 
              power_left, power_left * 100))
  
  # 결과를 리스트로 반환
  results <- list(
    n1 = n1, n2 = n2, s1 = s1, s2 = s2,
    df_pooled = df_pooled,
    pooled_var = pooled_var, pooled_sd = pooled_sd,
    se1 = se1, se2 = se2, se_pooled = se_pooled,
    delta = delta,
    ncp = ncp,
    t_critical_two_sided = mean_diff_critical_upper,
    t_critical_one_sided = mean_diff_one_sided,
    t_critical_std_two_sided = t_critical_std,
    t_critical_std_one_sided = t_one_sided_std,
    power_two_sided = total_power,
    power_one_sided_right = power_right,
    power_one_sided_left = power_left
  )
  
  invisible(results)
}
