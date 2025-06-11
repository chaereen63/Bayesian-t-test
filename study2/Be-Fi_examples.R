source("study2/Behrens-Fisher_distribution.R")

#### central Behrens-Fisher distribution examples (delta = 0) ####

# 예시 1: 모두 균등한 경우
plot_befi_central(n1=30, n2=30, s1=2, s2=2)

# 예시 2: 표본크기가 다른 경우
# 2-a: 집단1에 더 의존하는 경우(T1분포에 더 가까운 경우)
plot_befi_central(n1=10, n2=50, s1=2, s2=2)
# 2-b: 집단2에 더 의존하는 경우(T2분포에 더 가까운 경우)
plot_befi_central(n1=50, n2=10, s1=2, s2=2)

# 예시 3: 분산만 다른 경우
plot_befi_central(n1=30, n2=30, s1=4, s2=2)

# 예시 4: 분산과 표본크기가 모두 다른 경우
# 2-a: 집단1에 더 의존하는 경우(T1분포에 더 가까운 경우)
plot_befi_central(n1=10, n2=50, s1=4, s2=2)
# 2-b: 집단2에 더 의존하는 경우(T2분포에 더 가까운 경우)
plot_befi_central(n1=50, n2=10, s1=4, s2=2)

# 예시 5: 극단적인 경우
plot_befi_central(n1=5, n2=30, s1=2, s2=0.5)

#### location shifted Behrens-Fisher distribution examples (delta = 1.5로 설정) ####

# 표준화 효과크기(Cohen's d)를 0.5로 고정
target_cohens_d <- 0.5

# 각 조건에 맞는 delta 계산 함수
calculate_delta <- function(s1, s2, cohens_d = target_cohens_d) {
  variance <- (s1^2 + s2^2) / 2
  delta <- cohens_d * sqrt(variance)
  return(delta)
}

# 예시 1: 균등한 조건 (n1=n2, s1=s2)
cat("=== 예시 1: 균등한 조건 (n1=n2, s1=s2) ===\n")
delta1 <- calculate_delta(s1=2, s2=2)
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta1))
plot_befi_noncentral(n1=30, n2=30, s1=2, s2=2, delta=delta1)

# 예시 2: 표본크기만 다른 경우
delta2 <- calculate_delta(s1=2, s2=2)
cat("\n=== 예시 2-a: 집단1에 더 의존하는 경우 (n1 < n2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta2))
plot_befi_noncentral(n1=10, n2=50, s1=2, s2=2, delta=delta2)

cat("\n=== 예시 2-b: 집단2에 더 의존하는 경우 (n1 > n2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta2))
plot_befi_noncentral(n1=50, n2=10, s1=2, s2=2, delta=delta2)

# 예시 3: 분산만 다른 경우 (s1 >> s2)
cat("\n=== 예시 3: 분산만 다른 경우 (s1 >> s2) ===\n")
delta3 <- calculate_delta(s1=4, s2=2)
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta3))
plot_befi_noncentral(n1=30, n2=30, s1=4, s2=2, delta=delta3)

# 예시 4: 분산과 표본크기가 모두 다른 경우
delta4 <- calculate_delta(s1=4, s2=2)
cat("\n=== 예시 4-a: 집단1에 더 의존하는 경우 (n1 < n2, s1 > s2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta4))
plot_befi_noncentral(n1=10, n2=50, s1=4, s2=2, delta=delta4)

cat("\n=== 예시 4-b: 집단2에 더 의존하는 경우 (n1 > n2, s1 > s2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta4))
plot_befi_noncentral(n1=50, n2=10, s1=4, s2=2, delta=delta4)

# 각 조건별 delta 값 요약
{cat("\n=== Delta 값 요약 ===\n")
cat(sprintf("예시 1 (s1=2, s2=2): delta = %.3f\n", delta1))
cat(sprintf("예시 2 (s1=2, s2=2): delta = %.3f\n", delta2))
cat(sprintf("예시 3 (s1=16, s2=1): delta = %.3f\n", delta3))
cat(sprintf("예시 4 (s1=16, s2=1): delta = %.3f\n", delta4))
cat(sprintf("모든 경우에서 Cohen's d = %.2f로 동일\n", target_cohens_d))}

# 예시 5: 극단적인 경우
# case 4-a
plot_befi_noncentral(n1=30, n2=5, s1=2, s2=0.5, delta=1.5)
# case 4-b
plot_befi_noncentral(n1=5, n2=30, s1=2, s2=0.5, delta=1.5)

  # 분산비율이 고정되어 있을 경우, 자유도가 작은 분포를 따른다.

#### 추가 예시: 다양한 delta 값들의 비교(수정필요 2025.06.10) ####
cat("\n=== 추가 예시: 다양한 delta 값들의 비교 ===\n")

# delta 값들의 영향 비교 함수
compare_delta_effects <- function(n1, n2, s1, s2, deltas=c(0, 0.5, 1, 1.5, 2)) {
  
  # R 매개변수 계산
  R <- atan((s1/sqrt(n1))/(s2/sqrt(n2)))
  
  # x 값들 (scale free - delta에 따라 자동 조정)
  max_delta <- max(deltas)
  x_min <- -3
  x_max <- 3 + 1.5*max_delta  # delta가 클수록 오른쪽으로 확장
  x <- seq(x_min, x_max, length.out=200)
  
  # 색상 설정
  colors <- c("black", "blue", "red", "green", "purple")
  
  # 플롯 설정
  par(mfrow=c(1,2), mar=c(4,4,3,2), mgp=c(2.5,1,0))
  
  # PDF 플롯
  plot(x, rep(0, length(x)), type='n', 
       main=paste("Delta Effects on PDF\nn1=", n1, ", n2=", n2, ", s1=", s1, ", s2=", s2),
       xlab='x', ylab='Density',
       ylim=c(0, 0.5))
  
  for(i in 1:length(deltas)) {
    pdf_vals <- sapply(x, function(t) dbefi_noncentral(t, n1, n2, delta=deltas[i], s1=s1, s2=s2))
    lines(x, pdf_vals, col=colors[i], lwd=2, lty=i)
  }
  
  legend('topright', 
         legend=paste('δ =', deltas),
         col=colors, lty=1:length(deltas), lwd=2,
         cex=0.8)
  
  # CDF 플롯
  plot(x, rep(0, length(x)), type='n',
       main="Delta Effects on CDF",
       xlab='x', ylab='P(X ≤ x)', ylim=c(0,1))
  
  for(i in 1:length(deltas)) {
    cdf_vals <- pbefi_noncentral(x, n1, n2, delta=deltas[i], s1=s1, s2=s2)
    lines(x, cdf_vals, col=colors[i], lwd=2, lty=i)
  }
  
  legend('bottomright', 
         legend=paste('δ =', deltas),
         col=colors, lty=1:length(deltas), lwd=2,
         cex=0.8)
  
  # 분위수 비교
  cat("\n=== Delta에 따른 분위수 비교 ===\n")
  quantiles <- c(0.025, 0.05, 0.5, 0.95, 0.975)
  
  cat(sprintf("%-8s", "분위수"))
  for(d in deltas) {
    cat(sprintf(" %-10s", paste("δ=", d, sep="")))
  }
  cat("\n")
  
  cat(sprintf("%-8s", "------"))
  for(d in deltas) {
    cat(sprintf(" %-10s", "--------"))
  }
  cat("\n")
  
  for(q in quantiles) {
    cat(sprintf("%-8s", paste(q*100, "%", sep="")))
    for(d in deltas) {
      quant_val <- uniroot(function(x) pbefi_noncentral(x, n1, n2, delta=d, s1=s1, s2=s2) - q, 
                           interval=c(-10, 10))$root
      cat(sprintf(" %-10.3f", quant_val))
    }
    cat("\n")
  }
}

# 균등한 경우에서 delta 효과 비교
compare_delta_effects(n1=30, n2=30, s1=2, s2=2)

# 불균등한 경우에서 delta 효과 비교
cat("\n=== 불균등한 경우에서 Delta 효과 ===\n")
compare_delta_effects(n1=10, n2=50, s1=2, s2=2)

# 검정력 분석 함수
power_analysis_befi <- function(n1, n2, s1, s2, alpha=0.05, deltas=seq(0, 3, 0.1)) {
  
  # 임계값 계산 (중심 분포에서)
  critical_val <- uniroot(function(x) pbefi(x, n1, n2, s1=s1, s2=s2) - (1-alpha/2), 
                          interval=c(-10, 10))$root
  
  # 각 delta에 대한 검정력 계산
  power <- sapply(deltas, function(d) {
    1 - pbefi_noncentral(critical_val, n1, n2, delta=d, s1=s1, s2=s2) + 
      pbefi_noncentral(-critical_val, n1, n2, delta=d, s1=s1, s2=s2)
  })
  
  # 플롯
  plot(deltas, power, type='l', lwd=3, col='blue',
       main=paste("Power Analysis for Behrens-Fisher Test\nn1=", n1, ", n2=", n2, ", s1=", s1, ", s2=", s2, ", α=", alpha),
       xlab='Effect Size (δ = μ1 - μ2)', ylab='Power',
       ylim=c(0, 1))
  abline(h=0.8, col='red', lty=2, lwd=2)
  abline(h=alpha, col='gray', lty=3)
  
  # 80% 검정력을 위한 effect size 찾기
  if(max(power) > 0.8) {
    delta_80 <- deltas[min(which(power >= 0.8))]
    abline(v=delta_80, col='red', lty=2, lwd=2)
    text(delta_80, 0.5, paste("δ =", round(delta_80, 2), "\nfor 80% power"), pos=4)
  }
  
  legend('bottomright', 
         legend=c('Power Curve', '80% Power', paste('α =', alpha)),
         col=c('blue', 'red', 'gray'), lty=c(1,2,3), lwd=c(3,2,1),
         cex=0.8)
}

cat("\n=== 검정력 분석 ===\n")
power_analysis_befi(n1=30, n2=30, s1=2, s2=2)

cat("\n=== 불균등한 경우의 검정력 분석 ===\n")
power_analysis_befi(n1=10, n2=50, s1=2, s2=2)
