source("./study2/t-distribution.R")

#### centrality t-distribution ####
# 예시 1: 모두 균등한 경우
plot_t_distribution(n1=30, n2=30, s1=2, s2=2, xlim = c(-4, 4))

# 예시 2: 표본크기가 다른 경우
# 2-a: centrality parameter에 곱해지는 effective sample sizes(harmonic mean)이 작아지면서 SE 감소
plot_t_distribution(n1=10, n2=50, s1=2, s2=2, xlim = c(-4, 4))
# 2-b: same with 2-a
plot_t_distribution(n1=50, n2=10, s1=2, s2=2, xlim = c(-4, 4))

# 예시 3: 분산만 다른 경우
plot_t_distribution(n1=30, n2=30, s1=4, s2=2, xlim = c(-4, 4))

# 예시 4: 분산과 표본크기가 모두 다른 경우
# 2-a: 큰 분산에 작은 가중 -> 더 작은 SE
plot_t_distribution(n1=10, n2=50, s1=4, s2=2, xlim = c(-4, 4))
# 2-b: 큰 분산에 큰 가중 -> 더 큰 SE
plot_t_distribution(n1=50, n2=10, s1=4, s2=2, xlim = c(-4, 4))


#### noncentrality t-distribution ####
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
plot_noncentral_t(n1=30, n2=30, s1=2, s2=2, delta=delta1, xlim = c(-2, 5))

# 예시 2: 표본크기만 다른 경우
delta2 <- calculate_delta(s1=2, s2=2)
cat("\n=== 예시 2-a: small effective sample size ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta2))
plot_noncentral_t(n1=10, n2=50, s1=2, s2=2, delta=delta2, xlim = c(-2, 5))

cat("\n=== 예시 2-b: same with 2-b ===\n")

# 예시 3: 분산만 다른 경우 (s1 >> s2)
cat("\n=== 예시 3: 분산만 다른 경우 (s1 > s2) ===\n")
delta3 <- calculate_delta(s1=4, s2=2)
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta3))
plot_noncentral_t(n1=30, n2=30, s1=4, s2=2, delta=delta3, xlim = c(-2, 5))

# 예시 4: 분산과 표본크기가 모두 다른 경우
delta4 <- calculate_delta(s1=4, s2=2)
cat("\n=== 예시 4-a: 큰 분산에 작은 가중 (n1 < n2, s1 > s2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta4))
plot_noncentral_t(n1=10, n2=50, s1=4, s2=2, delta=delta4)

cat("\n=== 예시 4-b: 큰 분산에 큰 가중 (n1 > n2, s1 > s2) ===\n")
cat(sprintf("Cohen's d = %.2f, delta = %.3f\n", target_cohens_d, delta4))
plot_noncentral_t(n1=50, n2=10, s1=4, s2=2, delta=delta4)

# 각 조건별 delta 값 요약
{cat("\n=== Delta 값 요약 ===\n")
  cat(sprintf("예시 1 (s1=2, s2=2): delta = %.3f\n", delta1))
  cat(sprintf("예시 2 (s1=2, s2=2): delta = %.3f\n", delta2))
  cat(sprintf("예시 3 (s1=16, s2=1): delta = %.3f\n", delta3))
  cat(sprintf("예시 4 (s1=16, s2=1): delta = %.3f\n", delta4))
  cat(sprintf("모든 경우에서 Cohen's d = %.2f로 동일\n", target_cohens_d))}
