# 목표 통계량
target_mean1 <- 50
target_var1 <- 12
target_mean2 <- 55
target_var2 <- 40
n1 <- 20
n2 <- 12

# 시드별 결과를 저장할 데이터프레임
results <- data.frame()

# 많은 시드값을 시도
for(seed in 1:1000) {
  set.seed(seed)
  
  # 데이터 생성
  x1 <- rnorm(n1, mean = target_mean1, sd = sqrt(target_var1))
  x2 <- rnorm(n2, mean = target_mean2, sd = sqrt(target_var2))
  
  # 실제 통계량 계산
  mean1 <- mean(x1)
  var1 <- var(x1)
  mean2 <- mean(x2)
  var2 <- var(x2)
  
  # 목표값과의 차이 계산
  mean1_diff <- abs(mean1 - target_mean1)
  var1_diff <- abs(var1 - target_var1)
  mean2_diff <- abs(mean2 - target_mean2)
  var2_diff <- abs(var2 - target_var2)
  
  # 전체 오차 계산 (가중치 부여 가능)
  total_error <- mean1_diff + var1_diff + mean2_diff + var2_diff
  
  # 정규성 검정
  shapiro1_p <- shapiro.test(x1)$p.value
  shapiro2_p <- shapiro.test(x2)$p.value
  
  # 결과 저장
  results <- rbind(results, data.frame(
    seed = seed,
    mean1 = mean1,
    var1 = var1,
    mean2 = mean2,
    var2 = var2,
    total_error = total_error,
    shapiro1_p = shapiro1_p,
    shapiro2_p = shapiro2_p
  ))
}

# 정규성 가정을 만족하는 결과들 중에서 가장 좋은 시드 찾기
valid_results <- subset(results, shapiro1_p > 0.05 & shapiro2_p > 0.05)
best_seed <- valid_results[which.min(valid_results$total_error), ]

print("최적의 시드값의 결과:")
print(best_seed) # seed = 986

# 최적의 시드로 데이터 생성
set.seed(best_seed$seed)
x1_final <- rnorm(n1, mean = target_mean1, sd = sqrt(target_var1))
x2_final <- rnorm(n2, mean = target_mean2, sd = sqrt(target_var2))

cat("\n최종 생성된 데이터의 통계량:\n")
cat("집단 1 (목표: 평균=50, 분산=12):\n")
cat("  평균:", mean(x1_final), "\n")
cat("  분산:", var(x1_final), "\n")
cat("  정규성 검정 p값:", shapiro.test(x1_final)$p.value, "\n\n")

cat("집단 2 (목표: 평균=55, 분산=40):\n")
cat("  평균:", mean(x2_final), "\n")
cat("  분산:", var(x2_final), "\n")
cat("  정규성 검정 p값:", shapiro.test(x2_final)$p.value, "\n")

# 최종 데이터 출력
cat("\n최종 데이터:\n")
print("집단 1:")
print(x1_final)
print("\n집단 2:")
print(x2_final)

###################################
# 최적화된 시드값 사용
set.seed(986)

# 파라미터
n1 <- 20
mu1 <- 50
var1 <- 12
n2 <- 12
mu2 <- 55
var2 <- 40

# 데이터 생성
x1 <- rnorm(n1, mean = mu1, sd = sqrt(var1))
x2 <- rnorm(n2, mean = mu2, sd = sqrt(var2))

# 데이터 확인
cat("집단 1 요약통계:\n")
cat("평균:", mean(x1), "\n")
cat("분산:", var(x1), "\n")
cat("정규성 검정 p값:", shapiro.test(x1)$p.value, "\n\n")

cat("집단 2 요약통계:\n")
cat("평균:", mean(x2), "\n")
cat("분산:", var(x2), "\n")
cat("정규성 검정 p값:", shapiro.test(x2)$p.value, "\n\n")


# SD 함수 실행

library(rstan)
source("sd_ttest_stan.R", encoding = "UTF-8")
fit_wetzels <- SD_stan(x1, 
                  iters=5000, 
                  warmup=1000, 
                  chains=2, 
                  # sample=1,  # 두 집단 비교
                  # sig=2, 두 집단의 분산이 다름
                  prior='cauchy', 
                  plot=TRUE)
