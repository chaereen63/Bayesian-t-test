library(BayesFactor)
source(file = file.path("./New/functionsN.R"))

# n=25인 경우 데이터 생성
n_small <- 25
mu1 <- 0.2
mu2 <- -0.2
sd1 <- 2
sd2 <- 2

# 두 그룹 생성
group1_small <- rnorm(n_small, mean = mu1, sd = sd1)
group2_small <- rnorm(n_small, mean = mu2, sd = sd2)

# 결합된 데이터프레임 생성
data_small <- data.frame(
  value = c(group1_small, group2_small),
  group = factor(rep(c("Group1", "Group2"), each = n_small))
)

# t-test 실행
t_result_small <- t.test(value ~ group, data = data_small, var.equal = TRUE)
print(t_result_small)

#JZS
JZSraw <- ttestBF(group1_small, group2_small, rscale = 1)
print(JZSraw)
JZSstat <- ttest.tstat(t_result_small$statistic, n1=25, n2=25, rscale = 1, simple=T)
print(JZSstat) #ttestBF와 동일한 값을 생성하는지 확인
#Behrens-Fisher
BeFiraw <- BeFiBF(x1 = group1_small, x2 = group2_small)
print(BeFiraw$bf10)
# effect size
library(effectsize)
effect_size_small <- cohens_d(value ~ group, data = data_small)
print(effect_size_small)

log_JZS <- log10(JZSstat)
log_BeFi <- log10(BeFi$bf10)
print(c(log_JZS,log_BeFi))

# Rouder et al. (2009) critical t-value : N = 50, t = 2.68
Rouder_sample <- ttest.tstat(2.68, n1=50, rscale = 1, simple=T)
(R.sample_effect <- 2.68*sqrt(1/25 + 1/25))
 ## Rouder et al. (2009)이 논문에 제시한 t값과 표본크기로 보면, 효과크기는 0.76이어야 BF10 > 3

#### 표본크기를 이용한 power analysis ####
library(ggplot2)
library(BayesFactor)
# Girόn and Del Castillo (2021)
library(pracma)
gicaBF <- function(mean1, mean2, var1, var2, n1, n2) {
  
  # Calculate sample statistics
  xbar1 <- mean1
  xbar2 <- mean2
  var1 <- var1
  var2 <- var2
  n1 <- n1
  n2 <- n2
  
  # Calculate degrees of freedom
  df1 <- n1 - 1
  df2 <- n2 - 1
  
  # Location parameter
  d <- xbar2 - xbar1
  
  # Define the numerator integrand
  num_integrand <- function(delta) {
    term1 <- (1 + n1 * delta^2/(df1 * var1))^(-n1/2)
    term2 <- (1 + n2 * ((delta-d)^2)/(df2 * var2))^(-n2/2)
    return(term1 * term2)
  }
  
  # Define the denominator integrand
  den_integrand <- function(delta) {
    term1 <- (1 + (n1/(n1+1)) * delta^2/(df1 * var1))^(-n1/2)
    term2 <- (1 + (n2/(n2+1)) * ((delta-d)^2)/(df2 * var2))^(-n2/2)
    return(term1 * term2)
  }
  
  # Calculate reasonable bounds for the integral
  bound <- max(abs(d) + 10 * sqrt(max(var1, var2)), 100)
  
  # Compute the integrals
  num_integral <- integral(num_integrand, -bound, bound)
  den_integral <- integral(den_integrand, -bound, bound)
  
  # Calculate the product term
  prod_term <- sqrt(n1 + 1) * sqrt(n2 + 1)
  
  # Calculate final Bayes factor
  bf01 <- prod_term * num_integral/den_integral
  
  return(list(
    bf10 = 1/bf01
  ))
}

# data set
  # if total N = 800, groupsize = 400
groupsize <- 400
mu1 <- 0.2
mu2 <- -0.2
sd1 <- 2
sd2 <- 2

#calaculate t statistic
sd_p <- pooled_sd(2, 2, groupsize, groupsize)
tstat <- (mu1 - mu2)/(sd_p * sqrt(1/groupsize+1/groupsize))

#JZS : BayesFactor
JZSstat707 <- ttest.tstat(tstat, n1=groupsize, n2=groupsize, rscale = 1/sqrt(2), simple=T)
print(JZSstat707)
JZSstat1 <- ttest.tstat(tstat, n1=groupsize, n2=groupsize, rscale = 1, simple=T)
print(JZSstat1)
#BeFi
BeFi <- gicaBF(mean1 = 0.2, mean2 = -0.2, var1 = sd1^2, var2 = sd2^2, n1 = groupsize, n2=groupsize)

# log10
log_JZS707 <- log10(JZSstat707)
log_JZS1 <- log10(JZSstat1)
log_BeFi <- log10(BeFi$bf10)
print(c(log_JZS707,log_JZS1, log_BeFi))

#for문 실행
# 데이터 분석 파라미터 설정
mu1 <- 0.2
mu2 <- -0.2
sd1 <- 2
sd2 <- 2

# 표본 크기 범위 설정 (50부터 1000까지 10단위로)
sample_sizes <- seq(25, 500, by = 5)  # 각 그룹의 크기 (total N은 두 배)

# 결과를 저장할 데이터프레임 초기화
results <- data.frame(
  N = sample_sizes * 2,  # 총 표본 크기
  JZSr707 = numeric(length(sample_sizes)),
  JZSr1 = numeric(length(sample_sizes)),
  BeFi = numeric(length(sample_sizes)),
  log_JZSr707 = numeric(length(sample_sizes)),
  log_JZSr1 = numeric(length(sample_sizes)),
  log_BeFi = numeric(length(sample_sizes))
)

# 각 표본 크기에 대해 베이즈 인자 계산
for (i in 1:length(sample_sizes)) {
  groupsize <- sample_sizes[i]
  
  # 풀드 표준편차 계산
  sd_p <- pooled_sd(sd1, sd2, groupsize, groupsize)
  
  # t 통계량 계산
  tstat <- (mu1 - mu2)/(sd_p * sqrt(1/groupsize+1/groupsize))
  
  # JZS Bayes 인자 (r=1/sqrt(2)) 계산
  JZSr707 <- ttest.tstat(tstat, n1=groupsize, n2=groupsize, rscale = 1/sqrt(2), simple=TRUE)
  
  # JZS Bayes 인자 (r=1) 계산
  JZSr1 <- ttest.tstat(tstat, n1=groupsize, n2=groupsize, rscale = 1, simple=TRUE)
  
  # BeFi Bayes 인자 계산
  BeFi <- gicaBF(mean1 = mu1, mean2 = mu2, var1 = sd1^2, var2 = sd2^2, n1 = groupsize, n2 = groupsize)
  
  # 결과 저장
  results$JZSr707[i] <- JZSr707
  results$JZSr1[i] <- JZSr1
  results$BeFi[i] <- BeFi$bf10
  
  # log10 변환 값 저장
  results$log_JZSr707[i] <- log10(JZSr707)
  results$log_JZSr1[i] <- log10(JZSr1)
  results$log_BeFi[i] <- log10(BeFi$bf10)
}

# 플롯팅을 위한 데이터 변환
plot_data <- data.frame(
  N = rep(results$N, 3),
  log10BF = c(results$log_JZSr707, results$log_JZSr1, results$log_BeFi),
  Method = factor(rep(c("JZS (r=1/√2)", "JZS (r=1)", "BeFi"), each = nrow(results)))
)

# 각 방법이 log10(BF10) = 0을 넘는 지점 찾기
find_crossing <- function(df, y_val = 0) {
  result <- list()
  for (method in unique(df$Method)) {
    method_data <- df[df$Method == method, ]
    for (i in 1:(nrow(method_data)-1)) {
      # 0을 교차하는 지점 검사
      if ((method_data$log10BF[i] <= y_val && method_data$log10BF[i+1] >= y_val) ||
          (method_data$log10BF[i] >= y_val && method_data$log10BF[i+1] <= y_val)) {
        
        # 선형 보간법으로 교차점 계산
        x1 <- method_data$N[i]
        x2 <- method_data$N[i+1]
        y1 <- method_data$log10BF[i]
        y2 <- method_data$log10BF[i+1]
        
        x_cross <- x1 + (y_val - y1) * (x2 - x1) / (y2 - y1)
        result[[method]] <- data.frame(N = x_cross, log10BF = y_val, Method = method)
        break
      }
    }
  }
  # 결과 결합
  do.call(rbind, result)
}

# 교차점 찾기
crossing_points <- find_crossing(plot_data, y_val = 0)  # 0을 넘는 지점 찾기

# 텍스트 위치 오프셋 설정 (레이블이 겹치지 않도록)
crossing_points$y_offset <- 0.07  # 기본 y 오프셋
crossing_points$x_offset <- 0     # 기본 x 오프셋

# 각 방법별로 위치 조정 (레이블이 겹치지 않도록)
# BeFi(salmon) 레이블은 오른쪽 아래에 위치
crossing_points$y_offset[crossing_points$Method == "BeFi"] <- -0.07
crossing_points$x_offset[crossing_points$Method == "BeFi"] <- 40

# JZS(r=1)(#18392B) 레이블은 왼쪽 위에 위치
crossing_points$x_offset[crossing_points$Method == "JZS (r=1)"] <- -10
crossing_points$y_offset[crossing_points$Method == "JZS (r=1)"] <- 0.1

# JZS(r=1/√2)(#1B7D4F) 레이블은 왼쪽에 위치
crossing_points$x_offset[crossing_points$Method == "JZS (r=1/√2)"] <- -60
crossing_points$y_offset[crossing_points$Method == "JZS (r=1/√2)"] <- 0.03

# 빈도주의 t-test의 p < 0.05, 80% power 지점 계산 (효과 크기 d = 0.2)
# 공식: N = 2 * (1.96 + 0.84)^2 / d^2, 여기서 d = (mu1 - mu2) / sd = 0.4 / 2 = 0.2
p_05_power80_n <- 2 * (1.96 + 0.84)^2 * 4 / 0.4^2

# 그래프 생성
plot <- ggplot(plot_data, aes(x = N, y = log10BF, color = Method, linetype = Method)) +
  geom_line(size = 1) +
  labs(
    title = "Sample Size Power Analysis with Bayes Factors",
    x = "Sample Size (N = n1 + n2)",
    y = "log10(BF10)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("salmon", "#18392B", "#1B7D4F")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgray") +  # BF = 10 (substantial evidence)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +  # BF = 3.16 (moderate evidence)
  geom_hline(yintercept = 0, color = "black") +  # BF = 1 (no evidence)
  scale_x_continuous(breaks = seq(50, 1000, by = 100)) +  # x축 눈금을 100 단위로 설정
  
  # 교차점 표시
  geom_point(data = crossing_points, aes(x = N, y = log10BF), size = 3) +
  
  # 교차점 좌표 텍스트 추가 (겹치지 않게 위치 조정)
  geom_text(data = crossing_points, 
            aes(x = N + x_offset, y = log10BF + y_offset, 
                label = sprintf("N = %d", round(N))),
            hjust = 0.5, vjust = 0.5) +
  
  # p < 0.05, 80% power 지점 표시
  geom_vline(xintercept = p_05_power80_n, linetype = "dotted", color = "darkviolet", size = 0.7) +
  
  # p < 0.05, 80% power 텍스트 추가
  annotate("text", x = p_05_power80_n-1.2, y = 0.95, 
           label = paste0("t-test p<.05 (80% power): N = ", round(p_05_power80_n)),
           color = "darkviolet", angle = 90, hjust = 1, vjust = -0.2, size = 3.5) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()  # 더 깔끔한 그리드를 위해 작은 그리드 제거
  )

print(plot)

ggsave("sample_power_with_ttest.png", plot = plot, dpi = 600, width = 17, height = 15, units = "cm")

# 결과 테이블 출력
head(results)
write.csv(results, file = "samplesize_power.csv")

# 특정 표본 크기에서의 결과 출력 (예: N=100, 200, 500, 1000)
interesting_ns <- c(50, 100, 300, 500, 800, 1000)
subset_results <- results[results$N %in% interesting_ns, ]
print(subset_results)
