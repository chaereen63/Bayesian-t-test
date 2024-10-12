# Stan 모델 (이전과 동일)
stan_code <- "
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
"

# 데이터 준비 및 모델 실행
library(rstan)
set.seed(123)
y1 <- rnorm(30, 0, 1)
y2 <- rnorm(30, 0.5, 1.2)
data_list <- list(N1 = length(y1), N2 = length(y2), y1 = y1, y2 = y2)
fit <- stan(model_code = stan_code, data = data_list, iter = 5000, chains = 4)

# Savage-Dickey 밀도비 계산

# 사후 표본 추출
posterior_samples <- extract(fit)$delta

# 사후 밀도 추정
posterior_density <- density(posterior_samples)

# δ = 0에서의 사후 밀도 추정
posterior_at_zero <- approx(posterior_density$x, posterior_density$y, xout = 0)$y

# 사전 분포 (Cauchy(0, 1))에서 δ = 0일 때의 밀도
prior_at_zero <- dcauchy(0, 0, 1)

# 베이즈 팩터 계산 (BF10)
BF10 <- prior_at_zero / posterior_at_zero

print(paste("Bayes Factor (BF10):", BF10))

# 사후 분포 플롯
plot(posterior_density, main = "Posterior distribution of delta")
abline(v = 0, col = "red", lty = 2)

library(BayesFactor)
# BayesFactor 패키지를 사용한 베이즈 팩터 계산
bf_result <- ttestBF(y1, y2, rscale = 1)
BF10_BayesFactor <- exp(bf_result@bayesFactor$bf)
