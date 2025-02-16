# Box and Tiao (1973) 예제 데이터
x1 <- c(48.4, 48.0, 47.4, 47.8, 47.3, 46.8, 47.8, 48.3, 47.2, 46.7)
x2 <- c(47.2, 47.9, 47.5, 47.7, 47.5, 47.3, 47.2, 47.8, 47.8, 47.6)
set.seed(986)
x1 <- rnorm(20, mean = 50, sd = 12^0.5)
x2 <- rnorm(12, mean = 55, sd = 40^0.5)
samples <- c(x1,x2)
write.csv(samples, file = "samples.csv")
gicaBF(x1, x2)

fit_bf <- ttestBF(x = x1, y = x2)
# 사후분포 샘플 추출
posterior_jzs <- posterior(fit_bf, iterations = 1000)
# MCMC 객체를 데이터프레임으로 변환
posterior_df <- data.frame(delta = as.vector(posterior_jzs[, "delta"]))

# 시각화
ggplot(posterior_df, aes(x = delta)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "JZS Posterior Distribution",
       x = "Effect Size (δ)",
       y = "Density") +
  theme_minimal()

gicaMCMC <- function(x1, x2, posterior_samples = FALSE, n_iterations = 5000) {
  # Calculate sample statistics
  xbar1 <- mean(x1)
  xbar2 <- mean(x2)
  var1 <- var(x1)
  var2 <- var(x2)
  n1 <- length(x1)
  n2 <- length(x2)
  
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
  
  result <- list(
    bf01 = bf01,
    bf10 = 1/bf01,
    d = d,
    num_integral = num_integral,
    den_integral = den_integral,
    prod_term = prod_term
  )
  
  # MCMC sampling이 요청된 경우만 추가
  if(posterior_samples) {
    # Behrens-Fisher distribution parameters from Theorem 2
    m0 <- d                          # location: x̄2 - x̄1
    sigma0_sq <- var1/n1             # first component variance
    tau0_sq <- var2/n2              # second component variance
    total_scale <- sqrt(sigma0_sq + tau0_sq)  # total scale
    phi <- atan(sqrt(sigma0_sq/tau0_sq)) # angle parameter
    
    log_posterior <- function(delta) {
      # Prior 추가
      prior <- dnorm(delta, mean = m0, sd = total_scale, log = TRUE)
      
      # Likelihood 항
      term1 <- -(df1 + 1)/2 * log(1 + (delta^2)/(df1 * sigma0_sq))
      term2 <- -(df2 + 1)/2 * log(1 + ((delta - m0)^2)/(df2 * tau0_sq))
      
      return(prior + term1 + term2)
    }
    
    sample_posterior <- function(n_iter) {
      samples <- numeric(n_iter)
      current <- m0  # start from location parameter
      
      # 제안분포의 초기 스케일을 보수적으로 설정
      proposal_sd <- min(total_scale/10, abs(d)/10)
      
      # 더 긴 burn-in 기간
      burn_in <- floor(n_iter * 0.3)
      total_iter <- n_iter + burn_in
      recent_accepts <- numeric(50)
      
      for(i in 1:total_iter) {
        proposal <- rnorm(1, current, proposal_sd)
        log_ratio <- log_posterior(proposal) - log_posterior(current)
        
        if(log(runif(1)) < log_ratio) {
          current <- proposal
          recent_accepts[i %% 50 + 1] <- 1
        } else {
          recent_accepts[i %% 50 + 1] <- 0
        }
        
        if(i > burn_in) {
          samples[i - burn_in] <- current
          
          # 더 조심스러운 adaptive scaling
          if(i %% 50 == 0) {
            acc_rate <- mean(recent_accepts)
            if(acc_rate < 0.2) proposal_sd <- proposal_sd * 0.95
            if(acc_rate > 0.4) proposal_sd <- proposal_sd * 1.05
          }
        }
      }
      
      # 이상치 제거
      q1 <- quantile(samples, 0.025)
      q3 <- quantile(samples, 0.975)
      samples <- samples[samples >= q1 & samples <= q3]
      
      return(samples)
    }
    
    result$posterior_samples <- sample_posterior(n_iterations)
  }
  
  return(result)
}

# 함수 실행
fit <- gicaMCMC(x1, x2, posterior_samples = TRUE, n_iterations = 5000)
gicaRE <-  gicaBF(x1, x2)

# 사후분포 시각화
library(ggplot2)
ggplot(data.frame(delta = result_bf$posterior_samples), aes(x = delta)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = fit$d, linetype = "dashed", color = "red") +
  labs(title = "Behrens-Fisher Posterior Distribution",
       x = "Effect Size (δ)",
       y = "Density") +
  theme_minimal()


# 사전 분포
# Function to calculate Bayes factor for n=1
bf_n1 <- function(x1, x2, s1=1, s2=1) {
  # Using equation (1) from the paper with n1=n2=1
  d <- x2 - x1
  
  sqrt(2) * 
    integrate(function(delta) (1 + delta^2/s1^2)^(-1/2) * 
                (1 + (delta-d)^2/s2^2)^(-1/2), -Inf, Inf)$value /
    integrate(function(delta) (1 + 2*delta^2/s1^2)^(-1/2) * 
                (1 + 2*(delta-d)^2/s2^2)^(-1/2), -Inf, Inf)$value
}

# Calculate for different small differences
diffs <- c(0, 0.01, 0.05, 0.1)
bfs <- sapply(diffs, function(d) bf_n1(0, d))
posts <- 1/(1 + bfs)  # posterior prob of H1 with equal priors

# Print results
data.frame(
  difference = diffs,
  BF01 = bfs,
  post_H1 = posts
)

# JZS prior
library(BayesFactor)

# d=0인 경우
x1_0 <- c(0, 0.001)
x2_0 <- c(0, 0.001)
bf_0 <- ttestBF(x1_0, x2_0)

# d=0.01인 경우
x1_01 <- c(0, 0.001)
x2_01 <- c(0.01, 0.011)
bf_01 <- ttestBF(x1_01, x2_01)

# d=0.05인 경우
x1_05 <- c(0, 0.001)
x2_05 <- c(0.05, 0.051)
bf_05 <- ttestBF(x1_05, x2_05)

# d=0.1인 경우
x1_1 <- c(0, 0.001)
x2_1 <- c(0.1, 0.101)
bf_1 <- ttestBF(x1_1, x2_1)

# BF01로 변환하고 사후확률 계산
results <- data.frame(
  difference = c(0, 0.01, 0.05, 0.1),
  BF01 = c(1/exp(bf_0@bayesFactor$bf),
           1/exp(bf_01@bayesFactor$bf),
           1/exp(bf_05@bayesFactor$bf),
           1/exp(bf_1@bayesFactor$bf))
)
results$post_H1 <- 1/(1 + results$BF01)


# BayesFactor package
x1_0 <- c(0, 0.001)
x2_0 <- c(0, 0.001)
jzs_0 <- ttestBF(x1_0, x2_0)

x1_01 <- c(0, 0.001)
x2_01 <- c(0.01, 0.011)
jzs_01 <- ttestBF(x1_01, x2_01)

x1_05 <- c(0, 0.001)
x2_05 <- c(0.05, 0.051)
jzs_05 <- ttestBF(x1_05, x2_05)

x1_1 <- c(0, 0.001)
x2_1 <- c(0.1, 0.101)
jzs_1 <- ttestBF(x1_1, x2_1)

# Behrens-Fisher
bf_0 <- gicaBF(x1_0, x2_0)
bf_01 <- gicaBF(x1_01, x2_01)
bf_05 <- gicaBF(x1_05, x2_05)
bf_1 <- gicaBF(x1_1, x2_1)

# 결과 비교
results <- data.frame(
  difference = c(0, 0.01, 0.05, 0.1),
  JZS_BF01 = c(1/exp(jzs_0@bayesFactor$bf),
               1/exp(jzs_01@bayesFactor$bf),
               1/exp(jzs_05@bayesFactor$bf),
               1/exp(jzs_1@bayesFactor$bf)),
  BF_BF01 = c(bf_0$bf01, bf_01$bf01, bf_05$bf01, bf_1$bf01)
)

results$JZS_post_H1 <- 1/(1 + results$JZS_BF01)
results$BF_post_H1 <- 1/(1 + results$BF_BF01)
