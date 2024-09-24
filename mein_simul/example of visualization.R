library(tidyverse)

# 가상의 시뮬레이션 데이터 생성
set.seed(123)
n <- 1000

simulate_data <- function(n) {
  tibble(
    rho = runif(n, 0, 1),
    sdr = exp(rnorm(n, 0, 0.5)),  # log-normal distribution for SDR
    delta = sample(c(0, 0.5), n, replace = TRUE, prob = c(0.5, 0.5)),
    BF_M1 = exp(rnorm(n, 0, 1)),
    BF_M2 = exp(rnorm(n, 0, 1)),
    BF_M3 = exp(rnorm(n, 0, 1)),
    BF_M4 = exp(rnorm(n, 0, 1))
  ) %>%
    mutate(
      BF_M1 = ifelse(delta == 0, BF_M1 * 2, BF_M1 / 2),
      BF_M2 = ifelse(abs(rho - 0.5) > 0.2, BF_M2 * 2, BF_M2 / 2),
      BF_M3 = ifelse(delta == 0.5, BF_M3 * 2, BF_M3 / 2),
      BF_M4 = ifelse(delta == 0.5 & abs(rho - 0.5) > 0.2, BF_M4 * 2, BF_M4 / 2)
    )
}

results <- simulate_data(n)

# 데이터 정리
results_long <- results %>%
  pivot_longer(cols = starts_with("BF_"),
               names_to = "Model",
               values_to = "Bayes_Factor")

# rho에 따른 Bayes Factor 시각화
ggplot(results_long, aes(x = rho, y = Bayes_Factor, color = Model)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~delta, scales = "free_y", 
             labeller = labeller(delta = c("0" = "No Effect", "0.5" = "Medium Effect"))) +
  scale_y_log10() +
  labs(title = "Bayes Factors vs rho",
       x = "rho",
       y = "Bayes Factor (log scale)",
       color = "Model") +
  theme_minimal()

ggsave("bayes_factors_vs_rho.png", width = 10, height = 6)

# SDR에 따른 Bayes Factor 시각화
ggplot(results_long, aes(x = sdr, y = Bayes_Factor, color = Model)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  facet_wrap(~delta, scales = "free_y", 
             labeller = labeller(delta = c("0" = "No Effect", "0.5" = "Medium Effect"))) +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "Bayes Factors vs SDR",
       x = "SDR (log scale)",
       y = "Bayes Factor (log scale)",
       color = "Model") +
  theme_minimal()

ggsave("bayes_factors_vs_sdr.png", width = 10, height = 6)
