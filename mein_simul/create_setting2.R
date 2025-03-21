source("mein_simul/functions2.R")
library(tidyverse)

# set.seed(1234)

# 시나리오 정의
# scenarios <- list(
#   list(n1 = 60, n2 = 40, rho_params = c(10, 10)),    # rho 평균 약 0.5
#   list(n1 = 50, n2 = 50, rho_params = c(16, 8)),     # rho 평균 약 0.67
#   list(n1 = 40, n2 = 60, rho_params = c(16, 8)),     # rho 평균 약 0.67
#   list(n1 = 40, n2 = 60, rho_params = c(5, 10))      # rho 평균 약 0.33
# )

# 시나리오 정의2
scenarios <- list(
  list(n1 = 60, n2 = 40, rho_params = c(2, 2)),    # rho 평균 약 0.5 (SDR ≈ 1)
  list(n1 = 50, n2 = 50, rho_params = c(8, 2)),     # rho 평균 약 0.2 (SDR ≈ 2)
  list(n1 = 40, n2 = 60, rho_params = c(8, 2)),     # rho 평균 약 0.2 (SDR ≈ 2)
  list(n1 = 40, n2 = 60, rho_params = c(2, 8))      # rho 평균 약 0.8 (SDR ≈ 0.5)
)

# 효과 크기 설정
deltas <- c(0, 0.3, 0.5, 0.8)

# 설정 생성 함수
create_settings <- function(scenario, delta, replications) {
  tibble(
    scenario = scenario,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    rho_alpha = scenarios[[scenario]]$rho_params[1],
    rho_beta = scenarios[[scenario]]$rho_params[2],
    delta = delta,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 모든 설정 조합 생성
settings <- crossing(
  scenario = 1:4,
  delta = deltas
) %>%
  mutate(
    settings = map2(scenario, delta, ~create_settings(.x, .y, replications = 1000)) #pilot으로 replications를 줄임
  ) %>%
  select(-delta,-scenario) %>%  # crossing에서 생성된 delta 열 제거
  unnest(settings)

# rho 샘플링 및 표준편차 계산
grand_sd <- 1  # 전체 표준편차 설정

settings <- settings %>%
  rowwise() %>%
  mutate(
    rho = rbeta(1, shape1 = rho_alpha, shape2 = rho_beta),
    sd1 = grand_sd * sqrt(2 / (1 + sqrt(rho / (1 - rho)))),
    sd2 = grand_sd * sqrt(2 / (1 + sqrt((1 - rho) / rho))),
    sdr = ((1-rho)/rho)^(1/2)
  ) %>%
  ungroup()

# 그룹별 평균 계산
grand_mean <- 0
settings <- settings %>%
  mutate(
    pooled_sd = pooled_sd(sd1, sd2, n1, n2),
    mean1 = grand_mean + 0.5 * delta * pooled_sd,
    mean2 = grand_mean - 0.5 * delta * pooled_sd
  )

# 검증
settings <- settings %>%
  mutate(
    rho_check = (1/sd1^2) / (1/sd1^2 + 1/sd2^2),
    delta_check = cohens_d(mean1, mean2, sd1, sd2, n1, n2)
  )

summary_stats %>%  mutate( rho_check= 1 / (1 + mean_sdr^2)) -> check
all(abs(check$mean_rho - check$rho_check) < 1e-6)  # 허용 오차를 조금 더 크게 설정

# 검증 결과 확인
all(abs(settings$rho - settings$rho_check) < 1e-10)
all(abs(settings$delta - settings$delta_check) < 1e-10)

check %>%
  mutate(rho_check = 1 / (1 + mean_sdr^2),
         diff = abs(mean_rho - rho_check)) %>%
  select(scenario, delta, mean_rho, mean_sdr, rho_check, diff)

# 결과 저장
saveRDS(settings, file = "mein_simul/settings_delta.RDS")

# 요약 통계 출력
summary_stats <- settings %>%
  group_by(scenario, delta) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    mean_rho = mean(rho),
    sd_rho = sd(rho),
    mean_sdr = mean(sdr),
    sd_sdr = sd(sdr)
  )

print(summary_stats)

# rho 분포 시각화
ggplot(settings, aes(x = rho, fill = factor(scenario))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~scenario, scales = "free") +
  labs(title = "Distribution of rho across scenarios",
       x = "rho", y = "Density", fill = "Scenario") +
  theme_minimal()

# ggsave("rho_distributions.png", width = 10, height = 8)
