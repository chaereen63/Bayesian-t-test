source("simulation/functions.R")
library(tidyverse)

# 시나리오 정의
scenarios <- list(
  list(n1 = 60, n2 = 40, rho_params = c(10, 10)),    # rho 평균 약 0.5
  list(n1 = 50, n2 = 50, rho_params = c(16, 8)),     # rho 평균 약 0.67
  list(n1 = 40, n2 = 60, rho_params = c(16, 8)),     # rho 평균 약 0.67
  list(n1 = 40, n2 = 60, rho_params = c(5, 10))      # rho 평균 약 0.33
)

# 효과 크기 설정
deltas <- c(0, 0.5)

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
    seed = sample.int(1e6, replications)
  )
}

# 모든 설정 조합 생성
settings <- crossing(
  scenario = 1:4,
  delta = deltas
) %>%
  mutate(
    settings = map2(scenario, delta, ~create_settings(.x, .y, replications = 1000))
  ) %>%
  unnest(settings, names_repair = "unique")

# rho 샘플링 및 표준편차 계산
settings <- settings %>%
  rowwise() %>%
  mutate(
    rho = rbeta(1, shape1 = rho_alpha, shape2 = rho_beta),
    var1 = 1 / rho,
    var2 = 1 / (1 - rho),
    sd1 = sqrt(var1),
    sd2 = sqrt(var2),
    sdr = sd2 / sd1
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
    rho_check = 1 / var1 / (1 / var1 + 1 / var2),
    delta_check = cohens_d(mean1, mean2, sd1, sd2, n1, n2)
  )

# 검증 결과 확인
all(abs(settings$rho - settings$rho_check) < 1e-10)
all(abs(settings$delta - settings$delta_check) < 1e-10)

# 결과 저장
saveRDS(settings, file = "simulation/settings.RDS")

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

ggsave("rho_distributions.png", width = 10, height = 8)