source("study3/functions3.R")
library(tidyverse)

set.seed(2026)  # seed 목록 재현용

# 시나리오 정의 (frequentist 세팅과 통일, var 4:4 또는 4:2)
# A: 1:1, 1:1 / B: 1:1, 2:1 / C: 2:1, 1:1 / D: 2:1, 2:1 (positive pairing) / E: 1:2, 2:1 (negative pairing)
scenarios <- list(
  # total sample size = 60
  list(condition = "A", n1 = 30, n2 = 30, sd1 = 2, sd2 = 2),
  list(condition = "B", n1 = 30, n2 = 30, sd1 = 2, sd2 = sqrt(2)),
  list(condition = "C", n1 = 40, n2 = 20, sd1 = 2, sd2 = 2),
  list(condition = "D", n1 = 40, n2 = 20, sd1 = 2, sd2 = sqrt(2)),
  list(condition = "E", n1 = 20, n2 = 40, sd1 = 2, sd2 = sqrt(2)),

  # total sample size = 120
  list(condition = "A", n1 = 60, n2 = 60, sd1 = 2, sd2 = 2),
  list(condition = "B", n1 = 60, n2 = 60, sd1 = 2, sd2 = sqrt(2)),
  list(condition = "C", n1 = 80, n2 = 40, sd1 = 2, sd2 = 2),
  list(condition = "D", n1 = 80, n2 = 40, sd1 = 2, sd2 = sqrt(2)),
  list(condition = "E", n1 = 40, n2 = 80, sd1 = 2, sd2 = sqrt(2)),

  # total sample size = 240
  list(condition = "A", n1 = 120, n2 = 120, sd1 = 2, sd2 = 2),
  list(condition = "B", n1 = 120, n2 = 120, sd1 = 2, sd2 = sqrt(2)),
  list(condition = "C", n1 = 160, n2 = 80,  sd1 = 2, sd2 = 2),
  list(condition = "D", n1 = 160, n2 = 80,  sd1 = 2, sd2 = sqrt(2)),
  list(condition = "E", n1 = 80,  n2 = 160, sd1 = 2, sd2 = sqrt(2))
)

# 설정 생성 함수
create_settings <- function(scenario, replications) {
  tibble(
    scenario = scenario,
    condition = scenarios[[scenario]]$condition,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    sd1 = scenarios[[scenario]]$sd1,
    sd2 = scenarios[[scenario]]$sd2,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 모든 설정 조합 생성
settings <- tibble(scenario = 1:15) %>%
  crossing(delta = c(0, 0.2, 0.5, 0.8)) %>%  # 효과크기 조건
  mutate(
    settings = map(scenario, ~create_settings(.x, replications = 500))
  ) %>%
  select(-scenario) %>%
  unnest(settings) %>%
  mutate(
    mu_diff = delta * mean_sd(sd1, sd2),  # 평균 분산으로 표준화 (frequentist 세팅과 동일)
    mean1 = mu_diff/2,
    mean2 = -mu_diff/2,
    sdr = sd2/sd1
  )

# 시나리오별 통계 확인
scenario_stats <- settings %>%
  group_by(scenario, delta) %>%
  summarise(
    condition = first(condition),
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    sdr = first(sdr),
    mean1 = first(mean1),  # 평균 확인용
    mean2 = first(mean2),  # 평균 확인용
    n_rep = n()
  )
print(scenario_stats, n=90)

# 결과 저장
saveRDS(settings, file = "study3/settings3.RDS")
