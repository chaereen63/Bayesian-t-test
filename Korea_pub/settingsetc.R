source("Korea_pub/functionsK.R")
library(tidyverse)

# 시나리오 정의 (rho 값을 고정)
scenarios <- list(
  list(n1 = 15, n2 = 15, rho = 0.5),    # rho 고정값 0.5, SDR 1
  list(n1 = 15, n2 = 15, rho = 0.5),    # SDR 1
  list(n1 = 15, n2 = 15, rho = 0.2),    # rho 고정값 0.2, SDR 2
  list(n1 = 12, n2 = 18, rho = 0.2),    # rho 고정값 0.2, SDR 2
  list(n1 = 12, n2 = 18, rho = 0.8)     # rho 고정값 0.8, SDR 0.5
)

# 효과 크기 설정
deltas <- c(0, 0.2, 0.5, 0.8) #효과크기 좀 더 크게 변경 

# 설정 생성 함수
create_settings <- function(scenario, delta, replications) {
  tibble(
    scenario = scenario,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    rho = scenarios[[scenario]]$rho,
    delta = delta,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)  # 각 replication에 대한 고유한 시드 생성
  )
}

# 모든 설정 조합 생성
settings <- crossing(
  scenario = 1:5,
  delta = deltas
) %>%
  mutate(
    settings = map2(scenario, delta, ~create_settings(.x, .y, replications = 500)) #replication 필요하면 수정하기
  ) %>%
  select(-scenario, -delta) %>%  # 원래의 scenario와 delta 열 제거
  unnest(settings)

# 표준편차 계산
grand_sd <- 1  # 전체 표준편차 설정
settings <- settings %>%
  mutate(
    sd1 = grand_sd * sqrt(2 / (1 + sqrt(rho / (1 - rho)))),
    sd2 = grand_sd * sqrt(2 / (1 + sqrt((1 - rho) / rho))),
    sdr = ((1-rho)/rho)^(1/2)
  )

# 그룹별 평균 계산
grand_mean <- 0
settings <- settings %>%
  mutate(
    BeFi_sd = BeFi_sd(sd1, sd2, n1, n2),
    mean1 = grand_mean + 0.5 * delta * BeFi_sd,
    mean2 = grand_mean - 0.5 * delta * BeFi_sd
  )

# 시나리오별 rho 확인
rho_stats <- settings %>%
  group_by(scenario) %>%
  summarise(
    rho = first(rho)
  )

print(rho_stats)

# 시나리오별 sdr 확인
sdr_stats <- settings %>%
  group_by(scenario) %>%
  summarise(
    sdr = first(sdr)
  )

print(sdr_stats)

# 결과 저장
saveRDS(settings, file = "Korea_pub/settingsK30.RDS")
