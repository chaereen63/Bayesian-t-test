source("Korea_pub/functionsK.R")
library(tidyverse)

# 시나리오 정의 (rho 값을 고정)
scenarios <- list(
  list(n1 = 40, n2 = 60, rho = 0.5),    # rho 고정값 0.5, SDR 1
  list(n1 = 50, n2 = 50, rho = 0.2),    # rho 고정값 0.2, SDR 2
  list(n1 = 40, n2 = 60, rho = 0.2),    # rho 고정값 0.2, SDR 2
  list(n1 = 40, n2 = 60, rho = 0.8)     # rho 고정값 0.8, SDR 0.5
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
  scenario = 1:4,
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
    pooled_sd = pooled_sd(sd1, sd2, n1, n2),
    mean1 = grand_mean + 0.5 * delta * pooled_sd,
    mean2 = grand_mean - 0.5 * delta * pooled_sd
  )

# 검증
settings <- settings %>%
  mutate(
    sdr_check = ((1-rho)/rho)^(1/2),
    delta_check = cohens_d(mean1, mean2, sd1, sd2, n1, n2)
  )

# 검증 결과 확인
cat("SDR 검증 결과:", all(abs(settings$sdr - settings$sdr_check) < 1e-10), "\n")
cat("delta 검증 결과:", all(abs(settings$delta - settings$delta_check) < 1e-10), "\n")

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
saveRDS(settings, file = "Korea_pub/settingsK.RDS")

# 요약 통계 출력
summary_stats <- settings %>%
  group_by(scenario, delta) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    rho = first(rho),
    sdr = first(sdr),
    .groups = 'drop'
  )

# 모든 열이 표시되도록 출력 옵션 설정
options(tibble.width = Inf)

# 요약 통계 출력
print(summary_stats)

# 출력 옵션 원래대로 복원
options(tibble.width = NULL)


# rho 값 시각화 (시나리오 별)
ggplot(unique(settings[, c("scenario", "rho")]), aes(x = factor(scenario), y = rho, color = factor(scenario))) +
  geom_point(size = 5) +
  geom_line(group = 1) +
  labs(title = "Fixed rho values across scenarios",
       x = "Scenario", y = "rho", color = "Scenario") +
  theme_minimal() +
  ylim(0, 1)  # rho는 0에서 1 사이의 값이므로

# ggsave("rho_values.png", width = 10, height = 8)

# 각 시나리오별 replication 수 확인
replication_counts <- settings %>%
  group_by(scenario) %>%
  summarise(count = n())

print("Replication counts per scenario:")
print(replication_counts)
# ggsave("rho_values.png", width = 10, height = 8)

# 시드 확인을 위한 출력
print(settings %>% select(scenario, delta, replication, seed))
