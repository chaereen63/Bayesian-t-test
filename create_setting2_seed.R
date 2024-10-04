source("mein_simul/functions2.R")
library(tidyverse)

# 전체 시뮬레이션 seed 설정 옵션
use_global_seed <- TRUE  # TRUE 또는 FALSE로 설정
global_seed <- 603  # 원하는 숫자로 변경 가능

if (use_global_seed) {
  set.seed(global_seed)
}

# 시나리오 정의
  #replication =1로 설정하면서 표본크기 n_obs를 키움 ()
scenarios <- list(
  list(n1 = 400, n2 = 600, rho_params = c(2, 2)),    # rho 평균 약 0.5 (SDR ≈ 1)
  list(n1 = 500, n2 = 500, rho_params = c(8, 2)),     # rho 평균 약 0.2 (SDR ≈ 2)
  list(n1 = 400, n2 = 600, rho_params = c(8, 2)),     # rho 평균 약 0.2 (SDR ≈ 2)
  list(n1 = 400, n2 = 600, rho_params = c(2, 8))      # rho 평균 약 0.8 (SDR ≈ 0.5)
)

# 효과 크기 설정
deltas <- c(0, 0.8) #효과크기 좀 더 크게 변경 

# 설정 생성 함수
create_settings <- function(scenario, delta, replications) {
  tibble(
    scenario = scenario,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    rho_alpha = scenarios[[scenario]]$rho_params[1],
    rho_beta = scenarios[[scenario]]$rho_params[2],
    delta = delta,
    replication = 1:replications
  )
}

# 모든 설정 조합 생성
settings <- crossing(
  scenario = 1:4,
  delta = deltas
) %>%
  mutate(
    settings = map2(scenario, delta, ~create_settings(.x, .y, replications = 2)) #minimum size is 1, but set 2 for simulation loop check
  ) %>%
  select(-scenario, -delta
  ) %>%  # 원래의 scenario와 delta 열 제거
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
    sdr_check = ((1-rho)/rho)^(1/2),
    delta_check = cohens_d(mean1, mean2, sd1, sd2, n1, n2)
  )

# 검증 결과 확인
cat("SDR 검증 결과:", all(abs(settings$sdr - settings$sdr_check) < 1e-10), "\n")
cat("delta 검증 결과:", all(abs(settings$delta - settings$delta_check) < 1e-10), "\n")

#rho를 어덯게 검증할지는..? 아직 잘 모르겠음.

# 시나리오별 rho 평균 및 분산 확인
rho_stats <- settings %>%
  group_by(scenario) %>%
  summarise(
    mean_rho = mean(rho),
    var_rho = var(rho),
    expected_mean = rho_alpha / (rho_alpha + rho_beta),
    expected_var = (rho_alpha * rho_beta) / ((rho_alpha + rho_beta)^2 * (rho_alpha + rho_beta + 1))
  )

print(rho_stats)

# 결과 저장
# 결과 저장 시 설정 정보도 함께 저장
settings <- settings %>%
  mutate(
  use_global_seed = use_global_seed,
  global_seed = global_seed
)

saveRDS(settings, file = "mein_simul/settings2.RDS")

# 요약 통계 출력
summary_stats <- settings %>%
  group_by(scenario, delta) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    mean_rho = mean(rho),
    sd_rho = sd(rho),
    mean_sdr = mean(sdr),
    sd_sdr = sd(sdr),
    .groups = 'drop'
  )

# 모든 열이 표시되도록 출력 옵션 설정
options(tibble.width = Inf)

# 요약 통계 출력
print(summary_stats)

# 출력 옵션 원래대로 복원
options(tibble.width = NULL)

# rho 분포 시각화
ggplot(settings, aes(x = rho, fill = factor(scenario))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~scenario, scales = "free") +
  labs(title = "Distribution of rho across scenarios",
       x = "rho", y = "Density", fill = "Scenario") +
  theme_minimal()

# ggsave("rho_distributions.png", width = 10, height = 8)

  # create_setting2에서 sdr계산이 수치적으로 rho값과 비선형적인 관계를 가지고 있음. 
  # 나중에 필요할 경우 이 값을 어떻게 계산해서 사용할지 고려하기
