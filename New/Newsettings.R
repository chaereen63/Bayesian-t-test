source("New/functionsN.R")
library(tidyverse)

# 시나리오 정의
scenarios <- list(
  # total sample size = 30
  list(n1 = 15, n2 = 15, sd1 = 2, sd2 = 2),    
  list(n1 = 15, n2 = 15, sd1 = 2, sd2 = 1),    
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 2),    
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 1),    
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 4),
  list(n1 = 18, n2 = 12, sd1 = 2, sd2 = 2),
  list(n1 = 18, n2 = 12, sd1 = 2, sd2 = 1),
  list(n1 = 18, n2 = 12, sd1 = 2, sd2 = 4)
)

scenarios <- list(
  # total sample size = 100
  list(n1 = 50, n2 = 50, sd1 = 2, sd2 = 2),    
  list(n1 = 50, n2 = 50, sd1 = 2, sd2 = 1),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 2),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 1),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 4),
  list(n1 = 60, n2 = 40, sd1 = 2, sd2 = 2),
  list(n1 = 60, n2 = 40, sd1 = 2, sd2 = 1),
  list(n1 = 60, n2 = 40, sd1 = 2, sd2 = 4)
)
scenarios <- list(
  # total sample size = 200
  list(n1 = 100, n2 = 100, sd1 = 2, sd2 = 2),    
  list(n1 = 100, n2 = 100, sd1 = 2, sd2 = 1),    
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 2),    
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 1),    
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 4),
  list(n1 = 120, n2 = 80, sd1 = 2, sd2 = 2),
  list(n1 = 120, n2 = 80, sd1 = 2, sd2 = 1),
  list(n1 = 120, n2 = 80, sd1 = 2, sd2 = 4)
)

# 설정 생성 함수
create_settings <- function(scenario, replications) {
  tibble(
    scenario = scenario,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    sd1 = scenarios[[scenario]]$sd1,
    sd2 = scenarios[[scenario]]$sd2,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 모든 설정 조합 생성
# 여러 평균 차이 조건을 구현하려면 이런 방식으로 확장
settings <- tibble(scenario = 1:8) %>%
  crossing(mu_diff = 1) %>%  # 평균 차이 조건 추가 Moreno et al.(1999)
  mutate(
    settings = map(scenario, ~create_settings(.x, replications = 500))
  ) %>%
  select(-scenario) %>%
  unnest(settings) %>%
  mutate(
    mu1 = mu_diff/2,
    mu2 = -mu_diff/2,
    sdr = sd2/sd1,
    cohens_d = cohens_d(mu1, mu2, sd1, sd2),
    pooled_d = (mu1 - mu2) / pooled_sd(sd1, sd2, n1,  n2)
  )

# 시나리오별 통계 확인
scenario_stats <- settings %>%
  group_by(scenario, mu_diff) %>%  # mean_diff_condition도 그룹화에 추가
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    sdr = first(sdr),
    mu1 = first(mu1),  # 평균 확인용
    mu2 = first(mu2),  # 평균 확인용
    cohens_d = first(cohens_d),
    pooled_d = first(pooled_d)
  ) 
print(scenario_stats, n=90)

# 결과 저장
saveRDS(settings, file = "New/settingsm1_44.RDS")
