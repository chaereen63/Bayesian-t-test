source("Korea_pub/functionsK.R")
library(tidyverse)

# 시나리오 정의
scenarios <- list(
  # total sample size = 30
  list(n1 = 15, n2 = 15, sd1 = 2, sd2 = 2),    # 1:1, SDR=1
  list(n1 = 15, n2 = 15, sd1 = 2, sd2 = 1),    # 1:1, SDR=0.5
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 2),    # 2:3, SDR=1
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 1),    # 2:3, SDR=0.5
  list(n1 = 12, n2 = 18, sd1 = 2, sd2 = 4),    # 2:3, SDR=2
  
  # total sample size = 100
  list(n1 = 50, n2 = 50, sd1 = 2, sd2 = 2),    
  list(n1 = 50, n2 = 50, sd1 = 2, sd2 = 1),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 2),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 1),    
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 4),    
  
  # total sample size = 200
  list(n1 = 100, n2 = 100, sd1 = 2, sd2 = 2),  
  list(n1 = 100, n2 = 100, sd1 = 2, sd2 = 1),  
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 2),   
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 1),   
  list(n1 = 80, n2 = 120, sd1 = 2, sd2 = 4)    
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
settings <- tibble(scenario = 1:15) %>%
  mutate(
    settings = map(scenario, ~create_settings(.x, replications = 500))
  ) %>%
  select(-scenario) %>%
  unnest(settings) %>%
  mutate(
    mean1 = 0,
    mean2 = 1,    # 평균 차이를 1로 고정
    sdr = sd2/sd1,  # 참조용 SDR 계산
    cohens_d = cohens_d(mean1, mean2, sd1, sd2)  # 참조용 Cohen's d
  )

# 시나리오별 통계 확인
scenario_stats <- settings %>%
  group_by(scenario) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    sdr = first(sdr),
    cohens_d = first(cohens_d)
  ) 

print(scenario_stats)

# 결과 저장
saveRDS(settings, file = "Korea_pub/settingsNew.RDS")
