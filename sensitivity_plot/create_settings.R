library(tidyr)
library(dplyr)
home_dir <- "."

settings50 <- readRDS(file = file.path(home_dir, "/New/set_Fin50ES5.RDS")) 
settings100 <- readRDS(file = file.path(home_dir, "/New/set_Fin100ES5.RDS"))
settings200 <- readRDS(file = file.path(home_dir, "/New/set_Fin200ES5.RDS"))

# 시나리오별 통계 확인
scenario_stats50 <- settings50 %>%
  group_by(scenario, mu_diff) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    varr = first(varr),
    mu1 = first(mu1),
    mu2 = first(mu2)
  ) 

scenario_stats100 <- settings100 %>%
  group_by(scenario, mu_diff) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    varr = first(varr),
    mu1 = first(mu1),
    mu2 = first(mu2)
  ) 

scenario_stats200 <- settings200 %>%
  group_by(scenario, mu_diff) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    varr = first(varr),
    mu1 = first(mu1),
    mu2 = first(mu2)
  ) 
print(scenario_stats50);print(scenario_stats100);print(scenario_stats200)

# 각 데이터셋에서 시나리오 1만 필터링
scenario1_50 <- settings50 %>% filter(scenario == 1)
scenario1_100 <- settings100 %>% filter(scenario == 1)
scenario1_200 <- settings200 %>% filter(scenario == 1)

# 데이터셋에 샘플 크기 정보 추가
scenario1_50$sample_size <- 50
scenario1_100$sample_size <- 100
scenario1_200$sample_size <- 200

# 세 데이터셋 합치기
combined_scenario1 <- bind_rows(scenario1_50, scenario1_100, scenario1_200)

# 결과 확인
head(combined_scenario1)

# 각 통계 요약 데이터에서 시나리오 1만 필터링
scenario1_stats_50 <- scenario_stats50 %>% filter(scenario == 1)
scenario1_stats_100 <- scenario_stats100 %>% filter(scenario == 1)
scenario1_stats_200 <- scenario_stats200 %>% filter(scenario == 1)

# 샘플 크기 정보 추가
scenario1_stats_50$sample_size <- 50
scenario1_stats_100$sample_size <- 100
scenario1_stats_200$sample_size <- 200

# 세 데이터셋 합치기
combined_scenario1_stats <- bind_rows(scenario1_stats_50, scenario1_stats_100, scenario1_stats_200)

# 결과 확인
print(combined_scenario1_stats)
saveRDS(combined_scenario1, "sense5.RDS")
