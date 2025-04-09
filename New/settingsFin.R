library(tidyverse)

# 시나리오 설정 및 복제본 생성을 위한 통합 함수
create_settings <- function(var1, var2, effect_size, n1, n2, replications = 500) {
  # 표준편차 계산
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  
  # 평균 분산 계산
  av_var <- mean(c(var1, var2))
  
  # 효과 크기에 따른 평균 차이 계산
  mu_diff <- effect_size * sqrt(av_var)
  
  # 각 집단의 평균 계산
  mu1 <- mu_diff/2
  mu2 <- -mu_diff/2
  
  # 복제본 생성
  tibble(
    var1 = var1,
    var2 = var2,
    sd1 = sd1,
    sd2 = sd2,
    effect_size = effect_size,
    av_var = av_var,
    mu_diff = mu_diff,
    mu1 = mu1,
    mu2 = mu2,
    n1 = n1,
    n2 = n2,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 각 시나리오 설정 생성 effect size = {0, 0.2, 0.5, 0.8}
set1 <- create_settings(var1 = 4, var2 = 4, effect_size = 0.2, n1 = 300, n2 = 300, replications = 500)
set2 <- create_settings(var1 = 4, var2 = 4, effect_size = 0.2, n1 = 240, n2 = 360, replications = 500)
set3 <- create_settings(var1 = 4, var2 = 2, effect_size = 0.2, n1 = 300, n2 = 300, replications = 500)
set4 <- create_settings(var1 = 4, var2 = 2, effect_size = 0.2, n1 = 240, n2 = 360, replications = 500)
set5 <- create_settings(var1 = 4, var2 = 2, effect_size = 0.2, n1 = 360, n2 = 240, replications = 500)

# 시나리오 번호 추가하고 합치기
set1 <- set1 %>% mutate(scenario = 1)
set2 <- set2 %>% mutate(scenario = 2)
set3 <- set3 %>% mutate(scenario = 3)
set4 <- set4 %>% mutate(scenario = 4)
set5 <- set5 %>% mutate(scenario = 5)

# 모든 세트 합치기
all_settings <- bind_rows(set1, set2, set3, set4, set5)

# 표준편차 비율 추가 (run_simulation 함수에서 사용됨)
all_settings <- all_settings %>% 
  mutate(varr = var2 / var1)

# 확인
print(nrow(all_settings)) # 5 * 500 = 2500
head(all_settings)
scenario_stats <- all_settings %>%
  group_by(scenario, effect_size) %>%  # mean_diff_condition도 그룹화에 추가
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    total_n = first(n1) + first(n2),
    sd1 = first(sd1),
    sd2 = first(sd2),
    varr = first(varr),
    mu1 = first(mu1),  # 평균 확인용
    mu2 = first(mu2),  # 평균 확인용
    mu_diff = first(mu_diff)
  ) 
print(scenario_stats)

# 저장
saveRDS(all_settings, file = "New/addset600ES2.RDS")
