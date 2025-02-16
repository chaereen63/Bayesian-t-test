  #### 효과크기 별 JZS와 GICA의 크기 비교 ####
  # 전체적인 패턴을 보기 위한 함수 (ratio로 비교)
  calculate_proportions_ratio <- function(data) {
    data %>%
      mutate(
        effect_size_cat = case_when(
          abs(std_mean_diff) <= 0.2 ~ "no/small_effect",    # 0-0.2
          abs(std_mean_diff) < 0.8 ~ "medium_effect",       # 0.2-0.8 
          TRUE ~ "large_effect"                             # >=0.8
        ),
        bf_ratio = BF_gica/BF_jzs
      ) %>%
      group_by(scenario, effect_size_cat) %>%
      summarise(
        prop_gica_smaller = mean(bf_ratio < 1) * 100,
        median_ratio = median(bf_ratio),
        n = n()
      )
  }
  
  results_df2 %>%
    calculate_proportions_ratio() %>%
    arrange(scenario, effect_size_cat)
  
  #BF10 직접 비교
  calculate_proportions_bf <- function(data) {
    data %>%
      mutate(
        effect_size_cat = case_when(
          abs(std_mean_diff) <= 0.2 ~ "no/small_effect",    
          abs(std_mean_diff) < 0.8 ~ "medium_effect",       
          TRUE ~ "large_effect"                             
        ),
        gica_smaller = BF_gica < BF_jzs  # 직접 BF 비교
      ) %>%
      group_by(scenario, effect_size_cat) %>%
      summarise(
        prop_gica_smaller = mean(gica_smaller) * 100,
        median_gica = median(BF_gica),
        median_jzs = median(BF_jzs),
        n = n()
      )
  }
  
  results_df2 %>%
    calculate_proportions_bf() %>%
    arrange(scenario, effect_size_cat)
  
  #조건을 더 세분화
  calcul_per_effect <- function(data) {
    data %>%
      mutate(
        effect_size_cat = case_when(
          abs(std_mean_diff) < 0.1 ~ "0.0-0.1",
          abs(std_mean_diff) < 0.2 ~ "0.1-0.2",
          abs(std_mean_diff) < 0.3 ~ "0.2-0.3",
          abs(std_mean_diff) < 0.4 ~ "0.3-0.4",
          abs(std_mean_diff) < 0.5 ~ "0.4-0.5",
          abs(std_mean_diff) < 0.6 ~ "0.5-0.6", 
          abs(std_mean_diff) < 0.7 ~ "0.6-0.7",
          abs(std_mean_diff) < 0.8 ~ "0.7-0.8",
          TRUE ~ ">=0.8"
        ),
        gica_smaller = BF_gica < BF_jzs  # 직접 BF 비교
      ) %>%
      group_by(scenario, effect_size_cat) %>%
      summarise(
        prop_gica_smaller = mean(gica_smaller) * 100,
        median_gica = median(BF_gica),
        median_jzs = median(BF_jzs),
        n = n()
      )
  }
  
  results_df2 %>%
    calcul_per_effect() %>%
    arrange(scenario, effect_size_cat) %>%
  print(n = 36)
  
#### calculate mean of log BF ####
  calculate_means <- function(data) {
    data %>%
      group_by(scenario, model) %>%
      summarise(
        mean_BF = mean(log_BF),
        se_BF = sd(log_BF) / sqrt(n()),
        .groups = 'drop'
      ) %>%
      arrange(scenario, model)
  }
  
  # 네 가지 효과크기 조건에 대해 각각 계산
  means_no_effect <- result_long %>% 
    filter(effect_size == "no effect") %>%
    calculate_means()
  
  means_weak <- result_long %>%
    filter(effect_size == "weak") %>% 
    calculate_means()
  
  means_medium <- result_long %>%
    filter(effect_size == "medium") %>%
    calculate_means()
  
  means_strong <- result_long %>%
    filter(effect_size == "strong") %>%
    calculate_means()
  
  # 결과를 하나의 데이터프레임으로 합치기
  all_means <- bind_rows(
    mutate(means_no_effect, effect = "no effect"),
    mutate(means_weak, effect = "weak"),
    mutate(means_medium, effect = "medium"), 
    mutate(means_strong, effect = "strong")
  ) %>%
    arrange(effect, scenario, model) %>% print(n=32)  
  
#### 소수의 극단값 확인 ####
  #### 생성 데이터 효과도 저장한 결과 ####
  results_df3 <- readRDS("./Korea_pub/final_merged_resultsK_dsd.RDS")
  results_df3 %>% select(BF_jzs, BF_gica, std_mean_diff, sd_x, sd_y, sdr, delta, scenario) -> results_df3
  
  results_df3 %>%
    filter(
      #(scenario == 1 & BF_gica/BF_jzs < 0.5) |    # 시나리오1: GICA > JZS인 경우 찾기
        (scenario == 2 & abs(BF_gica/BF_jzs - 1) > 0.5) |  # 시나리오2: 1 중심에서 크게 벗어난 경우
        (scenario == 3 & BF_gica/BF_jzs > 1.5) |     # 시나리오3: GICA > JZS인 경우 찾기 
        (scenario == 4 & BF_gica/BF_jzs < 0.5)      # 시나리오4: GICA < JZS인 경우 찾기
    ) %>%
    arrange(scenario, BF_gica/BF_jzs) %>%
    select(scenario, std_mean_diff, sd_x, sd_y, sdr, delta, BF_jzs, BF_gica, ratio = BF_gica/BF_jzs) %>%
    print(n=30)
  
  results_df3 %>%
    filter(
      (scenario == 2 & abs(BF_gica/BF_jzs - 1) > 0.5) |  # 시나리오2: 1 중심에서 크게 벗어난 경우
        (scenario == 3 & BF_gica/BF_jzs > 1.5) |   # 시나리오3: GICA > JZS인 경우 찾기   
        (scenario == 4 & BF_gica/BF_jzs < 0.5)   # 시나리오4: GICA < JZS인 경우 찾기
    ) %>%
    arrange(scenario, BF_gica/BF_jzs) %>%
    mutate(sd_ratio = sd_x/sd_y) %>%
    select(scenario, std_mean_diff, sd_x, sd_y, sd_ratio, sdr, delta, BF_jzs, BF_gica, bf_ratio = BF_gica/BF_jzs) %>%
    print(n=30)
  results_df3 %>% mutate(sd_sdr = sd_x/sd_y) -> sd_sdr
  