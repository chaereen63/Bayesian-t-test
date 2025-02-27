source(file = file.path("./Korea_pub/functionsK.R"))

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
  
  results_df3 %>%
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
  
  results_df3 %>%
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
  
  results_df3 %>%
    calcul_per_effect() %>%
    arrange(scenario, effect_size_cat) %>%
  print(n = 45)
  
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
    arrange(effect, scenario, model) %>% print(n=40)  
  
#### 소수의 극단값 확인 ####
  #### 생성 데이터 효과도 저장한 결과 ####
  # results_df <- readRDS("./Korea_pub/final_merged_resultsKbig.RDS")
  results_df2 %>% 
    mutate(std_mean_diff=mean_diff/mean_sd(sd_x,sd_y)) %>%
    select(BF_jzs, BF_gica, std_mean_diff, sd_x, sd_y, sdr, delta, scenario) -> results_df4
  
  # 1단계: 기본 필터링
  results_df4 %>%
    filter(
      (scenario == 2 & abs(BF_gica/BF_jzs - 1) > 0.5) |
        (scenario == 3 & abs(BF_gica/BF_jzs - 1) > 0.5) |
        (scenario == 4 & BF_gica/BF_jzs > 1.5) |
        (scenario == 5 & BF_gica/BF_jzs < 0.5)
    ) -> filtered_df
  
  # 2단계: 정렬
  filtered_df %>%
    arrange(scenario, BF_gica/BF_jzs) -> arranged_df
  
  # 3단계: BF_gica 확인
  head(arranged_df$BF_gica)
  
  # 4단계: 포맷팅 시도
  filtered_df %>%
    mutate(log_jzs = format(BF_jzs, scientific = TRUE),
           log_gica = format(BF_gica, scientific = TRUE),sd_ratio = sd_x/sd_y) %>%
    select(scenario, std_mean_diff, sd_x, sd_y, sd_ratio, sdr, delta, log_jzs, 
           log_gica, ratio = BF_gica/BF_jzs) %>%
    print(n = 40)
  
  
  #### sample size = 200 ####
  # 데이터가 results_df4에 저장되어 있다고 가정
    # 각 시나리오별 전체 케이스 수와 기준을 벗어난 케이스 수 계산
  summary_stats <- results_df4 %>%
    group_by(scenario) %>%
    summarise(
      total_cases = n(),
      outlier_cases = sum(case_when(
        scenario == 2 ~ abs(BF_gica/BF_jzs - 1) > 0.5, #시나리오 추가
        scenario == 3 ~ abs(BF_gica/BF_jzs - 1) > 0.5,
        scenario == 4 ~ BF_gica/BF_jzs > 1.5,
        scenario == 5 ~ BF_gica/BF_jzs < 0.5,
        TRUE ~ FALSE
      )),
      outlier_percentage = round(outlier_cases / total_cases * 100, 2)
    )
  
  # 결과 출력
  print(summary_stats)
  
  # 추가적으로 이상치의 특성을 파악하기 위한 요약 통계
  detailed_stats <- results_df4 %>%
    filter(
      case_when(
        scenario == 2 ~ abs(BF_gica/BF_jzs - 1) > 0.5,
        scenario == 3 ~ abs(BF_gica/BF_jzs - 1) > 0.5,
        scenario == 4 ~ BF_gica/BF_jzs > 1.5,
        scenario == 5 ~ BF_gica/BF_jzs < 0.5,
        TRUE ~ FALSE
      )
    ) %>%
    group_by(scenario) %>%
    summarise(
      mean_ratio = mean(BF_gica/BF_jzs),
      median_ratio = median(BF_gica/BF_jzs),
      sd_ratio = sd(BF_gica/BF_jzs),
      min_ratio = min(BF_gica/BF_jzs),
      max_ratio = max(BF_gica/BF_jzs)
    )
  
  print(detailed_stats)
  