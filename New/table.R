library(dplyr)
library(tidyr)
library(BayesFactor)
library(knitr)
library(kableExtra)

# 모든 결과 파일 로드 함수
load_data <- function(expected_n_total, effect_sizes = c("0", "2", "5", "8")) {
  combined_data <- data.frame()
  
  for (es in effect_sizes) {
    # expected_n_total을 파일명에 사용 (50, 100, 200 등)
    file_path <- paste0("./New/mergedFin", expected_n_total, "ES", es, "_r1", ifelse(es == "5", "v2", ""), ".RDS")
    
    # 파일 읽기 시도
    tryCatch({
      results_df <- readRDS(file_path)
      
      # 데이터 처리
      results_long <- results_df %>%
        # 실제 n1+n2 합계로 표본크기 계산
        mutate(actual_n_total = n1 + n2) %>%
        # 로그 변환 적용
        mutate(
          log_BF_jzs = log10(BF_jzs),
          log_BF_gica = log10(BF_gica)
        ) %>%
        # long format으로 변환
        pivot_longer(
          cols = c(log_BF_jzs, log_BF_gica),
          names_to = "method",
          values_to = "log_BF"
        ) %>%
        # method 이름 간소화
        mutate(
          method = case_when(
            method == "log_BF_gica" ~ "BFGC",
            method == "log_BF_jzs" ~ "JZS",
            TRUE ~ method
          ),
          # 시나리오를 factor로 변환 및 설명적인 레이블 추가
          scenario = factor(as.character(scenario), 
                            levels = c("1", "2", "3", "4", "5"),
                            labels = c("등분산-균형", "등분산-불균형", "이분산-균형", 
                                       "이분산-불균형A", "이분산-불균형B")),
          # 효과크기 정보 추가
          effect_size = paste0("d=", ifelse(es == "0", "0.0", 
                                            ifelse(es == "2", "0.2", 
                                                   ifelse(es == "5", "0.5", "0.8")))),
          # 표본크기 정보 추가 - 실제 n1+n2 값 사용
          sample_size = paste0("N=", actual_n_total)
        ) %>%
        # 필요한 열만 선택
        select(scenario, method, log_BF, effect_size, sample_size, n1, n2)
      
      combined_data <- bind_rows(combined_data, results_long)
      
    }, error = function(e) {
      warning(paste("파일을 불러오는 중 오류 발생:", file_path, "\n", e))
    })
  }
  
  return(combined_data)
}

# 모든 표본크기 데이터 로드 및 결합
all_data <- bind_rows(
  load_data(50),   # 파일명에는 50이 있지만, 실제 n1+n2가 표본크기로 사용됨
  load_data(100),  # 파일명에는 100이 있지만, 실제 n1+n2가 표본크기로 사용됨
  load_data(200)   # 파일명에는 200이 있지만, 실제 n1+n2가 표본크기로 사용됨
)

# 효과크기와 표본크기를 factor로 변환하고 순서 지정
all_data <- all_data %>%
  mutate(
    effect_size = factor(effect_size, 
                         levels = c("d=0.0", "d=0.2", "d=0.5", "d=0.8")),
    sample_size = factor(sample_size, 
                         levels = c("N=50", "N=100", "N=200"))
  )

# JZS와 BFGC의 차이 계산
# 먼저 데이터에 고유 식별자를 추가하여 중복을 처리
all_data_with_id <- all_data %>%
  group_by(scenario, effect_size, sample_size, method) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

# 먼저 각 방법별 평균값 계산
avg_data <- all_data_with_id %>%
  group_by(scenario, effect_size, sample_size, method) %>%
  summarise(log_BF = mean(log_BF, na.rm = TRUE), .groups = "drop")

# 그 후 pivot_wider와 difference 계산
diff_data <- avg_data %>%
  pivot_wider(
    id_cols = c(scenario, effect_size, sample_size),
    names_from = method,
    values_from = log_BF
  ) %>%
  mutate(difference = BFGC - JZS) %>%
  pivot_longer(
    cols = c(JZS, BFGC, difference),
    names_to = "method",
    values_to = "log_BF"
  )

# 1. 시나리오별 요약 통계
scenario_summary <- diff_data %>%
  group_by(scenario, method) %>%
  summarise(
    mean = mean(log_BF, na.rm = TRUE),
    sd = sd(log_BF, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    id_cols = scenario,
    names_from = method,
    values_from = c(mean, sd),
    names_sep = "_"
  )

# 2. 효과크기별 요약 통계
effect_summary <- diff_data %>%
  group_by(effect_size, method) %>%
  summarise(
    mean = mean(log_BF, na.rm = TRUE),
    sd = sd(log_BF, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    id_cols = effect_size,
    names_from = method,
    values_from = c(mean, sd),
    names_sep = "_"
  )

# 3. 표본크기별 요약 통계
sample_summary <- diff_data %>%
  group_by(sample_size, method) %>%
  summarise(
    mean = mean(log_BF, na.rm = TRUE),
    sd = sd(log_BF, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    id_cols = sample_size,
    names_from = method,
    values_from = c(mean, sd),
    names_sep = "_"
  )

# 요약 통계를 결합하여 최종 표 구성
final_summary <- bind_rows(
  scenario_summary %>% mutate(category = "조건") %>% rename(factor = scenario),
  effect_summary %>% mutate(category = "효과크기") %>% rename(factor = effect_size),
  sample_summary %>% mutate(category = "표본크기") %>% rename(factor = sample_size)
) %>%
  select(category, factor, mean_JZS, sd_JZS, mean_BFGC, sd_BFGC, mean_difference)

# 최종 표 출력
kable(final_summary, digits = 3)

#### 베이지안 t-검정 함수 ####
run_bayes_ttest_improved <- function(data, grouping_var) {
  results <- data.frame(
    group = character(),
    bf10 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 각 그룹에 대해 t-검정 수행
  for (grp in unique(data[[grouping_var]])) {
    # 해당 그룹의 데이터 필터링
    group_data <- data %>% 
      filter(!!sym(grouping_var) == grp)
    
    # 각 방법에 대한 원시 데이터 추출 (개별 관측치)
    jzs_data <- group_data %>% 
      filter(method == "JZS") %>% 
      pull(log_BF)
    
    bfgc_data <- group_data %>% 
      filter(method == "BFGC") %>% 
      pull(log_BF)
    
    # 데이터가 충분한지 확인
    if (length(jzs_data) > 0 && length(bfgc_data) > 0) {
      # 베이지안 t-검정 수행 (paired=TRUE는 두 방법이 같은 데이터에 적용되었음을 의미)
      bf <- ttestBF(x = bfgc_data, y = jzs_data, paired = TRUE)
      
      # 결과 저장
      results <- rbind(results, data.frame(
        group = grp,
        bf10 = extractBF(bf)$bf,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}

# 시나리오별 베이지안 t-검정 (원본 데이터 사용)
scenario_ttests <- run_bayes_ttest_improved(all_data, "scenario")

# 효과크기별 베이지안 t-검정 (원본 데이터 사용)
effect_ttests <- run_bayes_ttest_improved(all_data, "effect_size")

# 표본크기별 베이지안 t-검정 (원본 데이터 사용)
sample_ttests <- run_bayes_ttest_improved(all_data, "sample_size")

# 모든 베이지안 t-검정 결과를 하나의 표로 통합
all_ttests <- bind_rows(
  scenario_ttests %>% mutate(category = "조건") %>% rename(factor = group),
  effect_ttests %>% mutate(category = "효과크기") %>% rename(factor = group),
  sample_ttests %>% mutate(category = "표본크기") %>% rename(factor = group)
)

# BF10 결과를 로그 변환하여 방향과 크기를 더 쉽게 해석할 수 있게 함
all_ttests <- all_ttests %>%
  mutate(
    log_bf10 = log10(bf10)
  )

print("종합 베이지안 t-검정 결과:")
print(all_ttests)

# 최종 결과표와 t-검정 결과를 통합 - difference의 sd 제외
final_result <- final_summary %>%
  left_join(
    all_ttests %>% select(category, factor, log_bf10),
    by = c("category", "factor")
  )

print("최종 통합 결과:")
print(final_result)
kable(final_result, digits = 3)

#### 표본크기별 효과크기 별 베이지안 t 검정 ####

# 데이터 로드 - 이미 로드되어 있다고 가정
results_es8 <- readRDS("./New/mergedFin50ES8_r1.RDS")  # 효과크기 0.8 파일

# 로그 베이즈 인자 계산
results_log <- results_es8 %>%
  mutate(
    log_BF_jzs = log10(BF_jzs),
    log_BF_gica = log10(BF_gica)
  )

# 각 시나리오별 분석
unique_scenarios <- unique(results_log$scenario)
cat("분석할 시나리오 개수:", length(unique_scenarios), "\n")

# 시나리오별 분석 함수
analyze_by_scenario <- function(data) {
  # 시나리오별 결과 저장할 데이터프레임
  scenario_results <- data.frame()
  
  # 각 시나리오별로 분석
  for (scen in unique(data$scenario)) {
    # 시나리오 데이터 필터링
    scenario_data <- data %>%
      filter(scenario == scen)
    
    # 로그 베이즈 인자 평균 계산
    mean_jzs <- mean(scenario_data$log_BF_jzs)
    mean_gica <- mean(scenario_data$log_BF_gica)
    mean_diff <- mean_gica - mean_jzs
    
    # 베이지안 paired t-test 수행
    t_result <- ttestBF(
      x = scenario_data$log_BF_gica,
      y = scenario_data$log_BF_jzs,
      paired = TRUE
    )
    
    # BF10 값 추출
    bf10 <- extractBF(t_result)$bf[1]
    
    # 결과 저장
    scenario_results <- rbind(
      scenario_results,
      data.frame(
        scenario = scen,
        mean_log_jzs = mean_jzs,
        mean_log_gica = mean_gica,
        mean_diff = mean_diff,
        bf10 = bf10,
        n_samples = nrow(scenario_data)
      )
    )
  }
  
  # NA나 무한값 처리하고 로그 변환만 수행
  scenario_results <- scenario_results %>%
    mutate(
      bf10_safe = ifelse(is.na(bf10) | is.infinite(bf10), .Machine$double.xmax, bf10),
      log_bf10 = log10(bf10_safe)
    )
  
  return(scenario_results)
}

# 각 시나리오별 분석 결과
for (scen in unique_scenarios) {
  cat("\n==== 시나리오", scen, "분석 결과 ====\n")
  
  # 시나리오 데이터 필터링
  scenario_data <- results_log %>%
    filter(scenario == scen)
  
  # 평균 계산
  mean_jzs <- mean(scenario_data$log_BF_jzs)
  mean_gica <- mean(scenario_data$log_BF_gica)
  mean_diff <- mean_gica - mean_jzs
  
  # paired t-test
  t_result <- ttestBF(
    x = scenario_data$log_BF_gica,
    y = scenario_data$log_BF_jzs,
    paired = TRUE
  )
  
  # BF10 값 추출
  bf10 <- extractBF(t_result)$bf[1]
  log_bf10 <- log10(bf10)
  
  # 결과 출력
  cat("로그 베이즈 인자 평균 (JZS):", round(mean_jzs, 3), "\n")
  cat("로그 베이즈 인자 평균 (GICA):", round(mean_gica, 3), "\n")
  cat("평균 차이 (GICA - JZS):", round(mean_diff, 3), "\n")
  cat("BF10:", bf10, "\n")
  cat("log10(BF10):", log_bf10, "\n")
  cat("시나리오 샘플 수:", nrow(scenario_data), "\n")
}

# 전체 시나리오 결과 테이블 생성
es8_results <- analyze_by_scenario(results_log)


# 전체 데이터에 대한 통합 분석
cat("\n==== 전체 데이터 통합 분석 결과 ====\n")
all_data_ttest <- ttestBF(
  x = results_log$log_BF_gica,
  y = results_log$log_BF_jzs,
  paired = TRUE
)

all_bf10 <- extractBF(all_data_ttest)$bf[1]
all_log_bf10 <- log10(all_bf10)

# 전체 결과 출력
cat("로그 베이즈 인자 평균 (JZS):", round(mean(results_log$log_BF_jzs), 3), "\n")
cat("로그 베이즈 인자 평균 (GICA):", round(mean(results_log$log_BF_gica), 3), "\n")
cat("평균 차이 (GICA - JZS):", round(mean(results_log$log_BF_gica) - mean(results_log$log_BF_jzs), 3), "\n")
cat("BF10:", all_bf10, "\n")
cat("log10(BF10):", all_log_bf10, "\n")
cat("전체 샘플 수:", nrow(results_log), "\n")

# 여러 효과크기 파일 분석 함수
analyze_multiple_es_files <- function(file_paths, es_labels) {
  # 결과 저장할 리스트
  all_results <- list()
  
  # 각 파일 분석
  for (i in 1:length(file_paths)) {
    # 파일 로드
    data <- readRDS(file_paths[i])
    
    # 로그 변환
    data_log <- data %>%
      mutate(
        log_BF_jzs = log10(BF_jzs),
        log_BF_gica = log10(BF_gica)
      )
    
    # 시나리오별 분석
    results <- analyze_by_scenario(data_log)
    
    # 효과크기 라벨 추가
    results$effect_size <- es_labels[i]
    
    # 결과 저장
    all_results[[i]] <- results
  }
  
  # 모든 결과 합치기
  combined_results <- bind_rows(all_results)
  
  return(combined_results)
}

file_paths <- c(
  "./New/mergedFin200ES0_r1.RDS",  # 효과크기 0.0
  "./New/mergedFin200ES2_r1.RDS",  # 효과크기 0.2
  "./New/mergedFin200ES5_r1v2.RDS",  # 효과크기 0.5
  "./New/mergedFin200ES8_r1.RDS"   # 효과크기 0.8
 )
 es_labels <- c("d=0.0", "d=0.2", "d=0.5", "d=0.8")
 
all_es_results <- analyze_multiple_es_files(file_paths, es_labels)
print(all_es_results)
