library(dplyr)
library(tidyr)
library(effectsize)  # effectsize 패키지 로드

## 효과크기별 결과 파일 로드
#표본크기 별로 바꿔보기
results_es8 <- readRDS("./New/mergedFin200ES8_r1.RDS")  # 효과크기 0.8
results_es5 <- readRDS("./New/mergedFin200ES5_r1v2.RDS")  # 효과크기 0.5
results_es2 <- readRDS("./New/mergedFin200ES2_r1.RDS")  # 효과크기 0.2
results_es0 <- readRDS("./New/mergedFin200ES0_r1.RDS")  # 효과크기 0.0

# 1. 데이터 처리 함수 - long format으로 변환
process_bf_data <- function(results_df, effect_size_label) {
  results_long <- results_df %>%
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
      # method를 factor로 변환하고 JZS를 reference로 설정
      method = factor(method, levels = c("JZS", "BFGC")),
      # 시나리오를 factor로 변환
      scenario = factor(as.character(scenario), levels = c("1", "2", "3", "4", "5")),
      # 효과크기 정보 추가
      effect_size = effect_size_label
    ) %>%
    # 필요한 열만 선택
    select(scenario, method, log_BF, effect_size)
  
  return(results_long)
}

# 2. 대비(contrast) 설정 함수
set_contrasts <- function(data) {
  # 시나리오 레벨이 1-5까지 모두 있는지 확인
  if(all(levels(data$scenario) %in% c("1", "2", "3", "4", "5")) && 
     length(levels(data$scenario)) >= 5) {
    
    contrasts(data$scenario) <- cbind(
      "Comp1" = c(1, -1, 0, 0, 0),
      "Comp2" = c(1, 1, -2, 0, 0),
      "Comp3" = c(-2, -2, -2, 3, 3),
      "Comp4" = c(0, 0, 0, 1, -1)
      
    )
    
    return(data)
  } else {
    warning("시나리오 레벨이 1-5까지 모두 존재하지 않습니다.")
    return(data)
  }
}

# 3. 개별 효과크기 데이터 분석 함수 - effectsize 패키지 활용
analyze_data <- function(data, effect_size_label) {
  cat("\n==================================================\n")
  cat("효과크기", effect_size_label, "에 대한 분석 결과\n")
  cat("==================================================\n")
  
  # ANOVA 모델
  model <- aov(log_BF ~ scenario * method, data = data)
  cat("ANOVA 결과:\n")
  print(summary(model))
  
  # effectsize 패키지를 사용한 부분 에타 제곱 계산
  cat("\n부분 에타 제곱(Partial Eta Squared) 효과크기:\n")
  eta_squared_results <- eta_squared(model, partial = TRUE)
  print(eta_squared_results)
  
  # 계획된 대비(contrast) 분석
  cat("\n대비(contrast) 결과:\n")
  contrast_results <- summary.lm(model)
  print(contrast_results)
  
  # 방법별 기술통계
  stats <- data %>%
    group_by(scenario, method) %>%
    summarise(
      mean_log_BF = mean(log_BF, na.rm = TRUE),
      sd_log_BF = sd(log_BF, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  cat("\n방법별 기술통계:\n")
  print(stats)
  
  return(model)  # 모델 객체 반환 (필요시 사용)
}

# 4. 데이터 처리 및 변환
long_data0 <- process_bf_data(results_es0, "0.0")
long_data2 <- process_bf_data(results_es2, "0.2")
long_data5 <- process_bf_data(results_es5, "0.5")
long_data8 <- process_bf_data(results_es8, "0.8")

# 5. 대비(contrast) 설정
long_data0 <- set_contrasts(long_data0)
long_data2 <- set_contrasts(long_data2)
long_data5 <- set_contrasts(long_data5)
long_data8 <- set_contrasts(long_data8)
# check
contrasts(long_data8$scenario)
# 6. 개별 효과크기 데이터 분석
model0 <- analyze_data(long_data0, "0.0")
model2 <- analyze_data(long_data2, "0.2")
model5 <- analyze_data(long_data5, "0.5")
model8 <- analyze_data(long_data8, "0.8")

# 7. 모든 데이터셋 통합
combined_data <- bind_rows(long_data0, long_data2, long_data5, long_data8) %>%
  mutate(
    effect_size = factor(effect_size, levels = c("0.0", "0.2", "0.5", "0.8"))
  )
combined_data <- set_contrasts(combined_data)
# 8. 통합 데이터 분석 - effectsize 패키지 활용
full_model <- aov(log_BF ~ effect_size * scenario * method, data = combined_data)
cat("모든 효과크기를 통합한 3원 ANOVA 분석 결과\n")
print(summary(full_model))

# effectsize 패키지를 사용한 부분 에타 제곱 계산
eta_squared_full <- eta_squared(full_model, partial = TRUE)
cat("\n통합 데이터에 대한 부분 에타 제곱(Partial Eta Squared) 효과크기:\n")
print(eta_squared_full)

# 계획된 대비(contrast) 분석
cat("\n대비(contrast) 결과:\n")
contrast_results_full <- summary.lm(full_model)
print(contrast_results_full)

# 9. 효과크기별 방법별 기술통계
overall_stats <- combined_data %>%
  group_by(effect_size, method) %>%
  summarise(
    mean_log_BF = mean(log_BF, na.rm = TRUE),
    sd_log_BF = sd(log_BF, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

cat("\n효과크기별 방법별 기술통계:\n")
print(overall_stats)

# 10. 조건별 방법별 기술통계
overall_stats2 <- combined_data %>%
  group_by(scenario, method) %>%
  summarise(
    mean_log_BF = mean(log_BF, na.rm = TRUE),
    sd_log_BF = sd(log_BF, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

cat("\n조건별 방법별 기술통계:\n")
print(overall_stats2)
