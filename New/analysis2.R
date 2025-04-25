library(dplyr)
library(tidyr)
library(effectsize)  # effectsize 패키지 로드

## upload file
data <- read.csv("./New/analysisData.csv")
# 0. method 별 분리
data %>% filter(method == "JZS") -> data_jzs
data %>% filter(method =="BFGC") -> data_bfgc
data_jzs %>% mutate(scenario = factor(scenario, levels = 1:5), 
                    effect_size = factor(effect_size, levels = c("0", "0.2", "0.5", "0.8")),
                    sample_size = factor(sample_size, levels = c("50", "100","200"))) -> data_jzs
data_bfgc %>% mutate(scenario = factor(scenario, levels = 1:5), 
                     effect_size = factor(effect_size, levels = c("0", "0.2", "0.5", "0.8")),
                     sample_size = factor(sample_size, levels = c("50", "100","200"))) -> data_bfgc
summary(data_jzs);summary(data_bfgc)


# 1. 대비(contrast) 설정 함수
set_contrasts <- function(data) {
  # 시나리오 레벨이 1-5까지 모두 있는지 확인
  if(all(levels(data$scenario) %in% c("1", "2", "3", "4", "5")) && 
     length(levels(data$scenario)) >= 5) {
    
    contrasts(data$scenario) <- cbind(
      "Comp1" = c(3, 3, -2, -2, -2),
      "Comp2" = c(1, -1, 0, 0, 0),
      "Comp3" = c(0, 0, 2, -1, -1),
      "Comp4" = c(0, 0, 0, 1, -1)
      
    )
    
    return(data)
  } else {
    warning("시나리오 레벨이 1-5까지 모두 존재하지 않습니다.")
    return(data)
  }
}

# 3. 개별 효과크기 데이터 분석 함수 - effectsize 패키지 활용
analyze_data <- function(data, method) {
  cat("\n==================================================\n")
  cat(method, "에 대한 분석 결과\n")
  cat("==================================================\n")
  
  # ANOVA 모델
  model <- aov(log_BF ~ scenario * effect_size * sample_size, data = data)
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
    group_by(scenario, effect_size, sample_size) %>%
    summarise(
      mean_log_BF = mean(log_BF, na.rm = TRUE),
      sd_log_BF = sd(log_BF, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  cat("\n방법별 기술통계:\n")
  print(stats, n = 60)
  
  return(model)  # 모델 객체 반환 (필요시 사용)
}

# 4. 대비(contrast) 설정
data_jzs <- set_contrasts(data_jzs)
data_bfgc <- set_contrasts(data_bfgc)
# check
contrasts(data_jzs$scenario)
# 6. 개별 효과크기 데이터 분석
analyze_data(data_jzs, "JZS")
analyze_data(data_bfgc, "BFGC")


#### Change score model ####
# 데이터를 wide format으로 변환
# 각 그룹 내에서 번호 부여
indexed_data <- data %>%
  group_by(sample_size, effect_size, scenario, method) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

# wide format으로 변환
wide_data <- indexed_data %>%
  pivot_wider(
    id_cols = c(sample_size, effect_size, scenario, rep_id),
    names_from = method,
    values_from = log_BF
  )
wide_data %>% mutate(scenario = factor(scenario, levels = 1:5), 
            effect_size = factor(effect_size, levels = c("0", "0.2", "0.5", "0.8")),
            sample_size = factor(sample_size, levels = c("50", "100","200"))) -> wide_data

# 두 방법 간의 차이 계산
wide_data$diff_JZS_BFGC <- wide_data$JZS - wide_data$BFGC
# 대비 설정
wide_data <- set_contrasts(wide_data)
contrasts(wide_data$scenario)
# 효과크기와 시나리오에 따른 방법 간 차이 분석
diff_model <- aov(diff_JZS_BFGC ~ sample_size * effect_size * scenario, data = wide_data)
cat("\n방법 간 차이의 효과크기와 시나리오에 따른 분석:\n")
print(summary(diff_model))

# 계획된 대비 분석
cat("\n대비(contrast) 결과:\n")
contrast_results <- summary.lm(diff_model)
print(contrast_results)
