library(effectsize)  # effectsize 패키지 로드
library(afex);library(emmeans)
library(tidyverse)
d50ES0 <- readRDS("./New/mergedFin50ES0_r1.RDS")
d50ES2 <- readRDS("./New/mergedFin50ES2_r1.RDS")
d50ES5 <- readRDS("./New/mergedFin50ES5_r1.RDS")
d50ES8 <- readRDS("./New/mergedFin50ES8_r1.RDS")
d100ES0 <- readRDS("./New/mergedFin100ES0_r1.RDS")
d100ES2 <- readRDS("./New/mergedFin100ES2_r1.RDS")
d100ES5 <- readRDS("./New/mergedFin100ES5_r1.RDS")
d100ES8 <- readRDS("./New/mergedFin100ES8_r1.RDS")
d200ES0 <- readRDS("./New/mergedFin200ES0_r1.RDS")
d200ES2 <- readRDS("./New/mergedFin200ES2_r1.RDS")
d200ES5 <- readRDS("./New/mergedFin200ES5_r1.RDS")
d200ES8 <- readRDS("./New/mergedFin200ES8_r1.RDS")
#### 0. ANOVA 분석을 위한 데이터 처리 ####
# 데이터셋 이름을 담은 벡터 생성
datasets <- c("d50ES0", "d50ES2", "d50ES5", "d50ES8", 
              "d100ES0", "d100ES2", "d100ES5", "d100ES8", 
              "d200ES0", "d200ES2", "d200ES5", "d200ES8")

# 각 데이터셋에 변환 적용
for (ds_name in datasets) {
  # 현재 환경에서 데이터셋 가져오기
  ds <- get(ds_name)
  
  # scenario 값의 범위 확인
  scenario_values <- unique(ds$scenario)
  print(paste("Dataset:", ds_name, "Scenario values:", paste(scenario_values, collapse=", ")))
  
  # 파일 이름에서 sample_size와 effect_size 추출
  sample_size_val <- as.numeric(sub("d(\\d+)ES.*", "\\1", ds_name))
  effect_size_val <- as.numeric(sub(".*ES(\\d+).*", "\\1", ds_name)) / 10
  
  # 원본 데이터에 새 변수 추가
  ds <- ds %>% 
    mutate(
      sample_size = sample_size_val,
      effect_size = effect_size_val,
      id = row_number()  # 각 행에 고유 ID 부여
    )
  
  # Wide에서 Long 형태로 변환 (factor 변환 전에 처리)
  ds_long <- ds %>%
    pivot_longer(
      cols = c(BF_jzs, BF_gica),
      names_to = "method",
      values_to = "BF"
    ) %>%
    mutate(
      method = case_when(
        method == "BF_jzs" ~ "JZS",
        method == "BF_gica" ~ "BFGC",
        TRUE ~ method
      ),
      log_BF = log10(BF)  # 밑이 10인 로그 변환
    )
  
  # 실제 scenario 값에 따라 factor 변환 적용
  ds_long <- ds_long %>%
    mutate(
      scenario = factor(scenario, 
                        levels = sort(unique(scenario)),  # 실제 데이터의 고유 값 사용
                        labels = LETTERS[1:length(unique(scenario))]),  # A, B, C, ... 로 레이블 지정
      sample_size = factor(sample_size, 
                           levels = c(50, 100, 200)),
      effect_size = factor(effect_size, 
                           levels = c(0, 0.2, 0.5, 0.8))
    ) %>%
    select(scenario, method, log_BF, effect_size, sample_size, id)
  
  # 원본 데이터셋에도 같은 방식으로 factor 변환 적용
  ds <- ds %>% 
    mutate(
      scenario = factor(scenario, 
                        levels = sort(unique(scenario)),
                        labels = LETTERS[1:length(unique(scenario))]),
      sample_size = factor(sample_size, 
                           levels = c(50, 100, 200)),
      effect_size = factor(effect_size, 
                           levels = c(0, 0.2, 0.5, 0.8))
    )
  
  # 변환된 원본 데이터셋 저장 (wide form)
  assign(ds_name, ds)
  
  # 변환된 long form 데이터셋 저장
  assign(paste0(ds_name, "_long"), ds_long)
  
  cat("Processed dataset:", ds_name, "and created", paste0(ds_name, "_long"), "\n")
}
# d200ES8_long %>% filter(method == "JZS") -> data_jzs2008
# d200ES8_long %>% filter(method =="BFGC") -> data_bfgc2008

#### 1. ANOVA 분석 ####
model <- aov(log_BF ~ scenario * method, data = d200ES8_long)
summary(model)
# 부분 에타 제곱 효과크기 계산
eta_squared(model, partial = TRUE)

# 2. 방법별 단순 효과 분석
simple <- emmeans(model, ~ scenario | method)
joint_tests(simple, by = "method")

# 3. Contrast 분석 - 직교비교 적용
# 직교비교 계수 설정
my_contrasts <- list(
  "Comp1" = c(3, 3, -2, -2, -2),
  "Comp2" = c(1, -1, 0, 0, 0),
  "Comp3" = c(0, 0, 2, -1, -1),
  "Comp4" = c(0, 0, 0, 1, -1)
)

ss_contrast <- emmeans(model, ~ scenario | method) %>%
  contrast(method = my_contrasts)
print(ss_contrast)

run_anova_analysis <- function(data_name) {
  # 데이터 가져오기
  data <- get(data_name)
  
  # 1. ANOVA 분석
  cat("\n\n================================================\n")
  cat("데이터셋:", data_name, "\n")
  cat("================================================\n\n")
  
  cat("1. ANOVA 분석 결과:\n")
  model <- aov(log_BF ~ scenario * method, data = data)
  model_summary <- summary(model)
  print(model_summary)
  
  # 부분 에타 제곱 효과크기 계산
  cat("\n2. 부분 에타 제곱 효과크기:\n")
  effect_size <- effectsize::eta_squared(model, partial = TRUE)
  print(effect_size)
  
  # 2. 방법별 단순 효과 분석
  cat("\n3. 방법별 단순 효과 분석:\n")
  simple <- emmeans(model, ~ scenario | method)
  simple_tests <- joint_tests(simple, by = "method")
  print(simple_tests)
  
  # 3. Contrast 분석 - 직교비교 적용
  # 직교비교 계수 설정
  my_contrasts <- list(
    "Comp1" = c(3, 3, -2, -2, -2),
    "Comp2" = c(1, -1, 0, 0, 0),
    "Comp3" = c(0, 0, 2, -1, -1),
    "Comp4" = c(0, 0, 0, 1, -1)
  )
  
  cat("\n4. Contrast 분석 결과:\n")
  ss_contrast <- emmeans(model, ~ scenario | method) %>%
    contrast(method = my_contrasts)
  print(ss_contrast)
  
  # 결과 객체 반환
  results <- list(
    model = model,
    model_summary = model_summary,
    effect_size = effect_size,
    simple_tests = simple_tests,
    contrasts = ss_contrast
  )
  
  return(results)
}

# 예시:
# results_d200ES8 <- run_anova_analysis("d200ES8_long")

# 모든 데이터셋에 대해 실행하기
run_all_analyses <- function() {
  dataset_patterns <- c(
    "d50ES0_long", "d50ES2_long", "d50ES5_long", "d50ES8_long",
    "d100ES0_long", "d100ES2_long", "d100ES5_long", "d100ES8_long",
    "d200ES0_long", "d200ES2_long", "d200ES5_long", "d200ES8_long"
  )
  
  all_results <- list()
  
  for (ds_name in dataset_patterns) {
    if (exists(ds_name)) {
      results <- run_anova_analysis(ds_name)
      all_results[[ds_name]] <- results
    } else {
      cat("데이터셋", ds_name, "이 존재하지 않습니다.\n")
    }
  }
  
  return(all_results)
}

# 모든 분석 실행 예시:
all_anova_results <- run_all_analyses()

#### 접근 방식 2: ANOVA 분석 함수를 수정하여 각 방법에 대해 별도로 분석 ####

run_anova_analysis_separate <- function(data_name) {
  # 데이터 가져오기
  data_long <- get(data_name)
  
  # JZS와 BFGC 데이터 분리
  data_jzs <- data_long %>% filter(method == "JZS")
  data_bfgc <- data_long %>% filter(method == "BFGC")
  
  # 1. JZS ANOVA 분석
  cat("\n\n================================================\n")
  cat("데이터셋:", data_name, "- JZS 분석\n")
  cat("================================================\n\n")
  
  cat("1. JZS ANOVA 분석 결과:\n")
  model_jzs <- aov(log_BF ~ scenario, data = data_jzs)
  model_summary_jzs <- summary(model_jzs)
  print(model_summary_jzs)
  
  # Contrast 분석 - 직교비교 적용
  # 직교비교 계수 설정
  my_contrasts <- list(
    "Comp1" = c(3, 3, -2, -2, -2),
    "Comp2" = c(1, -1, 0, 0, 0),
    "Comp3" = c(0, 0, 2, -1, -1),
    "Comp4" = c(0, 0, 0, 1, -1)
  )
  
  cat("\n3. JZS Contrast 분석 결과:\n")
  jzs_contrast <- emmeans(model_jzs, ~ scenario) %>%
    contrast(method = my_contrasts)
  print(jzs_contrast)
  
  # 2. BFGC ANOVA 분석
  cat("\n\n================================================\n")
  cat("데이터셋:", data_name, "- BFGC 분석\n")
  cat("================================================\n\n")
  
  cat("1. BFGC ANOVA 분석 결과:\n")
  model_bfgc <- aov(log_BF ~ scenario, data = data_bfgc)
  model_summary_bfgc <- summary(model_bfgc)
  print(model_summary_bfgc)

  # Contrast 분석
  cat("\n3. BFGC Contrast 분석 결과:\n")
  bfgc_contrast <- emmeans(model_bfgc, ~ scenario) %>%
    contrast(method = my_contrasts)
  print(bfgc_contrast)
  
  # 결과 객체 반환
  results <- list(
    jzs = list(
      model = model_jzs,
      model_summary = model_summary_jzs,
      contrasts = jzs_contrast
    ),
    bfgc = list(
      model = model_bfgc,
      model_summary = model_summary_bfgc,
      contrasts = bfgc_contrast
    )
  )
  
  return(results)
}

# 모든 데이터셋에 대해 분석 실행하는 함수
run_all_analyses_separate <- function() {
  dataset_patterns <- c(
    "d50ES0_long", "d50ES2_long", "d50ES5_long", "d50ES8_long",
    "d100ES0_long", "d100ES2_long", "d100ES5_long", "d100ES8_long",
    "d200ES0_long", "d200ES2_long", "d200ES5_long", "d200ES8_long"
  )
  
  all_results <- list()
  
  for (ds_name in dataset_patterns) {
    if (exists(ds_name)) {
      results <- run_anova_analysis_separate(ds_name)
      all_results[[ds_name]] <- results
    } else {
      cat("데이터셋", ds_name, "이 존재하지 않습니다.\n")
    }
  }
  
  return(all_results)
}

# 예시: 특정 데이터셋에 대한 분석
# results_d200ES8 <- run_anova_analysis_separate("d200ES8_long")

# 모든 데이터셋에 대한 분석 실행
all_anova_results_separate <- run_all_analyses_separate()
