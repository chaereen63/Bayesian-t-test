library(dplyr)
library(tidyr)
library(gridExtra)
library(purrr)
library(ggplot2)

#### Analysis ####
results <- readRDS("./study2/results_all/S22results120_E5_2.RDS")
head(results)
names(results)

# analysis
# Rejection rate
analyze_simulation_results <- function(results) {
  
  # rejection rate 분석 (bfTest 추가)
  reject_summary <- results %>%
    group_by(scenario, n1, n2, var1, var2, sd1, sd2, pop_effect_size) %>%
    summarise(
      # Rejection rates (bfTest 추가)
      student_reject = mean(student_p < 0.05, na.rm = TRUE),
      welch_reject = mean(welch_p < 0.05, na.rm = TRUE),
      bf_reject = mean(bf_p < 0.05, na.rm = TRUE),  # bfTest rejection rate 추가
      
      # Rejection rate 신뢰구간 (bfTest 추가)
      student_reject_se = sqrt(student_reject * (1 - student_reject) / n()),
      welch_reject_se = sqrt(welch_reject * (1 - welch_reject) / n()),
      bf_reject_se = sqrt(bf_reject * (1 - bf_reject) / n()),  # bfTest SE 추가
      
      # 효과크기 요약 통계
      mean_cohens_d = mean(cohens_d, na.rm = TRUE),
      sd_cohens_d = sd(cohens_d, na.rm = TRUE),
      mean_hedges_g = mean(hedges_g, na.rm = TRUE),
      sd_hedges_g = sd(hedges_g, na.rm = TRUE),
      
      # t-통계량 요약 (bfTest 추가)
      mean_student_t = mean(student_t, na.rm = TRUE),
      sd_student_t = sd(student_t, na.rm = TRUE),
      mean_welch_t = mean(welch_t, na.rm = TRUE),
      sd_welch_t = sd(welch_t, na.rm = TRUE),
      mean_bf_t = mean(bf_t, na.rm = TRUE),  # bfTest t 통계량 평균
      sd_bf_t = sd(bf_t, na.rm = TRUE),      # bfTest t 통계량 표준편차
      
      # 자유도 및 R 값
      student_df = first(student_df),
      welch_df = first(welch_df),
      mean_bf_R = mean(bf_R, na.rm = TRUE),  # bfTest R 값 평균
      sd_bf_R = sd(bf_R, na.rm = TRUE),      # bfTest R 값 표준편차
      
      # 표본 크기
      n_sims = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      # 신뢰구간 계산 (bfTest 추가)
      student_reject_ci_lower = pmax(0, student_reject - 1.96 * student_reject_se),
      student_reject_ci_upper = pmin(1, student_reject + 1.96 * student_reject_se),
      welch_reject_ci_lower = pmax(0, welch_reject - 1.96 * welch_reject_se),
      welch_reject_ci_upper = pmin(1, welch_reject + 1.96 * welch_reject_se),
      bf_reject_ci_lower = pmax(0, bf_reject - 1.96 * bf_reject_se),      # bfTest CI
      bf_reject_ci_upper = pmin(1, bf_reject + 1.96 * bf_reject_se)       # bfTest CI
    )
  
  return(reject_summary)
}

# effect size
analyze_effect_size_distribution <- function(results) {
  effect_size_summary <- results %>%
    group_by(scenario, pop_effect_size) %>%
    summarise(
      # Cohen's d 분포
      cohens_d_mean = mean(cohens_d, na.rm = TRUE),
      cohens_d_median = median(cohens_d, na.rm = TRUE),
      cohens_d_sd = sd(cohens_d, na.rm = TRUE),
      cohens_d_q25 = quantile(cohens_d, 0.25, na.rm = TRUE),
      cohens_d_q75 = quantile(cohens_d, 0.75, na.rm = TRUE),
      
      # Hedges' g 분포  
      hedges_g_mean = mean(hedges_g, na.rm = TRUE),
      hedges_g_median = median(hedges_g, na.rm = TRUE),
      hedges_g_sd = sd(hedges_g, na.rm = TRUE),
      
      # 효과크기의 정확도 (bias)
      cohens_d_bias = cohens_d_mean - first(pop_effect_size),
      hedges_g_bias = hedges_g_mean - first(pop_effect_size),
      
      # 효과크기 추정의 정밀도 (precision)  
      cohens_d_rmse = sqrt(mean((cohens_d - first(pop_effect_size))^2, na.rm = TRUE)),
      hedges_g_rmse = sqrt(mean((hedges_g - first(pop_effect_size))^2, na.rm = TRUE)),
      
      .groups = 'drop'
    )
  
  return(effect_size_summary)
}

reject_analysis <- analyze_simulation_results(results)

# 결과 전치
transpose_results <- function(reject_analysis) {
  
  # 시나리오별로 열을 만들고 변수명을 행으로
  transposed <- reject_analysis %>%
    # 시나리오를 기준으로 열 만들기
    pivot_longer(cols = -scenario, names_to = "variable", values_to = "value") %>%
    pivot_wider(names_from = scenario, values_from = value, names_prefix = "Scenario_") %>%
    # 변수명을 첫 번째 열로
    relocate(variable, .before = everything())
  
  return(transposed)
}

reject_analysis_t <- transpose_results(reject_analysis)
reject_analysis_t <- reject_analysis_t %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
print(reject_analysis_t, n = Inf)

# 효과크기 분포 분석
effect_size_analysis <- analyze_effect_size_distribution(results)
effect_size_t <- transpose_results(effect_size_analysis)

print(effect_size_t, n=Inf)


#### 이론적 분포 vs 표집분포 비교 시각화 ####
# 필요한 함수들 로드 (제공된 함수들)
source("./study2/t-distribution.R")  # 이론적 분포 함수들이 포함된 파일
source("./study2/Behrens-Fisher_distribution.R")

# 시뮬레이션 결과에서 각 시나리오별 모수 정보 추출
scenario_params <- results %>%
  group_by(scenario) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    s1 = first(sd1),  # 모집단 표준편차
    s2 = first(sd2),  # 모집단 표준편차
    var1 = first(var1),  # 모집단 분산
    var2 = first(var2),  # 모집단 분산
    pop_mean1 = first(pop_mean1),
    pop_mean2 = first(pop_mean2),
    pop_effect_size = first(pop_effect_size),  # 모집단 효과크기 사용
    equal_var = abs(first(var1) - first(var2)) < 0.001,  # 등분산 여부 판단
    .groups = 'drop'
  ) %>%
  split(.$scenario) %>%
  map(~as.list(.))

# 모수 정보 출력
{cat("=== 시나리오별 모수 정보 ===\n")
for(i in 1:5) {
  params <- scenario_params[[i]]
  cat(sprintf("Condition %d: n1=%d, n2=%d, s1=%.1f, s2=%.1f, pop_effect_size=%.2f, equal_var=%s\n", 
              i, params$n1, params$n2, params$s1, params$s2, params$pop_effect_size, params$equal_var))
}
cat("\n")
}
# 이론적 분포 계산 함수 (Welch's t)
calculate_theoretical_distribution <- function(scenario, x_range) {
  params <- scenario_params[[scenario]]
  n1 <- params$n1
  n2 <- params$n2  
  s1 <- params$s1
  s2 <- params$s2
  pop_mean1 <- params$pop_mean1
  pop_mean2 <- params$pop_mean2
  equal_var <- params$equal_var
  pop_effect_size <- params$pop_effect_size
  
  # 평균차이 계산
  mean_diff <- pop_mean1 - pop_mean2
  
  if (equal_var) {
    # 시나리오 1, 2: 등분산 조건 (Student's t)
    df_pooled <- n1 + n2 - 2
    pooled_var <- ((n1-1)*s1^2 + (n2-1)*s2^2) / df_pooled
    se_pooled <- sqrt(pooled_var * (1/n1 + 1/n2))
    
    if (abs(pop_effect_size) < 0.001) {
      # Central t-distribution
      density_vals <- dt(x_range, df_pooled)
      dist_name <- paste0("Central t(", df_pooled, ")")
    } else {
      # Noncentral t-distribution
      ncp <- mean_diff / se_pooled
      density_vals <- dt(x_range, df_pooled, ncp=ncp)
      dist_name <- paste0("Noncentral t(", df_pooled, ", δ=", round(ncp,3), ")")
    }
  } else {
    # 시나리오 3, 4, 5: 이분산 조건 - Welch 방법 사용
    se_welch <- sqrt(s1^2/n1 + s2^2/n2)
    df_welch <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
    
    if (abs(pop_effect_size) < 0.001) {
      # Central t-distribution (Welch's df)
      density_vals <- dt(x_range, df_welch)
      dist_name <- paste0("Central t(", round(df_welch,1), ") - Welch's df")
    } else {
      # Noncentral t-distribution (Welch's df)
      ncp_welch <- mean_diff / se_welch
      density_vals <- dt(x_range, df_welch, ncp=ncp_welch)
      dist_name <- paste0("Noncentral t(", round(df_welch,1), ", d=", round(ncp_welch,3), ") - Welch's df")
    }
  }
  
  return(list(density = density_vals, name = dist_name))
}


# 각 시나리오별 비교 플롯 생성 함수
create_comparison_plot <- function(scenario_num) {
  
  # 해당 시나리오 데이터 필터링
  scenario_data <- results %>% 
    filter(scenario == scenario_num) %>%
    select(student_t, welch_t, bf_t) %>%
    pivot_longer(cols = everything(), names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t", "bf_t"),
                           labels = c("Student's t", "Welch's t", "Behrens-Fisher t")))
  
  # 해당 시나리오 모수 정보
  params <- scenario_params[[scenario_num]]
  
  # x축 범위 설정 (데이터 범위 기반)
  x_min <- min(scenario_data$t_value, na.rm = TRUE)
  x_max <- max(scenario_data$t_value, na.rm = TRUE)
  x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 300)
  
  # 이론적 분포 계산
  theoretical <- calculate_theoretical_distribution(scenario_num, x_range)
  
  # 이론적 분포 데이터프레임 생성
  theoretical_df <- data.frame(
    x = x_range,
    density = theoretical$density,
    type = "Theoretical"
  )
  
  # 기본 플롯 생성
  p <- ggplot() +
    # 시뮬레이션 결과 (표집분포)
    geom_density(data = scenario_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.6, size = 0.8) +
    # 이론적 분포
    geom_line(data = theoretical_df, 
              aes(x = x, y = density), 
              color = "black", size = 1.5, linetype = "dashed") +
    
    facet_wrap(~method, scales = "free", ncol = 3) +
    
    labs(title = paste0("Condition ", scenario_num, " - Theoretical vs Sampling Distribution\n", 
                        theoretical$name),
         subtitle = paste0("Parameters: n1=", params$n1, ", n2=", params$n2,
                           ", s1=", params$s1, ", s2=", params$s2,
                           ", δ=", params$pop_effect_size),
         x = "t value", 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    
    scale_fill_manual(values = c("#FFA500", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#FFA500", "#08519C", "#56B4E9"))
  
  return(p)
}

# 모든 시나리오에 대한 플롯 생성
for(i in 1:5) {
  p <- create_comparison_plot(i)
  print(p)
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
}

# 통합 비교 플롯
create_integrated_comparison <- function() {
  
  # 모든 시나리오 데이터 준비
  all_data <- results %>%
    select(scenario, student_t, welch_t, bf_t) %>%
    pivot_longer(cols = c(student_t, welch_t, bf_t), 
                 names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t", "bf_t"),
                           labels = c("Student's t", "Welch's t", "Behrens-Fisher")),
           scenario = paste("Condition", scenario))
  
  # 각 시나리오별 이론적 분포 계산
  theoretical_data <- data.frame()
  
  for(s in 1:5) {
    scenario_subset <- all_data %>% filter(scenario == paste("Condition", s))
    x_min <- min(scenario_subset$t_value, na.rm = TRUE)
    x_max <- max(scenario_subset$t_value, na.rm = TRUE)
    x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 200)
    
    theoretical <- calculate_theoretical_distribution(s, x_range)
    
    temp_df <- data.frame(
      x = x_range,
      density = theoretical$density,
      scenario = paste("Condition", s),
      dist_name = theoretical$name
    )
    
    theoretical_data <- rbind(theoretical_data, temp_df)
  }
  
  # 통합 플롯 생성
  p <- ggplot() +
    # 표집분포
    geom_density(data = all_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.5, size = 0.6) +
    
    # 이론적 분포
    geom_line(data = theoretical_data, 
              aes(x = x, y = density), 
              color = "black", size = 1, linetype = "dashed") +
    
    facet_grid(scenario ~ method, scales = "free") +
    
    labs(title = paste0("Total Conditions: Theoretical vs Empirical Distribution (N=", 
                        params$n1+params$n2,", d=", params$pop_effect_size,")"),
         subtitle = "Black dash: theoretical distribution, Color areas: empirical distribution",
         x = "t value", 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8)
    ) +
    
    scale_fill_manual(values = c("#FFA500", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#FFA500", "#08519C", "#56B4E9"))
  
  return(p)
}

# 통합 비교 플롯 생성
integrated_plot <- create_integrated_comparison()
print(integrated_plot)


####Theoretical vs. Empirical sampling distribution (within Methods)####
source("./study2/t-distribution.R")  # 이론적 분포 함수들이 포함된 파일
source("./study2/Behrens-Fisher_distribution.R")
results <- readRDS("./study2/results_all/S22results60_E8_2.RDS")

# 시뮬레이션 결과에서 각 시나리오별 모수 정보 추출
scenario_params <- results %>%
  group_by(scenario) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    s1 = first(sd1),  # 모집단 표준편차
    s2 = first(sd2),  # 모집단 표준편차
    var1 = first(var1),  # 모집단 분산
    var2 = first(var2),  # 모집단 분산
    pop_mean1 = first(pop_mean1),
    pop_mean2 = first(pop_mean2),
    pop_effect_size = first(pop_effect_size),  # 모집단 효과크기 사용
    equal_var = abs(first(var1) - first(var2)) < 0.001,  # 등분산 여부 판단
    .groups = 'drop'
  ) %>%
  split(.$scenario) %>%
  map(~as.list(.))

# 모수 정보 출력
{cat("=== 시나리오별 모수 정보 ===\n")
  for(i in 1:5) {
    params <- scenario_params[[i]]
    cat(sprintf("Condition %d: n1=%d, n2=%d, s1=%.1f, s2=%.1f, pop_effect_size=%.2f, equal_var=%s\n", 
                i, params$n1, params$n2, params$s1, params$s2, params$pop_effect_size, params$equal_var))
  }
  cat("\n")
  }

# 방법별 이론적 분포 계산 함수
calculate_method_specific_distribution <- function(scenario, method_name, x_range) {
  params <- scenario_params[[scenario]]
  n1 <- params$n1
  n2 <- params$n2  
  s1 <- params$s1
  s2 <- params$s2
  pop_mean1 <- params$pop_mean1
  pop_mean2 <- params$pop_mean2
  pop_effect_size <- params$pop_effect_size
  
  # 평균차이 계산
  mean_diff <- pop_mean1 - pop_mean2
  
  if (method_name == "Student's t") {
    # Student's t: Pooled variance 사용
    df_pooled <- n1 + n2 - 2
    pooled_var <- ((n1-1)*s1^2 + (n2-1)*s2^2) / df_pooled
    se_pooled <- sqrt(pooled_var * (1/n1 + 1/n2))
    
    if (abs(pop_effect_size) < 0.001) {
      density_vals <- dt(x_range, df_pooled)
      dist_name <- paste0("Central t(", df_pooled, ")")
    } else {
      ncp <- mean_diff / se_pooled
      density_vals <- dt(x_range, df_pooled, ncp=ncp)
      dist_name <- paste0("Noncentral t(", df_pooled, ", δ=", round(ncp,3), ")")
    }
    
  } else if (method_name %in% c("Welch's t", "Fisher's t")) {
    # Welch's t와 Behrens-Fisher: Welch의 자유도 사용
    se_welch <- sqrt(s1^2/n1 + s2^2/n2)
    df_welch <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
    
    if (abs(pop_effect_size) < 0.001) {
      density_vals <- dt(x_range, df_welch)
      dist_name <- paste0("Central t(", round(df_welch,1), ") - Welch")
    } else {
      ncp_welch <- mean_diff / se_welch
      density_vals <- dt(x_range, df_welch, ncp=ncp_welch)
      dist_name <- paste0("Noncentral t(", round(df_welch,1), ", δ=", round(ncp_welch,3), ") - Welch")
    }
  }
  
  return(list(density = density_vals, name = dist_name))
}


# 각 시나리오별 비교 플롯 생성 함수 (방법별 이론 분포)
create_comparison_plot <- function(scenario_num) {
  
  # 해당 시나리오 데이터 필터링
  scenario_data <- results %>% 
    filter(scenario == scenario_num) %>%
    select(student_t, welch_t, bf_t) %>%
    pivot_longer(cols = everything(), names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t", "bf_t"),
                           labels = c("Student's t", "Welch's t", "Fisher's t")))
  
  # 해당 시나리오 모수 정보
  params <- scenario_params[[scenario_num]]
  
  # 각 방법별로 이론적 분포 계산
  methods <- c("Student's t", "Welch's t", "Fisher's t")
  theoretical_data <- data.frame()
  
  for(m in methods) {
    method_data <- scenario_data %>% filter(method == m)
    x_min <- min(method_data$t_value, na.rm = TRUE)
    x_max <- max(method_data$t_value, na.rm = TRUE)
    x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 300)
    
    theoretical <- calculate_method_specific_distribution(scenario_num, m, x_range)
    
    temp_df <- data.frame(
      x = x_range,
      density = theoretical$density,
      method = m,
      dist_name = theoretical$name
    )
    
    theoretical_data <- rbind(theoretical_data, temp_df)
  }
  
  # 기본 플롯 생성
  p <- ggplot() +
    # 시뮬레이션 결과 (표집분포)
    geom_density(data = scenario_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.6, size = 0.8) +
    # 방법별 이론적 분포
    geom_line(data = theoretical_data, 
              aes(x = x, y = density), 
              color = "black", size = 1.5, linetype = "dashed") +
    
    facet_wrap(~method, scales = "free", ncol = 3) +
    
    labs(title = paste0("Condition ", scenario_num, " - Theoretical vs Sampling Distribution"),
         subtitle = paste0("Parameters: n1=", params$n1, ", n2=", params$n2,
                           ", s1=", params$s1, ", s2=", params$s2,
                           ", d=", params$pop_effect_size),
         x = "t value", 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    
    scale_fill_manual(values = c("#FFA500", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#FFA500", "#08519C", "#56B4E9"))
  
  return(p)
}

# 통합 비교 플롯 (방법별 이론 분포)
create_integrated_comparison <- function() {
  
  # 모든 시나리오 데이터 준비
  all_data <- results %>%
    select(scenario, student_t, welch_t, bf_t) %>%
    pivot_longer(cols = c(student_t, welch_t, bf_t), 
                 names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t", "bf_t"),
                           labels = c("Student's t", "Welch's t", "Fisher's t")),
           scenario = paste("Condition", LETTERS[scenario]))  # A to E
  
  # 각 시나리오 및 방법별 이론적 분포 계산
  theoretical_data <- data.frame()
  methods <- c("Student's t", "Welch's t", "Fisher's t")
  
  for(s in 1:5) {
    for(m in methods) {
      scenario_label <- paste("Condition", LETTERS[s])
      scenario_method_subset <- all_data %>% 
        filter(scenario == scenario_label, method == m)
      
      x_min <- min(scenario_method_subset$t_value, na.rm = TRUE)
      x_max <- max(scenario_method_subset$t_value, na.rm = TRUE)
      x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 200)
      
      theoretical <- calculate_method_specific_distribution(s, m, x_range)
      
      temp_df <- data.frame(
        x = x_range,
        density = theoretical$density,
        scenario = scenario_label,
        method = m,
        dist_name = theoretical$name
      )
      
      theoretical_data <- rbind(theoretical_data, temp_df)
    }
  }
  
  # method를 factor로 변환하여 순서 고정
  theoretical_data <- theoretical_data %>%
    mutate(method = factor(method, 
                           levels = c("Student's t", "Welch's t", "Fisher's t")))
  
  # scenario를 factor로 변환하여 순서 고정 (A, B, C, D, E)
  all_data <- all_data %>%
    mutate(scenario = factor(scenario, 
                             levels = paste("Condition", LETTERS[1:5])))
  
  theoretical_data <- theoretical_data %>%
    mutate(scenario = factor(scenario, 
                             levels = paste("Condition", LETTERS[1:5])))
  
  # 각 방법별로 이론적 분포 선 그리기
  theoretical_student <- theoretical_data %>% filter(method == "Student's t")
  theoretical_welch <- theoretical_data %>% filter(method == "Welch's t")
  theoretical_bf <- theoretical_data %>% filter(method == "Fisher's t")
  
  # 대표 파라미터 추출 (첫 번째 시나리오 기준)
  first_params <- scenario_params[[1]]
  total_n <- first_params$n1 + first_params$n2
  effect_size <- first_params$pop_effect_size
  
  # 통합 플롯 생성
  p <- ggplot() +
    # 표집분포
    geom_density(data = all_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.5, size = 0.6) +
    
    # 이론적 분포 - 각각 다른 색으로
    geom_line(data = theoretical_student, 
              aes(x = x, y = density), 
              color = "gray50", size = 1.2, linetype = "dashed") +
    geom_line(data = theoretical_welch, 
              aes(x = x, y = density), 
              color = "#0066CC", size = 1.2, linetype = "dashed") +
    geom_line(data = theoretical_bf, 
              aes(x = x, y = density), 
              color = "#0066CC", size = 1.2, linetype = "dashed") +
    
    facet_grid(scenario ~ method, scales = "free") +
    
    labs(title = bquote("All Conditions: Theoretical vs Empirical Sampling Distribution"),
         subtitle = bquote(atop("N=" * .(total_n) * ", " * italic(d) * "=" * .(effect_size) * 
                                  ", 30,000 replications | Dashed lines: theoretical, Shaded areas: empirical")), 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 8)
    ) +
    
    scale_fill_manual(values = c("Student's t" = "#FFA500", 
                                 "Welch's t" = "#08519C",
                                 "Fisher's t" = "#56B4E9")) +
    scale_color_manual(values = c("Student's t" = "#FFA500", 
                                  "Welch's t" = "#08519C",
                                  "Fisher's t" = "#56B4E9"))
  
  return(p)
}

# 통합 비교 플롯 생성
integrated_plot <- create_integrated_comparison()
print(integrated_plot)

ggsave("./study2/results_all/sup_plot_60E8_2.png", integrated_plot, width = 11, height = 7, dpi = 600)
