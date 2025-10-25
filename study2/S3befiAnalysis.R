#### Be-Fi simulation 통계량 시각화 ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(purrr)

results <- readRDS("./study2/results2/S2befi100_E5.RDS")

head(results)
names(results)

# 데이터를 long format으로 변환
results_long <- results %>%
  select(scenario, student_t, welch_t, bf_t) %>%
  pivot_longer(cols = c(student_t, welch_t, bf_t),
               names_to = "t_type",
               values_to = "t_value") %>%
  mutate(t_type = factor(t_type, 
                         levels = c("student_t", "welch_t", "bf_t"),
                         labels = c("Student's t", "Welch's t", "Behrens-Fisher t")))

# 시나리오별로 각각 개별 플롯 생성하기
for(i in 1:5) {
  p <- results_long %>%
    filter(scenario == i) %>%
    ggplot(aes(x = t_value, fill = t_type, color = t_type)) +
    geom_density(alpha = 0.6, size = 0.8) +
    facet_wrap(~t_type, scales = "free", ncol = 3) +
    labs(title = paste("Condition", i, "- t-sampling distribution"),
         x = "t value",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")) +
    scale_fill_manual(values = c("#BDBDBD", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#BDBDBD", "#08519C", "#56B4E9"))
  
  print(p)
}

# 또는 모든 시나리오를 하나의 플롯에 표시 (선택사항)
all_scenarios_plot <- results_long %>%
  ggplot(aes(x = t_value, fill = t_type, color = t_type)) +
  geom_density(alpha = 0.5, size = 0.8) +
  facet_grid(scenario ~ t_type, scales = "free",
             labeller = labeller(scenario = function(x) paste("Condition", x))) +
  labs(title = "t-sampling Distribution",
       x = "t value",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text = element_text(size = 9, face = "bold")) +
  scale_fill_manual(values = c("#BDBDBD", "#08519C","#56B4E9")) +
  scale_color_manual(values = c("#BDBDBD", "#08519C", "#56B4E9"))

print(all_scenarios_plot)

# 기본 통계량 출력
summary_stats <- results_long %>%
  group_by(scenario, t_type) %>%
  summarise(
    mean = mean(t_value, na.rm = TRUE),
    sd = sd(t_value, na.rm = TRUE),
    median = median(t_value, na.rm = TRUE),
    .groups = 'drop'
  )

print("각 시나리오별 t값 기본 통계량:")
print(summary_stats)

#### 이론적 분포 vs 표집분포 비교 시각화 ####

results <- readRDS("./study2/results2/S2befi100_E5.RDS")
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
      dist_name <- paste0("Noncentral t(", round(df_welch,1), ", δ=", round(ncp_welch,3), ") - Welch's df")
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
    
    scale_fill_manual(values = c("#BDBDBD", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#BDBDBD", "#08519C", "#56B4E9"))
  
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
                        params$n1+params$n2,", δ=", params$pop_effect_size,")"),
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
    
    scale_fill_manual(values = c("#BDBDBD", "#08519C","#56B4E9")) +
    scale_color_manual(values = c("#BDBDBD", "#08519C", "#56B4E9"))
  
  return(p)
}

# 통합 비교 플롯 생성
integrated_plot <- create_integrated_comparison()
print(integrated_plot)


cat("\n분석 완료!\n")
cat("- 검은 점선: 각 시나리오의 이론적 분포\n")
cat("- 색깔 영역: 시뮬레이션으로 얻은 실제 표집분포\n")
cat("- 시나리오 1,2: t분포가 이론적 기준\n")
cat("- 시나리오 3,4,5: Behrens-Fisher분포가 이론적 기준\n")