## 연습용 시각화
## analysis script 참고해서 새로 작성 (2024-10-07)

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)


# 데이터 로드
results_df <- readRDS("./mein_simul/merged_results.RDS")

# 데이터 변환
# bain BF가 역수라서 변환
results_df %>% 
  select("scenario", "true_model", "rho", "sdr", "delta", "BF_robtt_effect", "BF_bain_student", "BF_bain_welch", "BF_bayesfactor") %>%
  mutate(BF_bain_student = (BF_bain_student)^(-1), BF_bain_welch = (BF_bain_welch)^(-1)) ->
  result_inv

# 1. 생성된 BF 분포 그리기
# 데이터 변환
result_long <- result_inv %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "model",
    values_to = "BF"
  ) %>%
  mutate(
    model = gsub("BF_", "", model),  # str_remove 대신 gsub 사용
    log_BF = log10(BF)  # log10 변환 추가
  )

# 시각화 - Ridgeline plot
library(tidyverse)

ggplot(result_long, aes(x = log_BF, y = factor(scenario), fill = model)) +
  geom_density_ridges(
    alpha = 0.6, 
    scale = 0.9, 
    quantile_lines = TRUE, 
    quantiles = 2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.5) + #log BF = 0, BF=1인 지점
  facet_wrap(~ true_model, scales = "free") +  # scales = "free" 추가
  theme_minimal() +
  labs(
    title = "Distribution of log10(Bayes Factors) by Scenario, Model, and True Hypothesis",
    subtitle = "Scales are different for H0 and H1",
    x = "log10(Bayes Factor)",
    y = "Scenario",
    fill = "Model"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")

# 시각화 - Boxplot
ggplot(result_long, aes(x = factor(scenario), y = log_BF, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.9), outlier.alpha = 0.3) +
  facet_wrap(~ true_model, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Distribution of log10(Bayes Factors) by Scenario, Model, and True Hypothesis",
    x = "Scenario",
    y = "log10(Bayes Factor)",
    fill = "Model"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))


# 2. 막대그래프 그리기
# 데이터 가공 및 로그 변환된 BF 계산
plot_data <- result_inv %>%
  mutate(across(starts_with("BF_"), ~log10(.))) %>%
  pivot_longer(cols = starts_with("BF_"), 
               names_to = "method", 
               values_to = "log_BF") %>%
  group_by(scenario, true_model, method, sdr) %>%
  summarise(
    mean_log_BF = mean(log_BF),
    se_log_BF = sd(log_BF) / sqrt(n()),
    ci_lower = mean_log_BF - 1.96 * se_log_BF,
    ci_upper = mean_log_BF + 1.96 * se_log_BF,
    rmse = (se_log_BF / abs(mean_log_BF)) * 100
  )



ggplot(result_long, aes(x = sdr, y = log_BF, color = model, shape = factor(scenario))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ delta, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(x = "SDR (Standard Deviation Ratio)", 
       y = "Log Bayes Factor",
       title = "SDR vs Log Bayes Factor",
       subtitle = "Faceted by Effect Size (Delta), Scenarios as Shapes",
       color = "Model",
       shape = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = c(16, 17, 18, 15)) + # 각 시나리오에 대해 다른 모양 사용
  scale_x_continuous(trans = 'log', breaks = c(0.5, 1, 2, 5, 10)) + # x축을 로그 스케일로 변경
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") # y=0 선 추가

# 필요하다면 y축 범위를 조정할 수 있습니다
# + ylim(-2, 10)

# 그래프 생성 함수
create_plot <- function(data, scenario_num, delta_val) {
  ggplot(data, aes(x = sdr, y = log_BF, color = model)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = sprintf("Scenario %d, Effect Size (Delta) = %.1f", scenario_num, delta_val),
      x = "SDR (Standard Deviation Ratio)",
      y = "Log Bayes Factor",
      color = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# 각 시나리오와 delta에 대해 그래프 생성
plots <- result_long %>%
  group_by(scenario, delta) %>%
  group_map(~ create_plot(.x, .y$scenario, .y$delta))

# 모든 그래프를 하나의 PDF 파일로 저장
pdf("SDR_vs_LogBF_plots.pdf", width = 10, height = 8)
for (plot in plots) {
  print(plot)
}
dev.off()

# 콘솔에 첫 번째 그래프 출력 (확인용)
print(plots[[1]])

# RMSE 결과 출력
print(plot_data %>% 
        select(scenario, true_model, method, mean_log_BF, rmse, sdr) %>% 
        arrange(scenario, true_model, method))



######## ANOVA 분석 스크립트 (broom 패키지 없이) (2024-10-08) #######

library(dplyr)
library(tidyr)
library(knitr)

# 데이터 변환 및 준비
anova_data <- result_inv %>% 
  select(scenario, true_model, rho, sdr, delta, starts_with("BF_")) %>%
  pivot_longer(cols = starts_with("BF_"), 
               names_to = "method", 
               values_to = "BF") %>%
  mutate(
    log_BF = log10(BF),
    method = as.factor(method),
    true_model = as.factor(true_model),
    sample_size = case_when(
      scenario %in% c(1,3,4) ~ "n1=40, n2=60",
      scenario == 2 ~ "n1=50, n2=50"
    ),
    sample_size = as.factor(sample_size),
    sdr = as.factor(sdr)
  )

# H0 모델에 대한 ANOVA (모든 효과 강제 표시)
h0_model <- aov(log_BF ~ sample_size * sdr * method, data = anova_data %>% filter(true_model == "H0"))
print(summary.aov(h0_model, intercept = TRUE))

# H1 모델에 대한 ANOVA (모든 효과 강제 표시)
h1_model <- aov(log_BF ~ sample_size * sdr * method, data = anova_data %>% filter(true_model == "H1"))
print(summary.aov(h1_model, intercept = TRUE))


# 개선된 Tukey's HSD 테이블 출력
tukey_result_h0 <- TukeyHSD(h0_model, "method")
tukey_result_h1 <- TukeyHSD(h1_model, "method") 
print(tukey_result_h0); print(tukey_result_h1)

# 시각화: 방법별 평균 log BF 비교
plot_data <- anova_data %>%
  group_by(true_model, method, sample_size, sdr) %>%
  summarise(mean_log_BF = mean(log_BF), .groups = 'drop')

ggplot(plot_data, aes(x = method, y = mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(true_model ~ sdr + sample_size, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Comparison of Mean log10(Bayes Factors)",
       x = "Method", y = "Mean log10(Bayes Factor)")


# 1. SDR과 방법론의 상호작용 효과 시각화 (H1 상황)
interaction_plot <- anova_data %>%
  filter(true_model == "H1") %>%
  group_by(sdr, method) %>%
  summarise(mean_log_BF = mean(log_BF), .groups = "drop") %>%
  ggplot(aes(x = sdr, y = mean_log_BF, color = method, group = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Interaction between SDR and Method (H1)",
       x = "SDR", y = "Mean log10(Bayes Factor)")

print(interaction_plot)