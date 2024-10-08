## 연습용 시각화
## analysis script 참고해서 새로 작성 (2024-10-07)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # 그래프를 쉽게 결합할 수 있는 패키지입니다.

# 데이터 로드
results_df <- readRDS("./post/final_merged_results.RDS")

# 데이터 변환
# bain BF가 역수라서 변환
results_df %>% 
  select("scenario", "true_model", "rho", "sdr", "delta", "BF_robtt_effect", "BF_bain_student", "BF_bain_welch", "BF_bayesfactor") %>%
  mutate(BF_bain_student = (BF_bain_student)^(-1), BF_bain_welch = (BF_bain_welch)^(-1)) ->
  result_inv

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

# 시나리오에 대한 레이블 생성 (SDR 포함)
# 시나리오 정의 (rho 값을 고정)
scenarios <- list(
  list(n1 = 40, n2 = 60, rho = 0.5),    # rho 고정값 0.5
  list(n1 = 50, n2 = 50, rho = 0.2),    # rho 고정값 0.2
  list(n1 = 40, n2 = 60, rho = 0.2),    # rho 고정값 0.2
  list(n1 = 40, n2 = 60, rho = 0.8)     # rho 고정값 0.8
)

scenario_labels <- result_inv %>%
  group_by(scenario) %>%
  summarise(sdr = first(sdr), 
            n1 = ifelse(first(scenario) %in% c(1,3,4), 40, 50),
            n2 = ifelse(first(scenario) %in% c(1,3,4), 60, 50)) %>%
  mutate(label = sprintf("n1=%d, n2=%d, SDR=%.2f", n1, n2, sdr)) %>%
  pull(label)

# 일관된 색상 팔레트 정의
method_colors <- c(
  "BF_bain_student" = "#95CAE4",  # 파스텔 파란색 (등분산 가정)
  "BF_bayesfactor" = "#7AC5CD",   # 파스텔 청록색 (등분산 가정)
  "BF_bain_welch" = "#FFAEB9",    # 파스텔 분홍색 (이분산 가정)
  "BF_robtt_effect" = "#FFB347"   # 파스텔 주황색 (이분산 가정)
)

# H0 모델에 대한 바 플롯 생성
H0 <- plot_data %>%
  filter(true_model == "H0") %>%
  ggplot(aes(x = method, y = -mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_errorbar(aes(ymin = -ci_lower, ymax = -ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", abs(x))) +
  labs(title = "Comparison of log10(Bayes Factors) for H0 across Scenarios and Methods",
       subtitle = "Error bars represent 95% confidence intervals",
       x = "Method", y = "Absolute log10(Bayes Factor)", fill = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

# H1 모델에 대한 바 플롯 생성
H1 <- plot_data %>%
  filter(true_model == "H1") %>%
  ggplot(aes(x = method, y = mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 3.5)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Comparison of log10(Bayes Factors) for H1 across Scenarios and Methods",
       subtitle = "Error bars represent 95% confidence intervals",
       x = "Method", y = "log10(Bayes Factor)", fill = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10))

# RMSE 결과 출력
print(plot_data %>% 
        select(scenario, true_model, method, mean_log_BF, rmse, sdr) %>% 
        arrange(scenario, true_model, method))

# 그래프 출력 (H0와 H1 그래프를 위아래로 배치)
print(H0)
print(H1)

# 동일 방법의 시나리오 별 비교를 위한 데이터 준비
plot_data <- plot_data %>%
  mutate(scenario_label = factor(scenario, 
                                 levels = 1:4, 
                                 labels = c("n1=40, n2=60, SDR=1.00",
                                            "n1=50, n2=50, SDR=2.00",
                                            "n1=40, n2=60, SDR=2.00",
                                            "n1=40, n2=60, SDR=0.50")))

# H0에 대한 그래프 생성 (y축 반전)
plot_H0 <- ggplot(plot_data %>% filter(true_model == "H0"), 
                  aes(x = scenario_label, y = -mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 0.8)) +
  geom_errorbar(aes(ymin = -ci_upper, ymax = -ci_lower),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ method, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", abs(x))) +
  labs(title = "Comparison of log10(Bayes Factors) for H0 Across Scenarios",
       subtitle = "Error bars represent 95% confidence intervals\nY-axis inverted for consistency with H1 graph",
       x = "Scenario", y = "Absolute log10(Bayes Factor)", fill = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

# H1에 대한 그래프 생성
plot_H1 <- ggplot(plot_data %>% filter(true_model == "H1"), 
                  aes(x = scenario_label, y = mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 3.5)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  facet_wrap(~ method, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Comparison of log10(Bayes Factors) for H1 Across Scenarios",
       subtitle = "Error bars represent 95% confidence intervals",
       x = "Scenario", y = "log10(Bayes Factor)", fill = "Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"))

# 그래프 출력 및 저장
print(plot_H0)
print(plot_H1)
# ggsave("BF_comparison_H0.png", plot_H0, width = 15, height = 10)
# ggsave("BF_comparison_H1.png", plot_H1, width = 15, height = 10)



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


## 추가 분석을 위한 R 코드 (2024-10-08)

library(ggplot2)
library(dplyr)
library(tidyr)

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
# ggsave("sdr_method_interaction_h1.png", interaction_plot, width = 10, height = 6)