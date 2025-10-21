#analysis the results

library(tidyverse)
# install.packages("asht")
library(asht)
source(file = file.path("./New/functionsN.R"))

n100es5 <- readRDS("./study2/S2merged100ES5.RDS")

head(n100es5)

# 시나리오 레이블 추가
n100es5 <- n100es5 %>%
  mutate(scenario_label = case_when(
    scenario == 1 ~ "A",
    scenario == 2 ~ "B",
    scenario == 3 ~ "C",
    scenario == 4 ~ "D",
    scenario == 5 ~ "E",
    TRUE ~ paste("시나리오", scenario)
  ))
print(n100es5)

#### 1. 효과크기 추정 편향 분석####
bias_analysis <- n100es5 %>%
  mutate(
    bias_s = d_s - true_d,
    bias_w = d_w - true_d,
    bias_c = d_c - true_d #pooled sd를 썼으므로 사실상 d_student와 동일해야함.
  ) %>%
  group_by(scenario_label, varr, n1, n2) %>%
  summarize(
    mean_bias_s = mean(bias_s),
    mean_bias_w = mean(bias_w),
    mean_bias_c = mean(bias_c),
    se_bias_s = sd(bias_s) / sqrt(n()),
    se_bias_w = sd(bias_w) / sqrt(n()),
    se_bias_c = sd(bias_c) / sqrt(n()),
    .groups = "drop"
  )

print(bias_analysis)

#### 2. 유의율 분석####
power_analysis <- n100es5 %>%
  group_by(scenario_label, varr, n1, n2) %>%
  summarize(
    power_student = mean(student_p < 0.05),
    power_welch = mean(welch_p < 0.05),
    power_diff = power_welch - power_student,
    se_power_s = sd(student_p < 0.05) / sqrt(n()),
    se_power_w = sd(welch_p < 0.05) / sqrt(n()),
    .groups = "drop"
  )

print(power_analysis)

error_analysis <- n100es5 %>%
  group_by(scenario_label, varr, n1, n2) %>%
  summarize(
    error_student = mean(student_p < 0.05),
    error_welch = mean(welch_p < 0.05),
    error_diff = error_welch - error_student,
    se_error_s = sd(student_p < 0.05) / sqrt(n()),
    se_error_w = sd(welch_p < 0.05) / sqrt(n()),
    .groups = "drop"
  )

print(error_analysis)

#### 3. 분산비-표본크기비 상호작용 시각화 ####
# 시나리오별 효과크기 편향 그래프
  # data
bias_long <- bias_analysis %>%
  select(scenario_label, mean_bias_s, mean_bias_w, se_bias_s, se_bias_w) %>%
  pivot_longer(
    cols = c(mean_bias_s, mean_bias_w),
    names_to = "method",
    values_to = "bias"
  ) %>%
  mutate(
    se = ifelse(method == "mean_bias_s", se_bias_s, se_bias_w),
    method = ifelse(method == "mean_bias_s", "Student", "Welch")
  )

  # plot
bias_plot <- ggplot(bias_long, aes(x = scenario_label, y = bias, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_errorbar(
    aes(ymin = bias - se, ymax = bias + se), 
    width = 0.1, 
    position = position_dodge(0.9)
  ) +
  labs(
    title = "조합별 효과크기 추정 편향 비교",
    subtitle = "n=200, effect size=0.0",
    x = "",
    y = "bias",
    fill = "methods"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Student" = "#00366C", "Welch" = "#F2A900"))

print(bias_plot)
ggsave("./plot2/bias_plotn200ES0.png", bias_plot, bg = "white")

#### 4. D와 E 조합 비교 (시나리오 4와 5) ####
D_E_comparison <- n100es5 %>%
  filter(scenario %in% c(4, 5)) %>%  # D와 E 시나리오
  mutate(
    bias_s = d_s - true_d,
    bias_w = d_w - true_d
  ) %>%
  group_by(scenario_label) %>%
  summarize(
    n1_n2 = paste(first(n1), ":", first(n2)),
    mean_bias_s = mean(bias_s),
    mean_bias_w = mean(bias_w),
    power_s = mean(student_p < 0.05),
    power_w = mean(welch_p < 0.05),
    .groups = "drop"
  )

print(D_E_comparison)

D_E_comparison <- n100es5 %>%
  filter(scenario %in% c(4, 5)) %>%  # D와 E 시나리오
  mutate(
    bias_s = d_s - true_d,
    bias_w = d_w - true_d
  ) %>%
  group_by(scenario_label) %>%
  summarize(
    n1_n2 = paste(first(n1), ":", first(n2)),
    mean_bias_s = mean(bias_s),
    mean_bias_w = mean(bias_w),
    error_s = mean(student_p < 0.05),
    error_w = mean(welch_p < 0.05),
    .groups = "drop"
  )

print(D_E_comparison)
#### 5. 모든 시나리오의 유의율 비교 ####
power_long <- power_analysis %>%
  select(scenario_label, power_student, power_welch, se_power_s, se_power_w) %>%
  pivot_longer(
    cols = c(power_student, power_welch),
    names_to = "method",
    values_to = "power"
  ) %>%
  mutate(
    se = ifelse(method == "power_student", se_power_s, se_power_w),
    method = ifelse(method == "power_student", "Student", "Welch")) %>%
  select(scenario_label, method, power, se)

  #plot
power_plot <- ggplot(power_long, aes(x = scenario_label, y = power, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_errorbar(
    aes(ymin = power - se, ymax = power + se), 
    width = 0.1, 
    position = position_dodge(0.9)
  ) +
  labs(
    title = "조합별 유의율 비교",
    subtitle = "n=200, effect size=0.0",
    x = "",
    y = "유의율 (p < .05 비율)",
    fill = "methods"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Student" = "#00366C", "Welch" = "#F2A900")) +
  coord_cartesian(ylim = c(0, 0.1)) +  # 적절한 범위 설정
  scale_y_continuous(breaks = seq(0, 1.1, 0.05))  # 0.05 간격으로 눈금 설정 눈금 간격 설정

print(power_plot)
ggsave("./plot2/power_n200es0.png",plot = power_plot, bg = "white")
