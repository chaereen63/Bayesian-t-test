## 논문용 시각화
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
source(file = file.path("./Korea_pub/functionsK.R"))

# 데이터 로드
results_df <- readRDS("./Korea_pub/final_merged_results100.RDS")

# 데이터 변환
results_df %>% 
  select("scenario", "sdr", "delta", 
         "BF_jzs", "BF_gica") ->
  result_temp
str(result_temp)

# 1. 생성된 BF 분포 그리기
# 기본 데이터 변환
result_long <- result_temp %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "model",
    values_to = "BF"
  ) %>%
  mutate(
    model = gsub("BF_", "", model),
    log_BF = log10(BF),
    effect_size = case_when(
      delta == 0 ~ "no effect",
      delta == 0.2 ~ "weak",
      delta == 0.5 ~ "medium",
      delta == 0.8 ~ "strong"
    ),
    effect_size = factor(effect_size, 
                         levels = c("no effect", "weak", "medium", "strong"))
  )

# 시나리오 레이블 생성
scenarios <- tibble(
  scenario = 1:5,
  n1 = c(40, 50, 50, 40, 40),
  n2 = c(60, 50, 50, 60, 60),
  sdr = c(1.00, 1.00, 2.00, 2.00, 0.50)
) %>%
  mutate(label = sprintf("n1=%d, n2=%d, SDR=%.2f", n1, n2, sdr))

  # total sample size = 200
  # 시나리오 레이블 생성
  scenarios <- tibble(
    scenario = 1:5,
    n1 = c(80, 100, 100, 80, 80),
    n2 = c(120, 100, 100, 120, 120),
    sdr = c(1.00, 1.00, 2.00, 2.00, 0.50)
  ) %>%
    mutate(label = sprintf("n1=%d, n2=%d, SDR=%.2f", n1, n2, sdr))

# 시나리오 레이블을 벡터로
scenario_labels <- scenarios$label

# 색상 팔레트와 레이블 정의
method_colors <- c(
  "jzs" = "turquoise3",
  "gica" = "coral"
)

method_labels <- c(
 # "robtt_homo" = expression(BF[RH]),
  "jzs" = expression(BF[JZS]),
  "log(BF)" = expression(logBF[10])
)

# 공통 테마 설정
theme_paper <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Noto Sans KR", size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90")
    )
}

# 1. Ridgeline plots
create_ridgeline <- function(data, effect_size_val) {
  filtered_data <- data %>% 
    filter(effect_size == effect_size_val)
  
  ggplot(filtered_data, aes(x = log_BF, fill = model)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_wrap(~ scenario, nrow = 2, 
               labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
    scale_fill_manual(values = method_colors,
                      labels = method_labels) +
    xlim(-1.5, NA) +
    labs(
      title = paste("Distribution of log(BF) -", effect_size_val),
      x = "log(BF)",
      y = "Density"
    ) +
    theme_paper()
}

# 2. Box plots
create_boxplot <- function(data, effect_size_val) {
  filtered_data <- data %>% 
    filter(effect_size == effect_size_val)
  
  ggplot(filtered_data, aes(x = model, y = log_BF, fill = model)) +
    geom_boxplot(outlier.alpha = 0.3) +
    facet_wrap(~ scenario, nrow = 2,
               labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
    scale_fill_manual(values = method_colors,
                      labels = method_labels) +
    scale_x_discrete(labels = method_labels) +
    labs(
      title = paste("Distribution of log(BF) -", effect_size_val),
      x = "Method",
      y = expression(log[10](BF))
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 3. Method panel Bar plots
create_method_panel_bar <- function(data, effect_size_val) {
  # 데이터 가공
  plot_data <- data %>%
    filter(effect_size == effect_size_val) %>%
    group_by(scenario, model) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      ci_lower = mean_log_BF - 1.96 * se_log_BF,
      ci_upper = mean_log_BF + 1.96 * se_log_BF,
      .groups = 'drop'
    )
  
  ggplot(plot_data, aes(x = factor(scenario), y = mean_log_BF, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.9),
      width = 0.25,
      linewidth = 0.5
    ) +
    facet_wrap(~ model, nrow = 2,
               labeller = labeller(model = method_labels)) +
    scale_fill_manual(values = method_colors,
                      labels = method_labels) +
    scale_x_discrete(labels = scenario_labels) +
    labs(
      title = paste("Mean log(BF) by Method -", effect_size_val),
      x = "Scenario",
      y = expression(log[10](BF))
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 4. Scenario panel Bar plots
create_scenario_panel_bar <- function(data, effect_size_val) {
  # 데이터 가공
  plot_data <- data %>%
    filter(effect_size == effect_size_val) %>%
    group_by(scenario, model) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      ci_lower = mean_log_BF - 1.96 * se_log_BF,
      ci_upper = mean_log_BF + 1.96 * se_log_BF,
      .groups = 'drop'
    )
  
  ggplot(plot_data, aes(x = model, y = mean_log_BF, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.9),
      width = 0.25,
      linewidth = 0.5
    ) +
    facet_wrap(~ scenario, nrow = 2,
               labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
    scale_fill_manual(values = method_colors,
                      labels = method_labels) +
    scale_x_discrete(labels = method_labels) +
    labs(
      title = paste("Mean log(BF) by Scenario -", effect_size_val),
      x = "Method",
      y = expression(log[10](BF))
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# 데이터의 effect_size 레벨 확인
effect_size_levels <- levels(result_long$effect_size)
print(effect_size_levels)

# Plot 생성
# Ridgeline plots
rid0 <- create_ridgeline(result_long, effect_size_levels[1])  # "no effect"
rid2 <- create_ridgeline(result_long, effect_size_levels[2])  # "weak"
rid5 <- create_ridgeline(result_long, effect_size_levels[3])  # "medium"
rid8 <- create_ridgeline(result_long, effect_size_levels[4])  # "strong"
rid0;rid2;rid5;rid8
# Box plots
box0 <- create_boxplot(result_long, effect_size_levels[1])
box2 <- create_boxplot(result_long, effect_size_levels[2])
box5 <- create_boxplot(result_long, effect_size_levels[3])
box8 <- create_boxplot(result_long, effect_size_levels[4])

# Method panel Bar plots
barm0 <- create_method_panel_bar(result_long, effect_size_levels[1])
barm2 <- create_method_panel_bar(result_long, effect_size_levels[2])
barm5 <- create_method_panel_bar(result_long, effect_size_levels[3])
barm8 <- create_method_panel_bar(result_long, effect_size_levels[4])
barm0;barm2;barm5;barm8
# Scenario panel Bar plots
bars0 <- create_scenario_panel_bar(result_long, effect_size_levels[1])
bars2 <- create_scenario_panel_bar(result_long, effect_size_levels[2])
bars5 <- create_scenario_panel_bar(result_long, effect_size_levels[3])
bars8 <- create_scenario_panel_bar(result_long, effect_size_levels[4])

# Plot 저장 및 출력
# 각 효과크기별로 4개의 plot을 한번에 보여주기 위한 함수
show_plots <- function(plots, title) {
  # 4개의 plot을 2x2로 배치
  arranged_plot <- grid.arrange(
    plots[[1]], plots[[2]], 
    plots[[3]], plots[[4]], 
    ncol = 2, 
    top = title
  )
  
  # 배치된 plot 저장
  ggsave(
    paste0(title, ".png"), 
    arranged_plot, 
    width = 20, 
    height = 16
  )
  
  # RStudio에서 plot 보여주기
  print(arranged_plot)
}

# Ridgeline plots 저장 및 출력
ridgeline_plots <- list(rid0, rid2, rid5, rid8)
show_plots(ridgeline_plots, "ridgeline_all")

# Box plots 저장 및 출력
box_plots <- list(box0, box2, box5, box8)
show_plots(box_plots, "boxplot_all")

# Method panel Bar plots 저장 및 출력
method_bar_plots <- list(barm0, barm2, barm5, barm8)
show_plots(method_bar_plots, "method_panel_all")

# Scenario panel Bar plots 저장 및 출력
scenario_bar_plots <- list(bars0, bars2, bars5, bars8)
show_plots(scenario_bar_plots, "scenario_panel_all")

# 개별 파일 저장
# Ridgeline plots
ggsave("rid0.png", rid0, width = 12, height = 8)
ggsave("rid2.png", rid2, width = 12, height = 8)
ggsave("rid5.png", rid5, width = 12, height = 8)
ggsave("rid8.png", rid8, width = 12, height = 8)

# Box plots
ggsave("box0.png", box0, width = 12, height = 8)
ggsave("box2.png", box2, width = 12, height = 8)
ggsave("box5.png", box5, width = 12, height = 8)
ggsave("box8.png", box8, width = 12, height = 8)

# Method panel Bar plots
ggsave("barm0.png", barm0, width = 12, height = 8)
ggsave("barm2.png", barm2, width = 12, height = 8)
ggsave("barm5.png", barm5, width = 12, height = 8)
ggsave("barm8.png", barm8, width = 12, height = 8)

# Scenario panel Bar plots
ggsave("bars0.png", bars0, width = 12, height = 8)
ggsave("bars2.png", bars2, width = 12, height = 8)
ggsave("bars5.png", bars5, width = 12, height = 8)
ggsave("bars8.png", bars8, width = 12, height = 8)

# BF 범주 분류 함수 추가
result_categorized <- result_long %>%
  mutate(
    evidence_category = case_when(
      log_BF < 0 ~ "Favors H0",
      log_BF < 0.5 ~ "Anecdotal",
      log_BF < 1 ~ "Substantial",
      log_BF < 1.5 ~ "Strong",
      log_BF < 2 ~ "Very Strong",
      TRUE ~ "Decisive"
    ),
    # 순서를 위한 factor 변환
    evidence_category = factor(evidence_category, 
                               levels = c("Favors H0", "Anecdotal", "Substantial", 
                                          "Strong", "Very Strong", "Decisive"))
  )

# 각 조건별 비율 계산
evidence_proportions <- result_categorized %>%
  group_by(scenario, effect_size, model, evidence_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(scenario, effect_size, model) %>%
  mutate(
    percentage = count / sum(count) * 100
  ) %>%
  ungroup()

# 누적 막대 그래프 생성
ggplot(evidence_proportions, 
       aes(x = model, y = percentage, fill = evidence_category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(effect_size ~ scenario, 
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  scale_fill_brewer(palette = "Blues") +
  labs(
    x = "Model",
    y = "Percentage",
    fill = "Evidence Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# 퍼센트를 깔끔하게 정리한 테이블 만들기
evidence_summary <- evidence_proportions %>%
  # 시나리오 레이블 추가
  left_join(scenarios %>% select(scenario, label), by = "scenario") %>%
  # 특히 medium effect size에 집중
  filter(effect_size == "medium") %>%
  # 보기 좋게 정리
  select(label, model, evidence_category, percentage) %>%
  # 퍼센트를 소수점 1자리까지만
  mutate(percentage = round(percentage, 1)) %>%
  # 보기 좋게 피벗
  pivot_wider(
    names_from = evidence_category,
    values_from = percentage
  )

# 결과 출력
print(evidence_summary, n = Inf)

# 퍼센트를 깔끔하게 정리한 테이블 만들기
evidence_summary2 <- evidence_proportions %>%
  # 시나리오 레이블 추가
  left_join(scenarios %>% select(scenario, label), by = "scenario") %>%
  # 특히 strong effect size에 집중
  filter(effect_size == "strong") %>%
  # 보기 좋게 정리
  select(label, model, evidence_category, percentage) %>%
  # 퍼센트를 소수점 1자리까지만
  mutate(percentage = round(percentage, 1)) %>%
  # 보기 좋게 피벗
  pivot_wider(
    names_from = evidence_category,
    values_from = percentage
  )

# 결과 출력
print(evidence_summary2, n = Inf)

#### 생성 데이터 효과도 저장한 결과 ####
results_df2 <- readRDS("./Korea_pub/final_merged_results100.RDS")
results_df2 %>% 
  mutate(
    # 시나리오별 n1, n2 할당
    n1 = case_when(
      scenario == 1 ~ 40,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 40,
      scenario == 5 ~ 40
    ),
    n2 = case_when(
      scenario == 1 ~ 60,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 60,
      scenario == 5 ~ 60
    )
     ) %>%
  select(BF_jzs, BF_gica,mean_diff, sdr, delta, scenario) -> results_df3

# 데이터 변환 및 시각화를 위한 코드
result_with_d <- results_df2 %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "model",
    values_to = "BF"
  ) %>%
  mutate(
    model = gsub("BF_", "", model),
    log_BF = log10(BF)
  )
result_ratio <- results_df3 %>%
  mutate(bf_ratio = BF_gica/BF_jzs) %>%
  select(scenario, mean_diff, bf_ratio)

# BF_ratio scatter plot with alpha
ggplot(result_ratio, aes(x = mean_diff, y = bf_ratio)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +  # ratio = 1 기준선
  geom_point(alpha = 0.1) +
  facet_wrap(~scenario) +
  scale_y_log10() +  # log scale로 보면 1을 중심으로 대칭적으로 볼 수 있음
  labs(
    x = "Standardized Mean Difference",
    y = "BF Ratio (GICA/JZS)"
  )

# ratio의 분포를 더 자세히 보기
summary(result_ratio$bf_ratio)

# 특이한 패턴이 있는지 scatter plot으로도 확인_not log
ggplot(result_ratio, aes(x = std_mean_diff, y = bf_ratio)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~scenario)

# BF_scatter plot with alpha
ggplot(result_with_d, aes(x = mean_diff, y = BF)) +
  geom_point(alpha = 0.1, size = 0.5) +  # alpha와 size는 조정 가능
  facet_grid(scenario ~ model) +
  labs(
    x = "Standardized Mean Difference",
    y = "Bayes Factor"
  ) +
  coord_cartesian(ylim = c(-3, 30)) +
  theme_minimal()


# 효과크기별 평균 BF와 RMSE 기반 신뢰구간 계산
result_summary <- result_long %>%
  group_by(scenario, effect_size, model) %>%
  summarise(
    mean_bf = mean(log_BF),
    rmse = sqrt(mean((log_BF - mean(log_BF))^2)),
    se = rmse/sqrt(n()),
    ci_lower = mean_bf - 1.96*se,
    ci_upper = mean_bf + 1.96*se,
    .groups = 'drop'
  )

# facet별로 다른 y축 범위 설정
ggplot(result_summary, 
       aes(x = model, y = mean_bf, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(0.9), width = 0.2) +
  facet_grid(effect_size ~ scenario, scales = "free_y") +  # y축 범위를 자유롭게
  labs(
    x = "Model",
    y = "Mean log10(BF)",
    title = "Mean Bayes Factors by Effect Size and Scenario"
  ) +
  scale_fill_manual(values = c("coral", "turquoise3")) +  # 색상 유지
  theme_minimal()

# GICA vs JZS 산점도 생성
results_df %>% 
  mutate(
    effect_size = case_when(
      delta == 0 ~ "no effect",
      delta == 0.2 ~ "weak",
      delta == 0.5 ~ "medium",
      delta == 0.8 ~ "strong"
    ),
    effect_size = factor(effect_size, 
                         levels = c("no effect", "weak", "medium", "strong"))
  ) %>%
  ggplot(aes(x = BF_jzs, y = BF_gica, color = effect_size)) +
  # 1:1 라인 추가
  geom_abline(intercept = 0, slope = 1, color = "red", 
              linetype = "dashed", alpha = 0.5) +
  # 점 그리기
  geom_point(alpha = 0.3, size = 0.8) +
  # 시나리오별 패널 분할
  facet_wrap(~scenario, nrow = 2,
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  # 색상 설정
  scale_color_brewer(palette = "Set2") +
  # 로그 스케일 적용
  scale_x_log10() +
  scale_y_log10() +
  # 라벨 설정
  labs(
    x = expression(BF[JZS]),
    y = expression(BF[GICA]),
    color = "Effect Size"
  ) +
  # 테마 설정
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 11, face = "bold")
  )

# data driven
results_df3 %>% 
  mutate(
    effect_magnitude = abs(std_mean_diff),
    # log BF로 변환
    log_BF_jzs = log10(BF_jzs),
    log_BF_gica = log10(BF_gica)
  ) %>%
  ggplot(aes(x = log_BF_jzs, y = log_BF_gica, color = effect_magnitude)) +
  # 1:1 라인 추가
  geom_abline(intercept = 0, slope = 1, color = "red", 
              linetype = "dashed", alpha = 0.5) +
  # 점 그리기
  geom_point(alpha = 0.3, size = 0.8) +
  # 시나리오별 패널 분할
  facet_wrap(~scenario, nrow = 2,
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  # 색상 설정
  scale_color_viridis_c(limits = c(0, 1)) +  # SMD 범위 0-1로 제한
  # 라벨 설정
  labs(
    x = expression(log[10](BF[JZS])),
    y = expression(log[10](BF[GICA])),
    color = "|SMD|"
  ) +
  # 테마 설정
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 11, face = "bold")
  )

####Regression####
# 시나리오별 R² 계산
# lm()을 사용하되 회귀계수를 1로 고정하는 방법
r_squares <- results_df3 %>%
  group_by(scenario) %>%
  summarise(
    model = list(lm(I(log10(BF_gica) - log10(BF_jzs)) ~ 0)),  # 차이를 0에 회귀
    r_squared = 1 - sum(residuals(model[[1]])^2) / 
      sum((log10(BF_gica) - mean(log10(BF_gica)))^2)
  )
# 시나리오별 상관계수 계산
results_df3 %>%
  group_by(scenario) %>%
  summarise(
    correlation = cor(log10(BF_gica), log10(BF_jzs))^2  # R²를 바로 계산
  )

# 그래프에 R² 추가
results_df3 %>% 
  mutate(
    effect_magnitude = abs(std_mean_diff),
    log_BF_jzs = log10(BF_jzs),
    log_BF_gica = log10(BF_gica)
  ) %>%
  ggplot(aes(x = log_BF_jzs, y = log_BF_gica, color = effect_magnitude)) +
  geom_abline(intercept = 0, slope = 1, color = "red", 
              linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.3, size = 0.8) +
  # R² 값 텍스트 추가
  geom_text(data = r_squares,
            aes(x = -Inf, y = Inf, 
                label = sprintf("R² = %.3f", r_squared)),
            hjust = -0.2, vjust = 2,
            color = "black",
            inherit.aes = FALSE) +
  facet_wrap(~scenario, nrow = 2,
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  scale_color_viridis_c(limits = c(0, 1)) +
  labs(
    x = expression(log[10](BF[JZS])),
    y = expression(log[10](BF[GICA])),
    color = "|SMD|"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 11, face = "bold")
  )
