## 시각화
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(gridExtra)
home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_30 <- readRDS("./New/mergedFin30ES0.RDS")
results_100 <- readRDS("./New/mergedFin100ES0.RDS")
results_200 <- readRDS("./New/mergedFin200ES0.RDS")

# 시나리오 정의 - sample size = 30인 경우
scenarios_30 <- tibble(
  scenario = 1:5,
  n1 = c(15, 12, 15, 12, 18),
  n2 = c(15, 18, 15, 18, 12),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )

# 시나리오 정의 - sample size = 100인 경우
scenarios_100 <- tibble(
  scenario = 1:5,
  n1 = c(50, 40, 50, 40, 60),
  n2 = c(50, 60, 50, 60, 40),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )

# 시나리오 정의 - sample size = 200인 경우
scenarios_200 <- tibble(
  scenario = 1:5,
  n1 = c(100, 80, 100, 80, 120),
  n2 = c(100, 120, 100, 120, 80),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )

bf_method_labels <- c(
  "jzs" = expression(BF[JZS]),
  "gica" = expression(BF[GICA])
)

# 공통 테마 설정
theme_paper <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Noto Sans KR", size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),  # 패널 제목(시나리오 레이블) 크기 줄임
      plot.title = element_text(size = 14, face = "bold"),  # 전체 그래프 제목 크기 줄임
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 11, hjust = 1, face = "italic"),  # 각주 크기 및 스타일
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),  # 범례 크기 증가
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.2, "lines"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# 데이터 가공 함수 정의
process_bf_data <- function(results_df, scenarios_df) {
  results_tidy <- results_df %>%
    # Log 변환 적용
    mutate(
      log_BF_jzs = log10(BF_jzs),
      log_BF_gica = log10(BF_gica)
    ) %>%
    # 그래프용 긴 형태로 변환
    pivot_longer(
      cols = c(log_BF_jzs, log_BF_gica),
      names_to = "method",
      values_to = "log_BF"
    ) %>%
    # 방법 레이블 정리
    mutate(
      method_short = gsub("log_BF_", "", method)
    ) %>%
    # 시나리오 정보 추가
    left_join(scenarios_df, by = "scenario")
  
  return(results_tidy)
}

results_30_tidy <- process_bf_data(results_30, scenarios_30) # 데이터 준비
results_100_tidy <- process_bf_data(results_100, scenarios_100) # 데이터 준비
results_200_tidy <- process_bf_data(results_200, scenarios_200) # 데이터 준비

#### 누적그래프 그리기: BF 범주 분류 함수 추가####
result_categorized <- results_30_tidy %>%
  mutate(
    evidence_category = case_when(
      #log_BF < 0 ~ "Null hypothesis supported.",
      #0 <= log_BF & log_BF < 1/2 ~ "Not worth more than a bare mention.",
      #1/2 <= log_BF & log_BF < 1 ~ "substantial.",
      #1 <= log_BF & log_BF < 3/2 ~ "strong.",
      #3/2 <= log_BF & log_BF < 2 ~ "very strong.",
      #2 <= log_BF ~ "decisive."
      log_BF > 0 ~ "Alternative hypothesis supported.",
      0 >= log_BF & log_BF > -1/2 ~ "Not worth more than a bare mention.",
      -1/2 >= log_BF & log_BF > -1 ~ "substantial.",
      -1 >= log_BF & log_BF > -3/2 ~ "strong.",
      -3/2 >= log_BF & log_BF > -2 ~ "very strong.",
      -2 >= log_BF ~ "decisive."
    ),
    # 순서를 위한 factor 변환 - 범주 순서 수정
    #evidence_category = factor(evidence_category, 
    #                          levels = c("Null hypothesis supported.", 
    #                                      "Not worth more than a bare mention.", 
    #                                      "substantial.", 
    #                                      "strong.", 
    #                                      "very strong.", 
    #                                      "decisive."))
    # 순서를 위한 factor 변환 - 영가설 지지 강도에 맞춘 순서
    evidence_category = factor(evidence_category, 
                               levels = c("decisive.", 
                                          "very strong.", 
                                          "strong.", 
                                          "substantial.", 
                                          "Not worth more than a bare mention.", 
                                          "Alternative hypothesis supported."))
  )

# 각 조건별 비율 계산
evidence_proportions <- result_categorized %>%
  group_by(scenario, method_short, evidence_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(scenario, method_short) %>%
  mutate(
    percentage = count / sum(count) * 100
  ) %>%
  ungroup()

# 시나리오와 레이블 연결을 위한 매핑 데이터프레임 생성
scenario_labels <- scenarios_200 %>%
  select(scenario, label)

# 두 방법 모두 비교하려면 (JZS와 GICA)
plot_data_both <- result_categorized %>%
  group_by(scenario, method_short, evidence_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(scenario, method_short) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  # 시나리오 레이블 연결
  left_join(scenario_labels, by = "scenario")

# 두 방법을 비교하는 그래프
ggplot(plot_data_both, aes(x = factor(scenario), y = percentage, fill = evidence_category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  facet_wrap(~ method_short, labeller = labeller(method_short = c("jzs" = "JZS Method", "gica" = "GICA Method"))) +
  labs(
    title = "Distribution of Evidence Categories by Scenario and Method",
    x = "Scenario",
    y = "Percentage (%)",
    fill = "Evidence Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  # 시나리오 번호를 각 시나리오의 레이블로 매핑
  scale_x_discrete(labels = function(x) {
    sapply(x, function(scenario_num) {
      scenario_labels$label[scenario_labels$scenario == as.numeric(scenario_num)]
    })
  })
