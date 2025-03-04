## 시각화
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_30 <- readRDS("./New/merged_resultsm1_30.RDS")
results_100 <- readRDS("./New/merged_resultsm1_100.RDS")
results_200 <- readRDS("./New/merged_resultsm1_200.RDS")

# Combine them into a single data frame using rbind #여기부터 2025-03-05
# Assuming they have the same structure/columns
results_df <- rbind(results_df1, results_df2, results_df3)
str(results_df)

# 데이터 변환
results_df %>% 
  select("scenario", 
         "BF_jzs", "BF_gica", "mean_diff", "sd_x","sd_y") ->
  result_temp

# 1. 생성된 BF 분포 그리기
# 기본 데이터 변환
# BF 
bf_data <- results_df %>%
  select(starts_with("BF_"), scenario, mean_diff, sd_x, sd_y) %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "model",
    values_to = "BF"
  ) %>%
  mutate(
    model = gsub("BF_", "", model),
    log_BF = log10(BF)
  )

# 시나리오 레이블 생성 (평균 차이 포함)
scenarios <- tibble(
  scenario = 1:15,
  n1 = c(15, 15, 12, 12, 12, 50, 50, 40, 40, 40, 100, 100, 80, 80, 80),
  n2 = c(15, 15, 18, 18, 18, 50, 50, 60, 60, 60, 100, 100, 120, 120, 120),
  sd1 = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
  sd2 = c(2, 1, 2, 1, 4, 2, 1, 2, 1, 4, 2, 1, 2, 1, 4)
) %>%
  mutate(
    sdr = sd2/sd1,
    total_n = n1 + n2,
    label = sprintf("n=%d (n1=%d, n2=%d), SDR=%.1f", total_n, n1, n2, sdr)
  )

# 시나리오와 평균 차이를 결합한 레이블 생성
mean_diff_values <- c(0, 1, 2.2, 4.2, 5, 10)
scenario_mean_labels <- tibble(
  scenario = rep(1:15, each = length(mean_diff_values)),
  mean_diff = rep(mean_diff_values, times = 15)
) %>%
  left_join(scenarios, by = "scenario") %>%
  mutate(
    scenario_mean_id = paste0(scenario, "_", mean_diff),  # 고유 ID 생성
    mean1 = -mean_diff/2,
    mean2 = mean_diff/2,
    cohens_d = (mean2 - mean1) / sqrt((sd1^2 + sd2^2)/2),  # 간략한 코헨의 d 계산
    label_with_mean = sprintf("%s, Mean diff=%.1f (d=%.2f)", 
                              label, mean_diff, cohens_d)
  )

# 시나리오 레이블을 벡터로
scenario_labels <- scenario_mean_labels$label

# 색상 팔레트와 레이블 정의
method_colors <- c(
  "jzs" = "turquoise3",
  "gica" = "coral"
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
      strip.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      # 가독성을 위해 패널 간격 조정
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1, "lines")
    )
}

# 누락된 add_scenario_info 함수 추가
add_scenario_info <- function(data) {
  data %>%
    left_join(scenarios, by = "scenario") %>%
    mutate(
      sample_size_group = case_when(
        total_n <= 30 ~ "n = 30",
        total_n <= 100 ~ "n = 100",
        TRUE ~ "n = 200"
      ),
      # factor로 변환하여 순서 지정
      sample_size_group = factor(sample_size_group, 
                                 levels = c("n = 30", "n = 100", "n = 200")),
      allocation = ifelse(n1 == n2, "Equal (1:1)", "Unequal (2:3)"),
      ratio = paste(n1, ":", n2),
      # 시나리오 라벨 생성
      scenario_label = sprintf("n=%d (n1=%d, n2=%d), SDR=%.1f", total_n, n1, n2, sdr)
    )
}

# Prepare data for plotting
plot_data <- bf_data %>%
  filter(scenario %in% n100_scenarios) %>%
  left_join(scenario_mean_labels, by = c("scenario", "mean_diff")) %>%
  # Rename models to match the legend in the image
  mutate(model = case_when(
    model == "gica" ~ "BF_GCA",
    model == "jzs" ~ "BF_JZS",
    TRUE ~ model
  ))

# Function to create nice labels
label_function <- function(x) {
  gsub("n=([0-9]+) \\(n1=([0-9]+), n2=([0-9]+)\\), SDR=([0-9\\.]+)", 
       "n=\\1 (n1=\\2, n2=\\3), SDR=\\4", x)
}


# Create faceted density plots
ggplot(plot_data, aes(x = log_BF, fill = model, color = model)) +
  geom_density(alpha = 0.6) +
  # Add vertical line at log(BF) = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Facet by scenario and mean difference
  facet_grid(
    mean_diff ~ label,
    labeller = labeller(
      mean_diff = function(x) paste0("Mean diff = ", x),
      label = label_function
    )
  ) +
  scale_fill_manual(values = c("BF_GCA" = "#FF9966", "BF_JZS" = "#66CCCC")) +
  scale_color_manual(values = c("BF_GCA" = "#FF6633", "BF_JZS" = "#339999")) +
  labs(
    title = "Distribution of log(BF) for n=100 Scenarios",
    x = "log₁₀(BF)",
    y = "Density",
    fill = "Method",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1, "cm")
  ) +
  ylim(0, 0.45) +  # Match y-axis limits from image
  xlim(-0.5, 6)    # Match x-axis limits from image
