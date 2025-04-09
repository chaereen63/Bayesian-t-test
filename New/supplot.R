## Visualization
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(gridExtra)
home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_50 <- readRDS("./New/mergedFin50ES8_r1.RDS")
results_100 <- readRDS("./New/robust_merge100E5_r1.RDS")
results_200 <- readRDS("./New/mergedFin200ES8_r1.RDS")

# Define scenarios - sample size = 50
scenarios_50 <- tibble(
  scenario = 1:5,
  n1 = c(25, 20, 25, 20, 30),
  n2 = c(25, 30, 25, 30, 20),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    # 시나리오 번호만 표시
    label = as.character(scenario)
  )

# Define scenarios - sample size = 100
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
    # 시나리오 번호만 표시
    label = as.character(scenario)
  )

# Define scenarios - sample size = 200
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
    # 시나리오 번호만 표시
    label = as.character(scenario)
  )

# Define color palette and labels
method_colors <- c(
  "jzs" = "#00366C",
  "BFGC" = "#F2A900"
)
# 명도 확인
#library(colorspace)
#desaturate(c("#00366C", "#4682B4", "#F2A900"))

# Define border colors (slightly darker)
method_border_colors <- c(
  "jzs" = "#00355C",
  "BFGC" = "#F2A903"
)

bf_method_labels <- c(
  "jzs" = expression(BF[JZS]),
  "BFGC" = expression(BF[BFGC])
)

# Common theme settings - 논문에 맞게 폰트 사이즈 축소
theme_paper <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman", size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 7),
      plot.caption = element_text(size = 7, hjust = 1, face = "italic"),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(0.8, "lines"),
      panel.spacing.y = unit(0.8, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}

# Data processing function
process_bf_data <- function(results_df, scenarios_df) {
  results_tidy <- results_df %>%
    # Apply log transformation
    mutate(
      log_BF_jzs = log10(BF_jzs),
      log_BF_gica = log10(BF_gica)
    ) %>%
    # Transform to long format for plotting
    pivot_longer(
      cols = c(log_BF_jzs, log_BF_gica),
      names_to = "method",
      values_to = "log_BF"
    ) %>%
    # Clean up method labels - change gica to BeFi here
    mutate(
      method_short = case_when(
        method == "log_BF_gica" ~ "BFGC",
        method == "log_BF_jzs" ~ "jzs",
        TRUE ~ gsub("log_BF_", "", method)
      )
    ) %>%
    # Add scenario information
    left_join(scenarios_df, by = "scenario")
  
  return(results_tidy)
}

# 히스토그램 함수 - 코드 완성을 위해 추가
plot_bf_hist <- function(bf_tidy, title) {
  ggplot(bf_tidy, aes(x = log_BF, fill = method_short)) +
    # Lower histogram transparency to 0.5 and add borders
    geom_histogram(aes(y = after_stat(density)), 
                   alpha = 0.5, 
                   position = "identity",
                   bins = 30,
                   color = "white",
                   linewidth = 0.2) +
    # Increase density curve transparency to 0.7 and line thickness to 1.2
    geom_density(aes(color = method_short), 
                 alpha = 0.7,
                 linewidth = 0.8) +
    facet_wrap(~ scenario, 
               scales = "fixed",
               labeller = labeller(scenario = setNames(bf_tidy$label[!duplicated(bf_tidy$scenario)], 
                                                       unique(bf_tidy$scenario)))) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_color_manual(values = method_border_colors,
                       labels = bf_method_labels) +
    labs(
      title = paste("Log Bayes Factor Distribution", title),
      x = "log(BF)",
      y = "Density",
      fill = "Method",
      color = "Method"
    ) +
    theme_paper()
}

# Boxplot function
plot_bf_boxplot <- function(bf_tidy, title) {
  ggplot(bf_tidy, aes(x = factor(scenario), y = log_BF, fill = method_short)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, linewidth = 0.3, outlier.size = 0.5) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_x_discrete(labels = bf_tidy$label[!duplicated(bf_tidy$scenario)]) +
    labs(
      title = paste("Log BF Box plot by Scenario", title),
      x = "Scenario",
      y = "log(BF)",
      fill = "Method"
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8)
    )
}

# Log difference function
log_bf_diff <- function(results_df, scenarios_df, title) {
  # Calculate log Bayes factor differences and prepare data
  scatter_data <- results_df %>%
    mutate(
      log_BF_gica = log10(BF_gica),
      log_BF_jzs = log10(BF_jzs),
      log_bf_diff = log_BF_jzs - log_BF_gica
    ) %>%
    left_join(scenarios_df, by = "scenario")
  
  # Create scatter plot by scenario
  ggplot(scatter_data, aes(x = -mean_diff, y = log_bf_diff)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(scatter_data$label[!duplicated(scatter_data$scenario)], 
                                                       unique(scatter_data$scenario))),
               nrow = 2, ncol = 3) +
    labs(
      title = paste("Log Bayes Factor Difference (log(BF_JZS) - log(BF_BFGC))", title),
      x = "Mean Difference (mean1-mean2)",
      y = "Log(BF) Difference (JZS - BFGC)"
    ) +
    theme_paper() +
    theme(
      aspect.ratio = 0.6,
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 7, margin = margin(b = 3)),
      panel.spacing = unit(0.8, "lines")
    )
}

# Mean bar graph function
plot_mean_bar <- function(bf_tidy, title) {
  # Calculate average log BF values by scenario and method
  summary_stats <- bf_tidy %>%
    group_by(scenario, method_short, label) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create graph dividing scenarios into individual panels (fixed y-axis scale)
  ggplot(summary_stats, aes(x = method_short, y = mean_log_BF, fill = method_short)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8, linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_log_BF - se_log_BF, ymax = mean_log_BF + se_log_BF),
                  width = 0.2, linewidth = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(summary_stats$label[!duplicated(summary_stats$scenario)], 
                                                       unique(summary_stats$scenario))),
               nrow = 2, ncol = 3) +
    labs(
      title = paste("Mean Log BF Comparison by Scenario", title),
      subtitle = "Error bars represent standard error (SE)",
      x = "Method",
      y = "Mean log(BF)",
      fill = "Method",
      caption = "Note: log base 10 is used"
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
      strip.text = element_text(size = 7),
      panel.spacing = unit(0.8, "lines")
    )
}

# === 데이터 준비 ===
results_50_tidy <- process_bf_data(results_50, scenarios_50)
results_100_tidy <- process_bf_data(results_100, scenarios_100)
results_200_tidy <- process_bf_data(results_200, scenarios_200)

# === 그래프 생성 및 저장 (13cm x 8cm, 600dpi) ===

# 저장 함수 정의 - cm 단위를 inch로 변환 (1 inch = 2.54 cm)
save_plot <- function(filename, plot, width = 13, height = 8, dpi = 600, units = "cm") {
  ggsave(filename, plot, width= width, height = height, dpi = dpi, units = units)
}

# Boxplot 생성 및 저장
box_50 <- plot_bf_boxplot(results_50_tidy, "(Total N=50, Effect size = 0.8)")
box_100 <- plot_bf_boxplot(results_100_tidy, "(Total N=100, Effect size = 0.5)")
box_200 <- plot_bf_boxplot(results_200_tidy, "(Total N=200, Effect size = 0.8)")

save_plot("boxplot_50E8.png", box_50)
save_plot("robustboxplot_100E5.png", box_100)
save_plot("boxplot_200E8.png", box_200)

# 로그 차이 산점도 생성 및 저장
log_diff_50 <- log_bf_diff(results_50, scenarios_50, "(Total N=50, Effect size = 0.8)")
log_diff_100 <- log_bf_diff(results_100, scenarios_100, "(Total N=100, Effect size = 0.8)")
log_diff_200 <- log_bf_diff(results_200, scenarios_200, "(Total N=200, Effect size = 0.8)")

save_plot("log_diff_50_8.png", log_diff_50, 20, 10)
save_plot("log_diff_100_8.png", log_diff_100, 20, 10)
save_plot("log_diff_200_8.png", log_diff_200, 20, 10)

#### 필요시 ####
# 평균 막대 그래프 생성 및 저장
mean_bar_50 <- plot_mean_bar(results_50_tidy, "(Total N=50, Effect size = 0.5)")
mean_bar_100 <- plot_mean_bar(results_100_tidy, "(Total N=100, Effect size = 0.5)")
mean_bar_200 <- plot_mean_bar(results_200_tidy, "(Total N=200, Effect size = 0.5)")

save_plot("mean_bar_50_5.png", mean_bar_50)
save_plot("robmean_bar_100_5.png", mean_bar_100)
save_plot("mean_bar_200_0.5.png", mean_bar_200)

# 히스토그램도 생성 및 저장
hist_50 <- plot_bf_hist(results_50_tidy, "(Total N=50, Effect size = 0.5)")
hist_100 <- plot_bf_hist(results_100_tidy, "(Total N=100, Effect size = 0.5)")
hist_200 <- plot_bf_hist(results_200_tidy, "(Total N=200, Effect size = 0.5)")

save_plot("hist_50.png", hist_50)
save_plot("hist_100.png", hist_100)
save_plot("hist_200.png", hist_200)