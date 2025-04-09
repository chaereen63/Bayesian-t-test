## Visualization
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_50 <- readRDS("./New/robust_merge50E5_r1.RDS")
results_100 <- readRDS("./New/robust_merge100E8_r1v3.RDS")
results_200 <- readRDS("./New/robust_merge200E5_r1.RDS")
results_600 <- readRDS("./New/addresults600E2.RDS")
results_1000 <- readRDS("./New/addresults1000E2.RDS")

# Scenarios
scenarios_600 <- tibble(
  scenario = 1:5,
  n1 = c(300, 240, 300, 240, 360),
  n2 = c(300, 360, 300, 360, 240),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    # Concise format for labels
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )
scenarios_1000 <- tibble(
  scenario = 1:5,
  n1 = c(500, 400, 500, 400, 600),
  n2 = c(500, 600, 500, 600, 400),
  var1 = c(4, 4, 4, 4, 4),
  var2 = c(4, 4, 2, 2, 2)
) %>%
  mutate(
    varr = var2/var1,
    total_n = n1 + n2,
    # Concise format for labels
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )

# Define scenarios - sample size = 30
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
    # Concise format for labels
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
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
    # Include varr in label
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
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
    # Concise format for labels
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )
# Define color palette and labels
method_colors <- c(
  "jzs" = "#00366C",
  "BFGC" = "#F2A900"
)

# Define border colors (slightly darker)
method_border_colors <- c(
  "jzs" = "turquoise4",
  "BFGC" = "coral4"
)

bf_method_labels <- c(
  "jzs" = expression(BF[JZS]),
  "BFGC" = expression(BF[BFGC])
)

# Common theme settings
theme_paper <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman", size = 10),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),  # Reduced panel title size
      plot.title = element_text(size = 10, face = "bold"),  # Reduced overall graph title size
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 9, hjust = 1, face = "italic"),  # Footnote size and style
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),  # Increased legend size
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.2, "lines"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
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

# Frequency distribution histogram function
plot_bf_histogram <- function(bf_tidy, title) {
  # Create histogram separated by scenario
  ggplot(bf_tidy, aes(x = log_BF, fill = method_short)) +
    # Lower histogram transparency to 0.5 and add borders
    geom_histogram(aes(y = after_stat(density)), 
                   alpha = 0.5, 
                   position = "identity",
                   bins = 30,
                   color = "white",  # White boundaries between histograms
                   linewidth = 0.2) +
    # Increase density curve transparency to 0.7 and line thickness to 1.2
    geom_density(aes(color = method_short), 
                 alpha = 0.7,
                 linewidth = 1.2) +
    facet_wrap(~ scenario, 
               scales = "fixed",
               labeller = labeller(scenario = setNames(bf_tidy$label[!duplicated(bf_tidy$scenario)], 
                                                       unique(bf_tidy$scenario)))) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_color_manual(values = method_border_colors,  # Use darker border colors
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
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
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
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )
}

# === results_50 ===
results_50_tidy <- process_bf_data(results_50, scenarios_50) # Prepare data
hist_50 <- plot_bf_histogram(results_50_tidy, "(Total N=50, Effect size = 0.5)") # Create histogram
box_50 <- plot_bf_boxplot(results_50_tidy, "(Total N=50, Effect size = 0.5)") # Create boxplot
print(hist_50)
print(box_50)
# === results_100 ===
results_100_tidy <- process_bf_data(results_100, scenarios_100) # Prepare data
hist_100 <- plot_bf_histogram(results_100_tidy, "(Total N=100, Effect size = 0.5)") # Create histogram
box_100 <- plot_bf_boxplot(results_100_tidy, "(Total N=100, Effect size = 0.5)") # Create boxplot
print(hist_100)
print(box_100)
# === results_200 ===
results_200_tidy <- process_bf_data(results_200, scenarios_200)
hist_200 <- plot_bf_histogram(results_200_tidy, "(Total N=200, Effect size = 0.5)")
box_200 <- plot_bf_boxplot(results_200_tidy, "(Total N=200, Effect size = 0.5)")
print(hist_200)
print(box_200)
# === results_600 ===
results_600_tidy <- process_bf_data(results_600, scenarios_600)
box_600 <- plot_bf_boxplot(results_600_tidy, "(Total N=600, Effect size = 0.2)")
# === results_600 ===
results_1000_tidy <- process_bf_data(results_1000, scenarios_1000)
box_1000 <- plot_bf_boxplot(results_1000_tidy, "(Total N=1000, Effect size = 0.2)")
# save added scenarios
ggsave("boxTotal_600.png", plot = box_600, width = 14, height = 10, dpi = 600, units = "cm")
ggsave("boxTotal_1000.png", plot = box_1000, width = 14, height = 10, dpi = 600, units = "cm")
# === Additional Analysis: Comparison of BF_jzs and BF_gica ===

# Calculate summary statistics by method
summary_stats_100 <- results_100_tidy %>% # Change sample if needed
  group_by(scenario, method_short) %>%
  summarise(
    mean_log_BF = mean(log_BF),
    median_log_BF = median(log_BF),
    sd_log_BF = sd(log_BF),
    q25 = quantile(log_BF, 0.25),
    q75 = quantile(log_BF, 0.75),
    .groups = "drop"
  )

print(summary_stats_100)

# Alternative approach: Using log differences (log(BF_JZS) - log(BF_BeFi))
log_bf_diff <- function(results_df, scenarios_df, title) {
  # Calculate log Bayes factor differences and prepare data
  scatter_data <- results_df %>%
    mutate(
      log_BF_gica = log10(BF_gica),  # Log transformation
      log_BF_jzs = log10(BF_jzs),    # Log transformation
      log_bf_diff = log_BF_jzs - log_BF_gica  # Log BF difference
    ) %>%
    left_join(scenarios_df, by = "scenario")
  
  # Create scatter plot by scenario
  ggplot(scatter_data, aes(x = -mean_diff, y = log_bf_diff)) +
    geom_point(alpha = 0.3, size = 0.8) +  # Set point transparency and size
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +  # Reference line for difference = 0
    #geom_smooth(method = "loess", color = "blue", se = FALSE, linewidth = 0.8) +  # LOESS regression line
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(scatter_data$label[!duplicated(scatter_data$scenario)], 
                                                       unique(scatter_data$scenario))),
               nrow = 2, ncol = 3) +  # 2 rows, 3 columns layout
    labs(
      title = paste("Log Bayes Factor Difference (log(BF_JZS) - log(BF_BeFi))", title),
      x = "Mean Difference (mean1-mean2)",
      y = "Log(BF) Difference (JZS - BeFi)"
    ) +
    theme_paper() +
    theme(
      aspect.ratio = 0.6,  # Aspect ratio
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 9, margin = margin(b = 5)),
      panel.spacing = unit(1, "lines")
    )
}

# Create and print log difference scatter plots
log_diff_50 <- log_bf_diff(results_50, scenarios_50, "(Total N=50, Effect size = 0.5)")
log_diff_100 <- log_bf_diff(results_100, scenarios_100, "(Total N=100, Effect size = 0.5)")
log_diff_200 <- log_bf_diff(results_200, scenarios_200, "(Total N=200, Effect size = 0.5)")
print(log_diff_50);print(log_diff_100);print(log_diff_200)

# Define mean bar graph function by scenario (free y-axis scale)
plot_mean_bar <- function(bf_tidy, title) {
  # Calculate average log BF values by scenario and method
  summary_stats <- bf_tidy %>%
    group_by(scenario, method_short, label) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      .groups = "drop"
    )
  
  # Create graph dividing scenarios into individual panels (free y-axis scale)
  ggplot(summary_stats, aes(x = method_short, y = mean_log_BF, fill = method_short)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_log_BF - se_log_BF, ymax = mean_log_BF + se_log_BF),
                  width = 0.25, linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    facet_wrap(~ scenario, scales = "fixed",  # Adjust y-axis freely for each panel
               labeller = labeller(scenario = setNames(summary_stats$label[!duplicated(summary_stats$scenario)], 
                                                       unique(summary_stats$scenario))),
               nrow = 2, ncol = 3) +  # 2 rows, 3 columns layout
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
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),  # Adjust x-axis text alignment
      strip.text = element_text(size = 9),  # Adjust panel title size
      panel.spacing = unit(1, "lines"),  # Adjust panel spacing
      legend.key.size = unit(0.4, "cm"),       # 범례 키 크기 조절
      legend.spacing.x = unit(0.2, "cm"),      # 범례 아이템 간 수평 간격
      legend.margin = margin(t = 0, r = 5, b = 0, l = 5),  # 범례 전체 마진
      legend.text = element_text(size = 9)     # 범례 글자 크기 조절
    )
}
mean_bar_50 <- plot_mean_bar(results_50_tidy, "(Total N=50, Effect size = 0.5)")
mean_bar_100 <- plot_mean_bar(results_100_tidy, "(Total N=100, Effect size = 0.5)")
mean_bar_200 <- plot_mean_bar(results_200_tidy, "(Total N=200, Effect size = 0.5)")
mean_bar_600 <- plot_mean_bar(results_600_tidy, "(Total N=600, Effect size = 0.2)")
mean_bar_1000 <- plot_mean_bar(results_1000_tidy, "(Total N=1000, Effect size = 0.2)")
# Output graphs
print(mean_bar_50)
print(mean_bar_100)
print(mean_bar_200)

# save added scenario
ggsave("barTotal_600.png", plot = mean_bar_600, width = 14, height = 10, dpi = 600, units = "cm")
ggsave("barTotal_1000.png", plot = mean_bar_1000, width = 14, height = 10, dpi = 600, units = "cm")

#### Additional visualization experiments below ####
# Correlation analysis function by scenario for BF_jzs and BF_gica (blue-yellow-red color contrast)
analyze_correlation <- function(results_df, scenarios_df, title) {
  # Data preparation
  corr_data <- results_df %>%
    select(scenario, BF_jzs, BF_gica, mean_diff) %>%
    mutate(
      log_BF_jzs = log10(BF_jzs),
      log_BF_gica = log10(BF_gica),
      abs_mean_diff = abs(mean_diff)
    ) %>%
    left_join(scenarios_df, by = "scenario")
  
  # Check absolute value range
  max_abs_diff <- max(corr_data$abs_mean_diff, na.rm = TRUE)
  
  # Calculate correlation coefficients by scenario
  corr_stats <- corr_data %>%
    group_by(scenario) %>%
    summarise(
      pearson_r = cor(log_BF_jzs, log_BF_gica, method = "pearson"),
      r_squared = pearson_r^2,
      n_obs = n(),
      label = first(label),
      .groups = "drop"
    )
  
  # Create scatter plot - blue-yellow-red color contrast
  scatter_plot <- ggplot(corr_data, aes(x = log_BF_jzs, y = log_BF_gica)) +
    # Set point color according to abs_mean_diff
    geom_point(aes(color = abs_mean_diff), alpha = 0.7, size = 1.2) +
    # Show 1:1 line in light gray
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
    # Set color scale - blue-yellow-red contrast
    scale_color_gradientn(
      colors = c("darkblue", "blue", "royalblue", "skyblue", 
                 "yellow", 
                 "orange", "red", "darkred"),
      values = scales::rescale(c(0, max_abs_diff*0.2, max_abs_diff*0.3, max_abs_diff*0.4, 
                                 max_abs_diff*0.5, max_abs_diff*0.6, max_abs_diff*0.8, max_abs_diff)),
      name = "|Mean Diff (SMD)|"
    ) +
    # Separate by scenario in panels
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(
                 corr_stats$label, 
                 corr_stats$scenario
               )),
               nrow = 2, ncol = 3) +
    # Display correlation coefficient in each panel
    geom_text(data = corr_stats,
              aes(label = sprintf("R² = %.3f", r_squared),
                  x = -Inf, y = Inf),
              hjust = -0.1, vjust = 1.2, size = 3) +
    # Graph title and axis labels
    labs(
      title = paste("Correlation Between Log BF Methods", title),
      subtitle = "Dashed line: 1:1 line, Point color: Absolute mean difference (|SMD|)",
      x = expression(log(BF[JZS])),
      y = expression(log(BF[BeFi])),
      caption = "Note: R² = coefficient of determination"
    ) +
    theme_paper() +
    theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom"
    )
  
  return(scatter_plot)
}

# Run correlation analysis for each dataset
scatter_30 <- analyze_correlation(results_30, scenarios_30, "(N=30, Effect size = 0.8)")
scatter_100 <- analyze_correlation(results_100, scenarios_100, "(N=100, Effect size = 0.8)")
scatter_200 <- analyze_correlation(results_200, scenarios_200, "(N=200, Effect size = 0.8)")

# Output results
print(scatter_30)
print(scatter_100)
print(scatter_200)