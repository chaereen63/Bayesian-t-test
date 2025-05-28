data1000 <- readRDS("./New/power_samplesize/addresults1000E2.RDS")
data600 <- readRDS("./New/power_samplesize/addresults600E2.RDS")
head(data600);head(data1000)

library(ggplot2)
library(dplyr)
library(afex);library(emmeans)

# Assuming data600 is already loaded
# First, let's transform the data for plotting

# Create a long format for the BF values
data600 <- data600 %>%
  select(scenario, BF_jzs, BF_gica) %>%
  tidyr::pivot_longer(
    cols = c(BF_jzs, BF_gica),
    names_to = "method",
    values_to = "BF"
  ) %>%
  mutate(
    method = case_when(
      method == "BF_jzs" ~ "JZS",
      method == "BF_gica" ~ "BFGC",
      TRUE ~ method
    ),
    # Taking log10 of BF for better visualization
    log_BF = log10(BF)
  ) %>% mutate(scenario = factor(scenario, levels = 1:5), sample_size = "600")

data1000 <- data1000 %>%
  select(scenario, BF_jzs, BF_gica) %>%
  tidyr::pivot_longer(
    cols = c(BF_jzs, BF_gica),
    names_to = "method",
    values_to = "BF"
  ) %>%
  mutate(
    method = case_when(
      method == "BF_jzs" ~ "JZS",
      method == "BF_gica" ~ "BFGC",
      TRUE ~ method
    ),
    # Taking log10 of BF for better visualization
    log_BF = log10(BF)
  ) %>% mutate(scenario = factor(scenario, levels = 1:5), sample_size = "1000")

plot_data <- bind_rows(data600, data1000)
plot_data$id <- ceiling(seq_len(nrow(plot_data)) / 2)

# 데이터 그룹화해서 평균 계산
library(dplyr)

plot_means <- plot_data %>%
  group_by(scenario, method, sample_size) %>%
  summarise(mean_log_BF = mean(log_BF), .groups = 'drop')

# 그래프 그리기
ggplot(plot_means, aes(x = scenario, 
                       y = mean_log_BF, 
                       color = interaction(method, sample_size),
                       shape = interaction(method, sample_size),
                       group = interaction(method, sample_size))) +
  geom_line() +
  geom_point(size = 2.3) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Mean~log[10]~BF[10]),
       x = expression("variance ratio" %*% "sample size ratio")) +
  scale_color_manual(name = "Method & Sample Size",
                     values = c("JZS.1000" = "darkviolet", 
                                "BFGC.1000" = "purple",
                                "JZS.600" = "gold", 
                                "BFGC.600" = "goldenrod"),
                     labels = c("JZS (n=1000)", "BFGC (n=1000)", 
                                "JZS (n=600)", "BFGC (n=600)")) +
  scale_shape_manual(name = "Method & Sample Size",
                     values = c("JZS.1000" = 16, 
                                "BFGC.1000" = 17,
                                "JZS.600" = 15, 
                                "BFGC.600" = 18),
                     labels = c("JZS (n=1000)", "BFGC (n=1000)", 
                                "JZS (n=600)", "BFGC (n=600)")) +
  theme_minimal() +
  theme(panel.spacing = unit(0, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(face = "bold", size = 11, family = "Times New Roman"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman"),
        axis.text = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 14, color = "black"))

# 각 시나리오별로 mean과 median 둘 다 확인
summary_stats <- plot_data %>%
  group_by(scenario, method, sample_size) %>%
  summarise(
    mean_log_BF = mean(log_BF),
    median_log_BF = median(log_BF),
    .groups = 'drop'
  )

print(summary_stats)

# 평균값의 왜곡이 심해 중앙값 사용
plot_median <- plot_data %>%
  group_by(scenario, method, sample_size) %>%
  summarise(median_log_BF = median(log_BF), .groups = 'drop')
# method를 factor로 변환하고 레벨 순서 지정
plot_median$method <- factor(plot_median$method, levels = c("JZS", "BFGC"))
plot_median$scenario <- factor(plot_median$scenario, levels = 1:5, labels = c("A", "B", "C", "D", "E"))
# 그래프 코드는 동일하게 유지
# 그래프 그리기
added_sample <- ggplot(plot_median, aes(x = scenario, 
                        y = median_log_BF, 
                        color = method,
                        linetype = method,
                        shape = sample_size,
                        group = interaction(method, sample_size))) +
  # 먼저 가로 점선 추가 (기존 그래프 요소보다 먼저 그려서 뒤에 배치)
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50") +
  # 기존 요소들
  geom_line() +
  geom_point(size = 2.3) +
  labs(title = "",
       y = expression(Median~log[10]~BF[10]),
       x = expression("variance ratio" %*% "sample size ratio")) +
  scale_color_viridis_d(option = "turbo", begin = 0.31, end = .95, direction = -1) +
  scale_linetype_manual(
                        values = c("JZS" = "solid", 
                                   "BFGC" = "dashed")) +
  scale_shape_manual(name = "Sample Size",
                     values = c("1000" = 16, 
                                "600" = 15),
                     labels = c("n=1000", "n=600")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(face = "bold", size = 8, family = "Times New Roman"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        text = element_text(family = "Times New Roman", size = 8),
        axis.title = element_text(family = "Times New Roman", size = 8),
        axis.text = element_text(family = "Times New Roman", size = 7),
        legend.text = element_text(family = "Times New Roman", size = 8),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 11, color = "black"),
        legend.key.width = unit(0.1, "lines"),  # 기호 영역 너비 줄이기
        legend.key.height = unit(0.1, "lines"), # 기호 영역 높이 줄이기
        legend.spacing.x = unit(0.2, "cm"),     # 범례 요소 간 간격 줄이기
        legend.margin = margin(0, 0, 0, 0, "cm")) # 범례 주변 여백 줄이기)



ggsave("sample600_10002.png", added_sample, dpi = 1000, units = "mm", width = 130, height = 80)
