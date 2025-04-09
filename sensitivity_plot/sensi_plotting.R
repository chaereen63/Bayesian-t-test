## Visualization
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
library(tidyverse)

home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_5 <- readRDS("./sensitivity_plot/merged_seni5.RDS")
results_0 <- readRDS("./sensitivity_plot/merged_seni0.RDS")

head(results_0)
nrow(results_0)

# 데이터 재구성 (long format으로 변환)
results_long <- results_5 %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "scale",
    values_to = "BF"
  ) %>%
  mutate(
    # log 변환
    logBF = log10(BF),
    # 스케일 이름 변경 및 팩터 설정
    scale = case_when(
      scale == "BF_jzsM" ~ "Medium",
      scale == "BF_jzsW" ~ "Wide",
      scale == "BF_jzsU" ~ "Ultrawide",
      scale == "BF_gica" ~ "BFGC",
      TRUE ~ scale
    ),
    scale = factor(scale, levels = c("Medium", "Wide", "Ultrawide", "BFGC"))
  )

# 평균 계산
results_mean <- results_long %>%
  group_by(scale) %>%
  summarise(mean_logBF = mean(logBF))

# BeFi 평균 추출
befi_mean <- results_mean %>% filter(scale == "BFGC") %>% pull(mean_logBF)

# JZS 데이터만 추출 (BeFi 제외)
results_jzs <- results_long %>% filter(scale != "BFGC")
results_jzs_mean <- results_mean %>% filter(scale != "BFGC")

# 그래프 생성 - 색상 통일
total <- ggplot() +
  # JZS 막대 그래프 - 모두 같은 색상 사용
  # geom_bar(data = results_jzs_mean, 
  #          aes(x = scale, y = mean_logBF),
  #          fill = "turquoise4", stat = "identity", alpha = 0.7, width = 0.6) +
  # JZS 평균값을 이은 선
  geom_line(data = results_jzs_mean,
            aes(x = scale, y = mean_logBF, group = 1),
            color = "#00366C", size = 0.8) +
  
  # JZS 평균값 포인트
  geom_point(data = results_jzs_mean,
             aes(x = scale, y = mean_logBF),
             color = "#00366C", size = 1) +
  
  # BeFi 가로선
  geom_hline(yintercept = befi_mean, linetype = "dashed", 
             color = "#F2A900", size = 0.8) +
  
  # BeFi 선 레이블
  annotate("text", x = 3, y = befi_mean + 0.01, 
           label = "BFGC", color = "#F2A900", family = "Times New Roman", size = 3) +
  
  # 그래프 제목 및 축 레이블 - 영어로 변경
  labs(title = "Comparison of log Bayes Factor Means by Scale",
       subtitle = "JZS (line) vs BeFi (dashed line)",
       x = "r-scale",
       y = "log(BF) Mean",
       caption = "Note: log10 scale; Cohen's d = 0") +
  
  # 테마 설정
  theme_minimal() +
  theme(
    legend.position = "none",  # 범례 제거 (모두 같은 색이므로)
    plot.title = element_text(face = "bold", hjust = 0.5, family = "Times New Roman", size = 10),
    plot.subtitle = element_text(hjust = 0.5, family = "Times New Roman", size = 9),
    axis.title = element_text(face = "bold", family = "Times New Roman", size = 9),
    plot.caption = element_text(hjust = 1, face = "italic", family = "Times New Roman", size = 8)
  )
print(total)
ggsave("total_sensitivity0.png", plot = total, width = 12, height = 9, units = "cm")

#### 표본크기 별로 ####
# 데이터 재구성 (long format으로 변환) 및 표본 크기 합 계산
results_long <- results_0 %>%
  mutate(sample_size_sum = n1 + n2) %>%  # 표본 크기 합 계산
  # 표본 크기 합 50, 100, 200만 필터링
  filter(sample_size_sum %in% c(50, 100, 200)) %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "scale",
    values_to = "BF"
  ) %>%
  mutate(
    # log 변환
    logBF = log10(BF),
    # 스케일 이름 변경 및 팩터 설정
    scale = case_when(
      scale == "BF_jzsM" ~ "Medium",
      scale == "BF_jzsW" ~ "Wide",
      scale == "BF_jzsU" ~ "Ultrawide",
      scale == "BF_gica" ~ "BFGC",
      TRUE ~ scale
    ),
    scale = factor(scale, levels = c("Medium", "Wide", "Ultrawide", "BFGC")),
    # 표본 크기 합을 요인으로 변환하고 순서 지정
    sample_size_sum = factor(sample_size_sum, levels = c(50, 100, 200))
  )

# 표본 크기 합 및 스케일별 평균 계산
results_mean <- results_long %>%
  group_by(sample_size_sum, scale) %>%
  summarise(mean_logBF = mean(logBF), .groups = 'drop')

# BeFi 데이터 추출 (각 표본 크기 합별)
befi_means <- results_mean %>% 
  filter(scale == "BFGC") %>%
  select(sample_size_sum, mean_logBF)

# JZS 데이터만 추출 (BeFi 제외)
results_jzs_mean <- results_mean %>% 
  filter(scale != "BFGC")

# 색상 정의
befi_color <- "#F2A900"  # 살구색/코랄
jzs_color <- "#00366C"   # 청록색

# 그래프 생성
plot <- ggplot() +
  # 표본 크기 합별로 패널 분할
  facet_wrap(~ sample_size_sum, labeller = labeller(sample_size_sum = function(x) paste0("n1 + n2 = ", x))) +
  
  # JZS 선 그래프
  geom_line(data = results_jzs_mean,
            aes(x = scale, y = mean_logBF, group = 1),
            color = jzs_color, size = 1) +
  
  # JZS 평균값 포인트
  geom_point(data = results_jzs_mean,
             aes(x = scale, y = mean_logBF),
             color = jzs_color, size = 3) +
  
  # 각 패널별로 BeFi 가로선 추가
  geom_hline(data = befi_means,
             aes(yintercept = mean_logBF),
             linetype = "dashed", color = befi_color, size = 1) +
  
  # 각 패널별로 BeFi 레이블 추가
  geom_text(data = befi_means,
            aes(x = 3, y = mean_logBF + 0.03, label = "BFGC"), #effect size= 0인 경우 +0.03
            color = befi_color, family = "Times New Roman") +
  
  # 그래프 제목 및 축 레이블
  labs(title = "Comparison of log Bayes Factor Means by Scale and Sample Size",
       subtitle = "JZS (bars) vs BFGC (line)",
       x = "r-scale",
       y = "log(BF) Mean",
       caption = "Note: log base 10 is used") +
  
  # 테마 설정
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(face = "bold", hjust = 0.5, family = "Times New Roman", size = 10),
    plot.subtitle = element_text(hjust = 0.5, family = "Times New Roman", size = 9),
    axis.title = element_text(face = "bold", family = "Times New Roman", size = 8),
    axis.text = element_text(family = "Times New Roman", size = 9),
    plot.caption = element_text(hjust = 1, face = "italic", family = "Times New Roman", size = 9),
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(face = "bold", family = "Times New Roman", size = 9)
  )
print(plot)
ggsave("sensitivity0.png", plot = plot, dpi = 600, width = 20, height = 10, units = "cm")
