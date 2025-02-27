library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 기존 코드에서 데이터 로드 부분
# RDS 파일 로드
results_df1 <- readRDS("./New/merged_resultsV.RDS")
results_df2 <- readRDS("./New/merged_resultsV2.RDS")
results_df3 <- readRDS("./New/merged_resultsV3.RDS")

# 데이터 병합
results_df <- rbind(results_df1, results_df2, results_df3)

# BF 데이터 변환 (log 변환 포함)
bf_data <- results_df %>%
  select(starts_with("BF_"), scenario, sd_x, sd_y, mean_diff) %>%
  pivot_longer(
    cols = starts_with("BF_"),
    names_to = "model",
    values_to = "BF"
  ) %>%
  mutate(
    model = gsub("BF_", "", model),
    log_BF = log10(BF)
  )

# 시나리오 정보 정의 (이전 버전에서 확장)
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
    ratio = ifelse(n1 == n2, "Equal (1:1)", "Unequal (2:3)"),
    ratio_numeric = paste0(n1, ":", n2),
    sample_size_group = case_when(
      total_n <= 30 ~ "n=30",
      total_n <= 100 ~ "n=100",
      TRUE ~ "n=200"
    ),
    # SDR과 비율의 조합에 대한 열 레이블 생성
    condition_column = case_when(
      sdr == 1.0 & ratio == "Equal (1:1)" ~ "SDR=1.0, ratio=1:1",
      sdr == 0.5 & ratio == "Equal (1:1)" ~ "SDR=0.5, ratio=1:1",
      sdr == 1.0 & ratio == "Unequal (2:3)" ~ "SDR=1.0, ratio=2:3",
      sdr == 0.5 & ratio == "Unequal (2:3)" ~ "SDR=0.5, ratio=2:3",
      sdr == 2.0 & ratio == "Unequal (2:3)" ~ "SDR=2.0, ratio=2:3"
    ),
    # 패널 라벨 (간결하게)
    panel_label = sprintf("n₁=%d, n₂=%d, SDR=%.1f", n1, n2, sdr)
  )

# 평균 차이 값 정의
mean_diff_values <- c(0, 1, 2.2, 4.2, 5, 10)

# 데이터와 시나리오 정보 조인
bf_data_with_info <- bf_data %>%
  left_join(scenarios, by = "scenario") %>%
  mutate(
    # 평균 차이를 하나의 요인으로 변환 (행으로 사용될 것임)
    mean_diff_factor = factor(mean_diff, levels = mean_diff_values),
    # Cohen's d 계산 (대략적인 값)
    cohens_d = mean_diff / sqrt((sd1^2 + sd2^2)/2),
    # 열 조건 레이블의 요인화 (순서 지정)
    condition_column = factor(condition_column, levels = c(
      "SDR=1.0, ratio=1:1", 
      "SDR=0.5, ratio=1:1",
      "SDR=1.0, ratio=2:3",
      "SDR=0.5, ratio=2:3",
      "SDR=2.0, ratio=2:3"
    )),
    # 샘플 크기 그룹의 요인화
    sample_size_group = factor(sample_size_group, levels = c("n=30", "n=100", "n=200"))
  )

# 색상 및 방법 레이블 정의
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
      text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.spacing = unit(0.5, "lines")
    )
}

# 각 샘플 크기 그룹에 대한 밀도 플롯 생성 함수
create_density_by_sample_size <- function(data, size_group) {
  # 해당 샘플 크기 그룹 필터링
  plot_data <- data %>%
    filter(sample_size_group == size_group)
  
  # 예쁜 제목 생성
  nice_title <- paste0("Distribution of log(BF) for ", size_group, " scenarios")
  
  # 밀도 플롯 생성 (5x6 패널 그리드: 5개 열 조건 x 6개 평균 차이)
  ggplot(plot_data, aes(x = log_BF, fill = model)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_grid(
      mean_diff_factor ~ condition_column,
      scales = "free_y",
      labeller = labeller(
        mean_diff_factor = function(x) paste0("Mean diff = ", x, " (d ≈ ", round(as.numeric(x)/2, 2), ")")
      )
    ) +
    scale_fill_manual(values = method_colors, labels = bf_method_labels) +
    xlim(-1.5, 6) +
    labs(
      title = nice_title,
      x = expression(log[10](BF)),
      y = "Density",
      fill = "Method"
    ) +
    theme_paper() +
    theme(
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8)
    )
}

# 각 샘플 크기에 대한 그림 생성
plot_n30 <- create_density_by_sample_size(bf_data_with_info, "n=30")
plot_n100 <- create_density_by_sample_size(bf_data_with_info, "n=100")
plot_n200 <- create_density_by_sample_size(bf_data_with_info, "n=200")

# 그림 출력
print(plot_n30)
print(plot_n100)
print(plot_n200)

# PDF 또는 PNG로 저장하려면 아래 코드 사용
# ggsave("plot_n30.pdf", plot_n30, width = 12, height = 10)
# ggsave("plot_n100.pdf", plot_n100, width = 12, height = 10)
# ggsave("plot_n200.pdf", plot_n200, width = 12, height = 10)

# 평균 BF 값을 요약하는 함수 (각 조건에 대한 평균 log BF)
summarize_bf_by_condition <- function(data) {
  data %>%
    group_by(sample_size_group, condition_column, mean_diff_factor, model) %>%
    summarize(
      mean_log_BF = mean(log_BF, na.rm = TRUE),
      sd_log_BF = sd(log_BF, na.rm = TRUE),
      median_log_BF = median(log_BF, na.rm = TRUE),
      count = n(),
      .groups = "drop"
    )
}

# 요약 통계 생성
bf_summary <- summarize_bf_by_condition(bf_data_with_info)

# 각 샘플 크기에 대한 막대 그래프 생성 함수
create_barplot_by_sample_size <- function(summary_data, size_group) {
  # 해당 샘플 크기 그룹 필터링
  plot_data <- summary_data %>%
    filter(sample_size_group == size_group)
  
  # 예쁜 제목 생성
  nice_title <- paste0("Average log(BF) by condition for ", size_group, " scenarios")
  
  # 막대 그래프 생성
  ggplot(plot_data, aes(x = model, y = mean_log_BF, fill = model)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_log_BF - sd_log_BF/sqrt(count), 
          ymax = mean_log_BF + sd_log_BF/sqrt(count)),
      width = 0.2, position = position_dodge(0.7)
    ) +
    facet_grid(
      mean_diff_factor ~ condition_column,
      scales = "free_y",
      labeller = labeller(
        mean_diff_factor = function(x) paste0("Mean diff = ", x)
      )
    ) +
    scale_fill_manual(values = method_colors, labels = bf_method_labels) +
    labs(
      title = nice_title,
      x = "Bayes Factor Method",
      y = "Average log10(BF)",
      fill = "Method"
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "lightgray", color = NA),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8)
    )
}

# 각 샘플 크기에 대한 막대 그래프 생성
barplot_n30 <- create_barplot_by_sample_size(bf_summary, "n=30")
barplot_n100 <- create_barplot_by_sample_size(bf_summary, "n=100")
barplot_n200 <- create_barplot_by_sample_size(bf_summary, "n=200")

# 그래프 출력
print(barplot_n30)
print(barplot_n100)
print(barplot_n200)

# 저장하려면 아래 코드 사용
# ggsave("barplot_n30.pdf", barplot_n30, width = 12, height = 10)
# ggsave("barplot_n100.pdf", barplot_n100, width = 12, height = 10)
# ggsave("barplot_n200.pdf", barplot_n200, width = 12, height = 10)