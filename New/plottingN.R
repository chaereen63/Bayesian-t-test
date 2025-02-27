## 시각화
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_df1 <- readRDS("./New/merged_resultsV.RDS")
results_df2 <- readRDS("./New/merged_resultsV2.RDS")
results_df3 <- readRDS("./New/merged_resultsV3.RDS")

# Combine them into a single data frame using rbind
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
  select(starts_with("BF_"), scenario, sd_x, sd_y) %>%
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

# 시나리오별 개별 플롯 함수
create_individual_bf_density <- function(data) {
  plot_data <- add_scenario_info(data)
  
  # 시나리오 순서 정렬을 위한 팩터 생성
  plot_data <- plot_data %>%
    mutate(
      # 시나리오를 명확하게 정렬: 표본 크기, 할당 비율, SDR 순
      scenario_ordered = factor(
        scenario_label,
        levels = scenarios %>% 
          arrange(total_n, n1 == n2, sdr) %>% 
          pull(label)  # 이전에는 scenario_label이었지만 scenarios 테이블에서는 'label'로 정의됨
      )
    )
  
  # 5x3 그리드 레이아웃으로 표시
  ggplot(plot_data, aes(x = log_BF, fill = model)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_wrap(
      ~ scenario_ordered,
      scales = "free_y",
      ncol = 5  # 행당 5개 시나리오
    ) +
    scale_fill_manual(values = method_colors[c("jzs", "gica")],
                      labels = bf_method_labels) +
    xlim(-1.5, 12) +  # 긴 꼬리를 포착하기 위해 x축 확장
    labs(
      title = "Distribution of log(BF) by Individual Scenario",
      x = expression(log[10](BF)),
      y = "Density",
      fill = "Method"
    ) +
    theme_paper() +
    theme(
      strip.text = element_text(size = 8),  # 패싯 라벨용 작은 텍스트
      strip.background = element_rect(fill = "lightgrey"),
      panel.spacing = unit(0.3, "lines")    # 패싯 간 간격 줄이기
    )
}

# 출력
plot <- create_individual_bf_density(bf_data)

# 최적화된 행 기반 BF 분포 플롯 함수
create_optimized_bf_density <- function(data, size_group) {
  # size_group은 "n = 30", "n = 100", "n = 200" 중 하나여야 함
  
  plot_data <- add_scenario_info(data) %>%
    filter(sample_size_group == size_group)
  
  # 간결한 패싯 라벨 생성
  plot_data <- plot_data %>%
    mutate(
      # 각 조건의 간결한 라벨 생성
      condition_label = sprintf("n1=%d, n2=%d, SDR=%.1f", n1, n2, sdr)
    )
  
  # 조건 정렬 순서 생성
  ordered_conditions <- plot_data %>%
    select(condition_label, n1, n2, sdr) %>%
    distinct() %>%
    arrange(n1 == n2, sdr) %>%  # 먼저 할당 비율로 정렬, 그 다음 SDR로 정렬
    pull(condition_label)
  
  # 정렬된 팩터 생성
  plot_data <- plot_data %>%
    mutate(
      condition_ordered = factor(condition_label, levels = ordered_conditions)
    )
  
  # 조건을 행으로, 하나의 열로 배치하고 화면을 최대한 활용
  ggplot(plot_data, aes(x = log_BF, fill = model)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_grid(
      condition_ordered ~ .,  # 조건을 행으로 배치
      scales = "free_y",
      switch = "y"  # 행 라벨을 왼쪽으로 이동
    ) +
    scale_fill_manual(values = method_colors[c("jzs", "gica")],
                      labels = bf_method_labels) +
    xlim(-1.5, 6) +
    labs(
      title = paste("Distribution of log(BF) for", size_group, "Scenarios"),
      x = expression(log[10](BF)),
      y = "Density",
      fill = "Method"
    ) +
    theme_minimal() +  # 더 깔끔한 테마 사용
    theme(
      text = element_text(family = "Noto Sans KR", size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      
      # 행 라벨 스타일 최적화
      strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1),
      strip.placement = "outside",
      strip.background = element_blank(),  # 회색 박스 제거
      
      # 여백 및 비율 조정
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"),
      aspect.ratio = 0.25,  # 가로로 매우 넓게
      
      # Y축 그리드 제거
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

# 가장 넓은 버전의 플롯 함수
create_widest_bf_density <- function(data, size_group) {
  # size_group은 "n = 30", "n = 100", "n = 200" 중 하나여야 함
  
  plot_data <- add_scenario_info(data) %>%
    filter(sample_size_group == size_group)
  
  # 간결한 패싯 라벨 생성
  plot_data <- plot_data %>%
    mutate(
      # 각 조건의 간결한 라벨 생성
      condition_label = sprintf("n1=%d, n2=%d, SDR=%.1f", n1, n2, sdr)
    )
  
  # 조건 정렬 순서 생성
  ordered_conditions <- plot_data %>%
    select(condition_label, n1, n2, sdr) %>%
    distinct() %>%
    arrange(n1 == n2, sdr) %>%  # 먼저 할당 비율로 정렬, 그 다음 SDR로 정렬
    pull(condition_label)
  
  # 정렬된 팩터 생성
  plot_data <- plot_data %>%
    mutate(
      condition_ordered = factor(condition_label, levels = ordered_conditions)
    )
  
  # 조건을 행으로, 넓은 그래프로 설정
  p <- ggplot(plot_data, aes(x = log_BF, fill = model)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_grid(
      condition_ordered ~ .,  # 조건을 행으로 배치
      scales = "free_y",
      switch = "y"  # 행 라벨을 왼쪽으로 이동
    ) +
    scale_fill_manual(values = method_colors[c("jzs", "gica")],
                      labels = bf_method_labels) +
    xlim(-1.5, 6) +
    labs(
      title = paste("Distribution of log(BF) for", size_group, "Scenarios"),
      x = expression(log[10](BF)),
      y = "Density",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Noto Sans KR", size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      
      # 행 라벨 스타일 최적화
      strip.text.y.left = element_text(size = 10, angle = 0, hjust = 1),
      strip.placement = "outside",
      strip.background = element_blank(),
      
      # 여백 및 비율 조정
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
      
      # Y축 그리드 제거
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # 논문 용으로 저장할 때는 width를 크게 설정
  # ggsave 예시: ggsave("bf_n30.pdf", p, width = 10, height = 6)
  
  return(p)
}

# 사용 예시:
plot_n30 <- create_optimized_bf_density(bf_data, "n = 30")
plot_n100 <- create_optimized_bf_density(bf_data, "n = 100")
plot_n200 <- create_optimized_bf_density(bf_data, "n = 200")

# 또는 가장 넓은 버전:
plot_n30_wide <- create_widest_bf_density(bf_data, "n = 30")
plot_n100_wide <- create_widest_bf_density(bf_data, "n = 100")
plot_n200_wide <- create_widest_bf_density(bf_data, "n = 200")
plot_n30_wide;plot_n100_wide;plot_n200_wide


# 더 간단한 방식
bf_summary <- bf_data %>%
  left_join(scenarios, by = "scenario") %>%
  group_by(scenario, model, n1, n2, sdr) %>%
  summarize(
    mean_log_BF = mean(log_BF, na.rm = TRUE),
    .groups = "drop"
  )

# 새로운 열 추가
bf_summary <- bf_summary %>%
  mutate(
    # 총 샘플 크기 계산
    total_n = n1 + n2,
    
    # 행 그룹 생성
    row = case_when(
      between(scenario, 1, 5) ~ "(n=30)",
      between(scenario, 6, 10) ~ "(n=100)",
      between(scenario, 11, 15) ~ "(n=200)"
    ),
    
    # 열 그룹 생성 (각 행마다 1~5번 열에 해당)
    col = case_when(
      scenario %in% c(1, 6, 11) ~ "Column 1",
      scenario %in% c(2, 7, 12) ~ "Column 2",
      scenario %in% c(3, 8, 13) ~ "Column 3",
      scenario %in% c(4, 9, 14) ~ "Column 4",
      scenario %in% c(5, 10, 15) ~ "Column 5"
    ),
    
    # 패널 제목용 레이블
    panel_label = sprintf("n₁=%d, n₂=%d, SDR=%.1f", n1, n2, sdr)
  )

# 데이터 준비
bf_summary <- bf_data %>%
  left_join(scenarios, by = "scenario") %>%
  group_by(scenario, model, n1, n2, sdr) %>%
  summarize(
    mean_log_BF = mean(log_BF, na.rm = TRUE),
    .groups = "drop"
  )

# 행과 열 생성 (모두 = 기호 사용)
bf_summary <- bf_summary %>%
  mutate(
    # 행: 모두 등호(=) 사용
    row = case_when(
      scenario <= 5 ~ "n=30",  # ≈ 대신 = 사용
      scenario <= 10 ~ "n=100", 
      TRUE ~ "n=200"
    ),
    # 열: SDR 및 sample ratio 기준
    col = case_when(
      scenario %in% c(1, 6, 11) ~ "SDR=1.0, ratio=1:1, d=0.5",
      scenario %in% c(2, 7, 12) ~ "SDR=0.5, ratio=1:1, d=0.63", 
      scenario %in% c(3, 8, 13) ~ "SDR=1.0, ratio=2:3, d=0.5",
      scenario %in% c(4, 9, 14) ~ "SDR=0.5, ratio=2:3, d=0.63",
      scenario %in% c(5, 10, 15) ~ "SDR=2.0, ratio=2:3, d=032"
    )
  )

# 팩터 변환으로 순서 지정
bf_summary <- bf_summary %>%
  mutate(
    row = factor(row, levels = c("n=30", "n=100", "n=200")),
    col = factor(col, levels = c(
      "SDR=1.0, ratio=1:1, d=0.5", 
      "SDR=0.5, ratio=1:1, d=0.63", 
      "SDR=1.0, ratio=2:3, d=0.5", 
      "SDR=0.5, ratio=2:3, d=0.63", 
      "SDR=2.0, ratio=2:3, d=032"
    ))
  )

# 3×5 그래프 생성
ggplot(bf_summary, aes(x = model, y = mean_log_BF, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_grid(row ~ col, scales = "free_y") +  # scales = "free_y" 추가
  scale_fill_manual(values = c("jzs" = "turquoise3", "gica" = "coral"),
                    labels = bf_method_labels) +
  labs(
    title = "Average log10(Bayes Factor) Comparison: JZS vs GICA",
    x = "Bayes Factor Method",
    y = "Average log10(Bayes Factor)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = NA),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

#각 총표본 별로
# 각 표본 크기별로 별도의 그래프 생성
for (sample_size in c("n=30", "n=100", "n=200")) {
  p <- bf_summary %>%
    filter(row == sample_size) %>%
    ggplot(aes(x = model, y = mean_log_BF, fill = model)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    facet_wrap(~ col, nrow = 1) +
    scale_fill_manual(values = c("jzs" = "turquoise3", "gica" = "coral"),
                      labels = bf_method_labels) +
    labs(
      title = paste("Average log10(Bayes Factor) Comparison:", sample_size),
      x = "Bayes Factor Method",
      y = "Average log10(Bayes Factor)",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = NA),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )
  
  print(p)
  # 필요하다면 ggsave로 저장
  # ggsave(paste0("bf_comparison_", sample_size, ".png"), p, width = 10, height = 4)
}
