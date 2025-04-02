library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)
library(gridExtra)

# 효과크기별 결과 파일 로드
results_100_es8 <- readRDS("./New/mergedFin200ES8_r1.RDS")  # 효과크기 0.8
results_100_es5 <- readRDS("./New/mergedFin200ES5_r1.RDS")  # 효과크기 0.5
results_100_es2 <- readRDS("./New/mergedFin200ES2_r1.RDS")  # 효과크기 0.2
results_100_es0 <- readRDS("./New/mergedFin200ES0_r1.RDS")  # 효과크기 0.0

# 시나리오 정의 - sample size = 50인 경우
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
    # 레이블에 b값 포함
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
    # 레이블에 b값 포함
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
    # 레이블에 b값 포함
    label = sprintf("n1=%d, n2=%d, VarR=%.1f", 
                    n1, n2, varr)
  )
# 색상 팔레트와 레이블 정의
method_colors <- c(
  "jzs" = "#00366C",
  "gica" = "#F2A900"
)

bf_method_labels <- c(
  "jzs" = expression(BF[JZS]),
  "gica" = expression(BF[GICA])
)

# 공통 테마 설정
theme_paper <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman", size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "bottom",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1, "lines"),
      plot.margin = margin(t = 1, r = 5, b = 1, l = 5, unit = "pt")
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

# 시나리오별 평균 막대그래프 함수 정의
create_plot <- function(results_df, title, is_left_column = FALSE, add_y_axis = FALSE, y_limits = NULL) {
  # 데이터 준비
  bf_tidy <- process_bf_data(results_df, scenarios_200)
  
  # 시나리오별, 방법별 로그 BF 평균값 계산
  summary_stats <- bf_tidy %>%
    group_by(scenario, method_short, label) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      .groups = "drop"
    )
  
  # y축 라벨 설정 - 모든 그래프에서 제거하고 나중에 하나만 추가
  y_title <- NULL
  
  # 그래프 생성
  p <- ggplot(summary_stats, aes(x = method_short, y = mean_log_BF, fill = method_short)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_log_BF - se_log_BF, ymax = mean_log_BF + se_log_BF),
                  width = 0.25, linewidth = 0.1) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.2) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_x_discrete(labels = c("jzs" = "JZS", "gica" = "GICA")) +
    # 각 효과 크기 별로 적절한 y축 눈금 설정
    scale_y_continuous(breaks = function(limits) pretty(limits, n = 5)) +
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(
                 summary_stats$label[!duplicated(summary_stats$scenario)],
                 unique(summary_stats$scenario)
               ))) +
    labs(
      title = title,
      y = y_title,
      x = NULL
    ) +
    theme_paper() +
    theme(
      strip.text = element_text(size = 8),
      # 1. x축의 두 방법 글씨 크기 축소
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
      # 2. y축의 단위 글씨 크기 축소
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 5, t = 2)),
      legend.position = "none",
      # 모든 그래프의 y축 제목 제거
      axis.title.y = element_blank()
    )
  
  # 각 효과 크기 그룹별로 y축 범위 설정
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

# y축 범위 설정
# 효과 크기 0.8 그룹의 y축 범위
es_0_8_y_limits <- c(0, 7) # 50: 1.25, 100: 2.7, 200: 7
# 효과 크기 0.5 그룹의 y축 범위 - 값이 더 작아 패턴이 잘 보이도록 범위 좁힘
es_0_5_y_limits <- c(0, 2) # 50: 0.3, 100: 0.8, 200: 2
# 효과 크기 0.2 및 0.0 그룹의 y축 범위
negative_y_limits <- c(-0.8, 0) # 50: -0.6, 100: -0.7, 200: -0.8

# 4개 그래프 생성 - 효과 크기별로 적절한 y축 범위 적용
p1 <- create_plot(results_100_es0, "(A) Effect size = 0.0", FALSE, y_limits = negative_y_limits)
p2 <- create_plot(results_100_es2, "(B) Effect size = 0.2", TRUE, y_limits = negative_y_limits)
p3 <- create_plot(results_100_es5, "(C) Effect size = 0.5", FALSE, y_limits = es_0_5_y_limits)
p4 <- create_plot(results_100_es8, "(D) Effect size = 0.8", TRUE, y_limits = es_0_8_y_limits)

# 범례 생성
legend_plot <- ggplot(process_bf_data(results_100_es8, scenarios_200), 
                      aes(x = method_short, y = log_BF, fill = method_short)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = method_colors, labels = bf_method_labels) +
  labs(fill = "Method") +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 0, 0, 0),
    legend.spacing.x = unit(0.2, "cm"),
    legend.title = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman")
  ) +
  guides(fill = guide_legend(nrow = 1))

# 범례 추출
legend <- get_legend(legend_plot)

# cowplot으로 2x2 그리드 생성
combined_plot <- plot_grid(
  p1, p2, p3, p4, 
  ncol = 2, 
  nrow = 2,
  align = 'v'
)

# y축 라벨을 추가하기 위한 함수
add_shared_y_label <- function(plot_grid, y_label = "Mean log(BF)") {
  # 3. y축의 레이블 글씨 크기 축소
  y_grob <- textGrob(
    y_label, 
    rot = 90, 
    gp = gpar(
      fontfamily = "Times New Roman",
      fontsize = 10
    ),
    vjust = 0.5
  )
  
  # 그리드에 y축 라벨 추가
  grid.arrange(
    arrangeGrob(
      plot_grid,
      left = y_grob,
      padding = unit(1, "lines")
    )
  )
}

# 최종 그래프 (범례 포함) - 추가 여백 설정
final_plot_with_label <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(20, 1)
)

# 공유 y축 라벨 추가
final_plot_with_y_label <- add_shared_y_label(final_plot_with_label)

# 최종 그래프 출력
print(final_plot_with_y_label)

# 최종 그래프 저장
ggsave("combined_plot_n200r1_final.png", final_plot_with_y_label, width = 23, height = 14, dpi = 600, units = "cm")
