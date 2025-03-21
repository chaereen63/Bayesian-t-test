## 시각화
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr);library(gridExtra)
home_dir <- "."
source(file = file.path("./New/functionsN.R"))

# Load all three RDS files
results_30 <- readRDS("./New/mergedFin30ES8.RDS")
results_100 <- readRDS("./New/mergedFin100ES8.RDS")
results_200 <- readRDS("./New/mergedFin200ES8.RDS")

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
    # 레이블에 b값 포함 (간결한 형식)
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
    # 레이블에 b값 포함 (더 간결한 형식)
    label = sprintf("n1=%d, n2=%d", 
                    n1, n2, varr)
  )
# 색상 팔레트와 레이블 정의
method_colors <- c(
  "jzs" = "turquoise3",
  "gica" = "coral"
)

# 테두리 색상 정의 (약간 더 진한 색상)
method_border_colors <- c(
  "jzs" = "turquoise4",
  "gica" = "coral4"
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

# 빈도분포 히스토그램 함수 정의
plot_bf_histogram <- function(bf_tidy, title) {
  # 각 시나리오별로 분리하여 히스토그램 생성
  ggplot(bf_tidy, aes(x = log_BF, fill = method_short)) +
    # 히스토그램의 투명도를 0.5로 낮추고 테두리 추가
    geom_histogram(aes(y = after_stat(density)), 
                   alpha = 0.5, 
                   position = "identity",
                   bins = 30,
                   color = "white",  # 히스토그램 사이 경계선을 흰색으로
                   linewidth = 0.2) +
    # 밀도 곡선의 투명도를 0.7로 올리고 선 두께를 1.2로 증가
    geom_density(aes(color = method_short), 
                 alpha = 0.7,
                 linewidth = 1.2) +
    facet_wrap(~ scenario, 
               scales = "fixed",
               labeller = labeller(scenario = setNames(bf_tidy$label[!duplicated(bf_tidy$scenario)], 
                                                       unique(bf_tidy$scenario)))) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_color_manual(values = method_border_colors,  # 더 진한 테두리 색상 사용
                       labels = bf_method_labels) +
    labs(
      title = paste("로그 베이즈 인자 분포", title),
      x = "log(BF)",
      y = "밀도",
      fill = "방법",
      color = "방법"#,
      #caption = "Note: "  # 각주 추가
    ) +
    theme_paper()
}

# 박스플롯 함수 정의
plot_bf_boxplot <- function(bf_tidy, title) {
  ggplot(bf_tidy, aes(x = factor(scenario), y = log_BF, fill = method_short)) +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    scale_x_discrete(labels = bf_tidy$label[!duplicated(bf_tidy$scenario)]) +
    labs(
      title = paste("시나리오별 Log BF 비교", title),
      x = "시나리오",
      y = "log(BF)",
      fill = "방법"
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
    )
}

# === results_30 ===
results_30_tidy <- process_bf_data(results_30, scenarios_30) # 데이터 준비
hist_30 <- plot_bf_histogram(results_30_tidy, "(N=30, Effect size = 0.5)") # 히스토그램 생성
box_30 <- plot_bf_boxplot(results_30_tidy, "(N=30, Effect size = 0.8)") # 박스플롯 생성
print(hist_30)
print(box_30)
# === results_100 ===
results_100_tidy <- process_bf_data(results_100, scenarios_100) # 데이터 준비
hist_100 <- plot_bf_histogram(results_100_tidy, "(N=100, Effect size = 0.5)") # 히스토그램 생성
box_100 <- plot_bf_boxplot(results_100_tidy, "(N=100, Effect size = 0.8)") # 박스플롯 생성
print(hist_100)
print(box_100)

# === results_200 ===
results_200_tidy <- process_bf_data(results_200, scenarios_200)
hist_200 <- plot_bf_histogram(results_200_tidy, "(N=200, Effect size = 0.5)")
box_200 <- plot_bf_boxplot(results_200_tidy, "(N=200, Effect size = 0.8)")
print(hist_200)
print(box_200)

# === 추가 분석: BF_jzs와 BF_gica의 비교 분석 ===

# 방법별 요약 통계 계산
summary_stats_200 <- results_200_tidy %>% # 표본 바꾸기
  group_by(scenario, method_short) %>%
  summarise(
    mean_log_BF = mean(log_BF),
    median_log_BF = median(log_BF),
    sd_log_BF = sd(log_BF),
    q25 = quantile(log_BF, 0.25),
    q75 = quantile(log_BF, 0.75),
    .groups = "drop"
  )

print(summary_stats_200)

# 다른 방법: 로그 차이 사용 (log(BF_GICA) - log(BF_JZS))
log_bf_diff <- function(results_df, scenarios_df, title) {
  # 로그 베이즈 인자 차이 계산 및 데이터 준비
  scatter_data <- results_df %>%
    mutate(
      log_BF_gica = log10(BF_gica),  # 로그 변환
      log_BF_jzs = log10(BF_jzs),    # 로그 변환
      log_bf_diff = log_BF_gica - log_BF_jzs  # 로그 BF 차이
    ) %>%
    left_join(scenarios_df, by = "scenario")
  
  # 시나리오별 산점도 생성
  ggplot(scatter_data, aes(x = -mean_diff, y = log_bf_diff)) +
    geom_point(alpha = 0.3, size = 0.8) +  # 점 투명도와 크기 설정
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +  # 차이 = 0 기준선
    #geom_smooth(method = "loess", color = "blue", se = FALSE, linewidth = 0.8) +  # LOESS 회귀선 추가
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(scatter_data$label[!duplicated(scatter_data$scenario)], 
                                                       unique(scatter_data$scenario))),
               nrow = 2, ncol = 3) +  # 2행 4열 배치
    labs(
      title = paste("로그 베이즈 인자 차이 (log(BF_GICA) - log(BF_JZS))", title),
      x = "평균 차이 (mean1-mean2)",
      y = "log(BF) 차이 (GICA - JZS)"
    ) +
    theme_paper() +
    theme(
      aspect.ratio = 0.6,  # 가로세로 비율
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      strip.text = element_text(size = 9, margin = margin(b = 5)),
      panel.spacing = unit(1, "lines")
    )
}

# 로그 차이 산점도 생성 및 출력
log_diff_30 <- log_bf_diff(results_30, scenarios_30, "(N=30, Effect size = 0.8)")
log_diff_100 <- log_bf_diff(results_100, scenarios_100, "(N=100, Effect size = 0.8)")
log_diff_200 <- log_bf_diff(results_200, scenarios_200, "(N=200, Effect size = 0.8)")
print(log_diff_30);print(log_diff_100);print(log_diff_200)

# 시나리오별 평균 막대그래프 함수 정의 (y축 스케일 자유롭게)
plot_mean_bar <- function(bf_tidy, title) {
  # 시나리오별, 방법별 로그 BF 평균값 계산
  summary_stats <- bf_tidy %>%
    group_by(scenario, method_short, label) %>%
    summarise(
      mean_log_BF = mean(log_BF),
      se_log_BF = sd(log_BF) / sqrt(n()),
      .groups = "drop"
    )
  
  # 시나리오를 개별 패널로 나누어 그래프 생성 (y축 스케일 자유롭게)
  ggplot(summary_stats, aes(x = method_short, y = mean_log_BF, fill = method_short)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean_log_BF - se_log_BF, ymax = mean_log_BF + se_log_BF),
                  width = 0.25, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    scale_fill_manual(values = method_colors, 
                      labels = bf_method_labels) +
    facet_wrap(~ scenario, scales = "fixed",  # 각 패널마다 y축 자유롭게 조정
               labeller = labeller(scenario = setNames(summary_stats$label[!duplicated(summary_stats$scenario)], 
                                                       unique(summary_stats$scenario))),
               nrow = 2, ncol = 3) +  # 2행 4열 배치
    labs(
      title = paste("시나리오별 평균 Log BF 비교", title),
      subtitle = "오차 막대는 표준 오차(SE)를 나타냄",
      x = "방법",
      y = "평균 log(BF)",
      fill = "방법",
      caption = "Note : log_10을 취함"
    ) +
    theme_paper() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),  # x축 텍스트 정렬 조정
      strip.text = element_text(size = 9),  # 패널 제목 크기 조정
      panel.spacing = unit(1, "lines")  # 패널 간격 조정
    )
}
mean_bar_30 <- plot_mean_bar(results_30_tidy, "(N=30, Effect size = 0.8)")
mean_bar_100 <- plot_mean_bar(results_100_tidy, "(N=100, Effect size = 0.8)")
mean_bar_200 <- plot_mean_bar(results_200_tidy, "(N=200, Effect size = 0.8)")
# 그래프 출력
print(mean_bar_30)
print(mean_bar_100)
print(mean_bar_200)

#### 아래부터는 심심해서 해본 그림 그리기 ####
# 시나리오별 BF_jzs와 BF_gica 간의 상관관계 분석 함수 (파랑-노랑-빨강 색상 대비)
analyze_correlation <- function(results_df, scenarios_df, title) {
  # 데이터 준비
  corr_data <- results_df %>%
    select(scenario, BF_jzs, BF_gica, mean_diff) %>%
    mutate(
      log_BF_jzs = log10(BF_jzs),
      log_BF_gica = log10(BF_gica),
      abs_mean_diff = abs(mean_diff)
    ) %>%
    left_join(scenarios_df, by = "scenario")
  
  # 절댓값 범위 확인
  max_abs_diff <- max(corr_data$abs_mean_diff, na.rm = TRUE)
  
  # 시나리오별 상관계수 계산
  corr_stats <- corr_data %>%
    group_by(scenario) %>%
    summarise(
      pearson_r = cor(log_BF_jzs, log_BF_gica, method = "pearson"),
      r_squared = pearson_r^2,
      n_obs = n(),
      label = first(label),
      .groups = "drop"
    )
  
  # 산점도 생성 - 파랑-노랑-빨강 색상 대비
  scatter_plot <- ggplot(corr_data, aes(x = log_BF_jzs, y = log_BF_gica)) +
    # 점 색상을 abs_mean_diff에 따라 설정
    geom_point(aes(color = abs_mean_diff), alpha = 0.7, size = 1.2) +
    # 1:1 선 연한 회색으로 표시
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.5) +
    # 색상 스케일 설정 - 파랑-노랑-빨강 대비
    scale_color_gradientn(
      colors = c("darkblue", "blue", "royalblue", "skyblue", 
                 "yellow", 
                 "orange", "red", "darkred"),
      values = scales::rescale(c(0, max_abs_diff*0.2, max_abs_diff*0.3, max_abs_diff*0.4, 
                                 max_abs_diff*0.5, max_abs_diff*0.6, max_abs_diff*0.8, max_abs_diff)),
      name = "|평균 차이(SMD)|"
    ) +
    # 각 시나리오별 패널 분리
    facet_wrap(~ scenario, scales = "fixed",
               labeller = labeller(scenario = setNames(
                 corr_stats$label, 
                 corr_stats$scenario
               )),
               nrow = 2, ncol = 3) +
    # 각 패널에 상관계수 표시
    geom_text(data = corr_stats,
              aes(label = sprintf("R² = %.3f", r_squared),
                  x = -Inf, y = Inf),
              hjust = -0.1, vjust = 1.2, size = 3) +
    # 그래프 제목 및 축 레이블
    labs(
      title = paste("Log BF 방법 간 상관관계", title),
      subtitle = "점선: 1:1 일치선, 점 색상: 평균 차이의 절댓값(|SMD|)",
      x = expression(log(BF[JZS])),
      y = expression(log(BF[GICA])),
      caption = "Note: R² = 결정계수"
    ) +
    theme_paper() +
    theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom"
    )
  
  return(scatter_plot)
}

# 각 데이터셋에 대한 상관관계 분석 실행
scatter_30 <- analyze_correlation(results_30, scenarios_30, "(N=30, Effect size = 0.8)")
scatter_100 <- analyze_correlation(results_100, scenarios_100, "(N=100, Effect size = 0.8)")
scatter_200 <- analyze_correlation(results_200, scenarios_200, "(N=200, Effect size = 0.8)")

# 결과 출력
print(scatter_30)
print(scatter_100)
print(scatter_200)
