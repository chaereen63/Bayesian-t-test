source("New/functionsN.R")
library(tidyverse)
library(gridExtra)

# 시나리오 정의 - 4개로 줄임
scenarios <- list(
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 2),    # SDR=1
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 4),    # SDR=2
  list(n1 = 40, n2 = 60, sd1 = 2, sd2 = 1),    # SDR=0.5
  list(n1 = 50, n2 = 50, sd1 = 2, sd2 = 4)     # 동일 표본 크기
)

# 설정 생성 함수
create_settings <- function(scenario, replications) {
  tibble(
    scenario = scenario,
    n1 = scenarios[[scenario]]$n1,
    n2 = scenarios[[scenario]]$n2,
    sd1 = scenarios[[scenario]]$sd1,
    sd2 = scenarios[[scenario]]$sd2,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 재현성 분석용 설정 (평균차 = 0)
reproducibility_settings <- tibble(scenario = 1:4) %>%
  mutate(
    settings = map(scenario, ~create_settings(.x, replications = 1000000))  # 10000회로 증가
  ) %>%
  select(-scenario) %>%
  unnest(settings) %>%
  mutate(
    mean1 = 0,
    mean2 = 0,
    sdr = sd2/sd1,
    cohens_d = 0
  )

# 재현성 분석용 시뮬레이션 함수 (베이지안 요소 제거)
run_reproducibility_simulation <- function(current_settings) {
  # Generate data
  data <- simulate_data(
    mean1 = current_settings$mean1,
    mean2 = current_settings$mean2,
    sd1   = current_settings$sd1,
    sd2   = current_settings$sd2,
    n1    = current_settings$n1,
    n2    = current_settings$n2,
    seed  = current_settings$seed
  )
  
  # 빈도주의 테스트만 실행
  student_p <- t.test(data$x1, data$x2, var.equal=TRUE)$p.value
  welch_p <- t.test(data$x1, data$x2, var.equal=FALSE)$p.value
  
  # 결과 수집 (베이지안 요소 제거)
  results <- list(
    student_p = student_p,
    welch_p = welch_p,
    scenario = current_settings$scenario,
    n1 = current_settings$n1,
    n2 = current_settings$n2,
    sd1 = current_settings$sd1,
    sd2 = current_settings$sd2,
    sdr = current_settings$sdr
  )
  
  return(results)
}

# 시뮬레이션 실행
reproducibility_results <- reproducibility_settings %>%
  rowwise() %>%
  mutate(
    result = list(run_reproducibility_simulation(cur_data()))
  ) %>%
  ungroup()

# 결과 데이터 추출 및 정리
analysis_data <- reproducibility_results %>%
  mutate(
    student_p = map_dbl(result, ~.x$student_p),
    welch_p = map_dbl(result, ~.x$welch_p),
    n1 = map_dbl(result, ~.x$n1),
    n2 = map_dbl(result, ~.x$n2),
    sd1 = map_dbl(result, ~.x$sd1),
    sd2 = map_dbl(result, ~.x$sd2),
    sdr = map_dbl(result, ~.x$sdr)
  ) %>%
  select(-result)  # 원본 result 열 제거

# 결과 저장 (최종 분석 데이터)
saveRDS(analysis_data, file = "Delacre(re)/reproduct_analysis_data.RDS")

# 원문과 유사한 레이아웃의 그래프 생성 함수
create_delacre_style_plots <- function(data) {
  # 시나리오 정보
  scenario_info <- data %>%
    group_by(scenario, n1, n2, sd1, sd2, sdr) %>%
    summarise(.groups = 'drop')
  
  # 각 시나리오에 라벨 추가
  scenario_labels <- c(
    "1" = sprintf("n1=%d n2=%d SDR=1", scenario_info$n1[1], scenario_info$n2[1]),
    "2" = sprintf("n1=%d n2=%d SDR=2", scenario_info$n1[2], scenario_info$n2[2]),
    "3" = sprintf("n1=%d n2=%d SDR=0.5", scenario_info$n1[3], scenario_info$n2[3]),
    "4" = sprintf("n1=%d n2=%d SDR=2", scenario_info$n1[4], scenario_info$n2[4])
  )
  
  # 2x2 패널 레이아웃을 위한 시나리오 쌍 정의
  panel_pairs <- list(
    list("1", "4"),   # n1=40,n2=60과 n1=50,n2=50 비교
    list("2", "3")  # SDR=2와 SDR=0.5 비교 (원문 Fig 2c & 2d)
  )
  
  all_plots <- list()
  
  # 각 패널 쌍에 대해 그래프 생성
  for(i in seq_along(panel_pairs)) {
    pair <- panel_pairs[[i]]
    
    plots <- list()
    # 각 시나리오에 대한 Student's t-test와 Welch's t-test 히스토그램 생성
    for(scen in pair) {
      scen_data <- data %>% filter(scenario == as.numeric(scen))
      
      # Student's t-test
      p_student <- ggplot(scen_data, aes(x = student_p)) +
        geom_histogram(bins = 20, fill = "skyblue", color = "darkblue", alpha = 0.7) +
        geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
        labs(
          title = paste("Student's t-test:", scenario_labels[scen]),
          x = "Observed p-value",
          y = "Frequency"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          panel.grid.major = element_line(color = "lightgray"),
          panel.border = element_rect(color = "black", fill = NA)
        )
      
      # Welch's t-test
      p_welch <- ggplot(scen_data, aes(x = welch_p)) +
        geom_histogram(bins = 20, fill = "coral", color = "darkred", alpha = 0.7) +
        geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
        labs(
          title = paste("Welch's t-test:", scenario_labels[scen]),
          x = "Observed p-value",
          y = "Frequency"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          panel.grid.major = element_line(color = "lightgray"),
          panel.border = element_rect(color = "black", fill = NA)
        )
      
      plots[[paste0(scen, "_student")]] <- p_student
      plots[[paste0(scen, "_welch")]] <- p_welch
    }
    
    # 2x2 레이아웃 (왼쪽: SDR=2, 오른쪽: SDR=0.5)
    # 위: Student's t-test, 아래: Welch's t-test
    panel_plot <- gridExtra::grid.arrange(
      plots[[paste0(pair[1], "_student")]],
      plots[[paste0(pair[2], "_student")]],
      plots[[paste0(pair[1], "_welch")]],
      plots[[paste0(pair[2], "_welch")]],
      ncol = 2,
      top = grid::textGrob(
        paste("Type I Error Rate Comparison - Panel", i),
        gp = grid::gpar(fontsize = 14, fontface = "bold")
      )
    )
    
    all_plots[[i]] <- panel_plot
  }
  
  return(all_plots)
}

# 원문 스타일 그래프 생성 및 저장
delacre_style_plots <- create_delacre_style_plots(analysis_data)

# 그래프 저장
for(i in seq_along(delacre_style_plots)) {
  ggsave(
    filename = paste0("Delacre(re)/delacre_style_panel_", i, ".png"),
    plot = delacre_style_plots[[i]],
    width = 12, height = 10
  )
}

# Type I error rate 계산 및 저장
type1_error <- analysis_data %>%
  group_by(scenario, n1, n2, sd1, sd2, sdr) %>%
  summarise(
    student_type1 = mean(student_p < 0.05, na.rm = TRUE),
    welch_type1 = mean(welch_p < 0.05, na.rm = TRUE),
    n_sims = n(),
    .groups = 'drop'
  )

saveRDS(type1_error, file = "Delacre(re)/reproduct_type1_error.RDS")
print(type1_error)