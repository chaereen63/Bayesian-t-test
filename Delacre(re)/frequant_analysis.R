source("New/functionsN.R")
library(tidyverse)
library(gridExtra)

# 재현성 분석용 시뮬레이션 함수 (베이지안 요소 제거)
  # effect size = {0.2, 0.5, 0.8}
run_reproducibility_simulation <- function(current_settings) {
  # Generate data
  data <- simulate_data(
    mean1 = current_settings$mu1,  
    mean2 = current_settings$mu2,  
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
    varr = current_settings$varr
  )
  
  return(results)
}

# 시뮬레이션 실행
reproducibility_results <- all_settings %>%
  rowwise() %>%
  mutate(
    result = list(run_reproducibility_simulation(cur_data()))
  ) %>%
  ungroup()
# 필요한 열만 추출
analysis_data <- reproducibility_results %>%
  mutate(
    student_p = map_dbl(result, ~.x$student_p),
    welch_p = map_dbl(result, ~.x$welch_p)
    # 이미 있는 열은 추출하지 않음
  ) %>%
  select(-result)  # 원본 result 열 제거
# save
saveRDS(analysis_data, file = "Delacre(re)/frequant_analaysis_data8.RDS")

#### plot ####
create_delacre_style_plots <- function(data) {
  # 시나리오 정보
  scenario_info <- data %>%
    group_by(scenario, n1, n2, sd1, sd2, varr) %>%
    summarise(.groups = 'drop')

  # 각 시나리오에 라벨 추가
  scenario_labels <- c(
    "1" = sprintf("n1=%d n2=%d varr=1", scenario_info$n1[1], scenario_info$n2[1]),
    "2" = sprintf("n1=%d n2=%d varr=1", scenario_info$n1[2], scenario_info$n2[2]),
    "3" = sprintf("n1=%d n2=%d varr=0.5", scenario_info$n1[3], scenario_info$n2[3]),
    "4" = sprintf("n1=%d n2=%d varr=0.5", scenario_info$n1[4], scenario_info$n2[4]),
    "5" = sprintf("n1=%d n2=%d varr=0.5", scenario_info$n1[5], scenario_info$n2[5])
  )
  
  # 5개 시나리오를 모두 포함
  scenarios <- c("1", "2", "3", "4", "5")
  
  plots <- list()
  
  # 모든 시나리오에 대한 히스토그램 생성
  for(scen in scenarios) {
    scen_data <- data %>% filter(scenario == as.numeric(scen))
    
    # Student's t-test
    p_student <- ggplot(scen_data, aes(x = student_p)) +
      geom_histogram(bins = 30, fill = "skyblue", color = "darkblue", alpha = 0.7) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_y_continuous(limits = c(0, 2500)) +
      labs(
        title = paste("Student's t-test:", scenario_labels[scen]),
        x = "Observed p-value",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.grid.major = element_line(color = "lightgray"),
        panel.border = element_rect(color = "black", fill = NA)
      )
    
    # Welch's t-test
    p_welch <- ggplot(scen_data, aes(x = welch_p)) +
      geom_histogram(bins = 20, fill = "coral", color = "darkred", alpha = 0.7) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_y_continuous(limits = c(0, 2500)) +
      labs(
        title = paste("Welch's t-test:", scenario_labels[scen]),
        x = "Observed p-value",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.grid.major = element_line(color = "lightgray"),
        panel.border = element_rect(color = "black", fill = NA)
      )
    
    plots[[paste0(scen, "_student")]] <- p_student
    plots[[paste0(scen, "_welch")]] <- p_welch
  }
  
  # 2행 5열 레이아웃으로 배치
  # 상단 행: 모든 Student's t-test
  # 하단 행: 모든 Welch's t-test
  panel_plot <- gridExtra::grid.arrange(
    # 상단 행: Student's t-test
    plots[["1_student"]],
    plots[["2_student"]],
    plots[["3_student"]],
    plots[["4_student"]],
    plots[["5_student"]],
    # 하단 행: Welch's t-test
    plots[["1_welch"]],
    plots[["2_welch"]],
    plots[["3_welch"]],
    plots[["4_welch"]],
    plots[["5_welch"]],
    ncol = 5,  # 5개 열
    nrow = 2,  # 2개 행
    top = grid::textGrob(
      "Type I Error Rate Comparison Across All Scenarios",
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
  )
  
  return(panel_plot)
}
delacre_style_plots <- create_delacre_style_plots(analysis_data)
