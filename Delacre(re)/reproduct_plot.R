# 원문과 유사한 레이아웃의 그래프 생성 함수
create_delacre_style_plots <- function(data) {
  # 시나리오 정보
  scenario_info <- data %>%
    group_by(scenario, n1, n2, sd1, sd2, sdr) %>%
    summarise(.groups = 'drop')
  
  # 각 시나리오에 라벨 추가 (8개 시나리오용으로 확장)
  scenario_labels <- character(8)
  for (i in 1:8) {
    scenario_labels[i] <- sprintf("n1=%d n2=%d sd1=%.1f sd2=%.1f", 
                                  scenario_info$n1[i], 
                                  scenario_info$n2[i], 
                                  scenario_info$sd1[i], 
                                  scenario_info$sd2[i])
  }
  names(scenario_labels) <- as.character(1:8)
  
  # 시나리오를 그룹화하여 패널 쌍 정의
  # 비슷한 특성을 가진 시나리오를 쌍으로 묶음
  panel_pairs <- list(
    list("1", "3"),    # arith_0.5, 다른 표본크기
    list("3", "6"),
    list("2", "4"),    # arith_0.632
    list("4", "7"),    # arith_0.632
    list("5", "8")    # arith_0.316
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
        geom_histogram(bins = 30, fill = "skyblue", color = "darkblue", alpha = 0.7) +
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
        geom_histogram(bins = 30, fill = "coral", color = "darkred", alpha = 0.7) +
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
    
    # 2x2 레이아웃 구성
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

# 모든 시나리오에 대한 요약 그래프 생성
summary_plot <- type1_error %>%
  pivot_longer(
    cols = c(student_type1, welch_type1),
    names_to = "test_type",
    values_to = "type1_rate"
  ) %>%
  mutate(
    test_type = factor(test_type, 
                       levels = c("student_type1", "welch_type1"),
                       labels = c("Student's t-test", "Welch's t-test"))
  ) %>%
  ggplot(aes(x = factor(scenario), y = type1_rate, fill = test_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Type I Error Rate by Scenario and Test Type",
    x = "Scenario",
    y = "Type I Error Rate",
    fill = "Test Type"
  ) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.01)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "lightgray"),
    panel.border = element_rect(color = "black", fill = NA)
  )

# 요약 그래프 저장
ggsave(
  filename = "Delacre(re)/type1_error_summary.png",
  plot = summary_plot,
  width = 14, height = 8
)

# 시나리오 정보 테이블 생성 및 저장
scenario_table <- type1_error %>%
  select(scenario, n1, n2, sd1, sd2, sdr, student_type1, welch_type1) %>%
  arrange(scenario) %>%
  mutate(
    student_type1 = round(student_type1, 4),
    welch_type1 = round(welch_type1, 4),
    type1_diff = round(student_type1 - welch_type1, 4)
  )

# 테이블을 CSV로 저장
write.csv(scenario_table, "Delacre(re)/scenario_summary_table.csv", row.names = FALSE)