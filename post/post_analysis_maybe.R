## 결과 파일 분석
#plotting의 result_long 이용

# 영가설이 참일 떄의 log BF값의 최댓값
(result_long %>% filter(true_model=="H0") %>% 
  pull(log_BF) %>%
  max() -> max_H0)
# 대안가설이 참일 떄의 log BF값의 최솟값
(result_long %>% filter(true_model=="H1") %>% 
  pull(log_BF) %>%
  min() -> min_H1)


# 오분류율 계산

# log10(3)와 log10(1/3) 계산
log_bf_threshold_positive <- log10(3)
log_bf_threshold_negative <- log10(1/3)

misclassification_rate <- result_long %>%
  mutate(misclassified = case_when(
    true_model == "H0" & log_BF > log_bf_threshold_negative ~ TRUE,
    true_model == "H1" & log_BF < log_bf_threshold_positive ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  group_by(model, scenario) %>%
  summarise(misclass_rate = mean(misclassified), .groups = 'drop')

print(misclassification_rate)

# 상세 분류 및 오류율 계산
detailed_classification <- result_long %>%
  mutate(classification = case_when(
    true_model == "H0" & log_BF > log_bf_threshold_positive ~ "Misclassified",
    true_model == "H1" & log_BF < log_bf_threshold_negative ~ "Misclassified",
    log_BF > log_bf_threshold_negative & log_BF < log_bf_threshold_positive ~ "Inconclusive",
    TRUE ~ "Correct"
  )) %>%
  group_by(model, scenario, true_model, classification) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(model, scenario, true_model) %>%
  mutate(total = sum(count),
         rate = count / total)

# 결과 출력
print(detailed_classification)

# 요약 통계

misclassification_summary <- detailed_classification %>%
  filter(classification != "Correct") %>%
  group_by(model, true_model, classification) %>%
  summarise(
    mean_rate = mean(rate),
    median_rate = median(rate),
    sd_rate = sd(rate),
    min_rate = min(rate),
    max_rate = max(rate),
    .groups = 'drop'
  ) %>%
  arrange(model, true_model, classification)

print(misclassification_summary)

# 전체 오분류율 (Inconclusive + Type I Error + Type II Error)
overall_misclassification <- detailed_classification %>%
  filter(classification != "Correct") %>%
  group_by(model, scenario, true_model) %>%
  summarise(total_misclassification = sum(rate), .groups = 'drop') %>%
  group_by(model) %>%
  summarise(
    mean_misclassification = mean(total_misclassification),
    median_misclassification = median(total_misclassification),
    sd_misclassification = sd(total_misclassification),
    min_misclassification = min(total_misclassification),
    max_misclassification = max(total_misclassification),
    .groups = 'drop'
  )

print(overall_misclassification)

write.csv(overall_misclassification, "overall_misclassification.csv")

  #시각화
  library(ggplot2)
  library(ggridges)
  
  ggplot(result_long, aes(x = log_BF, y = factor(scenario), fill = model)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_vline(xintercept = log_bf_threshold_positive, linetype = "dotted", color = "blue", size = 0.5) +
    geom_vline(xintercept = log_bf_threshold_negative, linetype = "dotted", color = "blue", size = 0.5) + #파란선은 BF가 0.33에서 3 사이의 범위
    facet_wrap(~ true_model, scales = "free") +
    geom_rect(data = data.frame(true_model = c("H0", "H1")),
              aes(xmin = ifelse(true_model == "H0", log_bf_threshold_negative, -Inf),
                  xmax = ifelse(true_model == "H0", Inf, log_bf_threshold_positive),
                  ymin = -Inf, ymax = Inf),
              fill = "salmon", alpha = 0.2, inherit.aes = FALSE) +
    theme_minimal() +
    labs(title = "Distribution of log10(Bayes Factors) with Misclassification Regions",
         subtitle = "Pink areas: Misclassified regions (H0: log10(BF) > -0.477, H1: log10(BF) < 0.477)",
         x = "log10(Bayes Factor)", y = "Scenario", fill = "Model") +
    theme(legend.position = "bottom")
  
  # 시각화_detail
  ggplot(detailed_classification, aes(x = factor(scenario), y = rate, fill = classification)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(true_model ~ model) +
    theme_minimal() +
    labs(title = "Detailed Classification of BF Results",
         subtitle = "Inconclusive: -0.477 < log10(BF) < 0.477",
         x = "Scenario", y = "Rate", fill = "Classification") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_text(size = 12),
          strip.text = element_text(size = 12, face = "bold")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("Correct" = "#91C483", "Inconclusive" = "#FFE162", 
                                 "Misclassified" = "#FF6464"))
  
  
  # 모형 간, 모형 내 오분류율 원인 분석
  library(tidyverse)
  library(broom)
  library(emmeans)
  library(knitr)
  
  # 데이터 준비 (result_long이 있다고 가정)
  misclassification_data <- result_long %>%
    mutate(
      misclassified = case_when(
        true_model == "H0" & log_BF > log_bf_threshold_negative ~ TRUE,
        true_model == "H1" & log_BF < log_bf_threshold_positive ~ TRUE,
        TRUE ~ FALSE
      ),
      scenario = as.factor(scenario),
      true_model = as.factor(true_model)
    )
  

    # 로지스틱 회귀
  overall_model <- glm(misclassified ~ model * scenario*true_model, data = misclassification_data, family = binomial)
  overall_emm <- emmeans(overall_model, ~model|scenario | true_model, type = "response")
  overall_pairs <- pairs(overall_emm)
  
  cat("\nOverall Model Comparison:\n")
  print(summary(overall_model))
  cat("\nEstimated Marginal Means for Each Model and Scenario:\n")
  print(overall_emm)
  cat("\nPairwise Comparisons between Models for Each Scenario:\n")
  print(overall_pairs)
  
  # 데이터 준비 및 표 생성
  plot_data <- as.data.frame(overall_emm)
  table_data <- plot_data %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    rename("Model" = model,
           "Scenario" = scenario,
           "Predicted Rate" = prob,
           "Lower CI" = asymp.LCL,
           "Upper CI" = asymp.UCL) 
  
  print(table_data, row.names = FALSE)
  
  # 시각화
  ggplot(plot_data, aes(x = scenario, y = prob, fill = true_model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                  position = position_dodge(width = 0.9), width = 0.25) +
    facet_wrap(~ model, scales = "free_y") +
    scale_fill_manual(values = c("H0" = "skyblue", "H1" = "lightgreen")) +
    theme_minimal() +
    labs(title = "Predicted Misclassification Rate by Model, Scenario, and True Model",
         x = "Scenario",
         y = "Predicted Misclassification Rate",
         fill = "True Model") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  