results_df %>% 
  select("scenario", "true_model", "rho", "sdr", "delta", 
         "BF_robtt_effect", "BF_bain_student", "BF_bain_welch", 
         "BF_bayesfactor", "BF_wetzels") %>%
  mutate(
    BF_bain_student = ifelse(!is.na(BF_bain_student), (BF_bain_student)^(-1), NA),
    BF_bain_welch = ifelse(!is.na(BF_bain_welch), (BF_bain_welch)^(-1), NA)
  ) -> result_inv

# NA가 발생한 조건 확인
result_long %>%
  filter(model == "wetzels", is.na(BF)) %>%
  group_by(rho, delta) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# 또는 전체 조건별 NA 비율
result_long %>%
  filter(model == "wetzels") %>%
  group_by(scenario, delta) %>%
  summarise(
    total = n(),
    na_count = sum(is.na(BF)),
    na_ratio = na_count/total
  ) %>%
  arrange(desc(na_ratio))

# 일관된 색상 팔레트 다시 정의
method_colors <- setNames(
  c("#95CAE4", "#7AC5CD", "#FFAEB9", "#FFB347", "#FF6464"),
  c("FBF_equal", "BF_JZS", "FBF_unequal", "BF_BMA", "BF_MCMC")
)

# method 이름을 보기 좋게 변경하는 함수
format_method_names <- function(methods) {
  methods %>%
    str_remove("BF_") %>%
    str_replace("bain_student", "FBF_equal") %>%
    str_replace("bain_welch", "FBF_unequal") %>%
    str_replace("bayesfactor", "BF_JZS") %>%
    str_replace("robtt_effect", "BF_BMA") %>%
    str_replace("wetzels", "BF_MCMC")
}

# 색상 팔레트 재정의
method_colors <- setNames(
  c("#7AC5CD", "#FF6464", "#FFB347", "#95CAE4", "#FFAEB9"),
  c("BF_JZS", "BF_MCMC", "BF_BMA", "FBF_equal", "FBF_unequal")
)

# method 이름 변경 및 순서 지정 함수
format_method_names <- function(methods) {
  methods <- methods %>%
    str_replace("BF_bain_student", "FBF_equal") %>%
    str_replace("BF_bain_welch", "FBF_unequal") %>%
    str_replace("BF_bayesfactor", "BF_JZS") %>%
    str_replace("BF_robtt_effect", "BF_BMA") %>%
    str_replace("BF_wetzels", "BF_MCMC")
  
  factor(methods, 
         levels = c("BF_JZS", "BF_MCMC", "BF_BMA", "FBF_equal", "FBF_unequal"))
}

plot_data <- result_inv %>%
  mutate(across(starts_with("BF_"), ~log10(.))) %>%
  pivot_longer(cols = starts_with("BF_"), 
               names_to = "method", 
               values_to = "log_BF") %>%
  group_by(scenario, method, sdr, delta) %>%
  summarise(
    mean_log_BF = mean(log_BF),
    se_log_BF = sd(log_BF) / sqrt(n()),
    ci_lower = mean_log_BF - 1.96 * se_log_BF,
    ci_upper = mean_log_BF + 1.96 * se_log_BF,
    rmse = (se_log_BF / abs(mean_log_BF)) * 100
  ) %>%
  mutate(method = format_method_names(method))

# 데이터 가공시에 NA 처리 수정
plot_NA <- result_inv %>%
  mutate(across(starts_with("BF_"), ~log10(.))) %>%
  pivot_longer(cols = starts_with("BF_"), 
               names_to = "method", 
               values_to = "log_BF") %>%
  group_by(scenario, method, sdr, delta) %>%
  summarise(
    mean_log_BF = mean(log_BF, na.rm = TRUE),  # na.rm = TRUE 추가
    se_log_BF = sd(log_BF, na.rm = TRUE) / sqrt(sum(!is.na(log_BF))),  # NA 아닌 값의 개수로 나눔
    ci_lower = mean_log_BF - 1.96 * se_log_BF,
    ci_upper = mean_log_BF + 1.96 * se_log_BF,
    rmse = (se_log_BF / abs(mean_log_BF)) * 100
  ) %>%
  mutate(method = format_method_names(method))

# delta = 0 그래프

p0 <- plot_NA %>%
  filter(delta == 0) %>%
  ggplot(aes(x = method, y = -mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = -ci_lower, ymax = -ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Log10(Bayes Factors) (delta = 0.0)",
       x = "Method", y = "-log10(Bayes Factor)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey90"),
    panel.border = element_rect(fill = NA, colour = "grey90")
  )

# delta = 0.3 그래프
p1 <- plot_NA %>%
  filter(delta == 0.3) %>%
  filter(!is.na(mean_log_BF)) %>%  # NA 값 제외
  ggplot(aes(x = method, y = -mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = -ci_lower, ymax = -ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Log10(Bayes Factors) (delta = 0.3)",
       x = "Method", y = "-log10(Bayes Factor)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey90"),
    panel.border = element_rect(fill = NA, colour = "grey90")
  )

# delta = 0.5 그래프
p2 <- plot_NA %>%
  filter(delta == 0.5) %>%
  filter(!is.na(mean_log_BF)) %>%  # NA 값 제외
  ggplot(aes(x = method, y = mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Log10(Bayes Factors) (delta = 0.5)",
       x = "Method", y = "log10(Bayes Factor)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey90"),
    panel.border = element_rect(fill = NA, colour = "grey90")
  )

# delta = 0.8 그래프
p3 <- plot_NA %>%
  filter(delta == 0.8) %>%
  filter(!is.na(mean_log_BF)) %>%  # NA 값 제외
  ggplot(aes(x = method, y = mean_log_BF, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.9), width = 0.25) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgray") +
  facet_wrap(~ scenario, scales = "free_y", 
             labeller = labeller(scenario = setNames(scenario_labels, 1:4))) +
  scale_fill_manual(values = method_colors) +
  labs(title = "Log10(Bayes Factors) (delta = 0.8)",
       x = "Method", y = "log10(Bayes Factor)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_line(colour = "grey90"),
    panel.border = element_rect(fill = NA, colour = "grey90")
  )

print(p0)
print(p1)
print(p2)
print(p3)
