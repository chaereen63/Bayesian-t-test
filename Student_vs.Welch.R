results_df2 %>% 
  mutate(
    # 시나리오별 n1, n2 할당
    n1 = case_when(
      scenario == 1 ~ 40,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 40,
      scenario == 5 ~ 40
    ),
    n2 = case_when(
      scenario == 1 ~ 60,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 60,
      scenario == 5 ~ 60
    ),
    # Student's t-test용 표준화된 평균차이
    std_mean_diff = mean_diff/pooled_sd(sd_x, sd_y, n1, n2)
  ) %>%
  select(BF_jzs, BF_gica, std_mean_diff, sdr, delta, scenario) -> results_S

# Welch용도 같은 방식으로:
results_df2 %>% 
  mutate(
    n1 = case_when(
      scenario == 1 ~ 40,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 40,
      scenario == 5 ~ 40
    ),
    n2 = case_when(
      scenario == 1 ~ 60,
      scenario == 2 ~ 50,
      scenario == 3 ~ 50,
      scenario == 4 ~ 60,
      scenario == 5 ~ 60
    ),
    # Welch's t-test용 표준화된 평균차이
    std_mean_diff = (mean_diff/approx_sd(sd_x, sd_y, n1, n2))/sqrt(n1+n2)
  ) %>%
  select(BF_jzs, BF_gica, std_mean_diff, sdr, delta, scenario) -> results_W
head(results_S)
head(results_W)

## 1. Student
result_ratio <- results_S %>%
  mutate(bf_ratio = BF_gica/BF_jzs) %>%
  select(scenario, std_mean_diff, bf_ratio)

# BF_ratio scatter plot with alpha
ggplot(result_ratio, aes(x = std_mean_diff, y = bf_ratio)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +  # ratio = 1 기준선
  geom_point(alpha = 0.1) +
  facet_wrap(~scenario) +
  scale_y_log10() +  # log scale로 보면 1을 중심으로 대칭적으로 볼 수 있음
  labs(
    x = "Standardized Mean Difference_Student",
    y = "BF Ratio (GICA/JZS)"
  )
## 2. Welch
result_ratio <- results_W %>%
  mutate(bf_ratio = BF_gica/BF_jzs) %>%
  select(scenario, std_mean_diff, bf_ratio)

# BF_ratio scatter plot with alpha
ggplot(result_ratio, aes(x = std_mean_diff, y = bf_ratio)) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +  # ratio = 1 기준선
  geom_point(alpha = 0.1) +
  facet_wrap(~scenario) +
  scale_y_log10() +  # log scale로 보면 1을 중심으로 대칭적으로 볼 수 있음
  labs(
    x = "Standardized Mean Difference_Welch",
    y = "BF Ratio (GICA/JZS)"
  )
####Regression####
# 시나리오별 R² 계산
# 1. Student
r_squares <- results_S %>%
  group_by(scenario) %>%
  summarise(
    model = list(lm(I(log10(BF_gica) - log10(BF_jzs)) ~ 0)),  # 차이를 0에 회귀
    r_squared = 1 - sum(residuals(model[[1]])^2) / 
      sum((log10(BF_gica) - mean(log10(BF_gica)))^2)
  )
# 시나리오별 상관계수 계산
results_S %>%
  group_by(scenario) %>%
  summarise(
    correlation = cor(log10(BF_gica), log10(BF_jzs))^2  # R²를 바로 계산
  )

# 그래프에 R² 추가
results_S %>% 
  mutate(
    effect_magnitude = abs(std_mean_diff),
    log_BF_jzs = log10(BF_jzs),
    log_BF_gica = log10(BF_gica)
  ) %>%
  ggplot(aes(x = log_BF_jzs, y = log_BF_gica, color = effect_magnitude)) +
  geom_abline(intercept = 0, slope = 1, color = "red", 
              linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.3, size = 0.8) +
  # R² 값 텍스트 추가
  geom_text(data = r_squares,
            aes(x = -Inf, y = Inf, 
                label = sprintf("R² = %.3f", r_squared)),
            hjust = -0.2, vjust = 2,
            color = "black",
            inherit.aes = FALSE) +
  facet_wrap(~scenario, nrow = 2,
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  scale_color_viridis_c(limits = c(0, 1)) +
  labs(
    x = expression(log[10](BF[JZS])),
    y = expression(log[10](BF[GICA])),
    color = "|SMD|"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 11, face = "bold")
  )
# 2. Welch
r_squares <- results_W %>%
  group_by(scenario) %>%
  summarise(
    model = list(lm(I(log10(BF_gica) - log10(BF_jzs)) ~ 0)),  # 차이를 0에 회귀
    r_squared = 1 - sum(residuals(model[[1]])^2) / 
      sum((log10(BF_gica) - mean(log10(BF_gica)))^2)
  )
# 시나리오별 상관계수 계산
results_W %>%
  group_by(scenario) %>%
  summarise(
    correlation = cor(log10(BF_gica), log10(BF_jzs))^2  # R²를 바로 계산
  )

# 그래프에 R² 추가
results_W %>% 
  mutate(
    effect_magnitude = abs(std_mean_diff),
    log_BF_jzs = log10(BF_jzs),
    log_BF_gica = log10(BF_gica)
  ) %>%
  ggplot(aes(x = log_BF_jzs, y = log_BF_gica, color = effect_magnitude)) +
  geom_abline(intercept = 0, slope = 1, color = "red", 
              linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.3, size = 0.8) +
  # R² 값 텍스트 추가
  geom_text(data = r_squares,
            aes(x = -Inf, y = Inf, 
                label = sprintf("R² = %.3f", r_squared)),
            hjust = -0.2, vjust = 2,
            color = "black",
            inherit.aes = FALSE) +
  facet_wrap(~scenario, nrow = 2,
             labeller = labeller(scenario = setNames(scenario_labels, 1:5))) +
  scale_color_viridis_c(limits = c(0, 1)) +
  labs(
    x = expression(log[10](BF[JZS])),
    y = expression(log[10](BF[GICA])),
    color = "|SMD|"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.text = element_text(size = 11, face = "bold")
  )

#### 백분율 계산 ####
#BF10 직접 비교
results_S %>% #Student
  calculate_proportions_bf() %>%
  arrange(scenario, effect_size_cat)
results_W %>% #Welch
  calculate_proportions_bf() %>%
  arrange(scenario, effect_size_cat)
