library(effectsize)  # effectsize 패키지 로드
library(afex);library(emmeans)
library(tidyverse)

## upload file
data <- read.csv("./New/analysisData.csv")
# 0. method 별 분리
data %>% filter(method == "JZS") -> data_jzs
data %>% filter(method =="BFGC") -> data_bfgc
data_jzs %>% mutate(scenario = factor(scenario, levels = 1:5), 
                    effect_size = factor(effect_size, levels = c("0", "0.2", "0.5", "0.8")),
                    sample_size = factor(sample_size, levels = c("50", "100","200"))) -> data_jzs
data_bfgc %>% mutate(scenario = factor(scenario, levels = 1:5), 
                     effect_size = factor(effect_size, levels = c("0", "0.2", "0.5", "0.8")),
                     sample_size = factor(sample_size, levels = c("50", "100","200"))) -> data_bfgc
summary(data_jzs);summary(data_bfgc)
# 직교비교 계수 설정
my_contrasts <- list(
  "Comp1" = c(3, 3, -2, -2, -2),
  "Comp2" = c(1, -1, 0, 0, 0),
  "Comp3" = c(0, 0, 2, -1, -1),
  "Comp4" = c(0, 0, 0, 1, -1)
)

#### JZS ####
# 1. ANOVA 분석
model <- aov(log_BF ~ scenario * effect_size * sample_size, data = data_jzs)
summary(model)

# 부분 에타 제곱 효과크기 계산
eta_squared(model, partial = TRUE)

# 2. 표본크기 수준별 단순 2원 상호작용 분석
simple <- emmeans(model, ~ scenario * effect_size | sample_size)
joint_tests(simple, by = "sample_size")

# 3. 단순 2원 상호작용이 유의한 sample_size 수준에서 단순단순 주효과 분석
  # scenario의 단순단순 주효과(effect_size별, sample_size별)
simple_simple <- emmeans(model, specs = ~ scenario | effect_size * sample_size)
joint_tests(simple_simple, by = c("effect_size", "sample_size"))

 # effect_size의 단순단순 주효과 (scenario별, sample_size별)
simple_simpleE <- emmeans(model, specs = ~ effect_size | scenario * sample_size)
joint_tests(simple_simpleE, by = c("scenario", "sample_size"))

# sample size의 단순단순 주효과 (scenario별, effect_size별)
simple_simple <- emmeans(model, specs = ~ sample_size | scenario * effect_size)
joint_tests(simple_simple, by = c("scenario", "effect_size"))


# 4. 단순단순 주효과의 pair wise 
  # scenarios (effect_size별, sample_size별)
ss_scen_pair <- emmeans(model, ~ scenario | effect_size * sample_size) %>%
  pairs()
print(ss_scen_pair)

  # effect_size (scenario별, sample_size별)
ss_effect_pair <- emmeans(model, ~ effect_size | scenario * sample_size) %>%
  pairs()
print(ss_effect_pair)

# sample_size (scenario별, effect_size별)
ss_sample_pair <- emmeans(model, ~ sample_size | scenario * effect_size) %>%
  pairs()
print(ss_sample_pair)


# 5. Contrast 분석 - 직교비교 적용
# 전체 contrast 분석
overall_contrast <- emmeans(model, ~ scenario) %>%
  contrast(method = my_contrasts)
print(overall_contrast)

# 단순단순 주효과의 맥락에서 contrast 분석
# 각 effect_size와 sample_size 조합별 contrast 분석
ss_contrast <- emmeans(model, ~ scenario | effect_size * sample_size) %>%
  contrast(method = my_contrasts)
print(ss_contrast)

# 결과를 표로 정리하여 시각화
# 상호작용 패턴 시각화
library(ggplot2)

# 평균값으로 상호작용 플롯 생성
interaction_means <- emmeans(model, ~ scenario * effect_size * sample_size)
plot_data <- as.data.frame(interaction_means)

# 상호작용 시각화
ggplot(plot_data, aes(x = scenario, y = emmean, color = effect_size, group = effect_size)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_size) +
  labs(title = "Interaction Plot of JZS: Scenario x Effect Size x Sample Size",
       y = "Mean log_BF") +
  theme_minimal()

#### BFGC ####
# 1. ANOVA 분석
model2 <- aov(log_BF ~ scenario * effect_size * sample_size, data = data_bfgc)
summary(model)

# 부분 에타 제곱 효과크기 계산
eta_squared(model2, partial = TRUE)

# 2. 표본크기 수준별 단순 2원 상호작용 분석
simple2 <- emmeans(model2, ~ scenario * effect_size | sample_size)
joint_tests(simple2, by = "sample_size")

# 3. 단순 2원 상호작용이 유의한 sample_size 수준에서 단순단순 주효과 분석
# scenario의 단순단순 주효과(effect_size별, sample_size별)
simple_simple2 <- emmeans(model2, specs = ~ scenario | effect_size * sample_size)
joint_tests(simple_simple2, by = c("effect_size", "sample_size"))

# effect_size의 단순단순 주효과 (scenario별, sample_size별)
simple_simpleE2 <- emmeans(model2, specs = ~ effect_size | scenario * sample_size)
joint_tests(simple_simpleE2, by = c("scenario", "sample_size"))

# sample size의 단순단순 주효과 (scenario별, effect_size별)
simple_simple2 <- emmeans(model2, specs = ~ sample_size | scenario * effect_size)
joint_tests(simple_simple2, by = c("scenario", "effect_size"))


# 4. 단순단순 주효과의 pair wise 
# scenarios (effect_size별, sample_size별)
ss_scen_pair2 <- emmeans(model2, ~ scenario | effect_size * sample_size) %>%
  pairs()
print(ss_scen_pair2)

# effect_size (scenario별, sample_size별)
ss_effect_pair2 <- emmeans(model2, ~ effect_size | scenario * sample_size) %>%
  pairs()
print(ss_effect_pair2)

# sample_size (scenario별, effect_size별)
ss_sample_pair2 <- emmeans(model2, ~ sample_size | scenario * effect_size) %>%
  pairs()
print(ss_sample_pair2)


# 5. Contrast 분석 - 직교비교 적용

# 전체 contrast 분석
overall_contrast2 <- emmeans(model, ~ scenario) %>%
  contrast(method = my_contrasts)
print(overall_contrast2)

# 단순단순 주효과의 맥락에서 contrast 분석
# 각 effect_size와 sample_size 조합별 contrast 분석
ss_contrast2 <- emmeans(model, ~ scenario | effect_size * sample_size) %>%
  contrast(method = my_contrasts)
print(ss_contrast2)

# 결과를 표로 정리하여 시각화
# 상호작용 패턴 시각화
# 평균값으로 상호작용 플롯 생성
interaction_means2 <- emmeans(model2, ~ scenario * effect_size * sample_size)
plot_data2 <- as.data.frame(interaction_means2)

# 상호작용 시각화
ggplot(plot_data2, aes(x = scenario, y = emmean, color = effect_size, group = effect_size)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_size) +
  labs(title = "Interaction Plot of BFGC: Scenario x Effect Size x Sample Size",
       y = "Mean log_BF") +
  theme_minimal()


#### within factor 도전 ####
data$id <- ceiling(seq_len(nrow(data)) / 2)
# within-subject 디자인 분석
model_within <- aov_car(
  log_BF ~ scenario * effect_size * sample_size * method + 
    Error(id/method),
  data = data
)

# 결과 출력
summary(model_within)

# 부분 에타 제곱 효과크기 계산
eta_squared(model_within, partial = TRUE)

# 2. 효과크기 수준별 단순 3원 상호작용 분석
simpleT <- emmeans(model_within, ~ scenario * sample_size * method | effect_size)
joint_tests(simpleT, by = "effect_size")

# 3. 단순 3원 상호작용이 유의한 effect_size 수준에서 단순단순 2원 상호작용 분석
# scenario의 단순단순 주효과(method별, effect_size별, sample_size별)
simple_simpleScT <- emmeans(model_within, specs = ~ scenario * method | effect_size * sample_size)
joint_tests(simple_simpleScT, by = c("effect_size", "sample_size"))

  # effect_size의 단순단순 주효과 (method별, scenario별, sample_size별)
  simple_simpleET <- emmeans(model_within, specs = ~ effect_size * method | scenario * sample_size)
  joint_tests(simple_simpleET, by = c("scenario", "sample_size"))

  # method의 단순단순 주효과 (sample size별, scenario별, effect_size별)
  simple_simpleSAT <- emmeans(model_within, specs = ~ sample_size * method | scenario * effect_size)
  joint_tests(simple_simpleSAT, by = c("scenario", "effect_size"))

# 4. 단순단순 2원 상호작용 분석이 유의한 effect size 수준에서 단순단순단순 주효과 분석
  # 시나리오 4와 5에서 방법 간 차이
  method_effect <- emmeans(model_within, 
                           ~ method | scenario * effect_size * sample_size,
                           subset = c(effect_size = c("0.2", "0.5", "0.8"), 
                                      sample_size = c("50", "100", "200")))
  pairs(method_effect)
# 시나리오의 단순단순단순 주효과
scenario_effect <- emmeans(model_within, 
                          ~ scenario | method * effect_size * sample_size,
                          subset = c(effect_size = c("0.2", "0.5", "0.8"), 
                                     sample_size = c("50", "100", "200")))
joint_tests(scenario_effect, by = c("method", "effect_size", "sample_size"))
summary(scenario_effect)
# contrast
scenario_effect2 <- emmeans(model_within, 
                         ~ scenario | method * effect_size * sample_size,
                         subset = c(effect_size = c("0.2", "0.5", "0.8"), 
                                    sample_size = c("50", "100", "200"))) %>%
  contrast(method = my_contrasts)
print(scenario_effect2)

# 결합된 데이터셋으로 시각화
interaction_means_combined <- emmeans(model_within, ~ scenario * effect_size * method | sample_size)
plot_data_combined <- as.data.frame(interaction_means_combined) %>% 
  mutate(scenario = factor(scenario, levels = 1:5,
                                     labels = c("A","B","C","D","E")))
# 방법론과 표본 크기별 그리드 시각화
ggplot(plot_data_combined, aes(x = scenario, y = emmean, color = method, group = method)) +
  geom_line() +
  geom_point() +
  facet_grid(effect_size ~ sample_size) +
  labs(title = "Comparison of Methods: Scenario x Effect Size x Sample Size",
       y =  expression(Mean~log[10]~BF[10])) +  # 수학 표기법 사용
  theme_minimal() +
  # 패널 간격 확대
  theme(panel.spacing = unit(0, "cm"),
        # 테두리 추가
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        # 그리드 제목 명확하게
        strip.text = element_text(face = "bold", size = 12),
        # 배경색 추가
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"))

#### 저장용 그림 그리기 ####
# 1. 방법론과 표본 크기별 그리드 시각화
save2 <- ggplot(plot_data_combined, aes(x = scenario, y = emmean, color = method, group = method, 
                                       linetype = method, shape = method)) +
  geom_line() +
  geom_point(size = 2.3) +
  facet_grid(effect_size ~ sample_size) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
                # expression("Comparison of Methods: Scenario " %*% " Effect Size " %*% " Sample Size"),
       y = expression(Expected~log[10]~BF[10]),
       x = expression("분산 비율" %*% "표본크기 비율 조합"))+
  scale_color_viridis_d() +  # 색상 차이가 큰 팔레트
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  theme_minimal() +
  # 패널 간격 확대
  theme(panel.spacing = unit(0, "cm"),
        # 테두리 추가
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        # 그리드 제목 명확하게
        strip.text = element_text(face = "bold", size = 11, family = "Times New Roman"),
        # 배경색 추가
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        # 글꼴 설정
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman"),
        axis.text = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 14, color = "black"))
print(save2)
ggsave("Analysis_colorK.png", save2, dpi = 1000, units = "mm", width = 185, height = 222)

# 2. y scale free
save_free <- ggplot(plot_data_combined, aes(x = scenario, y = emmean, color = method, group = method, 
                                            shape = method)) +
  geom_line() +
  geom_hline(data = hline_data, 
             aes(yintercept = y_value), 
             color = "red4",
             linetype = "dotted", 
             size = 0.6) +
  geom_point(size = 2.3) +
  facet_grid(effect_size ~ sample_size, scales = "free_y") +  # free_y 추가
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Expected~log[10]~BF[10]),
       x = expression("variance ratio" %*% "sample size ratio")) +
  scale_color_viridis_d() +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  theme_minimal() +
  theme(panel.spacing = unit(0, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(face = "bold", size = 10, family = "Times New Roman"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman"),
        axis.text = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 14, color = "black"))
print(save_free)
ggsave("Analysis_picture_free_scale.png", save_free, dpi = 1000, units = "mm", width = 185, height = 222)

# 2. 3*4 grid
# B 조건에 대한 수평선 데이터 (각 패널별로)
b_segments5 <- data.frame(
  sample_size = c(50, 100, 200),
  effect_size = rep(0.5, 3),  # 효과 크기 0.5에 대해서만 
  x_start = rep("A", 3),      # 시작 위치 (각 패널의 첫 번째 지점)
  x_end = rep("E", 3),        # 끝 위치 (각 패널의 마지막 지점)
  y_value = c(0.8815355, 0.60070821, 1.64158063)  # B 조건의 JZS 값
)
b_segments <- data.frame(
  sample_size = c(50, 100, 200),
  effect_size = c(0.8, 0.8, 0.8),  # 효과 크기 0.8에 대해서만 
  x_start = rep("A", 3),      # 시작 위치 (각 패널의 첫 번째 지점)
  x_end = rep("E", 3),        # 끝 위치 (각 패널의 마지막 지점)
  y_value = c(0.98, 2.39, 5.34)  # B 조건의 JZS 값
)
plot_data_combined$effect_size <- factor(plot_data_combined$effect_size, 
                                         levels = c(0.8, 0.5, 0.2, 0))

save_free2 <- ggplot(plot_data_combined, aes(x = scenario, y = emmean, color = method, group = method, 
                                            shape = method, linetype = method)) +
  geom_line() +
  geom_hline(data = b_segments5, 
             aes(yintercept = y_value), 
             color = "red4",
             linetype = "dotted", 
             size = 0.6) +
  geom_point(size = 1.5) +
  facet_grid(effect_size ~ sample_size, scales = "free_y") +  # free_y 추가
  facetted_pos_scales(
    y = list(
      effect_size == 0  ~ scale_y_continuous(limits = c(-0.8, -0.4)), # 0.6
      effect_size == 0.2 ~ scale_y_continuous(limits = c(-0.6, -0.2)), # 0.3
      effect_size == 0.5 ~ scale_y_continuous(limits = c(0, 2.0)), # # 2.1
      effect_size == 0.8 ~ scale_y_continuous(limits = c(0, 6.6))
    )
  ) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Factor",
       y = expression(Expected~log[10]~BF[10]),
       x = expression("variance ratio" %*% "sample size ratio")) +
  scale_color_viridis_d(option = "turbo", begin = 0.31, end = .95, direction = -1) +
  scale_linetype_manual(values = c("JZS" = "solid", 
                                   "BFGC" = "dashed")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  theme_minimal() +
  theme(aspect.ratio = 1,  # 모든 패널에 동일한 비율 적용
        panel.spacing = unit(0, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(face = "bold", size = 9, family = "Times New Roman"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman", size = 9),
        axis.text = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 10, color = "black"))
print(save_free2)
ggsave("4by3ANVOAplot.png", save_free2, units = "mm", dpi = 1000, width = 144, height = 160)
 # 2. 3 grid
library(ggh4x)
# 수평선을 추가하기 위한 데이터프레임 생성
hline_data <- data.frame(
  sample_size = c(50, 100, 200),
  y_value = c(0.98130610, 2.38736375, 5.33593648)
)

save_1x3_uniform <- ggplot(plot_data_combined, aes(x = scenario, y = emmean, 
                                                   color = method,
                                                   linetype = method,
                                                   shape = factor(effect_size),
                                                   group = interaction(method, effect_size))) +
  # 선 굵기를 줄임
  geom_line(size = 0.3) +
  geom_point(size = 1.5) +
  # 수평선도 얇게
  geom_hline(data = hline_data, 
             aes(yintercept = y_value), 
             color = "red4",
             linetype = "dotted", 
             size = 0.6) +
  facet_wrap(~ sample_size, nrow = 1, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      sample_size == 50  ~ scale_y_continuous(limits = c(-0.6, 1.25)),
      sample_size == 100 ~ scale_y_continuous(limits = c(-0.6, 3.25)),
      sample_size == 200 ~ scale_y_continuous(limits = c(-0.6, 7.25))
    )
  ) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Expected~log[10]~BF[10]),
       x = expression("variance ratio" %*% "sample size ratio")) +
  scale_color_viridis_d(option = "turbo", begin = 0.31, end = .95, direction = -1) +
  scale_linetype_manual(values = c("JZS" = "solid", 
                                   "BFGC" = "dashed")) +
  # 효과크기 0과 0.2에 빈 원 사용
  scale_shape_manual(name = "effect size",
                     values = c("0" = 1,  # 빈 원
                                "0.2" = 2, # 빈 삼각형
                                "0.5" = 15, # 채워진 사각형
                                "0.8" = 18), # 채워진 다이아몬드
                     labels = c("0", "0.2", "0.5", "0.8")) +
  theme_minimal() +
  theme(panel.spacing = unit(0.5, "cm"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(face = "bold", size = 9, family = "Times New Roman"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        text = element_text(family = "Times New Roman"),
        axis.title = element_text(family = "Times New Roman", size = 9),
        axis.text = element_text(family = "Times New Roman"),
        legend.title = element_text(family = "Times New Roman", size = 9),
        legend.text = element_text(family = "Times New Roman", size = 8),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.margin = margin(0, 5, 0, 5),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 11, color = "black")) +
  guides(color = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1))
print(save_1x3_uniform)
ggsave("1by3figure.png", plot = save_1x3_uniform, dpi = 1000, units = "mm", width = 144, height = 85)