# 시나리오 정보 설정
scenarios <- tibble(
  scenario = 1:5,
  n1 = c(1, 2, 1, 2, 3),
  n2 = c(1, 3, 1, 3, 2),
  var1 = c(1, 1, 2, 2, 2),
  var2 = c(1, 1, 1, 1, 1)
) %>%
  mutate(
    varR = var1/var2,
    # 좀 더 간결한 라벨 만들기 - 1/1 같은 경우는 그냥 1로 표시
    label = ifelse(
      n1 == n2,
      sprintf("v₁/v₂=%d, n₁/n₂=1", varR),
      sprintf("v₁/v₂=%d, n₁/n₂=%d/%d", varR, n1, n2)
    )
  ) %>%
  mutate(scenario = as.factor(scenario))

# 더 간결한 버전을 원하는 경우:
scenarios <- scenarios %>%
  mutate(
    minimal_label = ifelse(
      n1 == n2,
      sprintf("%d, 1", varR),
      sprintf("%d, %d/%d", varR, n1, n2)
    )
  )

# plot_data_combined 데이터에 scenarios 정보 추가
plot_data_combined <- plot_data_combined %>%
  left_join(scenarios, by = "scenario")

# 그래프 코드 수정
save3 <- ggplot(plot_data_combined, aes(x = factor(scenario), y = emmean, color = effect_size, 
                                        group = effect_size, linetype = effect_size, shape = effect_size)) +
  geom_line() +
  geom_point(size = 2.3) +
  facet_grid(sample_size ~ method) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Expected~log[10]~BF[10]),
       x = "") +
  # 간결한 라벨 사용
  scale_x_discrete(labels = function(x) {
    sapply(as.numeric(x), function(i) scenarios$label[scenarios$scenario == i])
  }) +
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
        axis.text.x = element_text(family = "Times New Roman", angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 11, color = "black"))

# 더 간결한 버전을 원한다면 이 줄을 대신 사용:
# scale_x_discrete(labels = function(x) {
#   sapply(as.numeric(x), function(i) scenarios$minimal_label[scenarios$scenario == i])
# }) +

print(save3)
ggsave("Analysis_colorK3.png", save3, dpi = 1000, units = "mm", width = 185, height = 222)

# 시나리오 정보 설정
scenarios <- tibble(
  scenario = 1:5,
  n1 = c(1, 2, 1, 2, 3),
  n2 = c(1, 3, 1, 3, 2),
  var1 = c(1, 1, 2, 2, 2),
  var2 = c(1, 1, 1, 1, 1)
) %>%
  mutate(
    varR = var2/var1,  # v1/v2에서 v2/v1로 변경
    # 좀 더 간결한 라벨 만들기 - 1/1 같은 경우는 그냥 1로 표시
    label = ifelse(
      n1 == n2,
      sprintf("v₂/v₁=%d, n₂/n₁=1", varR),  # v1/v2에서 v2/v1로 변경
      sprintf("v₂/v₁=%d, n₂/n₁=%d/%d", varR, n2, n1)  # v1/v2, n1/n2에서 v2/v1, n2/n1로 변경
    )
  ) %>%
  mutate(scenario = as.factor(scenario))

# 더 간결한 버전을 원하는 경우:
scenarios <- scenarios %>%
  mutate(
    minimal_label = ifelse(
      n1 == n2,
      sprintf("%d, 1", varR),
      sprintf("%d, %d/%d", varR, n2, n1)  # n1/n2에서 n2/n1로 변경
    )
  )

# plot_data_combined 데이터에 scenarios 정보 추가
plot_data_combined <- plot_data_combined %>%
  left_join(scenarios, by = "scenario")

# 그래프 코드 수정
save3 <- ggplot(plot_data_combined, aes(x = factor(scenario), y = emmean, color = effect_size, 
                                        group = effect_size, linetype = effect_size, shape = effect_size)) +
  geom_line() +
  geom_point(size = 2.3) +
  facet_grid(sample_size ~ method) +
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Expected~log[10]~BF[10]),
       x = "") +
  # 간결한 라벨 사용
  scale_x_discrete(labels = function(x) {
    sapply(as.numeric(x), function(i) scenarios$label[scenarios$scenario == i])
  }) +
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
        axis.text.x = element_text(family = "Times New Roman", angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 11, color = "black"))

# 더 간결한 버전을 원한다면 이 줄을 대신 사용:
# scale_x_discrete(labels = function(x) {
#   sapply(as.numeric(x), function(i) scenarios$minimal_label[scenarios$scenario == i])
# }) +

print(save3)
ggsave("Analysis_colorK3.png", save3, dpi = 1000, units = "mm", width = 185, height = 222)# 시나리오 정보 설정


scenarios <-tibble(
  scenario = 1:5,
  n1 = c(1, 2, 1, 2, 3),
  n2 = c(1, 3, 1, 3, 2),
  var1 = c(1, 1, 2, 2, 2),
  var2 = c(1, 1, 1, 1, 1)
) %>%
  mutate(
    label = case_when(
      var1 == var2 & n1 == n2 ~ sprintf("v₂/v₁=1, n₂/n₁=1"),
      var1 == var2 ~ sprintf("v₂/v₁=1, n₂/n₁=%d/%d", n2, n1),
      n1 == n2 ~ sprintf("v₂/v₁=%d/%d, n₂/n₁=1", var2, var1),
      TRUE ~ sprintf("v₂/v₁=%d/%d, n₂/n₁=%d/%d", var2, var1, n2, n1)
    )
  )%>%
  mutate(scenario = as.factor(scenario))


# plot_data_combined 데이터에 scenarios 정보 추가
plot_data_combined <- plot_data_combined %>%
  left_join(scenarios, by = "scenario")

# 그래프 코드 수정
save3 <- ggplot(plot_data_combined, aes(x = factor(scenario), y = emmean, color = effect_size, 
                                        group = effect_size, linetype = effect_size, shape = effect_size)) +
  geom_line() +
  geom_point(size = 2.3) +
  facet_grid(sample_size ~ method) + #free scale 추가: scales = "free_y"
  labs(title = "JZS vs. BFGC: Performance Comparison Across Scenarios",
       y = expression(Expected~log[10]~BF[10]),
       x = "") +
  # 간결한 라벨 사용
  scale_x_discrete(labels = function(x) {
    sapply(as.numeric(x), function(i) scenarios$label[scenarios$scenario == i])
  }) +
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
        axis.text.x = element_text(family = "Times New Roman", angle = 50, hjust = 1, size = 10),
        axis.text.y = element_text(family = "Times New Roman"),
        legend.text = element_text(family = "Times New Roman"),
        plot.title = element_text(family = "Times New Roman", face = "bold", size = 11, color = "black"))

print(save3)
ggsave("ANOVA1.png", save3, dpi = 1000, units = "mm", width = 140, height = 180)
