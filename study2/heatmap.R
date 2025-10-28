library(ggplot2)
library(gridExtra)
library(grid)

# 기본 alpha 설정
alpha <- 0.05

# Welch's t-test Power function
welch_power <- function(d, var1, var2, n1, n2, alpha = 0.05) {
  kappa <- n2/n1
  N <- (n1 + n2)
  rho <- var2/var1
  
  # Welch-Satterthwaite degrees of freedom
  df <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))
  
  # weight in NCP
  w <- sqrt((N/2) * ((kappa*(1+rho)) / ((1+kappa)*(kappa+rho))))
  
  # NCP
  ncp <- d * w
  
  # Power calculation (two-tailed test)
  power <- 1 - pt(qt(1-alpha/2, df=df), df=df, ncp=ncp) + 
    pt(qt(alpha/2, df=df), df=df, ncp=ncp)
  
  return(power)
}

# 효과크기(d) 범위 설정
d_values <- seq(0, 1, by = 0.01)

# 고정 조건 설정
N <- 120  # 총 표본크기
var1 <- 4  # 첫 번째 그룹의 분산
var2 <- 2  # 두 번째 그룹의 분산 (rho = 0.5)
rho <- var2/var1

# kappa 값과 라벨을 직접 정의
kappa_config <- data.frame(
  kappa_label = c("1/4", "1/3", "1/2", "2/3", "1", "3/2", "2", "3", "4"),
  n1 = c(96, 90, 80, 72, 60, 48, 40, 30, 24),
  n2 = c(24, 30, 40, 48, 60, 72, 80, 90, 96),
  stringsAsFactors = FALSE
)

# Power 계산
power_list <- list()

for(i in 1:nrow(kappa_config)) {
  n1 <- kappa_config$n1[i]
  n2 <- kappa_config$n2[i]
  kappa_label <- kappa_config$kappa_label[i]
  
  # Welch's t-test power 계산
  power_values <- numeric(length(d_values))
  for(j in 1:length(d_values)) {
    power_values[j] <- welch_power(d_values[j], var1, var2, n1, n2, alpha)
  }
  
  # 데이터프레임 생성
  temp_df <- data.frame(
    d = d_values,
    power = power_values,
    kappa_label = kappa_label,
    kappa = n2/n1,
    n1 = n1,
    n2 = n2
  )
  
  power_list[[i]] <- temp_df
}

# 모든 데이터 합치기
power_df <- do.call(rbind, power_list)

# kappa_label을 factor로 변환하여 순서 지정
power_df$kappa_label <- factor(power_df$kappa_label, 
                               levels = c("1/4", "1/3", "1/2", "2/3", "1", "3/2", "2", "3", "4"))

# 더 구분되는 색상과 선 스타일 조합
colors <- c(
  "#1f77b4",  # 1/4 (파랑)
  "#ff7f0e",  # 1/3 (주황)
  "#2ca02c",  # 1/2 (초록)
  "#d62728",  # 2/3 (빨강)
  "#000000",  # 1 (검정, 굵게)
  "#9467bd",  # 3/2 (보라)
  "#8c564b",  # 2 (갈색)
  "#e377c2",  # 3 (분홍)
  "#7f7f7f"   # 4 (회색)
)

# 선 타입 다양화
linetypes <- c("solid", "dashed", "dotted", "dotdash", "solid", 
               "dashed", "dotted", "dotdash", "longdash")

# 선 굵기
linewidths <- c(0.8, 0.8, 0.8, 0.8, 1.5, 0.8, 0.8, 0.8, 0.8)

# 플롯 생성
p <- ggplot(power_df, aes(x = d, y = power, color = kappa_label, 
                          linetype = kappa_label, group = kappa_label)) +
  geom_line(aes(size = kappa_label)) +
  scale_color_manual(values = colors,
                     name = bquote(kappa == n[2]/n[1])) +
  scale_linetype_manual(values = linetypes,
                        name = bquote(kappa == n[2]/n[1])) +
  scale_size_manual(values = linewidths,
                    guide = "none") +
  labs(title = bquote("Welch's t-test Power Functions: " ~ 
                        N == .(N) ~ ", " ~ rho == .(round(rho, 2))),
       x = "Effect Size (Cohen's d)",
       y = "Power") +
  scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0.5, 0.8)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50", alpha = 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm"),
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90")
  )

print(p)

# 특정 effect size에서 각 kappa의 power 확인
test_d <- 0.5

cat("\nPower at d =", test_d, ":\n")
for(i in 1:nrow(kappa_config)) {
  n1 <- kappa_config$n1[i]
  n2 <- kappa_config$n2[i]
  kappa_label <- kappa_config$kappa_label[i]
  
  power_val <- welch_power(test_d, var1, var2, n1, n2, alpha)
  
  cat(sprintf("κ = %s (n1=%d, n2=%d): Power = %.4f\n", 
              kappa_label, n1, n2, power_val))
}

#### kappa, rho function ####
library(ggplot2)
library(viridis)

# 고정 조건
N <- 120
alpha <- 0.05
d <- 0.5  # 고정된 effect size

# κ와 ρ의 범위
kappa_range <- seq(0.2, 5, length.out = 100)
rho_range <- seq(0.2, 5, length.out = 100)

# Power 계산을 위한 grid 생성
power_grid <- expand.grid(kappa = kappa_range, rho = rho_range)

# 각 조합에 대해 power 계산
power_grid$power <- sapply(1:nrow(power_grid), function(i) {
  kappa <- power_grid$kappa[i]
  rho <- power_grid$rho[i]
  
  # n1, n2 계산
  n1 <- round(N / (1 + kappa))
  n2 <- N - n1
  
  # var1을 1로 고정하고 var2 = rho * var1
  var1 <- 1
  var2 <- rho * var1
  
  # power 계산
  if(n1 < 2 || n2 < 2) return(NA)
  
  power <- welch_power(d, var1, var2, n1, n2, alpha)
  return(power)
})

# Contour plot
p1 <- ggplot(power_grid, aes(x = kappa, y = rho, z = power)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "white", alpha = 0.3, bins = 20) +
  scale_fill_viridis_d(option = "plasma", name = "Power") +
  labs(title = bquote("Welch's t-test Power: N=" ~ .(N) ~ ", d=" ~ .(d)),
       x = bquote(kappa == n[2]/n[1]),
       y = bquote(rho == sigma[2]^2/sigma[1]^2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "white", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "white", alpha = 0.5) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

print(p1)

# special condition points
special_points <- data.frame(
  kappa = c(1, 1/2, 1/3, 2/3),
  rho = c(1, 1/2, 1/2, 1/2),
  label = c("Reference", "Condition 1", "Condition 2", "Condition 3"),
  color = c("red", "black", "black", "black")
)

# Contour plot
p2 <- ggplot(power_grid, aes(x = kappa, y = rho, z = power)) +
  geom_contour_filled(bins = 20) +
  geom_contour(color = "white", alpha = 0.3, bins = 20) +
  scale_fill_viridis_d(option = "plasma", name = "Power") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "white", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "white", alpha = 0.5) +
  # add special points
  geom_point(data = special_points, 
             aes(x = kappa, y = rho, color = color),
             size = 1.5, shape = 16, inherit.aes = FALSE) +
  geom_text(data = special_points,
            aes(x = kappa, y = rho, label = paste0(round(kappa, 1))),
            vjust = 2, size = 2.3, color = "black",
            inherit.aes = FALSE) +
  scale_color_identity() +  # color를 직접 지정
  labs(title = bquote("Welch's t-test Power: N=" ~ .(N) ~ ", d=" ~ .(d)),
       x = bquote(kappa == n[2]/n[1]),
       y = bquote(rho == sigma[2]^2/sigma[1]^2)) +
  guides(fill = guide_legend(reverse = TRUE)) +  # legend
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )
print(p2)

ggsave("./study2/heatmap.png", p2, width = 8,height = 6,dpi=600)
