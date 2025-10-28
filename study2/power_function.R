#### Power Function of Welch's t-test ####
var1 <- 4
var2 <- 2
n1 <- 40
n2 <- 60
d <- 0.5 # Cohen's d
kappa <- n2/n1 # sample size ratio
rho <- var2/var1 # variance ratio
N <- (n1 + n2) # total sample size


wdf <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1)) #SWS degrees of freedom

w <- sqrt((N/2) * ((kappa*(1+rho)) / ((1+kappa)*(kappa+rho)))) # weight in NCP

NCPw <- d*w # Welch's NCP


1-pt(qt(.975,df=88.2),df=88.2,ncp=2.5)+pt(qt(.025,df=88.2),df=88.2,ncp=2.5) 
1-pt(qt(.975,df=97.618),df=97.618,ncp=2.3717)+pt(qt(.025,df=97.618),df=97.618,ncp=2.3717) 
1-pt(qt(.975,df=64.589),df=64.589,ncp=2.534)+pt(qt(.025,df=64.589),df=64.589,ncp=2.534)
#condition D
1-pt(qt(.975,df=64.58947),df=64.58947,ncp=2.371708)+pt(qt(.025,df=64.58947),df=64.58947,ncp=2.371708) 
#with variables
1-pt(qt(.975,df=wdf),df=wdf,ncp=NCPw)+pt(qt(.025,df=wdf),df=wdf,ncp=NCPw) 


#### Power Function: 4-panel comparison ####
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# 기본 alpha 설정
alpha <- 0.05

# 통합 t-test Power function
t_power <- function(d, var1, var2, n1, n2, alpha = 0.05, var_equal = FALSE) {
  kappa <- n2/n1
  N <- (n1 + n2)
  
  if (var_equal) {
    # Student's t-test (equal variance)
    rho <- 1
    df <- n1 + n2 - 2
  } else {
    # Welch's t-test (unequal variance)
    rho <- var2/var1
    # Welch-Satterthwaite degrees of freedom
    df <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))
  }
  
  # weight in NCP
  w <- sqrt((N/2) * ((kappa*(1+rho)) / ((1+kappa)*(kappa+rho))))
  
  # NCP
  ncp <- d * w
  
  # Power calculation (two-tailed test)
  power <- 1 - pt(qt(1-alpha/2, df=df), df=df, ncp=ncp) + 
    pt(qt(alpha/2, df=df), df=df, ncp=ncp)
  
  return(power)
}

# 각 조건에 대한 power 계산 함수
calculate_power_condition <- function(var1, var2, kappa, N_values, d_values, alpha) {
  power_list <- list()
  
  for(N in N_values) {
    # n1, n2 계산
    n1 <- round(N / (1 + kappa))
    n2 <- N - n1
    
    # Welch's t-test power 계산
    welch_power_values <- numeric(length(d_values))
    for(i in 1:length(d_values)) {
      welch_power_values[i] <- t_power(d_values[i], var1, var2, n1, n2, alpha, var_equal = FALSE)
    }
    
    # Student's t-test power 계산
    student_power_values <- numeric(length(d_values))
    for(i in 1:length(d_values)) {
      student_power_values[i] <- t_power(d_values[i], var1, var2, n1, n2, alpha, var_equal = TRUE)
    }
    
    # 데이터프레임 생성
    temp_df <- data.frame(
      d = rep(d_values, 2),
      power = c(welch_power_values, student_power_values),
      test = rep(c("Welch's t-test", "Student's t-test"), each = length(d_values)),
      N = N
    )
    
    power_list[[as.character(N)]] <- temp_df
  }
  
  # 모든 데이터 합치기
  power_df <- do.call(rbind, power_list)
  return(power_df)
}

# 패널 생성 함수
create_panel <- function(power_df, condition_name, kappa, rho, label_positions) {
  p <- ggplot(power_df, aes(x = d, y = power, color = test, linetype = test, group = interaction(test, N))) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Welch's t-test" = "#0066CC", 
                                  "Student's t-test" = "gray50")) +
    scale_linetype_manual(values = c("Welch's t-test" = "solid", 
                                     "Student's t-test" = "dashed")) +
    # Total sample size 라벨
    annotate("text", x = label_positions$N240[1], y = label_positions$N240[2], 
             label = "N=240", size = 3, hjust = 0) +
    annotate("text", x = label_positions$N120[1], y = label_positions$N120[2], 
             label = "N=120", size = 3, hjust = 0) +
    annotate("text", x = label_positions$N60[1], y = label_positions$N60[2], 
             label = "N=60", size = 3, hjust = 0) +
    labs(title = bquote(.(condition_name) ~ ": " ~ n[2]/n[1] == .(kappa) ~ ", " ~ rho == .(rho)),
         x = "Effect Size (Cohen's d)",
         y = "Power") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_line(color = "gray95"),
      panel.grid.major = element_line(color = "gray90")
    )
  
  return(p)
}

# 효과크기(d) 범위 설정
d_values <- seq(0, 1, by = 0.01)

# 총 표본크기 설정
N_values <- c(60, 120, 240)

# Condition B: n2/n1 = 2/3, rho = 1
power_df_B <- calculate_power_condition(var1 = 4, var2 = 4, kappa = 1/3, N_values, d_values, alpha)
panel_B <- create_panel(power_df_B, "B", "1/2", "1", 
                        list(N240 = c(0.4, 1.00), N120 = c(0.45, 0.9), N60 = c(0.65, 0.8)))

# Condition C: n2/n1 = 1, rho = 0.49
power_df_C <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 1, N_values, d_values, alpha)
panel_C <- create_panel(power_df_C, "C", "1", "0.49", 
                        list(N240 = c(0.4, 1.00), N120 = c(0.45, 0.9), N60 = c(0.6, 0.8)))

# Condition D: n2/n1 = 1/2, rho = 0.49
power_df_D <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 1/3, N_values, d_values, alpha)
panel_D <- create_panel(power_df_D, "D", "1/2", "0.49", 
                        list(N240 = c(0.4, 1.00), N120 = c(0.45, 0.9), N60 = c(0.6, 0.8)))

# Condition E: n2/n1 = 2, rho = 0.49
power_df_E <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 3, N_values, d_values, alpha)
panel_E <- create_panel(power_df_E, "E", "2", "0.49", 
                        list(N240 = c(0.4, 1.00), N120 = c(0.49, 0.9), N60 = c(0.65, 0.8)))

# 범례 생성 (별도로)
legend_df <- data.frame(
  d = c(0, 0.5),
  power = c(0.5, 0.6),
  test = c("Welch's t-test", "Student's t-test")
)

legend_plot <- ggplot(legend_df, aes(x = d, y = power, color = test, linetype = test)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Welch's t-test" = "#0066CC", 
                                "Student's t-test" = "gray50")) +
  scale_linetype_manual(values = c("Welch's t-test" = "solid", 
                                   "Student's t-test" = "dashed")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10)) +
  labs(color = "Test Type", linetype = "Test Type")

# 범례만 추출
legend <- get_legend(legend_plot)

# 4개 패널 배치 (2×2) - 수학 좌표계로
# 1행 (위): B(왼쪽), C(오른쪽)  <- 2사분면, 1사분면
# 2행 (아래): D(왼쪽), E(오른쪽)  <- 3사분면, 4사분면
combined_plot <- grid.arrange(
  panel_B, panel_C,
  panel_D, panel_E,
  ncol = 2,
  top = textGrob("Power Functions across Different Conditions", 
                 gp = gpar(fontsize = 16, fontface = "bold")),
  bottom = legend
)

# 출력
print(combined_plot)

#### Power Function: 4-panel comparison with Simulation Results (Publication version) ####
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# 기본 alpha 설정
alpha <- 0.05

# 시뮬레이션 데이터 (제공된 데이터 기반)
sim_data <- data.frame(
  Condition = rep(c("B", "C", "D", "E"), each = 24),
  N = rep(rep(c(60, 120, 240), each = 8), times = 4),
  d = rep(rep(c(0.0, 0.2, 0.5, 0.8), each = 2), times = 12),
  test = rep(c("Student's t-test", "Welch's t-test"), times = 48),
  power = c(
       # Condition B (n2/n1 = 1/3, rho = 1)
    # N=60
    0.051, 0.053, 0.101, 0.098, 0.377, 0.364, 0.757, 0.736,
    # N=120
    0.049, 0.049, 0.157, 0.155, 0.648, 0.638, 0.964, 0.96,
    # N=240
    0.052, 0.052, 0.264, 0.263, 0.915, 0.911, 1, 1,
    
    # Condition C (n2/n1 = 1, rho = 2)
    # N=60
    0.05, 0.05, 0.119, 0.118, 0.479, 0.477, 0.858, 0.857,
    # N=120
    0.049, 0.048, 0.191, 0.19, 0.773, 0.773, 0.991, 0.991,
    # N=240
    0.053, 0.053, 0.339, 0.338, 0.97, 0.97, 1, 1,
    
    # Condition D (n2/n1 = 1/3, rho = 2)
    # N=60
    0.021, 0.05, 0.056, 0.112, 0.305, 0.434, 0.711, 0.813,
    # N=120
    0.021, 0.052, 0.099, 0.179, 0.598, 0.722, 0.964, 0.983,
    # N=240
    0.02, 0.049, 0.196, 0.307, 0.911, 0.954, 1, 1,
    
    # Condition E (n2/n1 = 3, rho = 2)
    # N=60
    0.097, 0.052, 0.159, 0.091, 0.459, 0.322, 0.785, 0.655,
    # N=120
    0.096, 0.05, 0.224, 0.138, 0.701, 0.571, 0.967, 0.929,
    # N=240
    0.097, 0.048, 0.341, 0.232, 0.922, 0.863, 0.999, 0.999
  )
)

# 통합 t-test Power function
t_power <- function(d, var1, var2, n1, n2, alpha = 0.05, var_equal = FALSE) {
  kappa <- n2/n1
  N <- (n1 + n2)
  
  if (var_equal) {
    # Student's t-test (equal variance)
    rho <- 1
    df <- n1 + n2 - 2
  } else {
    # Welch's t-test (unequal variance)
    rho <- var2/var1
    # Welch-Satterthwaite degrees of freedom
    df <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))
  }
  
  # weight in NCP
  w <- sqrt((N/2) * ((kappa*(1+rho)) / ((1+kappa)*(kappa+rho))))
  
  # NCP
  ncp <- d * w
  
  # Power calculation (two-tailed test)
  power <- 1 - pt(qt(1-alpha/2, df=df), df=df, ncp=ncp) + 
    pt(qt(alpha/2, df=df), df=df, ncp=ncp)
  
  return(power)
}

# 각 조건에 대한 power 계산 함수
calculate_power_condition <- function(var1, var2, kappa, N_values, d_values, alpha) {
  power_list <- list()
  
  for(N in N_values) {
    # n1, n2 계산
    n1 <- round(N / (1 + kappa))
    n2 <- N - n1
    
    # Welch's t-test power 계산
    welch_power_values <- numeric(length(d_values))
    for(i in 1:length(d_values)) {
      welch_power_values[i] <- t_power(d_values[i], var1, var2, n1, n2, alpha, var_equal = FALSE)
    }
    
    # Student's t-test power 계산
    student_power_values <- numeric(length(d_values))
    for(i in 1:length(d_values)) {
      student_power_values[i] <- t_power(d_values[i], var1, var2, n1, n2, alpha, var_equal = TRUE)
    }
    
    # 데이터프레임 생성
    temp_df <- data.frame(
      d = rep(d_values, 2),
      power = c(welch_power_values, student_power_values),
      test = rep(c("Welch's t-test", "Student's t-test"), each = length(d_values)),
      N = N
    )
    
    power_list[[as.character(N)]] <- temp_df
  }
  
  # 모든 데이터 합치기
  power_df <- do.call(rbind, power_list)
  return(power_df)
}

# 패널 생성 함수
create_panel <- function(power_df, sim_df, condition_name, kappa, rho, label_positions) {
  # 제목: 조건명 + 파라미터 정보
  title_text <- bquote(bold(.(condition_name)) ~ ": " ~ n[2]/n[1] == .(kappa) ~ ", " ~ rho == .(rho))
  
  p <- ggplot(power_df, aes(x = d, y = power, color = test, linetype = test, group = interaction(test, N))) +
    geom_line(linewidth = 1.2) +
    # 시뮬레이션 결과 추가
    geom_point(data = sim_df, aes(x = d, y = power, shape = test, group = test), 
               color = "#E63946", size = 2.2, stroke = 1.2) +
    scale_color_manual(values = c("Welch's t-test" = "#0066CC", 
                                  "Student's t-test" = "gray50")) +
    scale_linetype_manual(values = c("Welch's t-test" = "solid", 
                                     "Student's t-test" = "longdash")) +
    scale_shape_manual(values = c("Welch's t-test" = 16, 
                                  "Student's t-test" = 1)) +
    # Total sample size 라벨 (크기 증가)
    annotate("text", x = label_positions$N240[1], y = label_positions$N240[2], 
             label = "N=240", size = 4, hjust = 0) +
    annotate("text", x = label_positions$N120[1], y = label_positions$N120[2], 
             label = "N=120", size = 4, hjust = 0) +
    annotate("text", x = label_positions$N60[1], y = label_positions$N60[2], 
             label = "N=60", size = 4, hjust = 0) +
    ggtitle(title_text) +
    labs(x = "Effect Size (Cohen's d)",
         y = "Power") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, size = 13, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_line(color = "gray95"),
      panel.grid.major = element_line(color = "gray90")
    )
  
  return(p)
}

# 효과크기(d) 범위 설정
d_values <- seq(0, 1, by = 0.01)

# 총 표본크기 설정
N_values <- c(60, 120, 240)

# Condition B (n2/n1 = 1/2, rho = 1)
power_df_B <- calculate_power_condition(var1 = 4, var2 = 4, kappa = 1/3, N_values, d_values, alpha)
sim_df_B <- sim_data[sim_data$Condition == "B", ]
panel_B <- create_panel(power_df_B, sim_df_B, "B", "1/3", "1", 
                        list(N240 = c(0.3, 0.9), N120 = c(0.45, 0.82), N60 = c(0.65, 0.75)))

# Condition C (n2/n1 = 1, rho = 0.49)
power_df_C <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 1, N_values, d_values, alpha)
sim_df_C <- sim_data[sim_data$Condition == "C", ]
panel_C <- create_panel(power_df_C, sim_df_C, "C", "1", "0.5", 
                        list(N240 = c(0.30, 1.00), N120 = c(0.44, 0.9), N60 = c(0.58, 0.8)))

# Condition D (n2/n1 = 1/2, rho = 0.49)
power_df_D <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 1/3, N_values, d_values, alpha)
sim_df_D <- sim_data[sim_data$Condition == "D", ]
panel_D <- create_panel(power_df_D, sim_df_D, "D", "1/3", "0.5", 
                        list(N240 = c(0.3, 0.9), N120 = c(0.44, 0.85), N60 = c(0.63, 0.8)))

# Condition E (n2/n1 = 2, rho = 0.49)
power_df_E <- calculate_power_condition(var1 = 4, var2 = 1.96, kappa = 3, N_values, d_values, alpha)
sim_df_E <- sim_data[sim_data$Condition == "E", ]
panel_E <- create_panel(power_df_E, sim_df_E, "E", "3", "0.5", 
                        list(N240 = c(0.35, 1.00), N120 = c(0.55, 0.9), N60 = c(0.62, 0.7)))

# 범례 생성 (별도로)
legend_df_line <- data.frame(
  d = c(0, 0.5, 0, 0.5),
  power = c(0.5, 0.6, 0.45, 0.55),
  test = rep(c("Welch's t-test", "Student's t-test"), each = 2)
)

legend_df_point <- data.frame(
  d = c(0.25, 0.25),
  power = c(0.55, 0.475),
  test = c("Welch's t-test", "Student's t-test")
)

legend_plot <- ggplot() +
  geom_line(data = legend_df_line, aes(x = d, y = power, color = test, linetype = test), linewidth = 1) +
  geom_point(data = legend_df_point, aes(x = d, y = power, shape = test), 
             color = "#E63946", size = 1.8, stroke = 1) +
  scale_color_manual(values = c("Welch's t-test" = "#0066CC", 
                                "Student's t-test" = "gray50"),
                     name = "Theoretical") +
  scale_linetype_manual(values = c("Welch's t-test" = "solid", 
                                   "Student's t-test" = "longdash"),
                        name = "Theoretical") +
  scale_shape_manual(values = c("Welch's t-test" = 16, 
                                "Student's t-test" = 1),
                     name = "Simulation") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9))

# 범례만 추출
legend <- get_legend(legend_plot)

# 4개 패널 배치 (2×2)
combined_plot <- grid.arrange(
  panel_B, panel_C,
  panel_D, panel_E,
  ncol = 2,
  bottom = legend
)

print(combined_plot)
# save
ggsave("power_function_4panelre.png", combined_plot, width = 9, height = 8, dpi = 600)
