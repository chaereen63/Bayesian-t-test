#### Comparison of Theoretical Distribution vs Sampling Distribution ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(moments)  # For skewness and kurtosis calculation

results <- readRDS("./study2/results2/S2befi200_E8.RDS")
# Load necessary functions
source("./study2/t-distribution.R")  # File containing theoretical distribution functions
source("./study2/Behrens-Fisher_distribution.R")

# Extract population parameters for each condition from simulation results
condition_params <- results %>%
  group_by(scenario) %>%
  summarise(
    n1 = first(n1),
    n2 = first(n2),
    s1 = first(sd1),  # Population standard deviation
    s2 = first(sd2),  # Population standard deviation
    var1 = first(var1),  # Population variance
    var2 = first(var2),  # Population variance
    pop_mean1 = first(pop_mean1),
    pop_mean2 = first(pop_mean2),
    pop_effect_size = first(pop_effect_size),  # Population effect size
    equal_var = abs(first(var1) - first(var2)) < 0.001,  # Check for equal variance
    .groups = 'drop'
  ) %>%
  split(.$scenario) %>%
  map(~as.list(.))

# Print population parameters
{cat("=== Population Parameters by Condition ===\n")
  for(i in 1:5) {
    params <- condition_params[[i]]
    cat(sprintf("Condition %d: n1=%d, n2=%d, σ1=%.1f, σ2=%.3f, pop_effect_size=%.2f, equal_var=%s\n", 
                i, params$n1, params$n2, params$s1, params$s2, params$pop_effect_size, params$equal_var))
  }
  cat("\n")
  }

#### Calculate Descriptive Statistics for Sampling Distributions ####
cat("=== Descriptive Statistics of Sampling Distributions (Mean, Skewness, Kurtosis) ===\n\n")

descriptive_stats <- results %>%
  group_by(scenario) %>%
  summarise(
    # Student's t
    student_mean = mean(student_t, na.rm = TRUE),
    student_skewness = skewness(student_t, na.rm = TRUE),
    student_kurtosis = kurtosis(student_t, na.rm = TRUE),
    
    # Welch's t
    welch_mean = mean(welch_t, na.rm = TRUE),
    welch_skewness = skewness(welch_t, na.rm = TRUE),
    welch_kurtosis = kurtosis(welch_t, na.rm = TRUE),
    
    .groups = 'drop'
  )

# Print results
for(i in 1:5) {
  cat(sprintf("Condition %d:\n", i))
  cat(sprintf("  Student's t - Mean: %.4f, Skewness: %.4f, Kurtosis: %.4f\n",
              descriptive_stats$student_mean[i],
              descriptive_stats$student_skewness[i],
              descriptive_stats$student_kurtosis[i]))
  cat(sprintf("  Welch's t   - Mean: %.4f, Skewness: %.4f, Kurtosis: %.4f\n\n",
              descriptive_stats$welch_mean[i],
              descriptive_stats$welch_skewness[i],
              descriptive_stats$welch_kurtosis[i]))
}

# Save descriptive statistics
write.csv(descriptive_stats, "./study2/results2/descriptive_stats.csv", row.names = FALSE)
cat("Descriptive statistics saved to './study2/results2/descriptive_stats.csv'\n\n")

#### Function to Calculate Theoretical Distribution ####
calculate_theoretical_distribution <- function(condition, x_range) {
  params <- condition_params[[condition]]
  n1 <- params$n1
  n2 <- params$n2  
  s1 <- params$s1
  s2 <- params$s2
  pop_mean1 <- params$pop_mean1
  pop_mean2 <- params$pop_mean2
  equal_var <- params$equal_var
  pop_effect_size <- params$pop_effect_size
  
  # Calculate mean difference
  mean_diff <- pop_mean1 - pop_mean2
  
  if (equal_var) {
    # Conditions 1-2: Equal variance condition (Student's t)
    df_pooled <- n1 + n2 - 2
    pooled_var <- ((n1-1)*s1^2 + (n2-1)*s2^2) / df_pooled
    se_pooled <- sqrt(pooled_var * (1/n1 + 1/n2))
    
    if (abs(pop_effect_size) < 0.001) {
      # Central t-distribution
      density_vals <- dt(x_range, df_pooled)
      dist_name <- paste0("Central t(", df_pooled, ")")
    } else {
      # Noncentral t-distribution
      ncp <- mean_diff / se_pooled
      density_vals <- dt(x_range, df_pooled, ncp=ncp)
      dist_name <- paste0("Noncentral t(", df_pooled, ", δ=", round(ncp,3), ")")
    }
  } else {
    # Conditions 3-5: Unequal variance condition - Welch method
    se_welch <- sqrt(s1^2/n1 + s2^2/n2)
    df_welch <- (s1^2/n1 + s2^2/n2)^2 / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
    
    if (abs(pop_effect_size) < 0.001) {
      # Central t-distribution (Welch's df)
      density_vals <- dt(x_range, df_welch)
      dist_name <- paste0("Central t(", round(df_welch,1), ") - Welch's df")
    } else {
      # Noncentral t-distribution (Welch's df)
      ncp_welch <- mean_diff / se_welch
      density_vals <- dt(x_range, df_welch, ncp=ncp_welch)
      dist_name <- paste0("Noncentral t(", round(df_welch,1), ", δ=", round(ncp_welch,3), ") - Welch's df")
    }
  }
  
  return(list(density = density_vals, name = dist_name))
}


#### Function to Create Comparison Plot for Each Condition ####
create_comparison_plot <- function(condition_num) {
  
  # Filter data for the condition
  condition_data <- results %>% 
    filter(scenario == condition_num) %>%
    select(student_t, welch_t) %>%
    pivot_longer(cols = everything(), names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t"),
                           labels = c("Student's t", "Welch's t")))
  
  # Get parameters for the condition
  params <- condition_params[[condition_num]]
  
  # Set x-axis range based on data range
  x_min <- min(condition_data$t_value, na.rm = TRUE)
  x_max <- max(condition_data$t_value, na.rm = TRUE)
  x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 300)
  
  # Calculate theoretical distribution
  theoretical <- calculate_theoretical_distribution(condition_num, x_range)
  
  # Create theoretical distribution data frame
  theoretical_df <- data.frame(
    x = x_range,
    density = theoretical$density,
    type = "Theoretical"
  )
  
  # Create base plot
  p <- ggplot() +
    # Simulation results (sampling distribution)
    geom_density(data = condition_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.6, size = 0.8) +
    # Theoretical distribution
    geom_line(data = theoretical_df, 
              aes(x = x, y = density), 
              color = "black", size = 1.5, linetype = "dashed") +
    
    facet_wrap(~method, scales = "free", ncol = 2) +
    
    labs(title = paste0("Condition ", condition_num, " - Theoretical vs Sampling Distribution"),
         subtitle = paste0("Parameters: n1=", params$n1, ", n2=", params$n2,
                           ", σ1=", params$s1, ", σ2=", round(params$s2,3),
                           ", d=", params$pop_effect_size),
         x = "t value", 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold", family = "Times New Roman"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, family = "Times New Roman"),
      strip.text = element_text(size = 10, face = "bold", family = "Times New Roman"),
      axis.title = element_text(family = "Times New Roman"),
      axis.text = element_text(family = "Times New Roman"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )+
    
    scale_fill_manual(values = c("#08519C", "#BDBDBD")) +
    scale_color_manual(values = c("#08519C", "#BDBDBD"))
  
  return(p)
}

# Generate plots for all conditions
for(i in 1:5) {
  p <- create_comparison_plot(i)
  print(p)
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
}

#### Integrated Comparison Plot ####
create_integrated_comparison <- function() {
  
  # Prepare data for all conditions
  all_data <- results %>%
    select(scenario, student_t, welch_t) %>%
    pivot_longer(cols = c(student_t, welch_t), 
                 names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t"),
                           labels = c("Student's t", "Welch's t")),
           scenario = paste("Condition", scenario))
  
  # Calculate theoretical distribution for each condition
  theoretical_data <- data.frame()
  
  for(s in 1:5) {
    condition_subset <- all_data %>% filter(scenario == paste("Condition", s))
    x_min <- min(condition_subset$t_value, na.rm = TRUE)
    x_max <- max(condition_subset$t_value, na.rm = TRUE)
    x_range <- seq(x_min - 0.5, x_max + 0.5, length.out = 200)
    
    theoretical <- calculate_theoretical_distribution(s, x_range)
    
    temp_df <- data.frame(
      x = x_range,
      density = theoretical$density,
      scenario = paste("Condition", s),
      dist_name = theoretical$name
    )
    
    theoretical_data <- rbind(theoretical_data, temp_df)
  }
  
  # Create integrated plot
  p <- ggplot() +
    # Sampling distribution
    geom_density(data = all_data, 
                 aes(x = t_value, fill = method, color = method), 
                 alpha = 0.5, size = 0.6) +
    
    # Theoretical distribution
    geom_line(data = theoretical_data, 
              aes(x = x, y = density), 
              color = "black", size = 1, linetype = "dashed") +
    
    facet_grid(scenario ~ method, scales = "free") +
    
    labs(title = "All Conditions: Theoretical vs Sampling Distribution",
         subtitle = "Black dashed line: Theoretical distribution, Colored area: Simulated sampling distribution",
         x = "t value", 
         y = "Density") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, family = "Times New Roman"),
      strip.text = element_text(size = 9, face = "bold", family = "Times New Roman"),
      axis.text = element_text(size = 8, family = "Times New Roman"),
      axis.title = element_text(family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )+
    
    scale_fill_manual(values = c("#08519C", "#BDBDBD")) +
    scale_color_manual(values = c("#08519C", "#BDBDBD"))
  
  return(p)
}

# Generate integrated comparison plot
integrated_plot <- create_integrated_comparison()
print(integrated_plot)

ggsave(
  filename = "./study2/results2/allc_200_8.png",
  plot = integrated_plot,
  width = 11,
  height = 8,
  dpi = 600,
  device = "png",
  bg = "white"
)

cat("\nAnalysis completed!\n")
cat("- Black dashed line: Theoretical distribution for each condition\n")
cat("- Colored area: Actual sampling distribution from simulation\n")
cat("- Conditions 1-2: t-distribution as theoretical reference\n")
cat("- Conditions 3-5: Welch's t-distribution as theoretical reference\n")

create_direct_overlay_plot <- function(condition_num) {
  
  condition_data <- results %>% 
    filter(scenario == condition_num) %>%
    select(student_t, welch_t) %>%
    pivot_longer(cols = everything(), names_to = "method", values_to = "t_value") %>%
    mutate(method = factor(method, 
                           levels = c("student_t", "welch_t"),
                           labels = c("Student's t", "Welch's t")))
  
  params <- condition_params[[condition_num]]
  
  # Calculate critical value
  se_welch <- sqrt(params$s1^2/params$n1 + params$s2^2/params$n2)
  df_welch <- (params$s1^2/params$n1 + params$s2^2/params$n2)^2 / 
    ((params$s1^2/params$n1)^2/(params$n1-1) + 
       (params$s2^2/params$n2)^2/(params$n2-1))
  critical_value <- qt(0.975, df_welch)
  
  # Calculate means
  means <- condition_data %>%
    group_by(method) %>%
    summarise(mean_t = mean(t_value, na.rm = TRUE), .groups = 'drop')
  
  p <- ggplot(condition_data, aes(x = t_value, fill = method, color = method)) +
    geom_density(alpha = 0.5, size = 1) +
    # Mean lines
    geom_vline(data = means, aes(xintercept = mean_t, color = method),
               linetype = "solid", size = 1.2) +
    # Critical value lines
    geom_vline(xintercept = c(-critical_value, critical_value), 
               linetype = "dashed", color = "red", size = 1) +
    
    labs(title = paste0("Condition ", condition_num, " - Direct Comparison of Distributions"),
         subtitle = paste0("Solid lines: Mean of each distribution, Dashed red: Critical value (α=.05)\n",
                           "Parameters: n1=", params$n1, ", n2=", params$n2,
                           ", σ1=", params$s1, ", σ2=", round(params$s2,3),
                           ", d=", params$pop_effect_size),
         x = "t value", 
         y = "Density",
         fill = "Method",
         color = "Method") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold", family = "Times New Roman"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, family = "Times New Roman"),
      axis.title = element_text(family = "Times New Roman"),
      axis.text = element_text(family = "Times New Roman"),
      legend.text = element_text(family = "Times New Roman"),
      legend.title = element_text(family = "Times New Roman"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    
    scale_fill_manual(values = c("#08519C", "#BDBDBD")) +
    scale_color_manual(values = c("#08519C", "#BDBDBD"))
  
  return(p)
}

# Generate for conditions 4, 5
for(i in 4:5) {
  p <- create_direct_overlay_plot(i)
  print(p)
  ggsave(
    filename = paste0("./study2/results2/condition_", i, "_overlay.png"),
    plot = p,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

create_welch_comparison_plot <- function() {
  
  # Extract Welch's t data for Conditions 4, 5
  welch_data <- results %>% 
    filter(scenario %in% c(4, 5)) %>%
    select(scenario, welch_t) %>%
    mutate(scenario = factor(scenario, 
                             levels = c(4, 5),
                             labels = c("Condition 4", "Condition 5")))
  
  # Get parameters for each condition
  params4 <- condition_params[[4]]
  params5 <- condition_params[[5]]
  
  # Calculate critical values for each condition
  se_welch4 <- sqrt(params4$s1^2/params4$n1 + params4$s2^2/params4$n2)
  df_welch4 <- (params4$s1^2/params4$n1 + params4$s2^2/params4$n2)^2 / 
    ((params4$s1^2/params4$n1)^2/(params4$n1-1) + 
       (params4$s2^2/params4$n2)^2/(params4$n2-1))
  critical_value4 <- qt(0.975, df_welch4)
  
  se_welch5 <- sqrt(params5$s1^2/params5$n1 + params5$s2^2/params5$n2)
  df_welch5 <- (params5$s1^2/params5$n1 + params5$s2^2/params5$n2)^2 / 
    ((params5$s1^2/params5$n1)^2/(params5$n1-1) + 
       (params5$s2^2/params5$n2)^2/(params5$n2-1))
  critical_value5 <- qt(0.975, df_welch5)
  
  # Calculate means
  means <- welch_data %>%
    group_by(scenario) %>%
    summarise(mean_t = mean(welch_t, na.rm = TRUE), .groups = 'drop')
  
  # Critical value data frame
  critical_df <- data.frame(
    scenario = factor(c("Condition 4", "Condition 5"), levels = c("Condition 4", "Condition 5")),
    critical_pos = c(critical_value4, critical_value5),
    critical_neg = c(-critical_value4, -critical_value5)
  )
  
  # Calculate theoretical distribution
  x_range <- seq(-3, 8, length.out = 300)
  
  theoretical4 <- calculate_theoretical_distribution(4, x_range)
  theoretical5 <- calculate_theoretical_distribution(5, x_range)
  
  theoretical_df <- rbind(
    data.frame(x = x_range, density = theoretical4$density, scenario = "Condition 4"),
    data.frame(x = x_range, density = theoretical5$density, scenario = "Condition 5")
  ) %>%
    mutate(scenario = factor(scenario, levels = c("Condition 4", "Condition 5")))
  
  # Create plot
  p <- ggplot(welch_data, aes(x = welch_t, fill = scenario, color = scenario)) +
    geom_density(alpha = 0.5, size = 1) +
    # Theoretical distribution
    geom_line(data = theoretical_df, 
              aes(x = x, y = density, color = scenario), 
              size = 1.5, linetype = "dashed") +
    # Mean lines
    geom_vline(data = means, aes(xintercept = mean_t, color = scenario),
               linetype = "solid", size = 1.2, alpha = 0.8) +
    # Critical value lines
    geom_vline(data = critical_df, aes(xintercept = critical_pos),
               linetype = "dotted", color = "red", size = 1) +
    geom_vline(data = critical_df, aes(xintercept = critical_neg),
               linetype = "dotted", color = "red", size = 1) +
    
    labs(title = "Welch's t Distribution Comparison: Condition 4 vs 5",
         subtitle = paste0("Red dotted: Critical values (α=.05)\n",
                           "Black dashed: Theoretical distribution\n",
                           "Condition 4: n1=", params4$n1, ", n2=", params4$n2,
                           ", σ1=", params4$s1, ", σ2=", round(params4$s2,3), ", d=", params4$pop_effect_size, "\n",
                           "Condition 5: n1=", params5$n1, ", n2=", params5$n2,
                           ", σ1=", params5$s1, ", σ2=", round(params5$s2,3), ", d=", params5$pop_effect_size),
         x = "t value", 
         y = "Density",
         fill = "Condition",
         color = "Condition") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, family = "Times New Roman"),
      axis.title = element_text(size = 11, family = "Times New Roman"),
      axis.text = element_text(size = 10, family = "Times New Roman"),
      legend.text = element_text(size = 10, family = "Times New Roman"),
      legend.title = element_text(size = 11, face = "bold", family = "Times New Roman"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    
    scale_fill_manual(values = c("#FF7F00", "#4DAF4A")) +
    scale_color_manual(values = c("#FF7F00", "#4DAF4A"))
  
  return(p)
}

# Generate and save plot
welch_comparison <- create_welch_comparison_plot()
print(welch_comparison)

ggsave(
  filename = "./study2/results2/welch_condition4_vs_5_comparison.png",
  plot = welch_comparison,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("\nWelch's t comparison plot for Conditions 4 and 5 saved!\n")