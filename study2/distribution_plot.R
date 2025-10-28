#### Comparison of Theoretical Distribution vs Sampling Distribution ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(moments)  # For skewness and kurtosis calculation
library(writexl)  # For Excel output

# Load necessary functions
source("./study2/t-distribution.R")  # File containing theoretical distribution functions
source("./study2/Behrens-Fisher_distribution.R")

#### Get all RDS files ####
rds_files <- list.files("./study2/results5/", pattern = "^S24results.*\\.RDS$", full.names = TRUE)
cat(sprintf("Found %d RDS files\n", length(rds_files)))
cat("Files:\n")
print(basename(rds_files))
cat("\n")

#### Initialize list to store all descriptive statistics ####
all_descriptive_stats <- list()

#### Process each RDS file ####
for (file_path in rds_files) {
  
  file_name <- basename(file_path)
  file_id <- gsub("S24results|\\.RDS", "", file_name)  # Extract identifier (e.g., "60_E0")
  
  # Extract sample size and effect size from file_id
  # file_id format: "60_E0", "120_E2", "240_E5", "240_E8" etc.
  parts <- strsplit(file_id, "_E")[[1]]
  total_n <- as.numeric(parts[1])
  effect_size <- as.numeric(parts[2]) / 10  # E0 -> 0, E2 -> 0.2, E5 -> 0.5, E8 -> 0.8
  
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse="")))
  cat(sprintf("Processing: %s (N=%d, d=%.1f)\n", file_name, total_n, effect_size))
  cat(sprintf("%s\n", paste(rep("=", 80), collapse="")))
  
  # Load results
  results <- readRDS(file_path)
  
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
  cat("\n=== Population Parameters by Condition ===\n")
  for(i in 1:5) {
    params <- condition_params[[i]]
    cat(sprintf("Condition %d: n1=%d, n2=%d, σ1=%.1f, σ2=%.3f, pop_effect_size=%.2f, equal_var=%s\n", 
                i, params$n1, params$n2, params$s1, params$s2, params$pop_effect_size, params$equal_var))
  }
  cat("\n")
  
  #### Calculate Descriptive Statistics for Sampling Distributions ####
  cat("=== Calculating Descriptive Statistics ===\n")
  
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
    ) %>%
    mutate(
      file_id = file_id,
      .before = 1
    )
  
  # Print results for this file
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
  
  # Add to the list
  all_descriptive_stats[[file_id]] <- descriptive_stats
  
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
  
  #### Function to Create Integrated Comparison Plot ####
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
    
    # Create plot title
    plot_title <- sprintf("All Conditions: Theoretical vs Sampling Distribution (N = %d, d = %.1f)", 
                          total_n, effect_size)
    
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
      
      labs(title = plot_title,
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
      
      scale_fill_manual(values = c("#FFA500", "#0066CC")) +
      scale_color_manual(values = c("#FFA500", "#0066CC"))
    
    return(p)
  }
  
  # Generate integrated comparison plot
  cat("=== Generating Plot ===\n")
  integrated_plot <- create_integrated_comparison()
  
  # Save plot
  output_filename <- paste0("./study2/results5/allcon_", file_id, ".png")
  ggsave(
    filename = output_filename,
    plot = integrated_plot,
    width = 11,
    height = 8,
    dpi = 600,
    device = "png",
    bg = "white"
  )
  cat(sprintf("Plot saved: %s\n", output_filename))
}

#### Combine all descriptive statistics and save to Excel ####
cat(sprintf("\n%s\n", paste(rep("=", 80), collapse="")))
cat("=== Combining All Descriptive Statistics ===\n")

combined_stats <- bind_rows(all_descriptive_stats)

# Round to 3 decimal places
combined_stats <- combined_stats %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Save to Excel
output_excel <- "./study2/results4/all_descriptive_stats.xlsx"
write_xlsx(combined_stats, output_excel)

cat(sprintf("All descriptive statistics saved to: %s\n", output_excel))
cat(sprintf("Total rows: %d (should be %d files × 5 conditions = %d)\n", 
            nrow(combined_stats), 
            length(rds_files), 
            length(rds_files) * 5))

# Also save as CSV as backup
output_csv <- "./study2/results4/all_descriptive_stats.csv"
write.csv(combined_stats, output_csv, row.names = FALSE)
cat(sprintf("CSV backup saved to: %s\n", output_csv))

cat(sprintf("\n%s\n", paste(rep("=", 80), collapse="")))
cat("=== Processing Complete ===\n")
cat(sprintf("%s\n", paste(rep("=", 80), collapse="")))
