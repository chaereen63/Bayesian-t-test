# Set the computational node
.libPaths(c(""))
home_dir   <- ""
output_dir <- ""
start_dir  <- getwd()

# Load required libraries
library(tidyverse)

# Helper function to safely extract values
safe_extract <- function(x, default = NA) {
  tryCatch(x, error = function(e) default)
}

# Function to safely extract and ensure single value
safe_extract_single <- function(x, default = NA) {
  value <- safe_extract(x)
  if (length(value) > 1) {
    warning("Multiple values found, using the first one")
    return(value[1])
  } else if (length(value) == 0) {
    warning("No value found, using default")
    return(default)
  } else {
    return(value)
  }
}

# Function to check if a simulation result is complete
is_simulation_complete <- function(result) {
  required_fields <- c("scenario", "delta", "rho", "sdr", "true_model", "robtt", "bain_student", "bain_welch", "bayes_factor")
  all(required_fields %in% names(result)) &&
    !is.null(result$robtt$fit_summary) &&
    !is.null(result$bain_student$fit) &&
    !is.null(result$bain_welch$fit) &&
    !is.null(result$bayes_factor@bayesFactor)
}

# Function to process a single result file
process_result_file <- function(f) {
  tryCatch({
    cat("Processing file:", f, "\n")
    result <- readRDS(f)
    
    # Check if the simulation is complete
    complete <- is_simulation_complete(result)
    cat("Simulation complete:", complete, "\n")
    
    if (!complete) {
      cat("Incomplete simulation. Skipping this file.\n")
      return(list(df = data.frame(), complete = FALSE, error = "Incomplete simulation"))
    }
    
    # Extract results individually with detailed debugging
    extracted <- list()
    
    # Basic information
    extracted$scenario <- safe_extract_single(result$scenario)
    extracted$delta <- safe_extract_single(result$delta)
    extracted$rho <- safe_extract_single(result$rho)
    extracted$sdr <- safe_extract_single(result$sdr)
    extracted$true_model <- safe_extract_single(result$true_model)
    
    # RoBTT result
    cat("Extracting RoBTT result...\n")
    if (!is.null(result$robtt) && !is.null(result$robtt$fit_summary)) {
      cat("RoBTT fit summary structure:\n")
      print(str(result$robtt$fit_summary))
      extracted$BF_robtt <- safe_extract_single(result$robtt$fit_summary$inclusionBF["Effect"])
    } else {
      cat("RoBTT result is missing or incomplete\n")
      extracted$BF_robtt <- NA
    }
    
    # Bain Student result
    cat("Extracting Bain Student result...\n")
    if (!is.null(result$bain_student) && !is.null(result$bain_student$fit)) {
      cat("Bain Student fit structure:\n")
      print(str(result$bain_student$fit))
      extracted$BF_bain_student <- safe_extract_single(result$bain_student$fit$BF[1, "BF.u"])
    } else {
      cat("Bain Student result is missing or incomplete\n")
      extracted$BF_bain_student <- NA
    }
    
    # Bain Welch result
    cat("Extracting Bain Welch result...\n")
    if (!is.null(result$bain_welch) && !is.null(result$bain_welch$fit)) {
      cat("Bain Welch fit structure:\n")
      print(str(result$bain_welch$fit))
      extracted$BF_bain_welch <- safe_extract_single(result$bain_welch$fit$BF[1, "BF.u"])
    } else {
      cat("Bain Welch result is missing or incomplete\n")
      extracted$BF_bain_welch <- NA
    }
    
    # BayesFactor result
    cat("Extracting BayesFactor result...\n")
    if (!is.null(result$bayes_factor) && !is.null(result$bayes_factor@bayesFactor)) {
      cat("BayesFactor structure:\n")
      print(str(result$bayes_factor@bayesFactor))
      extracted$BF_bf <- safe_extract_single(exp(result$bayes_factor@bayesFactor$bf))
    } else {
      cat("BayesFactor result is missing or incomplete\n")
      extracted$BF_bf <- NA
    }
    
    extracted$seed_used <- as.numeric(gsub("results_(\\d+)\\.RDS", "\\1", basename(f)))
    
    # Print extracted values and their lengths
    cat("Extracted values:\n")
    for (name in names(extracted)) {
      cat(name, ": ", extracted[[name]], " (length: ", length(extracted[[name]]), ")\n", sep = "")
    }
    
    # Create data frame
    df <- as.data.frame(extracted, stringsAsFactors = FALSE)
    
    cat("Data frame created with", nrow(df), "rows and", ncol(df), "columns\n")
    cat("Columns:", paste(names(df), collapse=", "), "\n")
    
    return(list(df = df, complete = TRUE, error = NULL))
  }, error = function(e) {
    cat("Error processing file:", f, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    return(list(df = data.frame(), complete = FALSE, error = conditionMessage(e)))
  })
}

# Function to merge results from different models
merge_results <- function(output_dir) {
  # List all result files
  files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  cat("Found", length(files), "result files\n")
  
  # Process each file
  results <- lapply(files, process_result_file)
  
  # Count complete and incomplete simulations
  complete_simulations <- sum(sapply(results, function(r) r$complete))
  incomplete_simulations <- length(results) - complete_simulations
  
  # Print summary before merging
  cat("\nSummary before merging:\n")
  cat("Complete simulations:", complete_simulations, "\n")
  cat("Incomplete simulations:", incomplete_simulations, "\n")
  
  # Ask user if they want to proceed with merging
  proceed <- readline(prompt="Do you want to proceed with merging the complete results? (y/n): ")
  
  if (tolower(proceed) != "y") {
    cat("Merging cancelled by user.\n")
    return(NULL)
  }
  
  # Combine all complete data frames
  merged_df <- do.call(rbind, lapply(results[sapply(results, function(r) r$complete)], function(r) r$df))
  cat("\nMerged data frame has", nrow(merged_df), "rows and", ncol(merged_df), "columns\n")
  
  # Print errors if any
  errors <- sapply(results, function(r) r$error)
  error_files <- files[!sapply(errors, is.null)]
  if (length(error_files) > 0) {
    cat("\nFiles with errors:\n")
    for (i in seq_along(error_files)) {
      cat(basename(error_files[i]), ":", errors[[i]], "\n")
    }
  }
  
  return(merged_df)
}

# Main execution
cat("Starting merge_results function...\n")
merged_results <- merge_results(output_dir)
cat("merge_results function completed\n")

# Check the structure of merged_results
if (!is.null(merged_results)) {
  cat("Structure of merged_results:\n")
  print(str(merged_results))
  cat("Column names of merged_results:\n")
  print(names(merged_results))
  cat("Summary of merged_results:\n")
  print(summary(merged_results))
  
  # Save merged results
  saveRDS(merged_results, file = file.path(home_dir, "merged_results.RDS"))
  cat("Merged results saved to", file.path(home_dir, "merged_results.RDS"), "\n")
  
  # Optional: Create summary statistics or visualizations
  if (nrow(merged_results) > 0) {
    summary_stats <- merged_results %>%
      group_by(scenario, delta, rho, sdr) %>%
      summarise(across(starts_with("BF_"), 
                       list(mean = ~mean(., na.rm = TRUE),
                            median = ~median(., na.rm = TRUE),
                            sd = ~sd(., na.rm = TRUE)),
                       .names = "{.col}_{.fn}"),
                .groups = "drop")
    
    # Save summary statistics
    write.csv(summary_stats, file = file.path(home_dir, "summary_statistics.csv"), row.names = FALSE)
    cat("Summary statistics saved to", file.path(home_dir, "summary_statistics.csv"), "\n")
    
    # Create visualization
    p <- ggplot(merged_results, aes(x = rho)) +
      geom_point(aes(y = log(BF_robtt), color = "RoBTT"), alpha = 0.5) +
      geom_point(aes(y = log(BF_bain_student), color = "Bain Student"), alpha = 0.5) +
      geom_point(aes(y = log(BF_bain_welch), color = "Bain Welch"), alpha = 0.5) +
      geom_point(aes(y = log(BF_bf), color = "BayesFactor"), alpha = 0.5) +
      geom_smooth(aes(y = log(BF_robtt), color = "RoBTT"), method = "loess", se = FALSE) +
      geom_smooth(aes(y = log(BF_bain_student), color = "Bain Student"), method = "loess", se = FALSE) +
      geom_smooth(aes(y = log(BF_bain_welch), color = "Bain Welch"), method = "loess", se = FALSE) +
      geom_smooth(aes(y = log(BF_bf), color = "BayesFactor"), method = "loess", se = FALSE) +
      facet_grid(delta ~ scenario) +
      labs(title = "Log Bayes Factor vs Rho",
           x = "Rho", y = "Log Bayes Factor", color = "Method") +
      theme_minimal() +
      scale_y_continuous(limits = c(NA, NA), oob = scales::squish)
    
    # Save the plot
    ggsave(file.path(home_dir, "BF_comparison.png"), plot = p, width = 12, height = 8)
    cat("Visualization saved to", file.path(home_dir, "BF_comparison.png"), "\n")
  } else {
    cat("No data available for summary statistics or visualization.\n")
  }
} else {
  cat("No merged results available for analysis.\n")
}

# Print a summary of missing values
if (!is.null(merged_results)) {
  cat("Summary of missing values:\n")
  print(colSums(is.na(merged_results)))
}

cat("Script execution completed.\n")