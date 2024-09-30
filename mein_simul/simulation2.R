# Modified simulation.R

# Set the computational node
.libPaths(c(""))
home_dir   <- ""
output_dir <- ""
start_dir  <- getwd()

# Load required libraries
library(rstan)
library(RoBTT)
library(bain)
library(BayesFactor)
library(tidyverse)

# Load the functions and settings
source(file = file.path(home_dir, "functions2.R"))
settings <- readRDS(file = file.path(home_dir, "settings2.RDS"))
tracker  <- "sim_loop"
max_time <- 23.5  # Set maximum runtime in hours


# Modified run_simulation function
run_simulation <- function(current_settings, use_common_seed = TRUE) {
  if (use_common_seed) {
    set.seed(current_settings$seed)
  } else {
    current_seed <- sample.int(1e6, 1)
    set.seed(current_seed)
    current_settings$seed <- current_seed
  }
  
  # Generate data using simulate_data function
  data <- simulate_data(
    mean1 = current_settings$mean1,
    mean2 = current_settings$mean2,
    sd1   = current_settings$sd1,
    sd2   = current_settings$sd2,
    n1    = current_settings$n1,
    n2    = current_settings$n2
  )
  
  # 1. RoBTT
  fit_robtt <- RoBTT(
    x1 = data$x1,
    x2 = data$x2,
    prior_d  = prior("cauchy", list(0, 1/sqrt(2))),
    prior_r  = prior("beta",   list(1.5, 1.5)),
    prior_nu = prior("exp",    list(1)),
    likelihood = c("normal", "t"),
    chains = 2, warmup = 500, iter = 2000,
    parallel = FALSE, seed = current_settings$seed
  )
  
  # 2. Bain - Student t-test
  t_result_student <- t.test(data$x1, data$x2, var.equal = TRUE)
  fit_bain_student <- bain(t_result_student, "x1 = x2")
  
  # 3. Bain - Welch t-test
  t_result_welch <- t.test(data$x1, data$x2, var.equal = FALSE)
  fit_bain_welch <- bain(t_result_welch, "x1 = x2")
  
  # 4. BayesFactor
  fit_bf <- ttestBF(x = data$x1, y = data$x2)
  
  # Collect results
  list(
    robtt = list(
      fit_summary = summary(fit_robtt, conditional = TRUE),
      diagnostics = summary(fit_robtt, diagnostics = TRUE),
      individual_models = summary(fit_robtt, type = "individual")
    ),
    bain_student = fit_bain_student,
    bain_welch = fit_bain_welch,
    bayes_factor = fit_bf,
    true_model = ifelse(current_settings$delta == 0, "H0", "H1"),
    rho = current_settings$rho,
    sdr = current_settings$sdr,
    delta = current_settings$delta,
    scenario = current_settings$scenario,
    seed_used = current_settings$seed
  )
}

# Run simulations
time_start <- Sys.time()
use_common_seed <- TRUE  # Set to FALSE if you want different seeds for each simulation

results <- list()
while(difftime(Sys.time(), time_start, units = "hours") < max_time) {
  # Get the current setting to compute
  loop <- get_loop(home_dir, tracker)
  if(loop > nrow(settings))
    break
  
  # Extract parameters for the current iteration
  current_settings <- settings[loop, ]
  
  # Run simulation
  result <- run_simulation(current_settings, use_common_seed)
  
  # Save the results
  saveRDS(result, file = file.path(output_dir, paste0("results_", current_settings$seed, ".RDS")))
  
  # Store result in the results list
  results[[loop]] <- result
  
  # Update the tracker
  update_tracker(home_dir, tracker)
}

# After simulations, you can add code here to collect and analyze the results
# For example:
collect_and_analyze_results <- function(results) {
  results_df <- map_dfr(results, ~data.frame(
    scenario = .x$scenario,
    delta = .x$delta,
    rho = .x$rho,
    sdr = .x$sdr,
    true_model = .x$true_model,
    BF_robtt = .x$robtt$fit_summary$BF10,
    BF_bain_student = .x$bain_student$fit$BF[1, "bf"],
    BF_bain_welch = .x$bain_welch$fit$BF[1, "bf"],
    BF_bf = exp(.x$bayes_factor@bayesFactor$bf),
    seed_used = .x$seed_used,
    stringsAsFactors = FALSE
  ))
  
  return(results_df)
}

# Analyze results
results_df <- collect_and_analyze_results(results)

# Example visualization
ggplot(results_df, aes(x = sdr)) +
  geom_point(aes(y = log(BF_robtt), color = "RoBTT"), alpha = 0.5) +
  geom_point(aes(y = log(BF_bain_student), color = "Bain Student"), alpha = 0.5) +
  geom_point(aes(y = log(BF_bain_welch), color = "Bain Welch"), alpha = 0.5) +
  geom_point(aes(y = log(BF_bf), color = "BayesFactor"), alpha = 0.5) +
  geom_smooth(aes(y = log(BF_robtt), color = "RoBTT"), method = "loess") +
  geom_smooth(aes(y = log(BF_bain_student), color = "Bain Student"), method = "loess") +
  geom_smooth(aes(y = log(BF_bain_welch), color = "Bain Welch"), method = "loess") +
  geom_smooth(aes(y = log(BF_bf), color = "BayesFactor"), method = "loess") +
  facet_grid(delta ~ scenario) +
  labs(title = "Log Bayes Factor vs SDR",
       x = "SDR", y = "Log Bayes Factor", color = "Method") +
  theme_minimal()

ggsave("BF_comparison.png", width = 12, height = 8)