#create setting2_seed 변수명에 대응되는 코드. 나중에 일반화 필요
# Load required libraries
library(rstan)
library(RoBTT)
library(bain)
library(BayesFactor)
library(tidyverse)

# set.seed(1234)

# Set the computational node
home_dir   <- "."
output_dir <- file.path("D:/post_results") #caution!! 경로 임시로 수정
start_dir  <- getwd()

# Load the functions and settings
source(file = file.path(home_dir, "./mein_simul/functions2.R"))
settings <- readRDS(file = file.path(home_dir, "/post/post_settings3.RDS")) #2 -> 3(more scenario)

print(head(settings))
print(nrow(settings))

tracker  <- "sim_loop"
max_time <- 24.0  # Set maximum runtime in hours

# Modified run_simulation function
run_simulation <- function(current_settings) {
  # Generate data using simulate_data function
  data <- simulate_data(
    mean1 = current_settings$mean1,
    mean2 = current_settings$mean2,
    sd1   = current_settings$sd1,
    sd2   = current_settings$sd2,
    n1    = current_settings$n1,
    n2    = current_settings$n2,
    rho   = current_settings$rho 
  )
  
  # 1. RoBTT
  fit_robtt <- RoBTT(
    x1 = data$x1,
    x2 = data$x2,
    prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
    prior_rho  = prior("beta",   list(1.5, 1.5)),
    prior_nu = prior("exp",    list(1)),
    chains = 2, warmup = 1000, iter = 5000,
    parallel = FALSE
  )
  
  robtt_summary <- summary(fit_robtt, conditional = TRUE)
  
  # 2. Bain - Student t-test
  t_result_student <- t_test(data$x1, data$x2, var.equal = TRUE)
  fit_bain_student <- bain(t_result_student, "x = y")
  
  # 3. Bain - Welch t-test
  t_result_welch <- t_test(data$x1, data$x2, var.equal = FALSE)
  fit_bain_welch <- bain(t_result_welch, "x = y")
  
  # 4. BayesFactor
  fit_bf <- ttestBF(x = data$x1, y = data$x2)
  
  # 5. Wetzels' MCMC t-test
  fit_wetzels <- wetzels_ttest(data$x1, data$x2)
  
  # Collect results
  results <- list(
    robtt = list(
      fit_summary = summary(fit_robtt, diagnostics = TRUE),
      individual_models = summary(fit_robtt, type = "individual")
    ),
    bain_student = fit_bain_student,
    bain_welch = fit_bain_welch,
    bayes_factor = fit_bf,
    wetzels = fit_wetzels,
    true_model = ifelse(current_settings$delta == 0, "H0", "H1"),
    rho = current_settings$rho,
    sdr = current_settings$sdr,
    delta = current_settings$delta,
    scenario = current_settings$scenario
  )
  
  return(results)
}

# Run simulations
time_start <- Sys.time()
print("Starting simulations...")

results <- list()
max_loop <- nrow(settings)

# 기존 결과 파일에서 가장 최근의 체크포인트 생성
create_checkpoint_from_results <- function(output_dir) {
  result_files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
  if (length(result_files) == 0) {
    print("No result files found.")
    return(1)  # 결과 파일이 없으면 1부터 시작
  }
  file_numbers <- as.numeric(gsub(".*results_(\\d+)\\.RDS", "\\1", result_files))
  last_number <- max(file_numbers)
  print(paste("Last completed simulation:", last_number))
  return(last_number + 1)  # 다음 번호부터 시작
}

# 시뮬레이션 시작 전 처리
start_loop <- create_checkpoint_from_results(output_dir)
print(paste("Starting simulation from loop:", start_loop))

# 기존 시뮬레이션 루프 수정
for (loop in start_loop:max_loop) {
  if(difftime(Sys.time(), time_start, units = "hours") >= max_time) break
  
  print(paste("Processing loop:", loop))
  
  # Extract parameters for the current iteration
  current_settings <- settings[loop, ]
  print(paste("Current settings:", toString(current_settings)))
  
  # Set the seed for this iteration
  set.seed(current_settings$seed)               #seed가 setting에 저장되어 있는 경우인 create_settings2에 맞게 쓰여 있음. 나중에 state를 만들 경우 수정 필요
  
  # Run simulation
  print("Running simulation...")
  tryCatch({
    result <- run_simulation(current_settings)
    
    # Save the results for each iteration
    result_file <- file.path(output_dir, paste0("results_", loop, ".RDS"))
    saveRDS(result, file = result_file)
    
    if(file.exists(result_file)) {
      print(paste("Results saved successfully for loop:", loop))
    } else {
      print(paste("Failed to save results for loop:", loop))
    }
    
    # Update the tracker
    update_tracker(home_dir, tracker)
    
    
  }, error = function(e) {
    print(paste("Error in simulation for loop:", loop, "- Error:", e$message))
    # Optionally, you can still save the error information
    saveRDS(list(error = e$message), file = file.path(output_dir, paste0("error_", loop, ".RDS")))
  })
  
  # Use get_loop to update the loop counter in the tracker file
  get_loop(home_dir, tracker, max_loop)
}

end_time <- Sys.time()
print(paste("Total execution time:", difftime(end_time, time_start, units = "hours"), "hours"))

print("Simulations completed.")
