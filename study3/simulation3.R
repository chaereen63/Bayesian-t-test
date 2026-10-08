# Study 3: JZS, GICA(BFGC), Bain, RoBTT 통합 시뮬레이션
# - 데이터 생성/GICA 함수: Korea_pub/functionsK.R
# - RoBTT, Bain 부분: post/post_simul.R 구조를 가져옴
# Load required libraries
library(rstan)
library(BayesFactor)
library(RoBTT)
library(bain)
library(tidyverse)

# Set the computational node
home_dir   <- "."
output_dir <- file.path(home_dir, "study3/results") # 필요하면 외장 경로로 수정 (예: "D:/resultsS3")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
start_dir  <- getwd()

# Load the functions and settings
source(file = file.path(home_dir, "./Korea_pub/functionsK.R"))
settings <- readRDS(file = file.path(home_dir, "./Korea_pub/settingsK30.RDS")) # study3용 settings로 교체하기

print(head(settings))
print(nrow(settings))

tracker  <- "sim_loop"
max_time <- 24.0  # Set maximum runtime in hours

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
    chains = 2, warmup = 1000, iter = 5000,
    parallel = FALSE
  )
  fit_robtt_ind <- summary(fit_robtt, type = "individual")
  
  # 2. Bain - Student t-test
  t_result_student <- t_test(data$x1, data$x2, var.equal = TRUE)
  fit_bain_student <- bain(t_result_student, "x = y")
  
  # 3. Bain - Welch t-test
  t_result_welch <- t_test(data$x1, data$x2, var.equal = FALSE)
  fit_bain_welch <- bain(t_result_welch, "x = y")
  
  # 4. BayesFactor_JZS
  fit_bf <- ttestBF(x = data$x1, y = data$x2)
  
  # 5. Giron & Del Castillo (GICA)
  fit_gica <- gicaBF(x1 = data$x1, x2 = data$x2)
  
  # Student's t와 Welch's t의 p-value만 저장
  student_p <- t.test(data$x1, data$x2, var.equal = TRUE)$p.value
  welch_p   <- t.test(data$x1, data$x2, var.equal = FALSE)$p.value
  
  # Collect results
  results <- list(
    robtt = list(
      fit_summary = summary(fit_robtt),  # components의 Effect inclusion BF (post_merge.R 방식)
      bf_effect = list(                  # 등분산/이분산 조건부 효과 BF (mergeK.R 방식)
        homoBF = as.numeric(fit_robtt_ind$models[[5]]$summary[4,2]) / as.numeric(fit_robtt_ind$models[[1]]$summary[4,2]),
        heteBF = as.numeric(fit_robtt_ind$models[[7]]$summary[4,2]) / as.numeric(fit_robtt_ind$models[[3]]$summary[4,2])
      )
    ),
    bain_student = fit_bain_student,
    bain_welch   = fit_bain_welch,
    bayes_factor = fit_bf,    # 전체 결과 저장
    gica         = fit_gica,  # 전체 결과 저장
    student_p    = student_p,
    welch_p      = welch_p,
    rho      = current_settings$rho,
    sdr      = current_settings$sdr,
    delta    = current_settings$delta,
    scenario = current_settings$scenario
  )
  
  return(results)
}

# Run simulations
time_start <- Sys.time()
print("Starting simulations...")

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

if (start_loop > max_loop) {
  print("All simulations already completed.")
} else {
  for (loop in start_loop:max_loop) {
    if (difftime(Sys.time(), time_start, units = "hours") >= max_time) break
    
    print(paste("Processing loop:", loop))
    
    # Extract parameters for the current iteration
    current_settings <- settings[loop, ]
    print(paste("Current settings:", toString(current_settings)))
    
    # Set the seed for this iteration (settings에 seed 열이 있어야 함)
    set.seed(current_settings$seed)
    
    # Run simulation
    print("Running simulation...")
    tryCatch({
      result <- run_simulation(current_settings)
      
      # Save the results for each iteration
      result_file <- file.path(output_dir, paste0("results_", loop, ".RDS"))
      saveRDS(result, file = result_file)
      
      if (file.exists(result_file)) {
        print(paste("Results saved successfully for loop:", loop))
      } else {
        print(paste("Failed to save results for loop:", loop))
      }
      
      # Update the tracker
      update_tracker(home_dir, tracker)
      
    }, error = function(e) {
      print(paste("Error in simulation for loop:", loop, "- Error:", e$message))
      saveRDS(list(error = e$message), file = file.path(output_dir, paste0("error_", loop, ".RDS")))
    })
    
    # Use get_loop to update the loop counter in the tracker file
    get_loop(home_dir, tracker, max_loop)
  }
}

end_time <- Sys.time()
print(paste("Total execution time:", difftime(end_time, time_start, units = "hours"), "hours"))

print("Simulations completed.")
