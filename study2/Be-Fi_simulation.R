source("New/functionsN.R")
library(tidyverse)
library(asht)  # bfTest

output_dir <- "./study2/results5"
intermediate_dir <- file.path("D:S24results") #효과크기 별로 폴더 구분
# 폴더 존재 확인
if (!dir.exists("./study2/results4")) dir.create("./study2/results5", recursive = TRUE)
if (!dir.exists("D:/S24results")) dir.create("D:/S24results", recursive = TRUE)
dir.exists("./study2/results5")
dir.exists("D:/S24results")

# 시나리오 설정 및 복제본 생성을 위한 통합 함수
create_settings <- function(var1, var2, effect_size, n1, n2, replications = 50000) {
  # 표준편차 계산
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  
  # 평균 분산 계산
  av_var <- mean(c(var1, var2))
  
  # 효과 크기에 따른 평균 차이 계산
  mu_diff <- effect_size * sqrt(av_var)
  
  # 각 집단의 평균 계산
  mu1 <- mu_diff/2
  mu2 <- -mu_diff/2
  
  # 복제본 생성
  tibble(
    var1 = var1,
    var2 = var2,
    sd1 = sd1,
    sd2 = sd2,
    effect_size = effect_size,
    av_var = av_var,
    mu_diff = mu_diff,
    mu1 = mu1,
    mu2 = mu2,
    n1 = n1,
    n2 = n2,
    varr = var2 / var1,
    replication = 1:replications,
    seed = sample.int(.Machine$integer.max, replications)
  )
}

# 각 시나리오 설정 생성(60, 120, 240)
scenarios_config <- list(
  list(var1 = 4, var2 = 4, n1 = 120, n2 = 120, scenario = 1),
  list(var1 = 4, var2 = 4, n1 = 144, n2 = 96, scenario = 2),
  list(var1 = 4, var2 = 2, n1 = 120, n2 = 120, scenario = 3),
  list(var1 = 4, var2 = 2, n1 = 144, n2 = 96, scenario = 4),
  list(var1 = 4, var2 = 2, n1 = 96, n2 = 144, scenario = 5)
)

# 시뮬레이션 함수 (메모리 효율적) - bfTest 추가로 수정됨
run_simulation_batch <- function(scenario_config, effect_size = 0, replications = 50000) {
  
  cat("Processing scenario", scenario_config$scenario, 
      "- n1:", scenario_config$n1, "n2:", scenario_config$n2, 
      "var1:", scenario_config$var1, "var2:", scenario_config$var2, "\n")
  
  # 설정 생성
  settings <- create_settings(
    var1 = scenario_config$var1,
    var2 = scenario_config$var2,
    effect_size = effect_size,
    n1 = scenario_config$n1,
    n2 = scenario_config$n2,
    replications = replications
  ) %>%
    mutate(scenario = scenario_config$scenario)
  
  # 결과 저장용 벡터 미리 할당 (메모리 효율성) - bfTest 벡터들 추가
  student_p_values <- numeric(replications)
  welch_p_values <- numeric(replications)
  bf_p_values <- numeric(replications)  # bfTest p-value 추가
  
  student_t_values <- numeric(replications)
  welch_t_values <- numeric(replications)
  bf_t_values <- numeric(replications)  # bfTest t-value 추가
  
  student_df <- numeric(replications)  # Student t-test 자유도
  welch_df <- numeric(replications)    # Welch t-test 자유도
  bf_R_values <- numeric(replications) # bfTest R 값 추가
  
  # 표본 통계량 저장용 벡터들
  sample_mean1 <- numeric(replications)
  sample_mean2 <- numeric(replications)
  sample_sd1 <- numeric(replications)
  sample_sd2 <- numeric(replications)
  
  # 효과크기 저장용 벡터
  pooled_d <- numeric(replications)
  hete_d <- numeric(replications) # 이분산 가정
  
  # 배치 처리 (메모리 관리를 위해)
  batch_size <- 1000  # 한 번에 처리할 시뮬레이션 수
  n_batches <- ceiling(replications / batch_size)
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, replications)
    current_batch_size <- end_idx - start_idx + 1
    
    if (batch %% 10 == 0) {
      cat("  Batch", batch, "of", n_batches, "completed\n")
    }
    
    # 현재 배치에 대한 시뮬레이션 실행
    for (i in 1:current_batch_size) {
      idx <- start_idx + i - 1
      current_row <- settings[idx, ]
      
      # 시드 설정
      set.seed(current_row$seed)
      
      # 데이터 생성
      data <- simulate_data(
        mean1 = current_row$mu1,
        mean2 = current_row$mu2,
        sd1   = current_row$sd1,
        sd2   = current_row$sd2,
        n1    = current_row$n1,
        n2    = current_row$n2,
        seed  = current_row$seed
      )
      
      # t-test 실행 및 bfTest 추가
      student_test <- t.test(data$x1, data$x2, var.equal = TRUE)
      welch_test <- t.test(data$x1, data$x2, var.equal = FALSE)
      bf_test <- bfTest(data$x1, data$x2)  # bfTest 추가
      
      # p-values와 t-statistics 저장 (bfTest 결과 추가)
      student_p_values[idx] <- student_test$p.value
      welch_p_values[idx] <- welch_test$p.value
      bf_p_values[idx] <- bf_test$p.value  # bfTest p-value
      
      student_t_values[idx] <- student_test$statistic
      welch_t_values[idx] <- welch_test$statistic
      bf_t_values[idx] <- bf_test$statistic  # bfTest t-value
      
      student_df[idx] <- student_test$parameter
      welch_df[idx] <- welch_test$parameter
      bf_R_values[idx] <- bf_test$parameter  # bfTest R 값
      
      # 표본 통계량 저장
      sample_mean1[idx] <- mean(data$x1)
      sample_mean2[idx] <- mean(data$x2)
      sample_sd1[idx] <- sd(data$x1)
      sample_sd2[idx] <- sd(data$x2)
      
      # 효과 크기 계산
      n1 <- current_row$n1
      n2 <- current_row$n2
      
      # Cohen's d (pooled standard deviation 사용)
      pooled_sd <- sqrt(((n1-1) * sample_sd1[idx]^2 + (n2-1) * sample_sd2[idx]^2) / (n1 + n2 - 2))
      pooled_d[idx] <- (sample_mean1[idx] - sample_mean2[idx]) / pooled_sd
      arith_sd <- sqrt((sample_sd1[idx]^2+sample_sd2[idx]^2)/2)
      hete_d[idx] <- (sample_mean1[idx] - sample_mean2[idx]) / arith_sd
    }
  }
  
  # 최종 결과 생성 (bfTest 결과 포함)
  results <- tibble(
    # 시나리오 정보
    scenario = scenario_config$scenario,
    replication = 1:replications,
    
    # 설계 정보 (모집단 모수)
    n1 = scenario_config$n1,
    n2 = scenario_config$n2,
    var1 = scenario_config$var1,
    var2 = scenario_config$var2,
    sd1 = sqrt(scenario_config$var1),
    sd2 = sqrt(scenario_config$var2),
    
    # 모집단 모수 (설계된 값)
    pop_mean1 = settings$mu1,
    pop_mean2 = settings$mu2,
    pop_effect_size = effect_size,  # 모집단 효과크기
    
    # 시드 정보 (재현성)
    seed = settings$seed,
    
    # 표본 통계량 (실제 생성된 데이터)
    sample_mean1 = sample_mean1,
    sample_mean2 = sample_mean2,
    sample_sd1 = sample_sd1,
    sample_sd2 = sample_sd2,
    sample_mean_diff = sample_mean1 - sample_mean2,
    
    # 검정 통계량 (bfTest 추가)
    student_t = student_t_values,
    welch_t = welch_t_values,
    bf_t = bf_t_values,  # bfTest t-value 추가
    
    student_p = student_p_values,
    welch_p = welch_p_values,
    bf_p = bf_p_values,  # bfTest p-value 추가
    
    student_df = student_df,
    welch_df = welch_df,
    bf_R = bf_R_values,  # bfTest R 값 추가
    
    # 표본 효과크기들
    cohens_d = pooled_d,           # 표준적인 Cohen's d (변수명 수정)
    hete_d = hete_d,              # 이분산 가정 효과크기
    
    # 표본 효과크기 추가 계산
    glass_delta1 = sample_mean_diff / sample_sd1,  # Glass's Δ (control group SD 사용)
    glass_delta2 = sample_mean_diff / sample_sd2,  # Glass's Δ (experimental group SD 사용)
    hedges_g = cohens_d * (1 - 3/(4*(n1 + n2 - 2) - 1))  # Hedges' g (small sample correction)
  )
  
  return(results)
}

# 메인 시뮬레이션 실행 함수
run_full_simulation <- function(effect_size = 0, replications = 50000, 
                                output_file = NULL, save_intermediate = TRUE, 
                                intermediate_dir = "intermediate_results") {
  
  start_time <- Sys.time()
  cat("Starting simulation with effect size =", effect_size, 
      "and", replications, "replications per scenario\n")
  
  # 중간 결과 저장 디렉토리 생성
  if (save_intermediate && !dir.exists(intermediate_dir)) {
    dir.create(intermediate_dir, recursive = TRUE)
    cat("Created intermediate directory:", intermediate_dir, "\n")
  }
  
  # 모든 시나리오 결과를 저장할 리스트
  all_results <- list()
  
  # 각 시나리오별로 시뮬레이션 실행
  for (i in seq_along(scenarios_config)) {
    scenario_start <- Sys.time()
    
    results <- run_simulation_batch(
      scenario_config = scenarios_config[[i]],
      effect_size = effect_size,
      replications = replications
    )
    
    all_results[[i]] <- results
    
    scenario_end <- Sys.time()
    cat("Scenario", i, "completed in", 
        round(difftime(scenario_end, scenario_start, units = "mins"), 2), "minutes\n")
    
    # 중간 저장 옵션
    if (save_intermediate) {
      if (!is.null(output_file)) {
        # output_file이 있는 경우: 해당 파일명 기반으로 저장
        base_name <- tools::file_path_sans_ext(basename(output_file))
        temp_file <- file.path(intermediate_dir, paste0(base_name, "_scenario_", i, ".RDS"))
      } else {
        # output_file이 없는 경우: 기본 파일명 사용
        temp_file <- file.path(intermediate_dir, paste0("simulation_scenario_", i, "_es_", effect_size, ".RDS"))
      }
      
      saveRDS(results, file = temp_file)
      cat("Intermediate results saved:", temp_file, "\n")
    }
    
    # 가비지 컬렉션으로 메모리 정리
    gc()
  }
  
  # 모든 결과 합치기
  cat("Combining all results...\n")
  final_results <- bind_rows(all_results)
  
  # 최종 저장
  if (!is.null(output_file)) {
    saveRDS(final_results, file = output_file)
    cat("Final results saved:", output_file, "\n")
  }
  
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "hours")
  cat("Total simulation completed in", round(total_time, 2), "hours\n")
  cat("Total simulations:", nrow(final_results), "\n")
  
  return(final_results)
}


# 실행
cat("simulation")
results <- run_full_simulation(
  effect_size = 0.0, # 효과크기 조건에 맞게 수정 
  replications = 30000,
  output_file = file.path(output_dir, "S24results240_E0.RDS"), #파일명 수정
  save_intermediate = TRUE,
  intermediate_dir = intermediate_dir
) 

 # 대략적인 내용 및 결측치 확인
print(results)
names(results)
sum(is.na(results))