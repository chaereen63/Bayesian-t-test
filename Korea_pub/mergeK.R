library(dplyr)
library(purrr)

home_dir <- "./Korea_pub"
output_dir <- "D:/results30" #저장경로 수정하기

# 안전한 추출 함수
safe_extract <- function(x, default = NA) {
  tryCatch(x, error = function(e) default)
}

# BayesFactor 객체에서 BF10 추출하는 함수
extract_bayes_factor <- function(bf_object) {
  tryCatch({
    if (inherits(bf_object, "BFBayesFactor")) {
      # numerator / denominator로 BF10 계산
      as.numeric(exp(bf_object@bayesFactor$bf))
    } else {
      NA
    }
  }, error = function(e) NA)
}

#1. 단일 결과를 처리하는 함수 (RoBTT포함)
process_result <- function(temp_result) {
  list(
    # RoBTT의 homo/hetero BF 추출
    BF_robtt_homo = safe_extract(temp_result$robtt$bf_effect$homoBF),
    BF_robtt_hete = safe_extract(temp_result$robtt$bf_effect$heteBF),
    
    # BayesFactor의 BF10 추출
    BF_jzs = extract_bayes_factor(temp_result$bayes_factor),
    
    # GICA의 BF10 추출
    BF_gica = safe_extract(temp_result$gica$bf10),
    
    # 시나리오 정보 추출
    rho = safe_extract(temp_result$rho),
    sdr = safe_extract(temp_result$sdr),
    delta = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario)
  )
}

# 모든 결과 파일 읽고 처리
files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
results <- map(files, ~{
  temp_result <- readRDS(.x)
  process_result(temp_result)
})

# 결과를 데이터 프레임으로 변환
results_df <- bind_rows(results)

# 행 이름 제거
rownames(results_df) <- NULL

# 결과 저장
saveRDS(results_df, file = file.path(home_dir, "final_merged_resultsK.RDS"))

# 시나리오별 요약 통계 계산
summary_by_scenario <- results_df %>%
  group_by(scenario) %>%
  summarise(
    mean_robtt_homo = mean(BF_robtt_homo, na.rm = TRUE),
    mean_robtt_hete = mean(BF_robtt_hete, na.rm = TRUE),
    mean_jzs = mean(BF_jzs, na.rm = TRUE),
    mean_gica = mean(BF_gica, na.rm = TRUE)
  )

# 결과 확인
print(str(results_df))
print(head(results_df))
print(summary_by_scenario)

#2. 표본 효과크기(표준화된 평균차 까지 저장)
process_resultd <- function(temp_result) {
  list(
    BF_jzs = extract_bayes_factor(temp_result$bayes_factor),
    BF_gica = safe_extract(temp_result$gica$bf10),
    std_mean_diff = safe_extract(temp_result$gica$d),  # 표준화된 평균차 추가
    rho = safe_extract(temp_result$rho),
    sdr = safe_extract(temp_result$sdr),
    delta = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario)
  )
}
# 모든 결과 파일 읽고 처리
filesd <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
results <- map(files, ~{
  temp_result <- readRDS(.x)
  process_resultd(temp_result)
})

# 결과를 데이터 프레임으로 변환
results_df_d <- bind_rows(results)

# 행 이름 제거
rownames(results_df_d) <- NULL

# 결과 저장
saveRDS(results_df_d, file = file.path(home_dir, "final_merged_resultsK_d.RDS"))

##3. 표준편차까지 계산
process_resultdsd <- function(temp_result) {
  # 각 그룹의 표준편차 계산
  data <- temp_result$bayes_factor@data
  sd_x <- sd(data$y[data$group == "x"])
  sd_y <- sd(data$y[data$group == "y"])
  
  list(
    BF_jzs = extract_bayes_factor(temp_result$bayes_factor),
    BF_gica = safe_extract(temp_result$gica$bf10),
    student_p = temp_result$student_p,
    welch_p = temp_result$welch_p,
    mean_diff = safe_extract(temp_result$gica$d),
    rho = safe_extract(temp_result$rho),
    sdr = safe_extract(temp_result$sdr),
    delta = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario),
    sd_x = sd_x,
    sd_y = sd_y
  )
}
# 모든 결과 파일 읽고 처리
filesd <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
results <- map(filesd, ~{
  temp_result <- readRDS(.x)
  process_resultdsd(temp_result)
})
# 결과를 데이터 프레임으로 변환
results_df_dsd <- bind_rows(results)

# 행 이름 제거
rownames(results_df_dsd) <- NULL

# 결과 저장
saveRDS(results_df_dsd, file = file.path(home_dir, "final_merged_results30.RDS"))
