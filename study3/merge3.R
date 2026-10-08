# Study 3: simulation3.R 결과 병합 (JZS, GICA, Bain, RoBTT)
# 모든 BF는 BF10(H1 지지) 방향으로 통일해서 저장
library(dplyr)
library(purrr)

home_dir   <- "./study3"
output_dir <- file.path(home_dir, "results") # simulation3.R의 output_dir과 같게 수정

# 안전한 추출 함수
safe_extract <- function(x, default = NA) {
  tryCatch({
    if (is.null(x) || length(x) == 0) default else as.numeric(x[1])
  }, error = function(e) default)
}

# BayesFactor 객체에서 BF10 추출
extract_bayes_factor <- function(bf_object) {
  tryCatch({
    if (inherits(bf_object, "BFBayesFactor")) {
      as.numeric(exp(bf_object@bayesFactor$bf))
    } else {
      NA
    }
  }, error = function(e) NA)
}

# RoBTT의 Effect inclusion BF 추출 (모델 평균)
extract_robtt_effect <- function(robtt_result) {
  tryCatch({
    components <- robtt_result$fit_summary$components
    if (is.data.frame(components) && "inclusion_BF" %in% colnames(components)) {
      as.numeric(components["Effect", "inclusion_BF"])
    } else {
      NA
    }
  }, error = function(e) NA)
}

# 표본 표준편차 계산 (JZS 객체에 저장된 데이터 사용)
extract_sample_sd <- function(bf_object, group) {
  tryCatch({
    data <- bf_object@data
    sd(data$y[data$group == group])
  }, error = function(e) NA)
}

# 단일 결과를 처리하는 함수
process_result <- function(temp_result) {
  list(
    # RoBTT
    BF_robtt_effect = extract_robtt_effect(temp_result$robtt),
    BF_robtt_homo   = safe_extract(temp_result$robtt$bf_effect$homoBF),
    BF_robtt_hete   = safe_extract(temp_result$robtt$bf_effect$heteBF),
    
    # Bain: BF.u는 H0(x = y) vs 비제약 가설의 BF0u → 역수로 BF10
    BF_bain_student = 1 / safe_extract(temp_result$bain_student$fit$BF.u),
    BF_bain_welch   = 1 / safe_extract(temp_result$bain_welch$fit$BF.u),
    
    # JZS (BayesFactor)
    BF_jzs  = extract_bayes_factor(temp_result$bayes_factor),
    
    # GICA
    BF_gica = safe_extract(temp_result$gica$bf10),
    
    # 빈도주의 p-value
    student_p = safe_extract(temp_result$student_p),
    welch_p   = safe_extract(temp_result$welch_p),
    
    # 표본 통계량
    mean_diff = safe_extract(temp_result$gica$d), # xbar2 - xbar1
    sd_x = extract_sample_sd(temp_result$bayes_factor, "x"),
    sd_y = extract_sample_sd(temp_result$bayes_factor, "y"),
    
    # 시나리오 정보
    rho      = safe_extract(temp_result$rho),
    sdr      = safe_extract(temp_result$sdr),
    delta    = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario)
  )
}

# 모든 결과 파일 읽고 처리
files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
loops <- as.numeric(gsub(".*results_(\\d+)\\.RDS", "\\1", files))
files <- files[order(loops)]
loops <- sort(loops)

results_df <- map2(files, loops, ~{
  temp_result <- readRDS(.x)
  c(loop = .y, process_result(temp_result))
}) %>%
  bind_rows()

rownames(results_df) <- NULL

# 에러 파일 확인
error_files <- list.files(output_dir, pattern = "^error_\\d+\\.RDS$", full.names = TRUE)
print(paste("Result files:", length(files), "/ Error files:", length(error_files)))
if (length(error_files) > 0) {
  error_df <- tibble(
    loop  = as.numeric(gsub(".*error_(\\d+)\\.RDS", "\\1", error_files)),
    error = map_chr(error_files, ~ as.character(readRDS(.x)$error))
  )
  print(error_df)
}

# 결과 저장
saveRDS(results_df, file = file.path(home_dir, "final_merged_results3.RDS"))

# 시나리오 x 효과크기별 요약 (BF는 log 척도의 중앙값)
summary_by_condition <- results_df %>%
  group_by(scenario, delta) %>%
  summarise(
    n_rep = n(),
    across(starts_with("BF_"), ~ median(log(.x), na.rm = TRUE), .names = "mdlog_{.col}"),
    .groups = "drop"
  )

# 결과 확인
str(results_df)
print(head(results_df))
print(summary_by_condition, width = Inf)
