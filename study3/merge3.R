library(dplyr)
library(purrr)

home_dir <- "./study3"
output_dir <- "D:/results_S3" #저장경로 수정하기 (simulation3.R의 output_dir과 동일하게)

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

# RoBTT의 Effect inclusion BF만 추출하는 함수
extract_robtt_effect <- function(robtt_result) {
  tryCatch({
    if (is.null(robtt_result) || is.null(robtt_result$fit_summary)) {
      return(NA)
    }

    components <- robtt_result$fit_summary$components
    if (is.data.frame(components) && "inclusion_BF" %in% colnames(components)) {
      effect_bf <- components["Effect", "inclusion_BF"]
      return(as.numeric(effect_bf))
    }

    return(NA)
  }, error = function(e) NA)
}

# 단일 결과를 처리하는 함수 (RoBTT, Bain, JZS, BFGC + 표준편차)
process_result <- function(temp_result) {
  # 각 그룹의 표준편차 계산
  data <- temp_result$bayes_factor@data
  sd_x <- sd(data$y[data$group == "x"])
  sd_y <- sd(data$y[data$group == "y"])

  list(
    # RoBTT: Effect inclusion BF (model-averaged) + homo/hetero BF
    BF_robtt_effect = extract_robtt_effect(temp_result$robtt),
    BF_robtt_homo = safe_extract(temp_result$robtt$bf_effect$homoBF),
    BF_robtt_hete = safe_extract(temp_result$robtt$bf_effect$heteBF),

    # Bain: "x = y" vs complement -> BF01 방향 (BF10으로 비교 시 1/x)
    BF_bain_student = safe_extract(temp_result$bain_student$fit$BF.c[1]),
    BF_bain_welch = safe_extract(temp_result$bain_welch$fit$BF.c[1]),

    # BayesFactor의 BF10 추출
    BF_jzs = extract_bayes_factor(temp_result$bayes_factor),

    # BFGC의 BF10 추출
    BF_bfgc = safe_extract(temp_result$bfgc$bf10),

    student_p = temp_result$student_p,
    welch_p = temp_result$welch_p,
    mean_diff = safe_extract(temp_result$bfgc$d),

    # 시나리오 정보 추출
    sdr = safe_extract(temp_result$sdr),
    delta = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario),
    condition = safe_extract(temp_result$condition),
    n1 = safe_extract(temp_result$n1),
    n2 = safe_extract(temp_result$n2),
    replication = safe_extract(temp_result$replication),
    sd_x = sd_x,
    sd_y = sd_y
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
saveRDS(results_df, file = file.path(home_dir, "final_merged_resultsS3.RDS"))

# 시나리오별 요약 통계 계산
summary_by_scenario <- results_df %>%
  group_by(scenario) %>%
  summarise(
    mean_robtt_effect = mean(BF_robtt_effect, na.rm = TRUE),
    mean_robtt_homo = mean(BF_robtt_homo, na.rm = TRUE),
    mean_robtt_hete = mean(BF_robtt_hete, na.rm = TRUE),
    mean_bain_student = mean(BF_bain_student, na.rm = TRUE),
    mean_bain_welch = mean(BF_bain_welch, na.rm = TRUE),
    mean_jzs = mean(BF_jzs, na.rm = TRUE),
    mean_bfgc = mean(BF_bfgc, na.rm = TRUE)
  )

# 결과 확인
print(str(results_df))
print(head(results_df))
print(summary_by_scenario)

# 에러 파일 확인
error_files <- list.files(output_dir, pattern = "^error_\\d+\\.RDS$", full.names = TRUE)
print(paste("Error files:", length(error_files)))
