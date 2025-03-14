library(dplyr)
library(purrr)

home_dir <- "./New"
output_dir <- "D:/resultsm2_30" #저장경로 수정하기

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



## 표준편차까지 계산
process_resultdsd <- function(temp_result) {
  # 각 그룹의 표준편차 계산
  data <- temp_result$bayes_factor@data
  sd_x <- sd(data$y[data$group == "x"])
  sd_y <- sd(data$y[data$group == "y"])
  
  list(
    scenario = safe_extract(temp_result$scenario),
    BF_jzs = extract_bayes_factor(temp_result$bayes_factor),
    BF_gica = safe_extract(temp_result$gica$bf10),
    mean_diff = safe_extract(temp_result$gica$d),
    sdr = safe_extract(temp_result$sdr),
    n1 = safe_extract(temp_result$n1),
    n2 = safe_extract(temp_result$n2),
    sd_x = sd_x,
    sd_y = sd_y,
    welch_p = safe_extract(temp_result$welch_p),
    student_p = safe_extract(temp_result$student_p)
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
saveRDS(results_df_dsd, file = file.path(home_dir, "merged_resultsm2_30.RDS"))
