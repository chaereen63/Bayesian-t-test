library(dplyr)
library(purrr)

home_dir <- "./study2"
output_dir <- "D:/S2_results50_E8" #저장경로 수정하기

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


## 결과 정리
process_result <- function(temp_result) {
  list(
    scenario = safe_extract(temp_result$scenario),
    varr = safe_extract(temp_result$varr),
    n1 = safe_extract(temp_result$n1),
    n2 = safe_extract(temp_result$n2),
    Mean1 = safe_extract(temp_result$Mean1),
    Mean2 = safe_extract(temp_result$Mean2),
    SD1 = safe_extract(temp_result$SD1),
    SD2 = safe_extract(temp_result$SD2),
    student_p = safe_extract(temp_result$student_p),
    welch_p = safe_extract(temp_result$welch_p),
    student_t = safe_extract(temp_result$student_t),
    welch_t = safe_extract(temp_result$welch_t),
    d_s = safe_extract(temp_result$d_s),
    d_w = safe_extract(temp_result$d_w),
    d_c = safe_extract(temp_result$d_c),
    true_d = safe_extract(temp_result$true_d)
  )
}
# 모든 결과 파일 읽고 처리
filesd <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
results <- map(filesd, ~{
  temp_result <- readRDS(.x)
  process_result(temp_result)
})
# 결과를 데이터 프레임으로 변환
results_df <- bind_rows(results)

# 행 이름 제거
rownames(results_df) <- NULL
print(results_df)
# 결과 저장
saveRDS(results_df, file = file.path(home_dir, "S2merged50ES8.RDS"))
