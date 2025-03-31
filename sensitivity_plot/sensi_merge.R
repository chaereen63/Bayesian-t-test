library(dplyr)
library(purrr)

home_dir <- "./sensitivity_plot"
output_dir <- "D:/sens5" #저장경로 수정하기

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
  data <- temp_result$jzs_medium@data
  sd_x <- sd(data$y[data$group == "x"])
  sd_y <- sd(data$y[data$group == "y"])
  
  list(
    scenario = safe_extract(temp_result$scenario),
    BF_jzsM = extract_bayes_factor(temp_result$jzs_medium),
    BF_jzsW = extract_bayes_factor(temp_result$jzs_wide),
    BF_jzsU = extract_bayes_factor(temp_result$jzs_ultrawide),
    BF_gica = safe_extract(temp_result$gica$bf10),
    mean_diff = safe_extract(temp_result$gica$d),
    varr = safe_extract(temp_result$varr),
    n1 = safe_extract(temp_result$n1),
    n2 = safe_extract(temp_result$n2),
    sd_x = sd_x,
    sd_y = sd_y #,
    # welch_p = safe_extract(temp_result$welch_p),
    # student_p = safe_extract(temp_result$student_p)
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
print(results_df_dsd)
# 결과 저장
saveRDS(results_df_dsd, file = file.path(home_dir, "merged_seni5.RDS"))
