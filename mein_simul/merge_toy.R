library(dplyr)
library(purrr)

## practice for practice 결과 merging ##

home_dir   <- "./mein_simul"
output_dir <- "./mein_simul/mein_results"

# 안전한 접근 함수 정의
safe_extract <- function(x, default = NA) {
  tryCatch(x, error = function(e) default)
}

# BayesFactor 객체에서 값을 안전하게 추출하는 함수
extract_bayes_factor <- function(bf_object) {
  tryCatch({
    if (inherits(bf_object, "BFBayesFactor")) {
      as.numeric(exp(bf_object@bayesFactor[, "bf"]))
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

# 단일 결과를 처리하는 함수
process_result <- function(temp_result) {
  list(
    BF_robtt_effect = extract_robtt_effect(temp_result$robtt),
    BF_bain_student = safe_extract(temp_result$bain_student$fit$BF.c[1]),
    BF_bain_welch = safe_extract(temp_result$bain_welch$fit$BF.c[1]),
    BF_bayesfactor = extract_bayes_factor(temp_result$bayes_factor),
    true_model = safe_extract(temp_result$true_model),
    rho = safe_extract(temp_result$rho),
    sdr = safe_extract(temp_result$sdr),
    delta = safe_extract(temp_result$delta),
    scenario = safe_extract(temp_result$scenario)
  )
}

# merge data from the simulations
files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)
results <- map(files, ~{
  temp_result <- readRDS(.x)
  process_result(temp_result)
})

# 결과를 데이터 프레임으로 변환
results_df <- bind_rows(results)

rownames(results_df) <- NULL
saveRDS(results_df, file = file.path(home_dir, "merged_results.RDS"))

# 결과 확인
print(str(results_df))
print(head(results_df))
