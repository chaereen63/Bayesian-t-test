# compare_robtt_structures <- function(n = 5) {
#   files <- list.files(output_dir, pattern = "^results_\\d+\\.RDS$", full.names = TRUE)[1:n]
#   structures <- lapply(files, function(f) {
#     temp_result <- readRDS(f)
#     capture.output(str(temp_result$robtt, max.level = 2))
#   })
  
  # 결과를 텍스트 파일로 저장
#   writeLines(unlist(structures), "robtt_structures_comparison.txt")
# }

# 처음 5개 파일의 RoBTT 구조 비교
# compare_robtt_structures(5)



library(ggplot2)
library(dplyr)
library(tidyr)

# 데이터 변환
  #bain BF가 역수라서 변환
  results_df %>% select("scenario", "true_model", "rho", "sdr", "delta", "BF_robtt_effect", "BF_bain_student","BF_bain_welch", "BF_bayesfactor" ) %>%
            mutate(BF_bain_student = (BF_bain_student)^(-1), BF_bain_welch=(BF_bain_welch)^(-1)) ->
            result_inv
  #시각화용 데이터
  results_re_long <- result_inv %>%
  pivot_longer(cols = starts_with("BF_"),
               names_to = "model",
               values_to = "BF")

#시각화
ggplot(results_re_long, aes(x = rho, y = BF, color = model, shape = factor(delta))) +
  geom_point() +
  facet_wrap(~ scenario) +
  scale_y_log10() + # BF는 종종 큰 값이므로 로그 스케일 사용
  labs(title = "모형의 베이지안 팩터(BF) 비교",
       x = "rho 값",
       y = "베이지안 팩터 (BF)",
       color = "모형",
       shape = "델타") +
  theme_minimal()

#시각화
ggplot(results_re_long, aes(x = sdr, y = BF, color = model, shape = factor(delta))) +
  geom_point() +
  scale_y_log10() + # BF는 종종 큰 값이므로 로그 스케일 사용
  labs(title = "모형의 베이지안 팩터(BF) 비교",
       x = "SDR 값",
       y = "베이지안 팩터 (BF)",
       color = "모형",
       shape = "델타") +
  theme_minimal()

# 시각화
ggplot(results_re_long, aes(x = sdr, y = BF, color = model)) +
  geom_point(aes(shape = factor(delta)), size = 3) +
  geom_smooth(aes(group = interaction(delta, model)), method = "lm", se = FALSE) + #델타와 모델 별로 그룹화
  scale_y_log10() + # BF는 종종 큰 값이므로 로그 스케일 사용
  labs(title = "모형의 베이지안 팩터(BF) 비교",
       x = "SDR 값",
       y = "베이지안 팩터 (BF)",
       color = "모형",
       shape = "델타") +
  theme_minimal()

# 델타가 0인 경우를 제외한 데이터 필터링
filtered_data <- results_re_long %>%
  filter(delta != 0)

# 필터링된 데이터로 시각화
ggplot(filtered_data, aes(x = sdr, y = BF, color = model)) +
  geom_point(aes(shape = factor(delta)), size = 3) +
  geom_smooth(aes(group = interaction(delta, model)), method = "lm", se = FALSE) +
  scale_y_log10() +
  labs(title = "모형의 베이지안 팩터(BF) 비교 (델타 0 제외)",
       x = "SDR 값",
       y = "베이지안 팩터 (BF)",
       color = "모형",
       shape = "델타") +
  theme_minimal()
