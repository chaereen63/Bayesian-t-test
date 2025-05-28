library(dplyr)
library(tidyr)
library(effectsize)  # effectsize 패키지 로드

## 효과크기별 결과 파일 로드
#표본크기 별로 바꿔보기
results_50es8 <- readRDS("./New/mergedFin50ES8_r1.RDS")  # 효과크기 0.8
results_50es5 <- readRDS("./New/mergedFin50ES5_r1.RDS")  # 효과크기 0.5
results_50es2 <- readRDS("./New/mergedFin50ES2_r1.RDS")  # 효과크기 0.2
results_50es0 <- readRDS("./New/mergedFin50ES0_r1.RDS")  # 효과크기 0.0
results_100es8 <- readRDS("./New/mergedFin100ES8_r1.RDS")  # 효과크기 0.8
results_100es5 <- readRDS("./New/mergedFin100ES5_r1.RDS")  # 효과크기 0.5
results_100es2 <- readRDS("./New/mergedFin100ES2_r1.RDS")  # 효과크기 0.2
results_100es0 <- readRDS("./New/mergedFin100ES0_r1.RDS")  # 효과크기 0.0
results_200es8 <- readRDS("./New/mergedFin200ES8_r1.RDS")  # 효과크기 0.8
results_200es5 <- readRDS("./New/mergedFin200ES5_r1.RDS")  # 효과크기 0.5
results_200es2 <- readRDS("./New/mergedFin200ES2_r1.RDS")  # 효과크기 0.2
results_200es0 <- readRDS("./New/mergedFin200ES0_r1.RDS")  # 효과크기 0.0

# 1. 데이터 처리 함수 - wide format으로 변환
process_bf_data <- function(results_df, effect_size_label, sample_size_label) {
  results_wide <- results_df %>%
    # 로그 변환 적용
    mutate(
      JZS = log10(BF_jzs),
      BFGC = log10(BF_gica)
    ) %>%
    # 필요한 열만 선택
    select(scenario, JZS, BFGC) %>%
    # 시나리오를 factor로 변환
    mutate(
      scenario = factor(as.character(scenario), levels = c("1", "2", "3", "4", "5"),
                        labels = c("A", "B", "C", "D", "E")),
      # 효과크기 정보 추가
      effect_size = effect_size_label,
      # 표본크기 정보 추가
      sample_size = sample_size_label
    )
  
  return(results_wide)
}

# 2. 데이터 처리 및 변환
wide_50data0 <- process_bf_data(results_50es0, "0.0", "50")
wide_50data2 <- process_bf_data(results_50es2, "0.2", "50")
wide_50data5 <- process_bf_data(results_50es5, "0.5", "50")
wide_50data8 <- process_bf_data(results_50es8, "0.8", "50")
wide_100data0 <- process_bf_data(results_100es0, "0.0", "100")
wide_100data2 <- process_bf_data(results_100es2, "0.2", "100")
wide_100data5 <- process_bf_data(results_100es5, "0.5", "100")
wide_100data8 <- process_bf_data(results_100es8, "0.8", "100")
wide_200data0 <- process_bf_data(results_200es0, "0.0", "200")
wide_200data2 <- process_bf_data(results_200es2, "0.2", "200")
wide_200data5 <- process_bf_data(results_200es5, "0.5", "200")
wide_200data8 <- process_bf_data(results_200es8, "0.8", "200")

# 3. 모든 데이터셋 통합
combined_data_wide50 <- bind_rows(wide_50data0, wide_50data2, wide_50data5, wide_50data8) %>%
  mutate(
    effect_size = factor(effect_size, levels = c("0.0", "0.2", "0.5", "0.8")))

combined_data_wide100 <- bind_rows(wide_100data0, wide_100data2, wide_100data5, wide_100data8) %>%
  mutate(
    effect_size = factor(effect_size, levels = c("0.0", "0.2", "0.5", "0.8")))

combined_data_wide200 <- bind_rows(wide_200data0, wide_200data2, wide_200data5, wide_200data8) %>%
  mutate(
    effect_size = factor(effect_size, levels = c("0.0", "0.2", "0.5", "0.8"))
  )

write.csv(combined_data_wide50, file = "./New/analysisData_wide50.csv", row.names = FALSE)
write.csv(combined_data_wide100, file = "./New/analysisData_wide100.csv", row.names = FALSE)
write.csv(combined_data_wide200, file = "./New/analysisData_wide200.csv", row.names = FALSE)
