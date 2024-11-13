library(bain)

# 파라미터 설정
n1 <- 20
n2 <- 12
mean1 <- 50
mean2 <- 55
var1 <- 12
var2 <- 40

# 방법 1의 t_test 객체 생성
set.seed(986)
group1 <- rnorm(n1, mean1, sqrt(var1))
group2 <- rnorm(n2, mean2, sqrt(var2))

# 생성된 데이터의 실제 평균과 분산 계산
real_mean1 <- mean(group1)
real_mean2 <- mean(group2)
real_var1 <- var(group1)
real_var2 <- var(group2)

cat("\n=== 생성된 데이터의 실제 통계량 ===\n")
cat("Group 1 - 평균:", real_mean1, "분산:", real_var1, "\n")
cat("Group 2 - 평균:", real_mean2, "분산:", real_var2, "\n")

data <- data.frame(
  value = c(group1, group2),
  group = factor(c(rep(1, n1), rep(2, n2)))
)

tresult <- t_test(value ~ group, data = data, paired = FALSE, var.equal = FALSE)

# t_test의 내부 구조 자세히 출력
cat("\n=== 방법 1: bain이 사용하는 t_test 객체의 Sigma 값 ===\n")
str(tresult$Sigma)
print(tresult$Sigma)

# 방법 2: 실제 생성된 데이터의 평균과 분산을 사용
create_t_test_object <- function(mean1, mean2, var1, var2, n1, n2) {
  t_obj <- list()
  
  t_obj$estimate <- c(mean1, mean2)
  names(t_obj$estimate) <- c("group1", "group2")
  
  t_obj$n <- c(n1, n2)
  t_obj$v <- c(var1, var2)
  
  t_obj$method <- "Two Sample t-test"
  
  t_obj$Sigma <- list(matrix(var1/n1), matrix(var2/n2))
  
  class(t_obj) <- "t_test"
  
  return(t_obj)
}

t_obj <- create_t_test_object(real_mean1, real_mean2, real_var1, real_var2, n1, n2)

cat("\n=== 방법 2: 직접 만든 t_test 객체의 Sigma 값 ===\n")
str(t_obj$Sigma)
print(t_obj$Sigma)

# bain 분석 실행 및 결과 비교
resultb <- bain(tresult, hypothesis = "group1 = group2")
result_direct <- bain(t_obj, hypothesis = "group1 = group2")

cat("\n=== Bain 결과 비교 ===\n")
cat("\n방법 1 결과:\n")
print(resultb)

cat("\n방법 2 결과:\n")
print(result_direct)

# t 통계량과 자유도 계산하여 비교
t_stat1 <- (tresult$estimate[1] - tresult$estimate[2]) / 
  sqrt(tresult$v[1]/n1 + tresult$v[2]/n2)
df1 <- (tresult$v[1]/n1 + tresult$v[2]/n2)^2 / 
  ((tresult$v[1]/n1)^2/(n1-1) + (tresult$v[2]/n2)^2/(n2-1))

t_stat2 <- (t_obj$estimate[1] - t_obj$estimate[2]) / 
  sqrt(t_obj$v[1]/n1 + t_obj$v[2]/n2)
df2 <- (t_obj$v[1]/n1 + t_obj$v[2]/n2)^2 / 
  ((t_obj$v[1]/n1)^2/(n1-1) + (t_obj$v[2]/n2)^2/(n2-1))

cat("\n=== t 통계량 비교 ===\n")
cat("방법 1 t-statistic:", t_stat1, "df:", df1, "\n")
cat("방법 2 t-statistic:", t_stat2, "df:", df2, "\n")
