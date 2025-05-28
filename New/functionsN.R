library(stats)
library(dplyr)

simulate_data <- function(n1, n2, mean1, mean2, sd1, sd2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # 데이터 생성
  x1 <- rnorm(n1, mean = mean1, sd = sd1)
  x2 <- rnorm(n2, mean = mean2, sd = sd2)
  
  return(list(
    x1 = x1,
    x2 = x2,
    mean1 = mean1,
    mean2 = mean2,
    sd1 = sd1,
    sd2 = sd2
  ))
}

mean_sd <- function(sd1, sd2) { #  모수 단위: 이분산 가정- 평균표준편차로 수정
  sqrt((sd1^2 + sd2^2) / 2)
}

pooled_sd <- function(sd1, sd2, n1, n2) { # Student
  sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
}

BeFi_sd <- function(sd1, sd2, n1, n2) { # Welch
  sqrt(sd1^2/n1 + sd2^2/n2)
}

cohen_d <- function(mean1, mean2, sd1, sd2, n1, n2) {
  sd_p <- pooled_sd(sd1, sd2, n1, n2)
  d = (mean1 - mean2) / sd_p
  
  return(d)
}

SES <- function(mean1, mean2, sd1, sd2) {
  (mean1 - mean2) / mean_sd(sd1, sd2)
}

calculate_sdr <- function(sd1, sd2) {
  return(sd2 / sd1)
}

calculate_rho <- function(var1, var2) {
  return(1 / var1 / (1 / var1 + 1 / var2))
}

# Helper functions for simulation tracking and visualization

get_loop <- function(dir, file, max_loop) {
  tracker_file <- file.path(dir, paste0(file, ".txt"))
  
  if (file.exists(tracker_file)) {
    loop <- as.numeric(readLines(tracker_file))
  } else {
    loop <- 1
  }
  
  # Ensure loop doesn't exceed max_loop
  loop <- min(loop, max_loop)
  
  # Increment and save the loop counter
  next_loop <- min(loop + 1, max_loop)
  writeLines(as.character(next_loop), tracker_file)
  
  return(loop)
}

update_tracker <- function(home_dir, tracker) {
  # This function can remain empty
}

empty_plot <- function() {
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', 
       xlim = c(0,1), ylim = c(0,1), mar = c(0,0,0,0))
}

# Giron and Castillo (2021)
library(pracma)

BeFiBF <- function(x1, x2) {
  
  # Calculate sample statistics
  xbar1 <- mean(x1)
  xbar2 <- mean(x2)
  var1 <- var(x1)
  var2 <- var(x2)
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Calculate degrees of freedom
  df1 <- n1 - 1
  df2 <- n2 - 1
  
  # Location parameter
  d <- xbar2 - xbar1
  
  # Define the numerator integrand
  num_integrand <- function(delta) {
    term1 <- (1 + n1 * delta^2/(df1 * var1))^(-n1/2)
    term2 <- (1 + n2 * ((delta-d)^2)/(df2 * var2))^(-n2/2)
    return(term1 * term2)
  }
  
  # Define the denominator integrand
  den_integrand <- function(delta) {
    term1 <- (1 + (n1/(n1+1)) * delta^2/(df1 * var1))^(-n1/2)
    term2 <- (1 + (n2/(n2+1)) * ((delta-d)^2)/(df2 * var2))^(-n2/2)
    return(term1 * term2)
  }
  
  # Calculate reasonable bounds for the integral
  bound <- max(abs(d) + 10 * sqrt(max(var1, var2)), 100)
  
  # Compute the integrals
  num_integral <- integral(num_integrand, -bound, bound)
  den_integral <- integral(den_integrand, -bound, bound)
  
  # Calculate the product term
  prod_term <- sqrt(n1 + 1) * sqrt(n2 + 1)
  
  # Calculate final Bayes factor
  bf01 <- prod_term * num_integral/den_integral
  
  return(list(
    bf01 = bf01,
    bf10 = 1/bf01,
    d = d,
    num_integral = num_integral,
    den_integral = den_integral,
    prod_term = prod_term
  ))
}
