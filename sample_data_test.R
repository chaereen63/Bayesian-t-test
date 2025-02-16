# set.seed(986)
# Function to run all Bayesian analyses for given parameters
run_bayes_analyses <- function(mean1, mean2, var1, var2, n1, n2) {
  library(BayesFactor)
  library(rjags)
  library(coda)
  library(RoBTT)
  library(bain)
  
  # Calculate standard deviations
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  
  # Calculate Welch's t-test statistics
  t_stat <- (mean1 - mean2) / sqrt(sd1^2/n1 + sd2^2/n2)
  df <- (sd1^2/n1 + sd2^2/n2)^2 / 
    ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
  
  # 1. JZS (already BF10)
  result_jzs <- ttest.tstat(t = t_stat, n1 = n1, n2 = n2, rscale = 1/sqrt(2))
  bf_jzs <- result_jzs$bf[1]
  
  # 2. BMA using RoBTT (already BF10)
  result_bma <- RoBTT(mean1 = mean1, mean2 = mean2,
                      N1 = n1, N2 = n2, sd1 = sd1, sd2 = sd2, 
                      prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
                      prior_rho  = prior("beta",   list(1.5, 1.5)),
                      prior_nu = prior_none(),
                      prior_delta_null = prior("spike", list(0)),
                      prior_rho_null = prior("spike", list(0.5)),
                      chains = 2, warmup = 2000, iter = 6000,
                      parallel = FALSE, thin = 1)
  
  bf_bma <- result_bma$RoBTT$inference$Effect$BF
  
  # 3. Bain (convert from BF01 to BF10)
  t_obj <- list()
  t_obj$estimate <- c(mean1, mean2)
  names(t_obj$estimate) <- c("group1", "group2")
  t_obj$n <- c(n1, n2)
  t_obj$v <- c(var1, var2)
  t_obj$method <- "Two Sample t-test"
  t_obj$Sigma <- list(matrix(var1/n1), matrix(var2/n2))
  class(t_obj) <- "t_test"
  
  result_bain <- bain(t_obj, hypothesis = "group1 = group2")
  bf_bain <- 1 / result_bain$fit$BF[1]  # Convert to BF10
  
  # 4. Wetzels (convert from BF01 to BF10)
  result_wetzels <- Wetzels(mean1 = mean1,
                            mean2 = mean2,
                            sd1 = sd1,
                            sd2 = sd2,
                            n1 = n1,
                            n2 = n2,
                            iters = 1000000,
                            burns = 800001,
                            chains = 1,
                            prior = 'cauchy')
  
  bf_wetzels <- 1 / result_wetzels$meanBF[1]  # Convert to BF10
  
  # Return results
  return(data.frame(
    mean_diff = mean2 - mean1,
    t_stat = t_stat,
    df = df,
    bf_jzs = bf_jzs,        # BF10
    bf_bma = bf_bma,        # BF10
    bf_bain = bf_bain,      # converted to BF10
    bf_wetzels = bf_wetzels # converted to BF10
  ))
}

# Function to run analyses for all cases
analyze_all_cases <- function(cases) {
  results <- do.call(rbind, lapply(1:nrow(cases), function(i) {
    cat(sprintf("\nAnalyzing case %d of %d (mean diff = %.2f)...\n", 
                i, nrow(cases), cases$mean_diff[i]))
    
    run_bayes_analyses(
      mean1 = cases$mean1[i],
      mean2 = cases$mean2[i],
      var1 = cases$var1[i],
      var2 = cases$var2[i],
      n1 = cases$n1[i],
      n2 = cases$n2[i]
    )
  }))
  
  # Add row names
  rownames(results) <- paste0("case_", 1:nrow(results))
  
  return(results)
}

# Run the analysis for your cases
cases <- data.frame(
  mean_diff = c(0.00, 2.20, 4.22, 5.00, 10.0),
  n1 = rep(20, 5),
  n2 = rep(12, 5),
  var1 = rep(12, 5),
  var2 = rep(40, 5)
)

# Calculate means
cases <- within(cases, {
  mean1 <- 50  # baseline
  mean2 <- mean1 + mean_diff
  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
})

# Run analyses
results <- analyze_all_cases(cases)

# Print results in a nice format
print(round(results, 4))

# Optional: Save results to CSV
write.csv(results, "bayes_analysis_results.csv")