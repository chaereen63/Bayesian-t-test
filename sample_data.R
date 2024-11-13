## Box & Tiao sample table ##
library(BayesFactor)  # JZS 
library(rjags)        # MCMC
library(RoBTT)        # BMA using RoBTT
library(bain)         # bain

# 연습 (pilot의 조건3의 구체적 예에 해당, 작은 표본이 더 큰 분산)
mean1 = 50
var1 = 12
sd1 = sqrt(12)
n1 = 20
mean2 = 55
var2 = 40
sd2 = sqrt(40)
n2 = 12

# 0. Welch's t-test
t_stat <- (mean1 - mean2) / sqrt(sd1^2/n1 + sd2^2/n2)
df <- (sd1^2/n1 + sd2^2/n2)^2 / 
  ((sd1^2/n1)^2/(n1-1) + (sd2^2/n2)^2/(n2-1))
p_value <- 2 * pt(-abs(t_stat), df)
print(c(t_stat, p_value))

# 1. JZS
bf <- ttest.tstat(t = t_stat, n1 = n1, n2 = n2, rscale = 1/sqrt(2))
bf$bf

# 2. BMA
resultM <- RoBTT(mean1 = mean1, mean2 = mean2,
                 N1 = n1, N2 = n2, sd1 = sd1, sd2 = sd2, 
                 prior_delta = prior("cauchy", list(0, 1/sqrt(2))),
                 prior_rho  = prior("beta",   list(1.5, 1.5)),
                 prior_nu = prior_none(),
                 prior_delta_null = prior("spike", list(0)),
                 prior_rho_null = prior("spike", list(0.5)),
                 chains = 2, warmup = 2000, iter = 6000,
                 parallel = FALSE, thin = 1)

resultM$RoBTT$inference$Effect$BF

# 3. Bain

# t_test 객체 생성
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
t_obj <- create_t_test_object(mean1, mean2, var1, var2, n1, n2)

# bain 분석 실행
resultb <- bain(t_obj, hypothesis = "group1 = group2")

# 결과 출력
print(resultb)
resultb$fit$BF[1]
# t 통계량과 자유도도 확인
t_stat <- (mean1 - mean2) / sqrt(var1/n1 + var2/n2)
df <- (var1/n1 + var2/n2)^2 / ((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1))

cat("\n=== t 통계량 ===\n")
cat("t-statistic:", t_stat, "df:", df, "\n")

# 4. Wetzels
# Modified JAGS model to match WinBUGS implementation
# Modified JAGS model for summary statistics
jags_model <- "
model {
  # Data level - using sufficient statistics
  y1_mean ~ dnorm(muData1, precision1/n1)  # Sample mean distribution
  y2_mean ~ dnorm(muData2, precision2/n2)  # Sample mean distribution
  
  # Parameter priors - matching WinBUGS Cauchy implementation
  delta ~ dnorm(0, lambdaDelta)
  lambdaDelta ~ dchisqr(1)
  
  # Sigma priors using same approach as WinBUGS
  sigma.1 ~ dnorm(0, sigmaL1)
  sigmaL1 ~ dchisqr(1)
  sigma1 <- abs(sigma.1)
  
  sigma.2 ~ dnorm(0, sigmaL2)
  sigmaL2 ~ dchisqr(1)
  sigma2 <- abs(sigma.2)
  
  # Variance calculations
  var1 <- pow(sigma1, 2)
  var2 <- pow(sigma2, 2)
  precision1 <- 1/var1
  precision2 <- 1/var2
  
  # Mean prior matching WinBUGS
  mu ~ dnorm(0, muL)
  muL ~ dchisqr(1)
  
  # Effect size calculation matching WinBUGS pooled variance approach
  pooled.var <- (var1*(n1-1) + var2*(n2-1))/(n1+n2-2)
  alpha <- delta * sqrt(pooled.var)
  muData1 <- mu + alpha*0.5
  muData2 <- mu - alpha*0.5
}
"

Wetzels <- function(mean1, mean2, sd1, sd2, n1, n2,
                     iters = 1000000,
                     burns = 800001,
                     chains = 1,
                     thins = 1,
                     prior = 'cauchy') {
  
  # Standardize summary statistics like WinBUGS version
  pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
  mean2_std <- (mean2 - mean1)/sd1
  mean1_std <- (mean1 - mean1)/sd1  # Will be 0
  sd1_std <- sd1/sd1  # Will be 1
  sd2_std <- sd2/sd1
  
  # Prepare data list for JAGS
  data <- list(
    y1_mean = mean1_std,
    y2_mean = mean2_std,
    n1 = n1,
    n2 = n2
  )
  
  # Initial values function matching WinBUGS approach
  inits <- function() {
    list(
      delta = rnorm(1, 0, 1),
      mu = rnorm(1, 0, 1),
      sigma.1 = runif(1, 0, 5),
      sigma.2 = runif(1, 0, 5)
    )
  }
  
  # Set up JAGS model
  model <- jags.model(textConnection(jags_model),
                      data = data,
                      inits = inits,
                      n.chains = chains)
  
  # Burn-in
  update(model, n.iter = burns)
  
  # Sample from posterior
  samples <- coda.samples(model,
                          variable.names = c('delta', 'mu', 'sigma1', 'sigma2', 'alpha'),
                          n.iter = (iters - burns),
                          thin = thins)
  
  # Calculate Bayes Factors like WinBUGS version
  dat1 <- matrix(nrow = chains, ncol = 3)
  rownames(dat1) <- paste('chain', 1:chains)
  colnames(dat1) <- c('no OR', 'OR1', 'OR2')
  
  # Prior densities
  if (prior == 'cauchy') {
    d0_prior1 <- dcauchy(0)
    d0_prior2 <- dcauchy(0) * 2
  } else {
    d0_prior1 <- dnorm(0, 0, 1)
    d0_prior2 <- dnorm(0, 0, 1) * 2
  }
  
  # Calculate BFs for each chain
  for (i in 1:chains) {
    post.dat <- as.matrix(samples[[i]])[, 'delta']
    m <- mean(post.dat)
    sd <- sd(post.dat)
    d0_post <- dnorm(0, m, sd)
    
    dat1[i, 1] <- d0_post/d0_prior1
    d0_post1 <- d0_post/pnorm(0, m, sd)
    dat1[i, 2] <- d0_post1/d0_prior2
    d0_post2 <- d0_post/(1 - pnorm(0, m, sd))
    dat1[i, 3] <- d0_post2/d0_prior2
  }
  
  # Prepare return object
  dat.SD <- list()
  dat.SD$Post.delta <- as.matrix(samples)[, 'delta']
  dat.SD$Post.mu <- as.matrix(samples)[, 'mu']
  dat.SD$Post.sigma1 <- as.matrix(samples)[, 'sigma1']
  dat.SD$Post.sigma2 <- as.matrix(samples)[, 'sigma2']
  dat.SD$Post.alpha <- as.matrix(samples)[, 'alpha']
  dat.SD$BF <- dat1
  dat.SD$meanBF <- apply(dat1, 2, mean)
  
  # Calculate prep like WinBUGS version
  meanBF <- mean(dat1[, 1])
  PH0D <- meanBF/(meanBF + 1)
  PH1D <- 1 - PH0D
  med <- mean1 - mean2
  
  dcount <- 0
  sa <- min(10000, (iters - burns)/thins)
  ind <- sample(1:((iters - burns)/thins), sa)
  
  for (i in 1:sa) {
    mu. <- as.matrix(samples)[ind[i], 'mu']
    si1. <- as.matrix(samples)[ind[i], 'sigma1']
    si2. <- as.matrix(samples)[ind[i], 'sigma2']
    al. <- as.matrix(samples)[ind[i], 'alpha']
    dd1_mean <- mu. + al./2
    dd2_mean <- mu. - al./2
    delta <- (dd1_mean - dd2_mean)/sqrt((si1.^2/n1 + si2.^2/n2)/2)
    if (delta * med > 0) dcount <- dcount + 1
  }
  
  prep <- PH0D * 0.5 + PH1D * (dcount/sa)
  dat.SD$prep <- prep
  
  # Print results
  cat('Results of a 2 sample SD Bayesian t-test.\n')
  if (chains > 1) {
    cat('The mean Bayes factor, calculated from', chains, 'chains =', apply(dat1, 2, mean), '.\n')
    cat('Bayesian p-rep =', prep, '.\n')
  }
  cat('All Bayes Factors are displayed below:\n\n')
  print(dat1)
  
  return(dat.SD)
}

result <- Wetzels(mean1 = mean1,    # 그룹1 평균
                   mean2 = mean2,     # 그룹2 평균
                   sd1 = sd1,        # 그룹1 표준편차
                   sd2 = sd2,      # 그룹2 표준편차
                   n1 = n1,        # 그룹1 표본크기
                   n2 = n2,        # 그룹2 표본크기
                   iters = 5000,
                   burns = 1000,
                   chains = 4,
                   prior = 'cauchy')
result$meanBF