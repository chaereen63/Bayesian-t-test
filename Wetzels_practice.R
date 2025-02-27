### Wetzels 데이터 생성하여 비교하기
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

cat("\n=== 생성된 데이터의 실제 통계량 ===\n")
cat("Group 1 - 평균:", real_mean1, "분산:", real_var1, "\n")
cat("Group 2 - 평균:", real_mean2, "분산:", real_var2, "\n")

data <- data.frame(
  value = c(group1, group2),
  group = factor(c(rep(1, n1), rep(2, n2)))
)

## 1.raw data function
jags_modelR <- "
model {
  # Data level - using raw data
  for(i in 1:n1) {
    group1[i] ~ dnorm(muData1, precision1)  
  }
  for(i in 1:n2) {
    group2[i] ~ dnorm(muData2, precision2)
  }
  
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

Wetzels_raw <- function(group1, group2,
                    iters = 1000000,
                    burns = 800001,
                    chains = 1,
                    thins = 1,
                    prior = 'cauchy') {
  
  # Get sample sizes
  n1 <- length(group1)
  n2 <- length(group2)
  
  # Standardize data like WinBUGS version
  group2 <- group2 - mean(group1)  # Center group2 by group1's mean
  group2 <- group2/sd(group1)      # Scale group2 by group1's sd
  group1 <- as.vector(scale(group1))  # Standardize group1
  
  # Prepare data list for JAGS
  data <- list(
    group1 = group1,
    group2 = group2,
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
  model <- jags.model(textConnection(jags_modelR),
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
  med <- mean(group1) - mean(group2)
  
  dcount <- 0
  sa <- min(10000, (iters - burns)/thins)
  ind <- sample(1:((iters - burns)/thins), sa)
  
  for (i in 1:sa) {
    mu. <- as.matrix(samples)[ind[i], 'mu']
    si1. <- as.matrix(samples)[ind[i], 'sigma1']
    si2. <- as.matrix(samples)[ind[i], 'sigma2']
    al. <- as.matrix(samples)[ind[i], 'alpha']
    dd1 <- rnorm(n1, mu. + al./2, si1.)
    dd2 <- rnorm(n2, mu. - al./2, si2.)
    delta <- (mean(dd1) - mean(dd2))/sqrt(((sd(dd1))^2*(n1-1) + (sd(dd2))^2*(n2-1))/(n1+n2-2))
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

# Wetzels 함수 실행 (chains=3으로 하여 수렴 진단)
result_raw <- Wetzels_raw(group1, group2, 
                  iters = 1000000,
                  burns = 800001,
                  chains = 1,
                  thins = 1,
                  prior = 'cauchy')

## 2. summary statistic

# 생성된 데이터의 실제 평균과 분산 계산
real_mean1 <- mean(group1)
real_mean2 <- mean(group2)
real_var1 <- var(group1)
real_var2 <- var(group2)
real_sd1 <- sqrt(real_var1)
real_sd2 <- sqrt(real_var2)


jags_model <- "
model {
  # Data level - using sufficient statistics
  y1_mean ~ dnorm(muData1, precision1)  
  y2_mean ~ dnorm(muData2, precision2) #for likelihood
  
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
  mean1_std <- 0
  sd1_std <- 1
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

result <- Wetzels(mean1 = real_mean1,    # 그룹1 평균
                  mean2 = real_mean2,     # 그룹2 평균
                  sd1 = real_sd1,        # 그룹1 표준편차
                  sd2 = real_sd2,      # 그룹2 표준편차
                  n1 = n1,        # 그룹1 표본크기
                  n2 = n2,        # 그룹2 표본크기
                  iters = 1000000,
                  burns = 800001,
                  chains = 1,
                  prior = 'cauchy')
result_raw$meanBF
result$meanBF
