# 필요한 패키지 설치 및 로드
install_required_packages <- function() {
  if (!require(rstan)) install.packages("rstan")
  if (!require(bridgesampling)) install.packages("bridgesampling")
  if (!require(ggplot2)) install.packages("ggplot2")
  
  library(rstan)
  library(bridgesampling)
  library(ggplot2)
  
  # Stan 설정
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

# Stan 모델 생성
create_stan_model <- function() {
  stan_code <- "
  functions {
    real cauchy_lpdf(real x, real mu, real sigma) {
      return -(log(pi()) + log(sigma) + log1p(pow((x - mu) / sigma, 2)));
    }
  }

  data {
    int<lower=0> n;
    vector[n] X;
    int<lower=0,upper=1> is_cauchy;
  }

  parameters {
    real delta;
    real<lower=0> sigma;
  }

  transformed parameters {
    real mu;
    mu = delta * sigma;
  }

  model {
    if (is_cauchy) {
      target += cauchy_lpdf(delta | 0, 1);
      target += cauchy_lpdf(sigma | 0, 1);
    } else {
      delta ~ normal(0, 1);
      sigma ~ normal(0, 1);
    }
    
    X ~ normal(mu, sigma);
  }

  generated quantities {
    vector[n] log_lik;
    real prior_delta;
    real post_delta;
    
    for (i in 1:n) {
      log_lik[i] = normal_lpdf(X[i] | mu, sigma);
    }
    
    if (is_cauchy) {
      prior_delta = cauchy_rng(0, 1);
    } else {
      prior_delta = normal_rng(0, 1);
    }
    post_delta = delta;
  }
  "
  
  writeLines(stan_code, "single_sample_ttest.stan")
  stan_model("single_sample_ttest.stan")
}

# 메인 함수
SD_stan <- function(
    group1,
    group2 = 0,
    iters = 10000,
    warmup = 5000,
    chains = 4,
    thin = 1,
    sample = 1,
    prior = 'cauchy',
    plot = TRUE,
    static_model = NULL
) {
  if (length(group2) == 1) {
    X <- group1
  } else {
    X <- group1 - group2
  }
  X <- scale(X)
  
  stan_data <- list(
    n = length(X),
    X = as.vector(X),
    is_cauchy = as.integer(prior == 'cauchy')
  )
  
  if (is.null(static_model)) {
    static_model <- create_stan_model()
  }
  
  fit <- sampling(static_model, 
                  data = stan_data,
                  iter = iters,
                  warmup = warmup,
                  chains = chains,
                  thin = thin)
  
  posterior <- as.data.frame(rstan::extract(fit))
  
  cat('\nResults of a single sample Bayesian t-test.\n')
  cat('Number of observations:', length(X), '\n')
  cat('Rhat values:', max(summary(fit)$summary[, "Rhat"]), '\n')
  
  if (plot) {
    posterior_plot <- ggplot(posterior, aes(x = delta)) +
      geom_density(fill = "lightblue", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(title = "Posterior Distribution of Effect Size",
           x = "Effect Size (delta)",
           y = "Density")
    print(posterior_plot)
  }
  
  result <- list(
    posterior = posterior,
    summary = summary(fit)$summary,
    model = fit,
    data = list(X = X)
  )
  
  return(result)
}