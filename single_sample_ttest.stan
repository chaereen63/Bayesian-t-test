
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
  
