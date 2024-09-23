source("simulation/functions.R")

# create the settings combinations
settings <- data.frame(expand.grid(
  n     = c(20,   50, 100),
  delta = c(0,   0.2, 0.5),
  rho   = c(1,   1.5, 2.0),
  nu    = c(Inf, 10,  5)
))

# choose parameters to fix
grand_mean <- 0
grand_sd   <- 1

# assign the shared parameters
settings$df1   <- settings$nu
settings$df2   <- settings$nu
settings$n1    <- settings$n
settings$n2    <- settings$n

# compute the population means and sds
settings$sd1   <- sqrt(2 * grand_sd^2 * settings$rho^2   / (settings$rho^2   + 1))
settings$sd2   <- sqrt(2 * grand_sd^2 * 1/settings$rho^2 / (1/settings$rho^2 + 1))

settings$mean1 <- grand_mean + 1/2 * settings$delta * pooled_sd(settings$sd1, settings$sd2, settings$n1, settings$n2) #pooled SD 사용!
settings$mean2 <- grand_mean - 1/2 * settings$delta * pooled_sd(settings$sd1, settings$sd2, settings$n1, settings$n2)

# verify that the delta and rho match
all(abs(settings$sd1 / settings$sd2 - settings$rho) < 1e-10)
all(abs(cohens_d(settings$mean1, settings$mean2, settings$sd1, settings$sd2, settings$n1, settings$n2) - settings$delta) < 1e-10)

# create the settings replications
settings      <- do.call(rbind, lapply(1:1000, function(x) settings))
settings$seed <- 1:nrow(settings)

saveRDS(settings, file = "simulation/settings.RDS")
