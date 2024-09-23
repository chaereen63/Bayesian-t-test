### a simulation evaluation script that for a distributed computational cluster ###
# set the computational node
.libPaths(c(""))
home_dir   <- ""
output_dir <- ""
start_dir  <- getwd()

library("rstan")
library("RoBTT")

# load the functions and settings
source(file = file.path(home_dir, "functions.R"))
settings <- readRDS(file = file.path(home_dir, "settings.RDS"))
tracker  <- "sim_loop"
max_time <- 23.5

# collect the simulation data
time_start <- Sys.time()
while(difftime(Sys.time(), time_start, units = "hours") < max_time){

  # get the current setting to compute
  loop <- get_loop(home_dir, tracker)
  if(loop > nrow(settings))
    stop("DONE!!!")
  
  # simulate data
  set.seed(settings$seed[loop])
  data <- simulate_data(
    mean1 = settings$mean1[loop],
    mean2 = settings$mean2[loop],
    sd1   = settings$sd1[loop],
    sd2   = settings$sd2[loop],
    df1   = settings$df1[loop],
    df2   = settings$df2[loop],
    n1    = settings$n1[loop],
    n2    = settings$n2[loop]
  )
  
  # fit the model
  fit <- RoBTT(
    x1 = data$x1,
    x2 = data$x2,
    prior_d  = prior("cauchy", list(0, 1/sqrt(2))),
    prior_r  = prior("beta",   list(1.5, 1.5)),
    prior_nu = prior("exp",    list(1)),
    likelihood = c("normal", "t"),
    chains = 2, warmup = 1000, iter = 5000,
    parallel = FALSE, seed = settings$seed[loop]
  )
  
  # collect the results
  saveRDS(summary(fit, conditional = TRUE),  file = file.path(output_dir, "ensemble",     paste0(settings$seed[loop], ".RDS")))
  saveRDS(summary(fit, diagnostics = TRUE),  file = file.path(output_dir, "diagnostics",  paste0(settings$seed[loop], ".RDS")))
  saveRDS(summary(fit, type = "individual"), file = file.path(output_dir, "models",       paste0(settings$seed[loop], ".RDS")))
}



