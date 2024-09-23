### a data merging script for a distributed computational cluster ###

# set the computational node
.libPaths(c(""))
home_dir   <- ""
output_dir <- ""
start_dir  <- getwd()

library("rstan")
library("RoBTT")

# merge data from the simulations
files   <- list.files(file.path(output_dir, "ensemble"))
results <- lapply(files, function(f){
  
  temp_ensemble <- readRDS(file = file.path(output_dir, "ensemble", f))
  temp_models   <- readRDS(file = file.path(output_dir, "models", f))
  
  return(data.frame(
    BF_effect        = temp_ensemble$overview["Effect", "Incl. BF"],
    BF_heterogeneity = temp_ensemble$overview["Heterogeneity", "Incl. BF"],
    BF_outliers      = temp_ensemble$overview["Outliers", "Incl. BF"],
    averaged_delta.est = temp_ensemble$parameters$averaged["d", "Mean"],
    averaged_delta.lCI = temp_ensemble$parameters$averaged["d", "0.025"],
    averaged_delta.uCI = temp_ensemble$parameters$averaged["d", "0.975"],
    averaged_rho.est   = temp_ensemble$parameters$averaged["r", "Mean"],
    averaged_rho.lCI   = temp_ensemble$parameters$averaged["r", "0.025"],
    averaged_rho.uCI   = temp_ensemble$parameters$averaged["r", "0.975"],
    averaged_nu.est    = temp_ensemble$parameters$averaged["nu", "Mean"],
    averaged_nu.lCI    = temp_ensemble$parameters$averaged["nu", "0.025"],
    averaged_nu.uCI    = temp_ensemble$parameters$averaged["nu", "0.975"],
    conditional_delta.est = temp_ensemble$parameters$conditional["d", "Mean"],
    conditional_delta.lCI = temp_ensemble$parameters$conditional["d", "0.025"],
    conditional_delta.uCI = temp_ensemble$parameters$conditional["d", "0.975"],
    conditional_rho.est   = temp_ensemble$parameters$conditional["r", "Mean"],
    conditional_rho.lCI   = temp_ensemble$parameters$conditional["r", "0.025"],
    conditional_rho.uCI   = temp_ensemble$parameters$conditional["r", "0.975"],
    conditional_nu.est    = temp_ensemble$parameters$conditional["nu", "Mean"],
    conditional_nu.lCI    = temp_ensemble$parameters$conditional["nu", "0.025"],
    conditional_nu.uCI    = temp_ensemble$parameters$conditional["nu", "0.975"],
    marg_lik   = t(sapply(temp_models$overview, function(m) m$marg_lik)),
    post_prob  = t(sapply(temp_models$overview, function(m) unname(m$posterior_prob))),
    delta      = t(sapply(temp_models$overview, function(m) {
      if(any("d[1]" == rownames(m$tab))){
        return(m$tab["d[1]", "Mean"])
      }else{
        return(0)
      }
    })),
    rho        = t(sapply(temp_models$overview, function(m) {
      if(any("r[1]" == rownames(m$tab))){
        return(m$tab["r[1]", "Mean"])
      }else{
        return(0.5)
      }
    })),
    seed = gsub(".RDS", "", f)
  ))
})
results <- do.call(rbind, results)

rownames(results) <- NULL
saveRDS(results, file = file.path(home_dir, "results.RDS"))