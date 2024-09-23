# computation
simulate_data <- function(mean1, mean2, sd1, sd2, df1, df2, n1, n2){
  
  # the first group
  if(is.infinite(df1)){
    x1 <- stats::rnorm(n1, mean = mean1, sd = sd1)
  }else{
    x1 <- extraDistr::rlst(n1, df = df1, mu = mean1, sigma = sd1 / sqrt(df1 / (df1 - 2)))
  }

  # the second group  
  if(is.infinite(df2)){
    x2 <- stats::rnorm(n2, mean = mean2, sd = sd2)
  }else{
    x2 <- extraDistr::rlst(n2, df = df2, mu = mean2, sigma = sd2 / sqrt(df2 / (df2 - 2)))
  }
  
  return(list(
    x1 = x1,
    x2 = x2
  ))
}
pooled_sd     <- function(sd1, sd2, n1, n2){
  sqrt( (sd1^2 * n1 + sd2^2 * n2) / (n1 + n2) )   #가중치가 들어간 SD를 사용함. 이거 근거 확인하기
}
cohens_d      <- function(mean1, mean2, sd1, sd2, n1, n2){
  (mean1 - mean2) / pooled_sd(sd1, sd2, n1, n2)
}

# helper functions
get_loop       <- function(dir, file){
  
  settings_file_loaded <<- FALSE
  
  while(settings_file_loaded == FALSE){
    tryCatch({
      loop <- read.table(file.path(dir, paste0(file, ".txt")))[1,1]
      settings_file_loaded <<- TRUE
    }, warning = function(w) {
    }, error   = function(e) {
    }, finally = {
    })
  }
  
  
  loop <- loop + 1
  
  
  settings_file_saved <<- FALSE
  
  while(settings_file_saved == FALSE){
    tryCatch({
      write.table(loop, file = file.path(dir, paste0(file, ".txt")), row.names = F, col.names = F)
      settings_file_saved <<- TRUE
    }, warning = function(w) {
    }, error   = function(e) {
    }, finally = {
    })
  }
  
  return(loop)
}
get_true_model <- function(delta, rho, nu){
  effect        <- delta != 0
  heterogeneity <- rho   != .5
  normal        <- is.infinite(nu)
  return(sapply(paste0(effect, heterogeneity, normal), function(M) switch(
    M,
    "FALSEFALSETRUE"  = 1,
    "FALSEFALSEFALSE" = 2,
    "FALSETRUETRUE"   = 3,
    "FALSETRUEFALSE"  = 4,
    "TRUEFALSETRUE"   = 5,
    "TRUEFALSEFALSE"  = 6,
    "TRUETRUETRUE"    = 7,
    "TRUETRUEFALSE"   = 8
  )))
}
empty_plot     <- function() plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0,1), mar = c(0,0,0,0))
get_RMSE       <- function(estimate, true){
  return(sqrt(mean((estimate - true)^2)))
}
