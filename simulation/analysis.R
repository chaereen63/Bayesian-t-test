### load the results and settings
source("simulation/functions.R")
results  <- readRDS("simulation/results.RDS")
settings <- readRDS("simulation/settings.RDS")

# the delta & rho estimates are in the opposite direction
results[,grep("delta", colnames(results))] <- - results[,grep("delta", colnames(results))]
results[,grep("rho",   colnames(results))] <- 1 - results[,grep("rho", colnames(results))]

# the rho in settings corresponds to the ratio, not the allocation estimate
settings$rho <- settings$rho^2 / (1 + settings$rho^2)

results  <- merge(results, settings, by = "seed")
deltas   <- sort(unique(results$delta)) 
rhos     <- sort(unique(results$rho))
nus      <- sort(unique(results$nu), decreasing = TRUE)
ns       <- sort(unique(results$n))

### check the model recovery
# M1 = no effect, no heterogeneity, normal
# M2 = no effect, no heterogeneity, t 
# M3 = no effect, heterogeneity, normal
# M4 = no effect, heterogeneity, t
# M5 = effect, no heterogeneity, normal
# M6 = effect, no heterogeneity, t
# M7 = effect, heterogeneity, normal
# M8 = effect, heterogeneity, t

results$model_best <- apply(results[, grep("marg_lik", colnames(results))], 1, which.max)
results$model_true <- get_true_model(results$delta, results$rho, results$nu)

# set the estimates from the true model
results$model_true_post_prob <- NA
results$model_true_delta     <- NA
results$model_true_rho       <- NA
results$model_true_nu        <- NA
results$model_true_BF        <- NA
results$model_best_delta     <- NA

for(m in unique(results$model_true)){
  results$model_true_post_prob[results$model_true == m] <- results[results$model_true == m, paste0("post_prob.", m)]
  results$model_true_delta[results$model_true == m]     <- results[results$model_true == m, paste0("delta.", m)]
  results$model_true_rho[results$model_true == m]       <- results[results$model_true == m, paste0("rho.", m)]
  results$model_true_nu[results$model_true == m]        <- results[results$model_true == m, paste0("nu.", m)]
  results$model_best_delta[results$model_best == m]     <- results[results$model_best == m, paste0("delta.", m)]
  results$model_true_BF[results$model_true == m]        <- exp(
    results[results$model_true == m, paste0("marg_lik.", switch(as.character(m), "1" = 5, "2" = 6, "3" = 7, "4" = 8, m))] - 
      results[results$model_true == m, paste0("marg_lik.", switch(as.character(m), "5" = 1, "6" = 2, "7" = 3, "8" = 4, m))])
}

temp_BF_57 <- exp(results$marg_lik.5 - results$marg_lik.7)
results$model_student_delta   <- results$delta.5
results$model_welch_delta     <- results$delta.7
results$model_SW_delta        <- results$delta.5 * (temp_BF_57 / (temp_BF_57 + 1)) + results$delta.7 * (1 - temp_BF_57 / (temp_BF_57 + 1))
results$model_student_BF      <- exp(results$marg_lik.5 - results$marg_lik.1)
results$model_welch_BF        <- exp(results$marg_lik.7 - results$marg_lik.3)
results$model_SW_BF           <- apply(exp(results[,c("marg_lik.5", "marg_lik.7")] - results$marg_lik.7) * 0.25, 1, sum) / apply(exp(results[,c("marg_lik.1", "marg_lik.3")] - results$marg_lik.7) * 0.25, 1, sum)

# compute full posterior based RMSE
model_estimates <- readRDS("simulation/model_estimates.RDS")
results  <- merge(results, model_estimates, by = "seed")

results[,c("mean_1", "mean_2", "mean_3", "mean_4", "mean_5", "mean_6", "mean_7", "mean_8")] <- 
  - results[,c("mean_1", "mean_2", "mean_3", "mean_4", "mean_5", "mean_6", "mean_7", "mean_8")]

results$model_student_pRMSE <- sqrt((results$mean_5 - results$delta)^2 + results$sd_5^2)
results$model_welch_pRMSE   <- sqrt((results$mean_7 - results$delta)^2 + results$sd_7^2)

results$model_SW_pBIAS      <- 
  (results$mean_5 - results$delta) * (temp_BF_57 / (temp_BF_57 + 1)) + 
  (results$mean_7 - results$delta) * (1 - temp_BF_57 / (temp_BF_57 + 1))
results$model_SW_pVAR       <- 
  (results$sd_5^2 + (results$mean_5 - results$model_SW_delta)^2) * (temp_BF_57 / (temp_BF_57 + 1)) + 
  (results$sd_7^2 + (results$mean_7 - results$model_SW_delta)^2) * (1 - temp_BF_57 / (temp_BF_57 + 1))
results$model_SW_pRMSE      <- sqrt(results$model_SW_pVAR + results$model_SW_pBIAS^2)

temp_weights_conditional        <- do.call(rbind, lapply(1:nrow(results), function(i){
  bridgesampling::post_prob(results$marg_lik.5[i], results$marg_lik.6[i], results$marg_lik.7[i], results$marg_lik.8[i])
}))
results$model_conditional_pBIAS <- 
  (results$mean_5 - results$delta) * temp_weights_conditional[,1] +
  (results$mean_6 - results$delta) * temp_weights_conditional[,2] +
  (results$mean_7 - results$delta) * temp_weights_conditional[,3] +
  (results$mean_8 - results$delta) * temp_weights_conditional[,4]
results$model_conditional_pVAR  <- 
  (results$sd_5^2 + (results$mean_5 - results$conditional_delta.est)^2) * temp_weights_conditional[,1] + 
  (results$sd_6^2 + (results$mean_6 - results$conditional_delta.est)^2) * temp_weights_conditional[,2] + 
  (results$sd_7^2 + (results$mean_7 - results$conditional_delta.est)^2) * temp_weights_conditional[,3] + 
  (results$sd_8^2 + (results$mean_8 - results$conditional_delta.est)^2) * temp_weights_conditional[,4]
results$model_conditional_pRMSE <- sqrt(results$model_conditional_pVAR + results$model_conditional_pBIAS^2)

rm(temp_weights_conditional)
rm(temp_BF_57)

saveRDS(results, file = "simulation/results_temp.RDS")
results <- readRDS(file = "simulation/results_temp.RDS")

# check recovery with increasing sample sizes for each of the conditions
plot(NA, type = "n", xaxt = "n", ylim = c(0, 1), xlim = c(0, 1), xlab = "N", ylab = "P(correct)", las = 1)
axis(1, seq(0, 1, length.out = length(ns)), ns)
for(delta in deltas){
  for(rho in rhos){
    for(nu in nus){
      temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
      temp_results <- sapply(ns, function(n)
        mean(temp_results$model_best[temp_results$n == n] == temp_results$model_true[temp_results$n == n]))
      lines(seq(0, 1, length.out = length(ns)), temp_results,
            col = ifelse(rho == 1, "blue", "red"),
            lty = ifelse(is.infinite(nu), 1, 2))
    }
  }
}

# check the probability of selecting the correct model
pdf(file = "simulation/model_selection.pdf", width = 7.5, height = 7.5)
{
  ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
  layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
  
  par(mar = c(0,0,0,0))
  empty_plot()
  
  for(delta in deltas){
    empty_plot()
    text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
  }
  for(rho in rhos){
    empty_plot()
    text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
  }
  
  par(mar = c(4, 4, 1, 1))
  par(mar = c(4, 4, 1, 1))
  for(rho in rhos){
    for(delta in deltas){
      
      temp_results <- results[results$delta == delta & results$rho == rho,]
      temp_results <- data.frame(
        est = temp_results[,"model_true_post_prob"],
        n   = factor(temp_results$n1, levels = ns),
        nu  = factor(temp_results$nu, levels = nus)
      )
      
      boxplot(est ~ nu + n, data = temp_results, xlab = "N", ylab = "P(correct)", las = 1, ylim = c(0, 1),
              at = c(1:3, 5:7, 9:11), col = c("green", "yellow", "red"),
              names = FALSE, xaxs = FALSE)
      axis(1, at = c(2, 6, 10), labels = ns)
      abline(h = 1/8, lty = 3)
    }
  }
  legend("topleft", fill = c("green", "yellow", "red"), legend = c(expression(nu==infinity), expression(nu==10), expression(nu==5)), bty = "n")
  dev.off()  
}

# check the inclusion BF
{
  for(BF in c("BF_effect", "BF_heterogeneity", "BF_outliers")){
    pdf(file = paste0("simulation/",BF,".pdf"), width = 7.5, height = 7.5)
    {
      ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
      layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
      
      par(mar = c(0,0,0,0))
      empty_plot()
      text(.5, .5, BF, cex = 1.75, adj = c(0.5, 0.5))  
      
      for(delta in deltas){
        empty_plot()
        text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
      }
      for(rho in rhos){
        empty_plot()
        text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
      }
      
      par(mar = c(4, 4, 1, 1))
      for(rho in rhos){
        for(delta in deltas){
          
          temp_results <- results[results$delta == delta & results$rho == rho,]
          temp_results <- data.frame(
            est = temp_results[,BF],
            n   = factor(temp_results$n1, levels = ns),
            nu  = factor(temp_results$nu, levels = nus)
          )
          
          boxplot(log10(est) ~ nu + n, data = temp_results, xlab = "N", ylab = "log10(BF)", las = 1, ylim = c(-10, 10),
                  at = c(1:3, 5:7, 9:11), col = c("green", "yellow", "red"),
                  names = FALSE, xaxs = FALSE)
          axis(1, at = c(2, 6, 10), labels = ns)
          abline(h = 0, lty = 3)
        }
      }
      legend("topleft", fill = c("green", "yellow", "red"), legend = c(expression(nu==infinity), expression(nu==10), expression(nu==5)), bty = "n")
      dev.off()   
    }
  }
}

# check the inclusion BF (zoomed)
{
  for(BF in c("BF_effect", "BF_heterogeneity", "BF_outliers")){
    pdf(file = paste0("simulation/",BF,"_zoomed.pdf"), width = 7.5, height = 7.5)
    {
      ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
      layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
      
      par(mar = c(0,0,0,0))
      empty_plot()
      text(.5, .5, BF, cex = 1.75, adj = c(0.5, 0.5))  
      
      for(delta in deltas){
        empty_plot()
        text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
      }
      for(rho in rhos){
        empty_plot()
        text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
      }
      
      par(mar = c(4, 4, 1, 1))
      for(rho in rhos){
        for(delta in deltas){
          
          temp_results <- results[results$delta == delta & results$rho == rho,]
          temp_results <- data.frame(
            est = temp_results[,BF],
            n   = factor(temp_results$n1, levels = ns),
            nu  = factor(temp_results$nu, levels = nus)
          )
          
          boxplot(log10(est) ~ nu + n, data = temp_results, xlab = "N", ylab = "log10(BF)", las = 1, ylim = c(-2, 2),
                  at = c(1:3, 5:7, 9:11), col = c("green", "yellow", "red"),
                  names = FALSE, xaxs = FALSE)
          axis(1, at = c(2, 6, 10), labels = ns)
          abline(h = 0, lty = 3)
        }
      }
      legend("topleft", fill = c("green", "yellow", "red"), legend = c(expression(nu==infinity), expression(nu==10), expression(nu==5)), bty = "n")
      dev.off()   
    }
  }
}

# check the correct model estimates
{
  for(par in c("delta", "rho", "nu")){
    pdf(file = paste0("simulation/bias_true_", par, ".pdf"), width = 7.5, height = 7.5)
    
    ylim <- switch(
      par,
      "delta" = c(-1, 1),
      "rho"   = c(-.5, .5),
      "nu"    = c(-10, 10),
    )
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, par, cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho,]
        temp_results <- data.frame(
          est = temp_results[,paste0("model_true_", par)] - temp_results[,par],
          n   = factor(temp_results$n1, levels = ns),
          nu  = factor(temp_results$nu, levels = nus)
        )
        
        
        boxplot(est ~ nu + n, data = temp_results, xlab = "N", ylab = "Bias", las = 1, ylim = ylim,
                at = c(1:3, 5:7, 9:11), col = c("green", "yellow", "red"),
                names = FALSE, xaxs = FALSE)
        axis(1, at = c(2, 6, 10), labels = ns)
        abline(h = 0, lty = 3)
      }
    }
    legend("topleft", fill = c("green", "yellow", "red"), legend = c(expression(nu==infinity), expression(nu==10), expression(nu==5)), bty = "n")
    dev.off()   
  }
}

# check estimates of different models
{
  for(nu in nus){
    
    pdf(file = paste0("simulation/estimates_delta_nu=", nu, ".pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, par, cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- rbind(
          data.frame(
            est  = temp_results[,"model_true_delta"] - temp_results[,"delta"],
            n    = factor(temp_results$n1, levels = ns),
            type = factor("true",   levels = c("true", "best", "ma", "ma-cond"))
          ),
          data.frame(
            est = temp_results[,"model_best_delta"] - temp_results[,"delta"],
            n   = factor(temp_results$n1, levels = ns),
            type= factor("best",    levels = c("true", "best", "ma", "ma-cond"))
          ),
          
          data.frame(
            est = temp_results[,"conditional_delta.est"] - temp_results[,"delta"],
            n   = factor(temp_results$n1, levels = ns),
            type= factor("ma-cond", levels = c("true", "best", "ma", "ma-cond"))
          ),
          data.frame(
            est = temp_results[,"averaged_delta.est"] - temp_results[,"delta"],
            n   = factor(temp_results$n1, levels = ns),
            type= factor("ma",      levels = c("true", "best", "ma", "ma-cond"))
          )
        )
        
        boxplot(est ~ type + n, data = temp_results, xlab = "N", ylab = "Bias", las = 1, ylim = c(-1, 1),
                at = c(1:4, 6:9, 11:14), col = c("green", "yellow", "red", "orange"), xaxt = "n",
                names = FALSE, xaxs = FALSE)
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("true", "best", "ma", "cond"), bty = "n")
    dev.off()   
  }
}

# check bias of different models
bias    <- function(x, y){
  mean(x - y)
}
bias.se <- function(x, y){
  sd(x - y)/sqrt(length(x))
}

{
  for(nu in nus){
    
    pdf(file = paste0("simulation/bias_delta_nu=", nu, ".pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = bias(this_results[,"model_true_delta"], this_results[,"delta"]),
              se   = bias.se(this_results[,"model_true_delta"], this_results[,"delta"]),
              n    = n,
              type = "true"
            ),
            data.frame(
              est  = bias(this_results[,"model_best_delta"], this_results[,"delta"]),
              se   = bias.se(this_results[,"model_best_delta"], this_results[,"delta"]),
              n    = n,
              type = "best"
            ),
            
            data.frame(
              est  = bias(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              se   = bias.se(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma-cond"
            ),
            data.frame(
              est  = bias(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              se   = bias.se(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma"
            )
          )
        }))
        
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(-.3, 0.3), las = 1, ylab = "bias", xlab = "N")
        i <- 1
        for(type in c("true", "best", "ma-cond", "ma")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("true", "best", "cond", "ma"), bty = "n")
    dev.off()   
  }
}

# check RMSE of different models
RMSE    <- function(x, y){
  sqrt(mean((x - y)^2))
}
RMSE.se <- function(x, y){
  sqrt( (length(x) - 1)/length(x) * sum( (sapply(1:length(x), function(i)RMSE(x[-i], y[-i])) - RMSE(x, y) )^2) )
}

{
  for(nu in nus){
    
    pdf(file = paste0("simulation/RMSE_delta_nu=", nu, ".pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = RMSE(this_results[,"model_true_delta"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"model_true_delta"], this_results[,"delta"]),
              n    = n,
              type = "true"
            ),
            data.frame(
              est  = RMSE(this_results[,"model_best_delta"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"model_best_delta"], this_results[,"delta"]),
              n    = n,
              type = "best"
            ),
            
            data.frame(
              est  = RMSE(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma-cond"
            ),
            data.frame(
              est  = RMSE(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma"
            )
          )
        }))
        
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(0, 0.5), las = 1, ylab = "RMSE", xlab = "N")
        i <- 1
        for(type in c("true", "best", "ma-cond", "ma")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("true", "best", "cond", "ma"), bty = "n")
    dev.off()   
  }
}


# check bias of different conditional estimates
bias    <- function(x, y){
  mean(x - y)
}
bias.se <- function(x, y){
  sd(x - y)/sqrt(length(x))
}

{
  for(nu in nus){
    
    pdf(file = paste0("simulation/bias_SWT_delta_nu=", nu, ".pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = bias(this_results[,"delta.5"], this_results[,"delta"]),
              se   = bias.se(this_results[,"delta.5"], this_results[,"delta"]),
              n    = n,
              type = "S"
            ),
            data.frame(
              est  = bias(this_results[,"delta.7"], this_results[,"delta"]),
              se   = bias.se(this_results[,"delta.7"], this_results[,"delta"]),
              n    = n,
              type = "W"
            ),
            data.frame(
              est  = bias(this_results[,"delta.6"], this_results[,"delta"]),
              se   = bias.se(this_results[,"delta.6"], this_results[,"delta"]),
              n    = n,
              type = "T"
            ),
            data.frame(
              est  = bias(this_results[,"delta.8"], this_results[,"delta"]),
              se   = bias.se(this_results[,"delta.8"], this_results[,"delta"]),
              n    = n,
              type = "WT"
            )
          )
        }))
        
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(-.3, 0.3), las = 1, ylab = "bias", xlab = "N")
        i <- 1
        for(type in c("S", "W", "T", "WT")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("S", "W", "T", "WT"), bty = "n")
    dev.off()   
  }
}

# check RMSE of different conditional estimates
RMSE    <- function(x, y){
  sqrt(mean((x - y)^2))
}
RMSE.se <- function(x, y){
  sqrt( (length(x) - 1)/length(x) * sum( (sapply(1:length(x), function(i)RMSE(x[-i], y[-i])) - RMSE(x, y) )^2) )
}

{
  for(nu in nus){
    
    pdf(file = paste0("simulation/RMSE_SWT_delta_nu=", nu, ".pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = RMSE(this_results[,"delta.5"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"delta.5"], this_results[,"delta"]),
              n    = n,
              type = "S"
            ),
            data.frame(
              est  = RMSE(this_results[,"delta.7"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"delta.7"], this_results[,"delta"]),
              n    = n,
              type = "W"
            ),
            data.frame(
              est  = RMSE(this_results[,"delta.6"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"delta.6"], this_results[,"delta"]),
              n    = n,
              type = "T"
            ),
            data.frame(
              est  = RMSE(this_results[,"delta.8"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"delta.8"], this_results[,"delta"]),
              n    = n,
              type = "WT"
            )
          )
        }))

        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(0, 0.5), las = 1, ylab = "RMSE", xlab = "N")
        i <- 1
        for(type in c("S", "W", "T", "WT")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("S", "W", "T", "WT"), bty = "n")
    dev.off()   
  }
}

# RMSE (nu in single figure)
pdf(file = "simulation/RMSE.pdf", width = 7.5, height = 7.5)
{
  ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
  layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
  
  par(mar = c(0,0,0,0))
  empty_plot()
  
  for(delta in deltas){
    empty_plot()
    text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
  }
  for(rho in rhos){
    empty_plot()
    text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
  }
  
  par(mar = c(4, 4, 1, 1))
  for(rho in rhos){
    for(delta in deltas){
      plot(NA, type = "n", xaxt = "n", ylim = c(0, .6), xlim = c(0, 1), xlab = "N", ylab = "RMSE", las = 1)
      axis(1, seq(0, 1, length.out = length(ns)), ns)
      for(nu in nus){
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- data.frame(do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n == n,]
          return(data.frame(
            "averaged" = get_RMSE(this_results$averaged_delta.est, this_results$delta),
            "true"     = get_RMSE(unlist(this_results[, paste0("delta.", this_results$model_true)]), this_results$delta),
            "best"     = get_RMSE(unlist(this_results[, paste0("delta.", this_results$model_best)]), this_results$delta),
            "complex"  = get_RMSE(this_results$delta.8, this_results$delta)
          ))
        })))

        lines(seq(0, 1, length.out = length(ns)), temp_results$averaged, col = "red",
              lty = switch(as.character(nu), "Inf" = 1, "10" = 2, "5" = 3))
        lines(seq(0, 1, length.out = length(ns)), temp_results$best, col = "green",
              lty = switch(as.character(nu), "Inf" = 1, "10" = 2, "5" = 3))
        lines(seq(0, 1, length.out = length(ns)), temp_results$true, col = "blue",
              lty = switch(as.character(nu), "Inf" = 1, "10" = 2, "5" = 3))
      }
    }
  }
  
  legend("topleft", c(expression(nu==infinity), expression(nu==10), expression(nu==5)), lty = c(1,2,3), bty = "n")
  legend("topright", c("MA", "best", "oracle"), col = c("red", "green", "blue"), bty = "n", lwd = 2)
  
  dev.off()  
}

# check RMSE increase from the correct model  (nu in single figure)
pdf(file = "simulation/RMSE_increase.pdf", width = 7.5, height = 7.5)
{
  ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
  layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
  
  par(mar = c(0,0,0,0))
  empty_plot()
  
  for(delta in deltas){
    empty_plot()
    text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
  }
  for(rho in rhos){
    empty_plot()
    text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
  }
  
  par(mar = c(4, 4, 1, 1))
  for(rho in rhos){
    for(delta in deltas){
      plot(NA, type = "n", xaxt = "n", ylim = c(-.2, .2), xlim = c(0, 1), xlab = "N", ylab = "RMSE", las = 1)
      axis(1, seq(0, 1, length.out = length(ns)), ns)
      for(nu in nus){
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- data.frame(do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n == n,]
          return(data.frame(
            "averaged" = get_RMSE(this_results$averaged_delta.est, this_results$delta),
            "true"     = get_RMSE(unlist(this_results[, paste0("delta.", this_results$model_true)]), this_results$delta),
            "best"     = get_RMSE(unlist(this_results[, paste0("delta.", this_results$model_best)]), this_results$delta),
            "complex"  = get_RMSE(this_results$delta.8, this_results$delta)
          ))
        })))
        
        lines(seq(0, 1, length.out = length(ns)), temp_results$averaged - temp_results$true, col = "red",
              lty = switch(as.character(nu), "Inf" = 1, "10" = 2, "5" = 3))
        lines(seq(0, 1, length.out = length(ns)), temp_results$best - temp_results$true, col = "green",
              lty = switch(as.character(nu), "Inf" = 1, "10" = 2, "5" = 3))
      }
    }
  }
  
  legend("topleft", c(expression(nu==infinity), expression(nu==10), expression(nu==5)), lty = c(1,2,3), bty = "n", lwd = 2)
  legend("topright", c("MA", "best"), col = c("red", "green"), bty = "n", lwd = 2)
  
  dev.off()  
}


# check bias RMSE in comparison to from simpler models
{
  for(nu in nus){
    
    pdf(file = paste0("simulation/bias_delta_nu=", nu, "_vs_simpler.pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = bias(this_results[,"delta.5"], this_results[,"delta"]),
              se   = bias.se(this_results[,"model_true_delta"], this_results[,"delta"]),
              n    = n,
              type = "simp"
            ),
            data.frame(
              est  = bias(this_results[,"delta.7"], this_results[,"delta"]),
              se   = bias.se(this_results[,"model_best_delta"], this_results[,"delta"]),
              n    = n,
              type = "heter"
            ),
            
            data.frame(
              est  = bias(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              se   = bias.se(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma-cond"
            ),
            data.frame(
              est  = bias(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              se   = bias.se(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma"
            )
          )
        }))
        
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(-.3, .3), las = 1, ylab = "bias", xlab = "N")
        i <- 1
        for(type in c("simp", "heter", "ma-cond", "ma")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("H1", "H1^r", "cond", "ma"), bty = "n")
    dev.off()   
  }
}
{
  for(nu in nus){
    
    pdf(file = paste0("simulation/RMSE_delta_nu=", nu, "_vs_simpler.pdf"), width = 7.5, height = 7.5)
    
    
    ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
    layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
    
    par(mar = c(0,0,0,0))
    empty_plot()
    text(.5, .5, "delta", cex = 1.75, adj = c(0.5, 0.5))  
    
    for(delta in deltas){
      empty_plot()
      text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
    }
    for(rho in rhos){
      empty_plot()
      text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
    }
    
    par(mar = c(4, 4, 1, 1))
    for(rho in rhos){
      for(delta in deltas){
        
        temp_results <- results[results$delta == delta & results$rho == rho & results$nu == nu,]
        temp_results <- do.call(rbind, lapply(ns, function(n){
          this_results <- temp_results[temp_results$n1 == n,]
          rbind(
            data.frame(
              est  = RMSE(this_results[,"delta.5"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"model_true_delta"], this_results[,"delta"]),
              n    = n,
              type = "simp"
            ),
            data.frame(
              est  = RMSE(this_results[,"delta.7"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"model_best_delta"], this_results[,"delta"]),
              n    = n,
              type = "heter"
            ),
            
            data.frame(
              est  = RMSE(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"conditional_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma-cond"
            ),
            data.frame(
              est  = RMSE(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              se   = RMSE.se(this_results[,"averaged_delta.est"], this_results[,"delta"]),
              n    = n,
              type = "ma"
            )
          )
        }))
        
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 14.5), ylim = c(0, 0.5), las = 1, ylab = "RMSE", xlab = "N")
        i <- 1
        for(type in c("simp", "heter", "ma-cond", "ma")){
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "red", "orange")[i])
          points(c(0,5,10) + i, temp_results[temp_results$type == type, "est"], pch = 0)
          
          arrows(
            x0 = c(0,5,10) + i,
            x1 = c(0,5,10) + i,
            y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
            y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
            length = .05, angle = 90, code = 3)
          i <- i + 1
          
        }
        axis(1, at = c(2.5, 7.5, 12.5), labels = ns)
        abline(h = 0, lty = 3)
        
      }
    }
    legend("topright", fill = c("green", "yellow", "red", "orange"), legend = c("H1", "H1^r", "cond", "ma"), bty = "n")
    dev.off()   
  }
}
# check evidence overestimation
{
 for(par in c("BF_effect", "model_student_BF", "model_welch_BF")){
   pdf(file = paste0("simulation/evidence_overestimation_",par,".pdf"), width = 7.5, height = 7.5)
   {
     ly <- rbind(1:4, cbind(5:7, matrix(8:(8+9-1), ncol = 3, byrow = T)))
     layout(ly, widths = c(0.5, 1, 1, 1), heights = c(0.25, 1, 1, 1))
     
     par(mar = c(0,0,0,0))
     empty_plot()
     text(.5, .5, par, cex = 1.75, adj = c(0.5, 0.5))  
     
     for(delta in deltas){
       empty_plot()
       text(.5, .5, bquote(delta==.(delta)), cex = 1.75, adj = c(0.5, 0.5))  
     }
     for(rho in rhos){
       empty_plot()
       text(.5, .5, bquote(rho==.(round(rho,2))), cex = 1.75, adj = c(0.5, 0.5))  
     }
     
     par(mar = c(4, 4, 1, 1))
     for(rho in rhos){
       for(delta in deltas){
         
         temp_results <- results[results$delta == delta & results$rho == rho,]
         temp_results <- data.frame(
           est = log10(temp_results[,par]) - log10(temp_results$model_true_BF),
           n   = factor(temp_results$n1, levels = ns),
           nu  = factor(temp_results$nu, levels = nus)
         )
         
         boxplot(est ~ nu + n, data = temp_results, xlab = "N", ylab = "log10(BF)", las = 1, ylim = c(-10, 10),
                 at = c(1:3, 5:7, 9:11), col = c("green", "yellow", "red"),
                 names = FALSE, xaxs = FALSE)
         axis(1, at = c(2, 6, 10), labels = ns)
         abline(h = 0, lty = 3)
       }
     }
     legend("topleft", fill = c("green", "yellow", "red"), legend = c(expression(nu==infinity), expression(nu==10), expression(nu==5)), bty = "n")
     dev.off()   
   }
 }
}


# RMSE comparison for paper
paper_conditions <- data.frame(
  delta = deltas[c(3,3,3)],
  rho   = rhos[c(1,3,3)],
  nu    = nus[c(1,1,3)],
  lab   = c("A) Effect", "B) Effect + Uneq. Variance", "C) Effect + Uneq. Variance + Outliers")
)
{
  pdf(file = paste0("simulation/paper_RMSE.pdf"), width = 7.5, height = 3.5)
  
  
  par(mfrow = c(1, 3))
 
  for(j in 1:nrow(paper_conditions)){
      
      temp_results <- results[results$delta == paper_conditions$delta[j] & results$rho == paper_conditions$rho[j] & results$nu == paper_conditions$nu[j],]
      temp_results <- do.call(rbind, lapply(ns, function(n){
        this_results <- temp_results[temp_results$n1 == n,]
        rbind(
          data.frame(
            est  = RMSE(this_results[,"model_student_delta"], this_results[,"delta"]),
            se   = RMSE.se(this_results[,"model_student_delta"], this_results[,"delta"]),
            n    = n,
            type = "Student"
          ),
          data.frame(
            est  = RMSE(this_results[,"model_welch_delta"], this_results[,"delta"]),
            se   = RMSE.se(this_results[,"model_welch_delta"], this_results[,"delta"]),
            n    = n,
            type = "Welch"
          ),
          data.frame(
            est  = RMSE(this_results[,"model_SW_delta"], this_results[,"delta"]),
            se   = RMSE.se(this_results[,"model_SW_delta"], this_results[,"delta"]),
            n    = n,
            type = "MB"
          ),
          data.frame(
            est  = RMSE(this_results[,"conditional_delta.est"], this_results[,"delta"]),
            se   = RMSE.se(this_results[,"conditional_delta.est"], this_results[,"delta"]),
            n    = n,
            type = "RoMB"
          )
        )
      }))
      
      if(j == 1){
        par(mar = c(4,4,1.5,0))
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 12.5), ylim = c(0, 0.4), las = 1, ylab = "RMSE", xlab = "N", bty = "n")
      }else{
        par(mar = c(4,2,1.5,1))
        plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 12.5), ylim = c(0, 0.4), las = 1, ylab = "", yaxt = "n", xlab = "N", bty = "n")  
      }
      
      i <- 1
      for(type in c("Student", "Welch", "MB", "RoMB")){
        points(c(0,4,8) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "blue", "red")[i])
        points(c(0,4,8) + i, temp_results[temp_results$type == type, "est"], pch = 0)
        
        arrows(
          x0 = c(0,4,8) + i,
          x1 = c(0,4,8) + i,
          y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
          y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
          length = .05, angle = 90, code = 3)
        i <- i + 1
        
      }
      mtext(text = paper_conditions$lab[j], side = 3, at = 0, line = .5, adj = 0, cex = .75)
      mtext(text = ns, side = 1, at = c(2.5, 6.5, 10.5), line = 1, cex = .75)
      abline(h = 0, lty = 3)
      
  }
  legend("topright", fill = c("green", "yellow", "blue", "red"), legend = c("Student", "Welch", "MB", "RoMB"), bty = "n")
  dev.off()   
}


# full posterior RMSE comparison for paper
paper_conditions <- data.frame( 
  delta = deltas[c(3,3,3)],
  rho   = rhos[c(1,3,3)],
  nu    = nus[c(1,1,3)],
  lab   = c("A) Effect", "B) Effect + Uneq. Variance", "C) Effect + Uneq. Variance + Outliers")
)
{
  pdf(file = paste0("simulation/paper_posterior_RMSE.pdf"), width = 7.5, height = 3.5)
  
  
  par(mfrow = c(1, 3))
  
  for(j in 1:nrow(paper_conditions)){
    
    temp_results <- results[results$delta == paper_conditions$delta[j] & results$rho == paper_conditions$rho[j] & results$nu == paper_conditions$nu[j],]
    temp_results <- do.call(rbind, lapply(ns, function(n){
      this_results <- temp_results[temp_results$n1 == n,]
      rbind(
        data.frame(
          est  = mean(this_results$model_student_pRMSE),
          se   = sd(this_results$model_student_pRMSE)/sqrt(length(this_results$model_student_pRMSE)),
          n    = n,
          type = "Student"
        ),
        data.frame(
          est  = mean(this_results$model_welch_pRMSE),
          se   = sd(this_results$model_welch_pRMSE)/sqrt(length(this_results$model_welch_pRMSE)),
          n    = n,
          type = "Welch"
        ),
        data.frame(
          est  = mean(this_results$model_SW_pRMSE),
          se   = sd(this_results$model_SW_pRMSE)/sqrt(length(this_results$model_SW_pRMSE)),
          n    = n,
          type = "MB"
        ),
        data.frame(
          est  = mean(this_results$model_conditional_pRMSE),
          se   = sd(this_results$model_conditional_pRMSE)/sqrt(length(this_results$model_conditional_pRMSE)),
          n    = n,
          type = "RoMB"
        )
      )
    }))
    
    if(j == 1){
      par(mar = c(4,4,1.5,0))
      plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 12.5), ylim = c(0, .6), las = 1, ylab = "Posterior RMSE", xlab = "N", bty = "n")
    }else{
      par(mar = c(4,2,1.5,1))
      plot(NA, type = "n", xaxt = "n", xlim = c(0.5, 12.5), ylim = c(0, .6), las = 1, ylab = "", yaxt = "n", xlab = "N", bty = "n")  
    }
    
    i <- 1
    for(type in c("Student", "Welch", "MB", "RoMB")){
      points(c(0,4,8) + i, temp_results[temp_results$type == type, "est"], pch = 15, col = c("green", "yellow", "blue", "red")[i])
      points(c(0,4,8) + i, temp_results[temp_results$type == type, "est"], pch = 0)
      
      arrows(
        x0 = c(0,4,8) + i,
        x1 = c(0,4,8) + i,
        y0 = temp_results[temp_results$type == type, "est"] + qnorm(.025) * temp_results[temp_results$type == type, "se"],
        y1 = temp_results[temp_results$type == type, "est"] + qnorm(.975) * temp_results[temp_results$type == type, "se"],
        length = .05, angle = 90, code = 3)
      i <- i + 1
      
    }
    mtext(text = paper_conditions$lab[j], side = 3, at = 0, line = .5, adj = 0, cex = .75)
    mtext(text = ns, side = 1, at = c(2.5, 6.5, 10.5), line = 1, cex = .75)
    abline(h = 0, lty = 3)
    
  }
  legend("topright", fill = c("green", "yellow", "blue", "red"), legend = c("Student", "Welch", "MB", "RoMB"), bty = "n")
  dev.off()   
}

### show the full distribution
pdf(file = paste0("simulation/paper_pRMSE.pdf"), width = 7.5, height = 3.5)
paper_conditions <- data.frame( 
  delta = deltas[c(3,3,3)],
  rho   = rhos[c(1,3,3)],
  nu    = nus[c(1,1,3)],
  lab   = c("A) Effect", "B) Effect + Uneq. Variance", "C) Effect + Uneq. Variance + Outliers")
)
par(mfrow = c(1, 3), mar = c(2, 4, 4, 1))
for(j in 1:nrow(paper_conditions)){
  
  temp_results <- results[results$delta == paper_conditions$delta[j] & results$rho == paper_conditions$rho[j] & results$nu == paper_conditions$nu[j],]
  
  this_results <- rbind(
    data.frame(
      est  = (temp_results$model_student_pRMSE),
      n    = temp_results$n1,
      type = "Student"
    ),
    data.frame(
      est  = (temp_results$model_welch_pRMSE),
      n    = temp_results$n1,
      type = "Welch"
    ),
    data.frame(
      est  = (temp_results$model_SW_pRMSE),
      n    = temp_results$n1,
      type = "MB"
    ),
    data.frame(
      est  = (temp_results$model_conditional_pRMSE),
      n    = temp_results$n1,
      type = "RoMB"
    )
  )
  this_results$type <- factor(this_results$type, levels = c("Student", "Welch", "MB", "RoMB"))
  boxplot(est ~ type * n, data = this_results, at = c(1:4, 6:9, 11:14),
          col = rep(alpha(c("green", "yellow", "blue", "red"), alpha = .50), 4),
          ylab = "Posterior RMSE", frame = FALSE, xaxt = "n", yaxt = "n", las = 1, names = NA, ylim = c(0, 1.3))
  axis(2, at = seq(0, 1.2, 0.2), las = 1)
  mtext(ns[1], side = 1, at = 2.5, line = 0, cex = .75)
  mtext(ns[2], side = 1, at = 7.5, line = 0, cex = .75)
  mtext(ns[3], side = 1, at = 12.5, line = 0, cex = .75)
  mtext(text = paper_conditions$lab[j], side = 3, at = 0, line = .5, adj = 0, cex = .75)
  abline(h = 0, lty = 3)
  # temp_x <- 0
  # for(n in unique(this_results$n)){
  #   for(m in c("Student", "Welch", "MB", "RoMB")){
  #     temp_x <- temp_x + 1
  #     points(temp_x, mean(this_results$est[this_results$n == n & this_results$type == m]), pch = 16)
  #   }
  #   temp_x <- temp_x + 1
  # }
  
}
legend("topright", fill = c("green", "yellow", "blue", "red"), legend = c("Student", "Welch", "MB", "RoMB"), bty = "n")
dev.off()


# evidence Distortion for paper
paper_conditions <- data.frame(
  delta = deltas[c(1,1,1,3,3,3)],
  rho   = rhos[c(1,3,3,1,3,3)],
  nu    = nus[c(1,1,3,1,1,3)],
  lab   = c("A1) No Effect", "B1) No Effect + Uneq. Variance", "C1) No Effect + Uneq. Variance + Outliers",
            "A2) Effect", "B2) Effect + Uneq. Variance", "C2) Effect + Uneq. Variance + Outliers")
)
{
  pdf(file = paste0("simulation/paper_evidence.pdf"), width = 7.5, height = 7.5)
  
  
  par(mfrow = c(2, 3))
  
  for(j in 1:nrow(paper_conditions)){
    
    temp_results <- results[results$delta == paper_conditions$delta[j] & results$rho == paper_conditions$rho[j] & results$nu == paper_conditions$nu[j],]
    temp_results <- do.call(rbind, lapply(ns, function(n){
      this_results <- temp_results[temp_results$n1 == n,]
      rbind(
        data.frame(
          est  = log10(this_results$model_student_BF) - log10(this_results$model_true_BF),
          n    = n,
          type = "Student"
        ),
        data.frame(
          est  = log10(this_results$model_welch_BF) - log10(this_results$model_true_BF),
          n    = n,
          type = "Welch"
        ),
        data.frame(
          est  = log10(this_results$model_SW_BF) - log10(this_results$model_true_BF),
          n    = n,
          type = "MB"
        ),
        data.frame(
          est  = log10(this_results$BF_effect) - log10(this_results$model_true_BF),
          n    = n,
          type = "RoMB"
        )
      )
    }))
    temp_results$type <- factor(temp_results$type, levels = c("Student", "Welch", "MB", "RoMB"))
    
    if(j %in% c(1,4)){
      par(mar = c(4.5,4,1.5,0))
      boxplot(est ~ type + n, frame = FALSE, data = temp_results, xaxt = "n", xlim = c(0.5, 14.5), ylim = c(-2, 2), las = 1, yaxt = "n", ylab = "Evidence Distortion", xlab = "N", bty = "n",
              at = c(1:4, 6:9, 11:14), col = c("green", "yellow", "blue", "red"))
      axis(2, at = seq(-2, 2), labels = c(expression(1/100), expression(1/10), 1, 10, 100), las = 1)
    }else{
      par(mar = c(4.5,2,1.5,1))
      boxplot(est ~ type + n,  frame = FALSE, data = temp_results, xaxt = "n", xlim = c(0.5, 14.5), ylim = c(-2, 2), las = 1, ylab = "", yaxt = "n", xlab = "N", bty = "n",
              at = c(1:4, 6:9, 11:14), col = c("green", "yellow", "blue", "red"))  
    }
    

    mtext(text = paper_conditions$lab[j], side = 3, at = 0, line = .5, adj = 0, cex = .75)
    mtext(text = ns, side = 1, at = c(2.5, 7.5, 12.5), line = 1, cex = .75)
    abline(h = 0, lty = 3)
    
  }
  legend("topleft", fill = c("green", "yellow", "blue", "red"), legend = c("Student", "Welch", "MB", "RoMB"), bty = "n")
  dev.off()   
}

# dichotomoized evidence Distortion for paper
paper_conditions <- data.frame(
  delta = deltas[c(1,1,1,3,3,3)],
  rho   = rhos[c(1,3,3,1,3,3)],
  nu    = nus[c(1,1,3,1,1,3)],
  lab   = c("A1) No Effect", "B1) No Effect + Heterogeneity", "C1) No Effect + Heterogeneity + Outliers",
            "A2) Effect", "B2) Effect + Heterogeneity", "C2) Effect + Heterogeneity + Outliers")
)
BF_lables <- function(bf){
  ifelse(bf < 1/10, -3, 
         ifelse(bf < 1/3, -2,
                ifelse(bf < 1, -1,
                       ifelse(bf < 3, 1, 
                              ifelse(bf < 10, 2, 3)))))
}
{
  pdf(file = paste0("simulation/paper_evidence2.pdf"), width = 7.5, height = 7.5)
  
  
  par(mfrow = c(2, 3))
  
  for(j in 1:nrow(paper_conditions)){
    
    temp_results <- results[results$delta == paper_conditions$delta[j] & results$rho == paper_conditions$rho[j] & results$nu == paper_conditions$nu[j],]
    temp_results <- do.call(rbind, lapply(ns, function(n){
      this_results <- temp_results[temp_results$n1 == n,]
      rbind(
        data.frame(
          bf   = BF_lables(this_results$model_student_BF),
          n    = n,
          type = "Student"
        ),
        data.frame(
          bf   = BF_lables(this_results$model_welch_BF),
          n    = n,
          type = "Welch"
        ),
        
        data.frame(
          bf   = BF_lables(this_results$BF_effect),
          n    = n,
          type = "RoMB"
        )
      )
    }))
    temp_results$type <- factor(paste0(temp_results$type, "-", temp_results$n), levels = 
                                  c("Student-10", "Welch-10", "RoMB-10",
                                    "Student-100", "Welch-100", "RoMB-100",
                                    "Student-1000", "Welch-1000", "RoMB-1000"))
    temp_results$bf   <- factor(temp_results$bf, levels = c(-3, -2, -1, 1, 2, 3)) 
    temp_table <- table(temp_results$bf, temp_results$type)/1000

    
    if(j %in% c(1,4)){
      par(mar = c(4.5,4,1.5,0))
      barplot(temp_table, las = 1, col = c(
        scales::alpha("blue", .75), scales::alpha("blue", .50),scales::alpha("blue", .25),
        scales::alpha("red",  .25), scales::alpha("red",  .50),scales::alpha("red",  .75)), xaxt = "n",
        xlab = "", ylab = "", main = "", cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, cex.sub = 1.3, cex = 1.3, 
        space = c(0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.5, 0.1, 0.1))
      mtext(c(10, 100, 1000), side = 1, line = 1, at = c(1.7, 5.4, 9.1), cex = 1)
      
    }else{
      par(mar = c(4.5,2,1.5,1))
      barplot(temp_table, las = 1, col = c(
        scales::alpha("blue", .75), scales::alpha("blue", .50),scales::alpha("blue", .25),
        scales::alpha("red",  .25), scales::alpha("red",  .50),scales::alpha("red",  .75)), xaxt = "n",
        xlab = "", ylab = "", main = "", cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, cex.sub = 1.3, cex = 1.3, yaxt = "n",
        space = c(0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.5, 0.1, 0.1))
      mtext(c(10, 100, 1000), side = 1, line = 1, at = c(1.7, 5.4, 9.1), cex = 1)
      
    }
    
    
    mtext(text = paper_conditions$lab[j], side = 3, at = 0, line = .5, adj = 0, cex = .70)
    
  }
  dev.off()   
}

