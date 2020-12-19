# Simulated Results (including the precision rate bar plot)

# Load Libraries ------------------------------------------------------------------------------------------

rm(list=ls())
library(caret)
library(tidyverse)
library(astsa)
library(ggplot2)
library(dplyr)
library(timeSeries)
library(forecast)
library(quantmod)
library(depmixS4)
source("SIMULATIONS.R")
source("rSTARS.R")

# Data Set 1 Shift in Mean ----------------------------------------------------------------------------------
# HMM -----------------------------------
sim_mean_hmm <- depmix(data.set1 ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=data.set1))
sim_mean_hmmfit <- fit(sim_mean_hmm, verbose = FALSE)
sim_mean_post_probs <- posterior(sim_mean_hmmfit)

# 2 5 10

par(mfrow = c(2, 1))
plot.ts(data.set1, main="Simulated Data with Shifts in Mean", xlab="", ylab="", cex = 1, xy.labels = T)
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(sim_mean_post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability', xlab = 'Index')

# Data Set 2 Shift in Variance ------------------------------------------------------------------------------
# HMM -----------------------------------
sim_var_hmm <- depmix(data.set2 ~ 1, family = gaussian(), nstates = 10, data=data.frame(returns=data.set2))
sim_var_hmmfit <- fit(sim_var_hmm, verbose = FALSE)
sim_var_post_probs <- posterior(sim_var_hmmfit)

par(mfrow = c(2, 1))
plot.ts(data.set2, main="Simulated Data with Shifts in Variance", xlab="", ylab="")
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(sim_var_post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability', xlab = 'Index')
# 2 5 10

# Data Set 3 Shift in Mean and Variance ---------------------------------------------------------------------
# HMM --------------------------------------
sim_meanvar_hmm <- depmix(data.set3 ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=data.set3))
sim_meanvar_hmmfit <- fit(sim_meanvar_hmm, verbose = FALSE)
sim_meanvar_post_probs <- posterior(sim_meanvar_hmmfit)
# 2 3 5
par(mfrow = c(2, 1))
plot.ts(data.set3, main="Simulated Data with Shifts in Mean and Variance", xlab="", ylab="")
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(sim_meanvar_post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability', xlab = 'Index')

# HMM Confusion Stats ----------------------------------------------------------------------------------------

predicted_shifts <- function(posterior.states){ 
  # outputs a vector indicating where shifts occurred
  # 1 indicates no shift, 0 indicates a shift
  
  N <- length(posterior.states)
  pred_shifts <- numeric(N)
  for(i in 2:N+1){
    ifelse(posterior.states[i] - posterior.states[i-1] == 0, pred_shifts[i] <- 1, pred_shifts[i] <- 0)
  }
  return(pred_shifts)
}

hmm_confusionMatrix <- function(shift.type, .data){

  # shift.type - regime shift type - i.e. mean, variance, mean and variance
  # .data is the simulated data set - i.e. data.set1, data.set2, data.set3
  
  set.seed(1)
  n_states_mean_or_var <- c(2, 5, 10)
  n_states_mean_and_var <- c(2, 3, 5)
  
  if (shift.type == "mean" | shift.type == "variance"){
    N = 100
    n_selected_states <- n_states_mean_and_var
  } else {
    N = 200
    n_selected_states <- n_states_mean_or_var
  }
  
  confusion_table <- matrix(NA, nrow = n_selected_states[length(n_selected_states)], ncol = length(n_selected_states))
  colnames(confusion_table) <- c('Precision', 'Misclass', 'Accuracy')
  
    for(i in n_selected_states){
       observed <- rep(c(rep(1, N-1), 0), length(.data)/N)
       observed[length(observed)] <- 1
       sim_hmm <- depmix(.data ~ 1, family = gaussian(), nstates = i, data=data.frame(returns=.data))
       sim_hmmfit <- fit(sim_hmm, verbose = FALSE)
       sim_post_probs <- posterior(sim_hmmfit)
       pred <- predicted_shifts(sim_hmmfit@posterior$state)
       pred[length(pred)] <- 1
       pred_indices <- which(pred == 0)
       rounded_indices <- round(pred_indices, -2)
       
        for(m in 1:length(rounded_indices)){
            if(is.na(rounded_indices[m] - pred_indices[m])){
              break 
            } else if(abs(rounded_indices[m] - pred_indices[m]) <= 5){
              pred[rounded_indices[m]] <- 0 
              pred[pred_indices[m]] <- 1
            } else {
              pred[pred_indices[m]] <- 0
            }
        }
          confusion_stats <- confusionMatrix(as.factor(as.character(pred)), as.factor(as.character(observed)))
          .precision <- confusion_stats$byClass[5]
          .accuracy <- confusion_stats$overall[1]
          .misclass <- 1 - .accuracy
          confusion_table[i,1:3] <- c(.precision, .misclass, .accuracy)
    }
  
        return(summary.stats.table = as.data.frame(na.omit(confusion_table)))
        
}
       
hmm_confusion_mean <- hmm_confusionMatrix("mean", data.set1)  
hmm_confusion_mean
hmm_confusion_variance <- hmm_confusionMatrix("variance", data.set2)
hmm_confusion_variance
hmm_confusion_mean_and_variance <- hmm_confusionMatrix("mean and variance", data.set3)
hmm_confusion_mean_and_variance
  

# STARS Results ---------------------------------------------------------------------------------------------

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "A1:H1001"))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

sim_stars_results <- read_excel_allsheets('SIMULATIONS RESULTS.xlsm')

# MEAN and VARIANCE
stars_confusion_stats <- function(.data = sim_stars_results, N, n){
  # .data is the complete data set from the excel file that has all the RSI values and time series data for all parameters
  # N is the interval of the simulated shifts, every N = 100 or N = 200 observations
  # n denotes the end of a specific shift type data. This value is then taken 8 steps back to get the entire data set. i.e.:
    # n = 9 means that the function will consider the first nine 
  
  
  
  confusion_stats <- list()
  confusion_table <- matrix(NA, nrow = n, ncol = 5)
  confusion_table[,1] <- c(rep(20, 3), rep(50, 3), rep(100, 3))
  confusion_table[,2] <- c(rep(c(0.01, 0.05, 0.1), 3))
  colnames(confusion_table) <- c('Cutoff', 'p-value', 'Precision', 'Misclass', 'Accuracy')

  for(i in (n-8):n){
  observed <- rep(c(rep(1, N-1), 0), length(.data[[i]][,3])/N)
  pred <- ifelse(.data[[i]][,3] == 0, 1, 0)
  pred[length(pred)] <- 1
  observed[length(observed)] <- 1
  pred_indices <- which(pred == 0)
  rounded_indices <- round(pred_indices, -2)
  
  for(m in 1:length(rounded_indices)){
    if(is.na(rounded_indices[m] - pred_indices[m])){
      break 
      } else if(abs(rounded_indices[m] - pred_indices[m]) <= 5){
        pred[rounded_indices[m]] <- 0 
        pred[pred_indices[m]] <- 1
     } else {
       pred[pred_indices[m]] <- 0
       }
  }
  
  confusion_stats[[i]] <- confusionMatrix(as.factor(as.character(pred)), as.factor(as.character(observed)))
   .precision <- confusion_stats[[i]]$byClass[5]
   .accuracy <- confusion_stats[[i]]$overall[1]
   .misclass <- 1 - .accuracy
   confusion_table[i,3:5] <- c(.precision, .misclass, .accuracy)
   
  }
  
  return(summary.stats.table = as.data.frame(confusion_table))
  
}

stars_confusion_mean <- stars_confusion_stats(sim_stars_results, N=100, n=9)
stars_confusion_mean
stars_confusion_var <- na.omit(stars_confusion_stats(sim_stars_results, N=100, n=18))
stars_confusion_var
stars_confusion_meanvar_rs_mean <- na.omit(stars_confusion_stats(sim_stars_results, N=200, n=27))
stars_confusion_meanvar_rs_mean
stars_confusion_meanvar_rs_var <- stars_confusion_stats(sim_stars_results, N=100, n=36)[28:36,]
stars_confusion_meanvar_rs_var


# Comparing HMM and STARS

data_1_STARS <- c(1, 0.007, 0.993)
data_2_STARS <- c(0.4, 0.01, 0.99)
data_3_STARS <- c(0.5, 0.007, 0.99)
data_stars_df <- data.frame(data_1_STARS, data_2_STARS, data_3_STARS)
data_stars_df

data_1_HMM <- c(1, 0.008, 0.992)
data_2_HMM <- c(0.5, 0.009, 0.991)
data_3_HMM <- c(0.667, 0.003, 0.997)
data_hmm_df <- data.frame(data_1_HMM, data_2_HMM, data_3_HMM)
data_hmm_df

data_mat <- matrix(NA, nrow = 6, ncol =5)
data_mat[1,3:5] <- data_1_STARS
data_mat[2,3:5] <- data_2_STARS
data_mat[3,3:5] <- data_3_STARS
data_mat[4,3:5] <- data_1_HMM
data_mat[5,3:5] <- data_2_HMM
data_mat[6,3:5] <- data_3_HMM
data_mat[,2] <- data_sets
data_mat[,1] <- methods_stars_hmm
d1d2d3 <- as.data.frame(data_mat)
d1d2d3[,3:5] <- as.numeric(data_mat[,3:5])
colnames(d1d2d3) <- c('Method', 'Data Set', 'Precision', 'Misclassification', 'Accuracy')
barplot_d1d2d3 <- t(d1d2d3[,3:5])
colnames(barplot_d1d2d3) <- rep(c("D1", "D2", "D3"), 2)
d1d2d3 

prec_rates <- matrix(NA, nrow = 3, ncol = 2) 
prec_rates[1,] <- c(d1d2d3[1, 3], d1d2d3[4, 3])
prec_rates[2,] <- c(d1d2d3[2, 3], d1d2d3[5, 3])
prec_rates[3,] <- c(d1d2d3[3, 3], d1d2d3[6, 3])
colnames(prec_rates) <- c('STARS', 'HMM')
rownames(prec_rates) <- c('Precision', 'Misclassification', 'Accuracy')
prec_rates

# Bar Graph
par(mfrow=c(1,1))
barplot(
  prec_rates,
  col = colors()[c(23, 100, 12)],
  beside = TRUE,
  main = 'Precision Rates by Methods',
  xlab = "Methods",
  args.legend = list(title = "Shift Types", x = "topright", cex = .7),
  border = "grey",
  legend = c('Mean', 'Variance', 'Mean and Variance')
)




