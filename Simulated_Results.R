# Simulated Results (including the precision rate and standard deviation bar plots)

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
library(data.table)
library(depmixS4)
library(readxl)
source("SIMULATIONS.R")
source("rSTARS.R")

# Data Set 1 Shift in Mean ----------------------------------------------------------------------------------
# HMM EXAMPLE -----------------------------------
set.seed(1)
df.eg <- data.frame(matrix(unlist(dataset_1), ncol=length(dataset_1), byrow=F))
sim_mean_hmm <- depmix(df.eg[,3] ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=df.eg[,3]))
sim_mean_hmmfit <- fit(sim_mean_hmm, verbose = FALSE)
sim_mean_post_probs <- posterior(sim_mean_hmmfit)

plot.ts(df.eg[,3], main="Simulated Data with Shifts in Mean", xlab="", ylab="", cex = 1, xy.labels = T)
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(sim_mean_post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability', xlab = 'Index')

# HMM Confusion Stats ----------------------------------------------------------------------------------------

predicted_shifts <- function(posterior.states){ 
  # outputs a vector indicating where shifts occurred
  # 1 indicates no shift, 0 indicates a shift
  
  N <- length(posterior.states)
  pred_shifts <- c()
  pred_shifts[1] <- 1
  for(i in 2:N){
    ifelse(posterior.states[i] - posterior.states[i-1] == 0, pred_shifts[i] <- 1, pred_shifts[i] <- 0)
  }
  return(pred_shifts)
}


hmm_confusionMatrix <- function(shift.type, .data){

  # shift.type (character) - regime shift type - i.e. mean, variance, mean and variance
  # .data (list) is the simulated data set - i.e. dataset_1, dataset_2, dataset_3
  
  dataset_df <- data.frame(matrix(unlist(.data), ncol=length(.data), byrow=F))
  set.seed(1)
  n_states_mean_or_var <- c(2, 5, 10)
  n_states_mean_and_var <- c(2, 5, 10)
  
  if (shift.type == "mean" | shift.type == "variance"){
    N = 100
    n_selected_states <- n_states_mean_or_var
  } else {
    N = 200
    n_selected_states <- n_states_mean_and_var
  }
  
  confusion_table <- matrix(NA, nrow = length(n_selected_states), ncol = 3)
  confusion_table[,1] <- c(as.character(n_selected_states[1]), as.character(n_selected_states[2]), as.character(n_selected_states[3]))
  colnames(confusion_table) <- c('No. of States', 'Avg Precision', 'Standard Deviation')
  
  confusion_storage_matrix <- matrix(NA, nrow = ncol(dataset_df), ncol = 3)
  colnames(confusion_storage_matrix) <- c(as.character(n_selected_states[1]), as.character(n_selected_states[2]), as.character(n_selected_states[3]))
  

  for (k in 1:100){
    
    confusion_list <- list()
    
    for(i in n_selected_states){
       observed <- rep(c(rep(1, N-1), 0), nrow(dataset_df)/N)
       observed[length(observed)] <- 1
       sim_hmm <- depmix(dataset_df[,k] ~ 1, family = gaussian(), nstates = i, data=data.frame(returns=dataset_df[,k]))
       sim_hmmfit <- fit(sim_hmm, verbose = FALSE)
       sim_post_probs <- posterior(sim_hmmfit)
       pred <- predicted_shifts(sim_post_probs$state)
       pred[length(pred)] <- 1
       pred_indices <- which(pred == 0)
       rounded_indices <- round(pred_indices, -2)
       
        for(m in 1:length(rounded_indices)){
            if(is.na(rounded_indices[m] - pred_indices[m])){
              break 
            } else if(abs(rounded_indices[m] - pred_indices[m]) == 0){
              break
            }
              else if(abs(rounded_indices[m] - pred_indices[m]) <= 5){
              pred[rounded_indices[m]] <- 0 
              pred[pred_indices[m]] <- 1
            } else {
              pred[pred_indices[m]] <- 0
            }
        }
          
       if(i == n_selected_states[1]){
         store.i = 1
       } else if(i == n_selected_states[2]){
         store.i = 2
       } else {
         store.i = 3
       }
       
          confusion_list[[store.i]] <- confusionMatrix(as.factor(as.character(pred)), as.factor(as.character(observed)))
          .precision <- confusion_list[[store.i]]$byClass[5]
          confusion_storage_matrix[k,store.i] <- .precision
          
    }
  
}
  
  confusion_table[,3] <- c(sd(confusion_storage_matrix[,1], na.rm = T), sd(confusion_storage_matrix[,2], na.rm = T), sd(confusion_storage_matrix[,3], na.rm = T))
  confusion_table[,2] <- c(mean(confusion_storage_matrix[,1], na.rm = T), mean(confusion_storage_matrix[,2], na.rm = T), mean(confusion_storage_matrix[,3], na.rm = T))
  return(list(summary.stats.table = confusion_table, iterations = confusion_storage_matrix))
  
}

hmm_mean_prec <- hmm_confusionMatrix("mean", dataset_1)  
hmm_mean_prec
hmm_var_prec <- hmm_confusionMatrix('variance', dataset_2)
hmm_var_prec
hmm_meanvar_prec <- hmm_confusionMatrix('mean and variance', dataset_3)
hmm_meanvar_prec

# Results Bar Plot ------------------------------------------------
# For Data Set 1 - Shift in Mean
barplot_hmm_d1 <- hmm_mean_prec$summary.stats.table[,2] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d1) <- c('2', '5', '10') 
barplot(barplot_hmm_d1,
  col = colors()[c(23)],
  beside = T,
  xlab = "Number of States",
  ylab = "Average Precision Rate",
  ylim = c(0, 1),
  border = "grey",
  space = c(0.2, 0.1, 0.1)
)

barplot_hmm_d1_error <- hmm_mean_prec$summary.stats.table[,3] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d1_error) <- c('2', '5', '10') 
barplot(barplot_hmm_d1_error,
        col = colors()[c(23)],
        beside = T,
        xlab = "Number of States",
        ylab = "Standard Deviation",
        border = "grey",
        space = c(0.2, 0.1, 0.1)
)

# For Data Set 2 - Shift in Variance
barplot_hmm_d2 <- hmm_var_prec$summary.stats.table[,2] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d2) <- c('2', '5', '10') 
barplot(barplot_hmm_d2,
        col = colors()[c(23)],
        beside = T,
        xlab = "Number of States",
        ylab = "Average Precision Rate",
        ylim = c(0, 1),
        border = "grey",
        space = c(0.2, 0.1, 0.1)
)

barplot_hmm_d2_error <- hmm_var_prec$summary.stats.table[,3] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d2_error) <- c('2', '5', '10') 
barplot(barplot_hmm_d2_error,
        col = colors()[c(23)],
        beside = T,
        xlab = "Number of States",
        ylab = "Standard Deviation",
        border = "grey",
        space = c(0.2, 0.1, 0.1)
)

# For Data Set 3 - Shift in Mean and Variance
barplot_hmm_d3 <- hmm_meanvar_prec$summary.stats.table[,2] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d3) <- c('2', '5', '10') 
barplot(barplot_hmm_d3,
        col = colors()[c(23)],
        beside = T,
        xlab = "Number of States",
        ylab = "Average Precision Rate",
        ylim = c(0, 1),
        border = "grey",
        space = c(0.2, 0.1, 0.1)
)

barplot_hmm_d3_error <- hmm_meanvar_prec$summary.stats.table[,3] %>% 
  as.numeric() %>%
  as.matrix() %>%
  t()
colnames(barplot_hmm_d3_error) <- c('2', '5', '10') 
barplot(barplot_hmm_d3_error,
        col = colors()[c(23)],
        beside = T,
        xlab = "Number of States",
        ylab = "Standard Deviation",
        border = "grey",
        space = c(0.2, 0.1, 0.1)
)

# STARS Results ---------------------------------------------------------------------------------------------
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "A1:H1001"))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

file.list <- list.files(pattern='Data Set') 
df.list <- lapply(file.list, read_excel_allsheets) # List containing the STARS results of the data sets
dataset_ones <- file.list[1:9]
dataset_twos <- file.list[10:18]
dataset_threes_mean <- file.list[19:27]
dataset_threes_var <- file.list[28:36]

# STARS Function
stars_confusion_stats <- function(.data = df.list, dataset_, N){
  # .data (list) is the complete data set from the excel file that has all the RSI values and time series data for all parameters
  # N (numerical) is the interval of the simulated shifts, i.e. if N = 100, then that means that there was a shift at every 100 shift..
  # dataset_ (character), look at the dataset_'' objects above. eg. 'dataset_ones'.  
  
  if(dataset_ == 'dataset_ones'){
    start_index = 1
    end_index = 9
  } else if(dataset_ == 'dataset_twos'){
    start_index = 10
    end_index = 18
  } else if (dataset_ == 'dataset_threes_mean'){
    start_index = 19
    end_index = 27
  } else {
    start_index = 28
    end_index =36
  }
  
  # Table for storing the average precision rate for the data sets 
  final_confusion_stats <- list()
  final_confusion_table <- matrix(NA, nrow = 9, ncol = 4)
  final_confusion_table[,1] <- c(rep(100, 3), rep(20, 3), rep(50, 3))
  final_confusion_table[,2] <- c(rep(c(0.01, 0.05, 0.1), 3))
  colnames(final_confusion_table) <- c('Cutoff', 'p-value', 'Avg Precision', 'Std Deviation')
  
  for (k in start_index:end_index){
    
    confusion_list <- list()
    confusion_vector <- c()
    
  for(i in 1:100){
    observed <- rep(c(rep(1, N-1), 0), length(.data[[k]][[i]][,3])/N)
    pred <- ifelse(.data[[k]][[i]][,3] == 0, 1, 0)
    pred[length(pred)] <- 1
    observed[length(observed)] <- 1
    pred_indices <- which(pred == 0)
    rounded_indices <- round(pred_indices, -2)
  
  for(m in 1:length(rounded_indices)){
    if(is.na(rounded_indices[m] - pred_indices[m])){
      break 
    } else if(abs(rounded_indices[m] - pred_indices[m]) == 0){
      break
    }
        else if(abs(rounded_indices[m] - pred_indices[m]) <= 5){
        pred[rounded_indices[m]] <- 0 
        pred[pred_indices[m]] <- 1
     } else {
       pred[pred_indices[m]] <- 0
       }
  }
  
  confusion_list[[i]] <- confusionMatrix(as.factor(as.character(pred)), as.factor(as.character(observed)))
   .precision <- confusion_list[[i]]$byClass[5]
   confusion_vector[i] <- .precision
   
  }
  
  if(dataset_ == 'dataset_ones'){
    store.k <- k
  }  else if(dataset_ == 'dataset_twos'){
    store.k <- k - 9
  } else if(dataset_ == 'dataset_threes_mean'){
    store.k <- k - 18
  } else {
    store.k <- k - 27
  }
    
  final_confusion_table[store.k, 3] <- mean(confusion_vector, na.rm = T)
  final_confusion_table[store.k, 4] <- sd(confusion_vector, na.rm = T)

  }
  return(summary.stats.table = as.data.frame(final_confusion_table))
}

avg_prec_stars_mean <- stars_confusion_stats(df.list, 'dataset_ones', 100)
avg_prec_stars_mean
avg_prec_stars_var <- stars_confusion_stats(df.list, 'dataset_twos', 100)
avg_prec_stars_var
avg_prec_stars_meanvar_mean <- stars_confusion_stats(df.list, 'dataset_threes_mean', 100)
avg_prec_stars_meanvar_mean
avg_prec_stars_meanvar_var <- stars_confusion_stats(df.list, 'dataset_threes_var', 100)
avg_prec_stars_meanvar_var

# Results Bar Plots -----------------------------------------------------------
# For Data Set 1 - Shift in Mean
barplot_stars_mean <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_mean) <- c('20', '50', '100')
rownames(barplot_stars_mean) <- c('0.01', '0.05', '0.1')
barplot_stars_mean[,1] <- avg_prec_stars_mean[4:6, 3]
barplot_stars_mean[,2] <- avg_prec_stars_mean[7:9, 3]
barplot_stars_mean[,3] <- avg_prec_stars_mean[1:3, 3]
colnames(barplot_stars_mean) <- c('20', '50', '100') 
barplot(barplot_stars_mean,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Average Precision Rate",
        ylim = c(0, 0.8),
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

barplot_stars_mean_error <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_mean_error) <- c('20', '50', '100')
rownames(barplot_stars_mean_error) <- c('0.01', '0.05', '0.1')
barplot_stars_mean_error[,1] <- avg_prec_stars_mean[4:6, 4]
barplot_stars_mean_error[,2] <- avg_prec_stars_mean[7:9, 4]
barplot_stars_mean_error[,3] <- avg_prec_stars_mean[1:3, 4]
colnames(barplot_stars_mean_error) <- c('20', '50', '100') 
barplot(barplot_stars_mean_error,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylim = c(0, 0.03),
        ylab = "Standard Deviation",
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

# For Data Set 2 - Shift in the Variance
barplot_stars_var <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_var) <- c('20', '50', '100')
rownames(barplot_stars_var) <- c('0.01', '0.05', '0.1')
barplot_stars_var[,1] <- avg_prec_stars_var[4:6, 3]
barplot_stars_var[,2] <- avg_prec_stars_var[7:9, 3]
barplot_stars_var[,3] <- avg_prec_stars_var[1:3, 3]
colnames(barplot_stars_var) <- c('20', '50', '100') 
barplot(barplot_stars_var,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Average Precision Rate",
        ylim = c(0, 0.8),
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

barplot_stars_var_error <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_var_error) <- c('20', '50', '100')
rownames(barplot_stars_var_error) <- c('0.01', '0.05', '0.1')
barplot_stars_var_error[,1] <- avg_prec_stars_var[4:6, 4]
barplot_stars_var_error[,2] <- avg_prec_stars_var[7:9, 4]
barplot_stars_var_error[,3] <- avg_prec_stars_var[1:3, 4]
colnames(barplot_stars_var_error) <- c('20', '50', '100') 
barplot(barplot_stars_var_error,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Standard Deviation",
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

# For Data Set 3  - Shift in the Mean
barplot_stars_meanvar_mean <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_meanvar_mean) <- c('20', '50', '100')
rownames(barplot_stars_meanvar_mean) <- c('0.01', '0.05', '0.1')
barplot_stars_meanvar_mean[,1] <- avg_prec_stars_meanvar_mean[4:6, 3]
barplot_stars_meanvar_mean[,2] <- avg_prec_stars_meanvar_mean[7:9, 3]
barplot_stars_meanvar_mean[,3] <- avg_prec_stars_meanvar_mean[1:3, 3]
colnames(barplot_stars_meanvar_mean) <- c('20', '50', '100') 
barplot(barplot_stars_meanvar_mean,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Average Precision Rate",
        ylim = c(0, 0.8),
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

barplot_stars_meanvar_mean_error <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_meanvar_mean_error) <- c('20', '50', '100')
rownames(barplot_stars_meanvar_mean_error) <- c('0.01', '0.05', '0.1')
barplot_stars_meanvar_mean_error[,1] <- avg_prec_stars_meanvar_mean[4:6, 4]
barplot_stars_meanvar_mean_error[,2] <- avg_prec_stars_meanvar_mean[7:9, 4]
barplot_stars_meanvar_mean_error[,3] <- avg_prec_stars_meanvar_mean[1:3, 4]
colnames(barplot_stars_meanvar_mean_error) <- c('20', '50', '100') 
barplot(barplot_stars_meanvar_mean_error,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Standard Deviation",
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

# For Data Set 3 - Shift in Variance
barplot_stars_meanvar_var <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_meanvar_var) <- c('20', '50', '100')
rownames(barplot_stars_meanvar_var) <- c('0.01', '0.05', '0.1')
barplot_stars_meanvar_var[,1] <- avg_prec_stars_meanvar_var[4:6, 3]
barplot_stars_meanvar_var[,2] <- avg_prec_stars_meanvar_var[7:9, 3]
barplot_stars_meanvar_var[,3] <- avg_prec_stars_meanvar_var[1:3, 3]
colnames(barplot_stars_meanvar_var) <- c('20', '50', '100') 
barplot(barplot_stars_meanvar_var,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Average Precision Rate",
        ylim = c(0, 0.8),
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)

barplot_stars_meanvar_var_error <- matrix(NA, nrow = 3, ncol = 3)
colnames(barplot_stars_meanvar_var_error) <- c('20', '50', '100')
rownames(barplot_stars_meanvar_var_error) <- c('0.01', '0.05', '0.1')
barplot_stars_meanvar_var_error[,1] <- avg_prec_stars_meanvar_var[4:6, 4]
barplot_stars_meanvar_var_error[,2] <- avg_prec_stars_meanvar_var[7:9, 4]
barplot_stars_meanvar_var_error[,3] <- avg_prec_stars_meanvar_var[1:3, 4]
colnames(barplot_stars_meanvar_var_error) <- c('20', '50', '100') 
barplot(barplot_stars_meanvar_var_error,
        col = colors()[c(23, 100, 12)],
        beside = T,
        xlab = "Cutoff Length",
        ylab = "Standard Deviation",
        args.legend = list(title = "p Value", x = "topright", cex = 1.1, box.col = 'grey'),
        border = "grey",
        legend = c('0.01', '0.05', '0.1')
)


# Comparing HMM and STARS -----------------------------------------------------------------------------------------------------------------
# For the Highest Average Precision Rate ---------------------------------------------------------------------------------------------------
# STARS
best_data_1_STARS <- avg_prec_stars_mean[which.max(avg_prec_stars_mean$`Avg Precision`),]
best_data_2_STARS <- avg_prec_stars_var[which.max(avg_prec_stars_var$`Avg Precision`),]
best_prec_d3_mean <- avg_prec_stars_var[which.max(avg_prec_stars_meanvar_mean$`Avg Precision`),]
best_prec_d3_var <- avg_prec_stars_var[which.max(avg_prec_stars_meanvar_var$`Avg Precision`),]
best_data_3_STARS <- (best_prec_d3_mean[3] + best_prec_d3_var[3])/2
best_stars <-data.frame(best_data_1_STARS[3], best_data_2_STARS[3], best_data_3_STARS)
best_stars
best_data_3_STARS_error <- (best_prec_d3_mean[4] + best_prec_d3_var[4])/2
best_stars_error <-data.frame(best_data_1_STARS[4], best_data_2_STARS[4], best_data_3_STARS_error)
best_stars_error

# HMM
best_data_1_HMM <- hmm_mean_prec$summary.stats.table[which.max(as.numeric(hmm_mean_prec$summary.stats.table[,2])),] %>% as.numeric()
best_data_2_HMM <- hmm_var_prec$summary.stats.table[which.max(as.numeric(hmm_var_prec$summary.stats.table[,2])),] %>% as.numeric()
best_data_3_HMM <- hmm_meanvar_prec$summary.stats.table[which.max(as.numeric(hmm_meanvar_prec$summary.stats.table[,2])),] %>% as.numeric()
best_data_hmm <- data.frame(best_data_1_HMM[2], best_data_2_HMM[2], best_data_3_HMM[2])
best_data_hmm
best_hmm_error <-data.frame(best_data_1_HMM[3], best_data_2_HMM[3], best_data_3_HMM[3])
best_hmm_error

# Average Precision Rate  Summary
data_sets <- rep(c(1, 2, 3), 2)
methods_stars_hmm <- c(rep('STARS', 3), rep('HMM', 3))
data_mat <- matrix(5, nrow = 6, ncol =3) %>% as.data.frame() 
data_mat[1,3] <- best_data_1_STARS[3]
data_mat[2,3] <- best_data_2_STARS[3]
data_mat[3,3] <- best_data_3_STARS
data_mat[4,3] <- best_data_1_HMM[2]
data_mat[5,3] <- best_data_2_HMM[2]
data_mat[6,3] <- best_data_3_HMM[2]
data_mat[,2] <- as.factor(data_sets) 
data_mat[,1] <- methods_stars_hmm
d1d2d3 <- as.data.frame(data_mat)
d1d2d3[,3] <- as.numeric(data_mat[,3])
colnames(d1d2d3) <- c('Method', 'Data Set', 'Precision')
barplot_d1d2d3 <- t(d1d2d3[,3])
colnames(barplot_d1d2d3) <- rep(c("D1", "D2", "D3"), 2)
d1d2d3 

prec_rates <- matrix(NA, nrow = 3, ncol = 2) 
prec_rates[1,] <- c(d1d2d3[1, 3], d1d2d3[4, 3])
prec_rates[2,] <- c(d1d2d3[2, 3], d1d2d3[5, 3])
prec_rates[3,] <- c(d1d2d3[3, 3], d1d2d3[6, 3])
colnames(prec_rates) <- c('STARS', 'HMM')
rownames(prec_rates) <- c('Data Set 1', 'Data Set 2', 'Data Set 3')
prec_rates

barplot(    # Plot
  prec_rates,
  col = colors()[c(23, 100, 12)],
  beside = TRUE,
  xlab = "Methods",
  ylab = "Average Precision Rate",
  ylim = c(0, 1),
  args.legend = list('topright', title = "Shift Types", cex = 0.8, box.col = 'grey'),
  border = "grey",
  legend = c('Mean', 'Variance', 'Mean and Variance')
)

# Std Deviation Rate Summary
data_sets_error <- rep(c(1, 2, 3), 2)
methods_stars_hmm_error <- c(rep('STARS', 3), rep('HMM', 3))
data_mat_error <- matrix(0, nrow = 6, ncol =3) %>% as.data.frame()
data_mat_error[1,3] <- best_data_1_STARS[4]
data_mat_error[2,3] <- best_data_2_STARS[4]
data_mat_error[3,3] <- best_data_3_STARS
data_mat_error[4,3] <- best_data_1_HMM[3]
data_mat_error[5,3] <- best_data_2_HMM[3]
data_mat_error[6,3] <- best_data_3_HMM[3]
data_mat_error[,2] <- as.factor(data_sets_error) 
data_mat_error[,1] <- methods_stars_hmm_error
d1d2d3_error <- as.data.frame(data_mat_error)
d1d2d3_error[,3] <- as.numeric(data_mat_error[,3])
colnames(d1d2d3_error) <- c('Method', 'Data Set', 'Std Error')
barplot_d1d2d3_error <- t(d1d2d3_error[,3])
colnames(barplot_d1d2d3_error) <- rep(c("D1", "D2", "D3"), 2)
d1d2d3_error 

se <- matrix(NA, nrow = 3, ncol = 2) 
se[1,] <- c(d1d2d3_error[1, 3], d1d2d3_error[4, 3])
se[2,] <- c(d1d2d3_error[2, 3], d1d2d3_error[5, 3])
se[3,] <- c(d1d2d3_error[3, 3], d1d2d3_error[6, 3])
colnames(se) <- c('STARS', 'HMM')
rownames(se) <- c('Data Set 1', 'Data Set 2', 'Data Set 3')
se

barplot(               # Plot
  se,
  col = colors()[c(23, 100, 12)],
  beside = TRUE,
  xlab = "Methods",
  ylab = "Standard Deviation",
  args.legend = list('topright', title = "Shift Types", cex = 0.8, box.col = 'grey'),
  border = "grey",
  legend = c('Mean', 'Variance', 'Mean and Variance')
)

# For the Lowest Standard Deviation ---------------------------------------------------------------------------------------------------------
# STARS
best_data_1_STARS <- avg_prec_stars_mean[which.min(avg_prec_stars_mean$`Std Error`),]
best_data_2_STARS <- avg_prec_stars_var[which.min(avg_prec_stars_var$`Std Error`),]
best_prec_d3_mean <- avg_prec_stars_var[which.min(avg_prec_stars_meanvar_mean$`Std Error`),]
best_prec_d3_var <- avg_prec_stars_var[which.min(avg_prec_stars_meanvar_var$`Std Error`),]
best_data_3_STARS <- (best_prec_d3_mean[3] + best_prec_d3_var[3])/2
best_stars <-data.frame(best_data_1_STARS[3], best_data_2_STARS[3], best_data_3_STARS)
best_stars
best_data_3_STARS_error <- (best_prec_d3_mean[4] + best_prec_d3_var[4])/2
best_stars_error <-data.frame(best_data_1_STARS[4], best_data_2_STARS[4], best_data_3_STARS_error)
best_stars_error

# HMM
best_data_1_HMM <- hmm_mean_prec$summary.stats.table[which.min(as.numeric(hmm_mean_prec$summary.stats.table[,3])),] %>% as.numeric()
best_data_2_HMM <- hmm_var_prec$summary.stats.table[which.min(as.numeric(hmm_var_prec$summary.stats.table[,3])),] %>% as.numeric()
best_data_3_HMM <- hmm_meanvar_prec$summary.stats.table[which.min(as.numeric(hmm_meanvar_prec$summary.stats.table[,3])),] %>% as.numeric()
best_data_hmm <- data.frame(best_data_1_HMM[2], best_data_2_HMM[2], best_data_3_HMM[2])
best_data_hmm
best_hmm_error <-data.frame(best_data_1_HMM[3], best_data_2_HMM[3], best_data_3_HMM[3])
best_hmm_error

# Average Precision Rate  Summary
data_sets <- rep(c(1, 2, 3), 2)
methods_stars_hmm <- c(rep('STARS', 3), rep('HMM', 3))
data_mat <- matrix(5, nrow = 6, ncol =3) %>% as.data.frame() 
data_mat[1,3] <- best_data_1_STARS[3]
data_mat[2,3] <- best_data_2_STARS[3]
data_mat[3,3] <- best_data_3_STARS
data_mat[4,3] <- best_data_1_HMM[2]
data_mat[5,3] <- best_data_2_HMM[2]
data_mat[6,3] <- best_data_3_HMM[2]
data_mat[,2] <- as.factor(data_sets) 
data_mat[,1] <- methods_stars_hmm
d1d2d3 <- as.data.frame(data_mat)
d1d2d3[,3] <- as.numeric(data_mat[,3])
colnames(d1d2d3) <- c('Method', 'Data Set', 'Precision')
barplot_d1d2d3 <- t(d1d2d3[,3])
colnames(barplot_d1d2d3) <- rep(c("D1", "D2", "D3"), 2)
d1d2d3 

prec_rates <- matrix(NA, nrow = 3, ncol = 2) 
prec_rates[1,] <- c(d1d2d3[1, 3], d1d2d3[4, 3])
prec_rates[2,] <- c(d1d2d3[2, 3], d1d2d3[5, 3])
prec_rates[3,] <- c(d1d2d3[3, 3], d1d2d3[6, 3])
colnames(prec_rates) <- c('STARS', 'HMM')
rownames(prec_rates) <- c('Data Set 1', 'Data Set 2', 'Data Set 3')
prec_rates

barplot(    # Plot
  prec_rates,
  col = colors()[c(23, 100, 12)],
  beside = TRUE,
  xlab = "Methods",
  ylab = "Average Precision Rate",
  ylim = c(0, 1),
  args.legend = list('topright', title = "Shift Types", cex = 0.8, box.col = 'grey'),
  border = "grey",
  legend = c('Mean', 'Variance', 'Mean and Variance')
)

# Std Deviation Rate Summary
data_sets_error <- rep(c(1, 2, 3), 2)
methods_stars_hmm_error <- c(rep('STARS', 3), rep('HMM', 3))
data_mat_error <- matrix(0, nrow = 6, ncol =3) %>% as.data.frame()
data_mat_error[1,3] <- best_data_1_STARS[4]
data_mat_error[2,3] <- best_data_2_STARS[4]
data_mat_error[3,3] <- best_data_3_STARS
data_mat_error[4,3] <- best_data_1_HMM[3]
data_mat_error[5,3] <- best_data_2_HMM[3]
data_mat_error[6,3] <- best_data_3_HMM[3]
data_mat_error[,2] <- as.factor(data_sets_error) 
data_mat_error[,1] <- methods_stars_hmm_error
d1d2d3_error <- as.data.frame(data_mat_error)
d1d2d3_error[,3] <- as.numeric(data_mat_error[,3])
colnames(d1d2d3_error) <- c('Method', 'Data Set', 'Std Error')
barplot_d1d2d3_error <- t(d1d2d3_error[,3])
colnames(barplot_d1d2d3_error) <- rep(c("D1", "D2", "D3"), 2)
d1d2d3_error 

se <- matrix(NA, nrow = 3, ncol = 2) 
se[1,] <- c(d1d2d3_error[1, 3], d1d2d3_error[4, 3])
se[2,] <- c(d1d2d3_error[2, 3], d1d2d3_error[5, 3])
se[3,] <- c(d1d2d3_error[3, 3], d1d2d3_error[6, 3])
colnames(se) <- c('STARS', 'HMM')
rownames(se) <- c('Data Set 1', 'Data Set 2', 'Data Set 3')
se

barplot(               # Plot
  se,
  col = colors()[c(23, 100, 12)],
  beside = TRUE,
  xlab = "Methods",
  ylab = "Standard Deviation",
  args.legend = list(x='topright', title = "Shift Types", cex = 0.8, box.col = 'grey'),
  border = "grey",
  legend = c('Mean', 'Variance', 'Mean and Variance')
)
