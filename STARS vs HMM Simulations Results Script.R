###### STARS VS HMM ########

rm(list=ls())
library(tidyverse)
library(reshape)
library(caret)
library(astsa)
library(dplyr)
library(timeSeries)
library(forecast)
library(quantmod)
library(data.table)
library(depmixS4)
library(readxl)
library(ggplot2)
library(ggthemes)


#### Run STARS ####

# Read in Excel Simulation Datasets 

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "B1:CW1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
stars_mean_data <- read_excel_allsheets("STARS Mean Results.xlsm")
stars_variance_data <- read_excel_allsheets("STARS Variance Results 2.xlsm")




### Performance Metrics ###

confusion_metrics <- function(.data, cutoff_length){
  
  # Prepare objects 
  
  
  total_detections_of_x <- matrix(NA, nrow = 100, ncol = 3)
  colnames(total_detections_of_x) <- c("0.01", "0.05", "0.1")
  
  total_detections_at_time_t <- matrix(NA, nrow = 1000, ncol = 3)
  colnames(total_detections_at_time_t) <- c("0.01", "0.05", "0.1")
  
  FP <- matrix(NA, nrow = 100, ncol = 3) 
  colnames(FP) <- c("0.01", "0.05", "0.1")
  
  FN <- matrix(NA, nrow = 100, ncol = 3)
  colnames(FN) <- c("0.01", "0.05", "0.1")
  
  TP <- matrix(NA, nrow = 100, ncol = 3)
  colnames(TP) <- c("0.01", "0.05", "0.1")
  
  TN <- matrix(NA, nrow = 100, ncol = 3)
  colnames(TN) <- c("0.01", "0.05", "0.1")
  
  FPR <- matrix(NA, nrow = 100, ncol = 3)
  colnames(FPR) <- c("0.01", "0.05", "0.1")
  
  FNR <- matrix(NA, nrow = 100, ncol = 3)
  colnames(FNR) <- c("0.01", "0.05", "0.1")
  
  TPR <- matrix(NA, nrow = 100, ncol = 3)
  colnames(TPR) <- c("0.01", "0.05", "0.1")
  
  TNR <- matrix(NA, nrow = 100, ncol = 3)
  colnames(TNR) <- c("0.01", "0.05", "0.1")
  
  LOC <- matrix(NA, nrow = 1000, ncol = 3)
  colnames(LOC) <- c("0.01", "0.05", "0.1")
  
  method_detections_list <- list()
  
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  if(cutoff_length == 100){
    beg <- 1
    end <- 3
  } else if(cutoff_length == 250){
    beg <- 4
    end <- 6
  } else if(cutoff_length == 500){
    beg <- 7
    end <- 9
  }
  
  for(k in beg:end){
    
    method_detections <- matrix(NA,nrow = 1000,ncol = 100)
    
    for(j in 1:100){
      
      # Turning RSIs into zeros and ones (for STARS method)
      for (t in 1:1000){
        if(.data[[k]][t,j] == 0){
          method_detections[t,j] <- 0
        } else {
          method_detections[t,j] <- 1
        }
      }
      
      # Where the method detected shifts
      method_detections_index <- which(method_detections[, j] == 1)
      
      # Where the simulation detected a shift , i.e. at time t = 500
      simulated_shifts <- rep(0, 1000)
      simulated_shifts[500] <- 1
      simulated_shift_index <- which(simulated_shifts == 1)
      
      # Accommodating early/late shifts (within 10 units)
      if(method_detections[simulated_shift_index, j] == 1){ # if the method has the shift at t = 500, leave it
        
        method_detections[simulated_shift_index, j] <- 1
        
      } else if(method_detections[simulated_shift_index, j] != 1 & any(method_detections_index <= abs(simulated_shift_index + 10) & method_detections_index >= abs(simulated_shift_index - 10))){ 
        
        index_of_shifts_within_10_units <- method_detections_index[which(method_detections_index <= (simulated_shift_index + 10)|method_detections_index >= (simulated_shift_index - 10) )] # choose the observations that are nearest to the simulated shift, within a 10 unit threshold. 
        index_of_nearest_shifts <- index_of_shifts_within_10_units[which(abs(index_of_shifts_within_10_units - simulated_shift_index) == min(abs(index_of_shifts_within_10_units - simulated_shift_index)))]
        index_of_nearest_one <- sample.vec(index_of_nearest_shifts, 1) # if there were two shifts that fulfill the above condition, for example, one at t=501 and 502, randomly select one of those shifts.
        method_detections[index_of_nearest_one, j] <- 0
        method_detections[simulated_shift_index, j] <- 1 # now consider the nearest shift within 10 units as a shift.
        
      } 
      
      # Total number of shifts made by the method for each time series simulation x
      if(beg == 1){
        total_detections_of_x[j, k] <- sum(method_detections[, j])
      } else if(beg == 4){
        total_detections_of_x[j, k-3] <- sum(method_detections[, j])
      } else if(beg == 7){
        total_detections_of_x[j, k-6] <- sum(method_detections[, j])
      }
      
      
      # Classification Metrics - evaluating each observation to see if it's a true positive, negative, etc.
      
      # Prepare for storing                    
      TP_count <- c()
      TN_count <- c()
      FP_count <- c()
      FN_count <- c()
      
      for(i in 1:1000){
        
        if(method_detections[i,j] == 1 & simulated_shifts[i] == 1){
          TP_count[i] <- 1
        } else if(method_detections[i,j] == 0 & simulated_shifts[i] == 0){
          TN_count[i] <- 1
        } else if(method_detections[i,j] == 1 & simulated_shifts[i] == 0){
          FP_count[i] <- 1
        } else if(method_detections[i,j] == 0 & simulated_shifts[i] == 1){
          FN_count[i] <- 1
        }
        
      }
      
      if(beg == 1){
        TP[j, k] <- sum(TP_count, na.rm = T)
        TN[j, k] <- sum(TN_count, na.rm = T)
        FP[j, k] <- sum(FP_count, na.rm = T)
        FN[j, k] <- sum(FN_count, na.rm = T)
      } else if(beg == 4){
        TP[j, k-3] <- sum(TP_count, na.rm = T)
        TN[j, k-3] <- sum(TN_count, na.rm = T)
        FP[j, k-3] <- sum(FP_count, na.rm = T)
        FN[j, k-3] <- sum(FN_count, na.rm = T)
      } else if(beg == 7){
        TP[j, k-6] <- sum(TP_count, na.rm = T)
        TN[j, k-6] <- sum(TN_count, na.rm = T)
        FP[j, k-6] <- sum(FP_count, na.rm = T)
        FN[j, k-6] <- sum(FN_count, na.rm = T)
      }
      
      # TPR[j, k] <- TP[j, k]/(TP[j, k] + FN[j, k])
      # TNR[j, k] <- TN[j, k]/(FP[j, k] + TN[j, k])
      # FPR[j, k] <- FP[j, k]/(FP[j, k] + TN[j, k]) 
      # FNR[j, k] <- FN[j, k]/(FN[j, k] + TP[j, k])
      if(beg == 1){
        total_detections_at_time_t[, k] <- rowSums(method_detections)
      } else if(beg == 4){
        total_detections_at_time_t[, k-3] <- rowSums(method_detections)
      } else if(beg == 7){
        total_detections_at_time_t[, k-6] <- rowSums(method_detections)
      }
      
    } # for j
    
    if(beg == 1){
      method_detections_list[[k]] <- as.data.frame(method_detections)
    } else if(beg == 4){
      method_detections_list[[k-3]] <- as.data.frame(method_detections)
    } else if(beg == 7){
      method_detections_list[[k-6]] <- as.data.frame(method_detections)
    }
    
  } # for k
  
  # Total number of detections made at time t, i.e. rowsums
  #total_detections_at_time_t[, k] <- rowSums(method_detections)
  
  # Prepare Objects for Output list
  total_detections_at_time_t_df <- data.frame(total_detections_at_time_t, row.names = as.character(seq(1,1000,1)))
  total_detections_of_x_df <- data.frame(total_detections_of_x)
  FP_df <- data.frame(FP)
  FN_df <- data.frame(FN)
  FNR_df <- data.frame(FNR)
  TNR_df <- data.frame(TNR)
  TPR_df <- data.frame(TPR)
  TP_df <- data.frame(TP)
  TN_df <- data.frame(TN)
  FPR_df <- data.frame(FPR)
  object_list <- list(total_detections_of_x = total_detections_of_x_df,
                      total_detections_at_time_t = total_detections_at_time_t_df,
                      FP = FP_df,
                      FN = FN_df,
                      TN = TN_df,
                      TP = TP_df,
                      # FPR = FPR_df,
                      # FNR = FNR_df,
                      # TNR = TNR_df,
                      # TPR = TPR_df,
                      Predictions = .data, 
                      Method_Predictions_Encoded = method_detections_list)
  return(object_list)
  
} 

stars_mean_confusion_metrics_100 <- confusion_metrics(stars_mean_data, cutoff_length = 100)
stars_var_confusion_metrics_100 <- confusion_metrics(stars_variance_data, cutoff_length = 100)
stars_mean_confusion_metrics_250 <- confusion_metrics(stars_mean_data, cutoff_length = 250)
stars_var_confusion_metrics_250 <- confusion_metrics(stars_variance_data, cutoff_length = 250)
stars_mean_confusion_metrics_500 <- confusion_metrics(stars_mean_data, cutoff_length = 500)
stars_var_confusion_metrics_500 <- confusion_metrics(stars_variance_data, cutoff_length = 500)


#### Run HMM ####

### Read in Excel Simulation Datasets ###

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "A1:CV1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
all_datasets <- read_excel_allsheets("Simulations.xlsx")


### Applying HMM ###

hmm_detections_func <- function(posterior.states){ 
  # outputs a vector indicating where shifts occurred
  # 1 indicates a shift and 0 not a shift
  
  n <- length(posterior.states)
  pred_shifts <- c()
  pred_shifts[1] <- 0
  pred_shifts[n] <- 0
  for(i in 2:(n-1)){
    ifelse(posterior.states[i] - posterior.states[i-1] == 0, pred_shifts[i] <- 0, pred_shifts[i] <- 1)
  }
  return(pred_shifts)
}

hmm_apply <- function(.data, type){
  
  hmm_detections <- matrix(data=NA, nrow = 1000, ncol = 100)
  detections_list <- list()
  states <- c(2)
  data_index <- c(1, 2, 3)
  if (type == "mean"){
    data_index = 1
  } else if (type == "variance") {
    data_index = 2
  } else if (type == "mean and variance") {
    data_index = 3
  }
  df <- data.frame(matrix(unlist(all_datasets[[data_index]]), ncol=length(all_datasets[[data_index]]), byrow=F)) %>%
    mutate_all( function(x) as.numeric(x))
  number_of_series <- dim(df)[2]
  
  for(i in 1:number_of_series){
    hmm <- depmix(df[,i] ~ 1, family = gaussian(), nstates = states, data=data.frame(returns=df[,i]))
    hmmfit <- fit(hmm, verbose = FALSE)
    post_probs <- posterior(hmmfit)
    hmm_detections[,i] <- hmm_detections_func(post_probs[,1])
  }
  
  hmm_detections <- data.frame(hmm_detections, check.names = TRUE)
  detections_list[[1]] <- hmm_detections
  
  return(detections_list)
  
}

hmm_apply_mean <- hmm_apply(all_datasets, "mean")
hmm_apply_var <- hmm_apply(all_datasets,"variance")


### Performance Metrics ###

confusion_metrics <- function(.data){
  
  # Prepare objects 
  unknown_shifts <- c()
  total_detections_of_x <- c()
  FP <- matrix(NA, nrow = 100, ncol = 1)
  FN <- matrix(NA, nrow = 100, ncol = 1)
  TP <- matrix(NA, nrow = 100, ncol = 1)
  TN <- matrix(NA, nrow = 100, ncol = 1)
  FPR <- matrix(NA, nrow = 100, ncol = 1)
  FNR <- matrix(NA, nrow = 100, ncol = 1)
  TPR <- matrix(NA, nrow = 100, ncol = 1)
  TNR <- matrix(NA, nrow = 100, ncol = 1)
  LOC <- matrix(NA, nrow = 1000, ncol = 1)
  method_detections <- matrix(NA,nrow = 1000,ncol = 100)
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  for(k in 1:1){
    
    for(j in 1:100){
      
      # Turnning RSIs into zeros and ones (for STARS method)
      # for (i in 1000) {
      #  ifelse(.data[[k]][i,j] == 0, method_detections[i, j] <- 0, method_detections[i, j] <- 1) 
      #}
      
      # Where the method detected shifts
      method_detections <- as.data.frame(.data)
      method_detections_index <- which(method_detections[, j] == 1)
      
      # Where the simulation detected a shift , i.e. at time t = 500
      simulated_shifts <- rep(0, 1000)
      simulated_shifts[500] <- 1
      simulated_shift_index <- which(simulated_shifts == 1)
      
      # Accommodating early/late shifts (within 10 units)
      if(method_detections[simulated_shift_index, j] == 1){ # if the method has the shift at t = 500, leave it
        
        method_detections[simulated_shift_index, j] <- 1
        
      } else if(method_detections[simulated_shift_index, j] != 1 & any(method_detections_index <= abs(simulated_shift_index + 10) & method_detections_index >= abs(simulated_shift_index - 10))){ 
        
        index_of_shifts_within_10_units <- method_detections_index[which(method_detections_index <= (simulated_shift_index + 10)|method_detections_index >= (simulated_shift_index - 10) )] # choose the observations that are nearest to the simulated shift, within a 10 unit threshold. 
        index_of_nearest_shifts <- index_of_shifts_within_10_units[which(abs(index_of_shifts_within_10_units - simulated_shift_index) == min(abs(index_of_shifts_within_10_units - simulated_shift_index)))]
        index_of_nearest_one <- sample.vec(index_of_nearest_shifts, 1) # if there were two shifts that fulfill the above condition, for example, one at t=501 and 502, randomly select one of those shifts.
        method_detections[index_of_nearest_one, j] <- 0
        method_detections[simulated_shift_index, j] <- 1 # now consider the nearest shift within 10 units as a shift.
        
      } 
      
      # Total number of shifts made by the method for each time series simulation x
      total_detections_of_x <- sum(method_detections[, j])
      
      
      # Classification Metrics - evaluating each observation to see if it's a true positive, negative, etc.
      
      # Prepare for storing 
      TP_count <- c()
      TN_count <- c()
      FP_count <- c()
      FN_count <- c()
      
      for(i in 1:1000){
        ifelse(method_detections[i,j] == 1 & simulated_shifts[i] == 1, TP_count[i] <- 1, TP_count[i] <- 0) 
        ifelse(method_detections[i,j] == 0 & simulated_shifts[i] == 0, TN_count[i] <- 1, TN_count[i] <- 0) 
        ifelse(method_detections[i,j] == 1 & simulated_shifts[i] == 0, FP_count[i] <- 1, FP_count[i] <- 0) 
        ifelse(method_detections[i,j] == 0 & simulated_shifts[i] == 1, FN_count[i] <- 1, FN_count[i] <- 0) 
      }
      
      TP[j] <- sum(TP_count)
      TN[j] <- sum(TN_count)
      FP[j] <- sum(FP_count)
      FN[j] <- sum(FN_count)
      
      TPR[j] <- TP[j]/(TP[j] + FN[j])
      TNR[j] <- TN[j]/(FP[j] + TN[j])
      FPR[j] <- FP[j]/(FP[j] + TN[j]) 
      FNR[j] <- FN[j]/(FN[j] + TP[j])
      
      
    } # for j
    
  } # for k
  
  # Total number of detections made at time t, i.e. rowsums
  total_detections_at_time_t <- rowSums(method_detections)
  
  # Prepare Objects for Output list
  total_detections_at_time_t_df <- data.frame(total_detections_at_time_t, row.names = as.character(seq(1,1000,1)))
  total_detections_of_x_df <- data.frame(total_detections_of_x)
  FP_df <- data.frame(FP)
  FN_df <- data.frame(FN)
  FNR_df <- data.frame(FNR)
  TNR_df <- data.frame(TNR)
  TPR_df <- data.frame(TPR)
  TP_df <- data.frame(TP)
  TN_df <- data.frame(TN)
  FPR_df <- data.frame(FPR)
  object_list <- list(total_detections_of_x = total_detections_of_x_df,
                      total_detections_at_time_t = total_detections_at_time_t_df,
                      FP = FP_df,
                      FN = FN_df,
                      TN = TN_df,
                      TP = TP_df,
                      FPR = FPR_df,
                      FNR = FNR_df,
                      TNR = TNR_df,
                      TPR = TPR_df,
                      Predictions = .data)
  return(object_list)
  
} 

hmm_mean_confusion_metrics <- confusion_metrics(hmm_apply_mean)
hmm_var_confusion_metrics <- confusion_metrics(hmm_apply_var)


mean_barplot_vect_stars <- c( sum(stars_mean_confusion_metrics_100$TP[,1])/100, length(which(stars_mean_confusion_metrics_100$FP[,1] >= 1))/100 )
var_barplot_vect_stars <- c( sum(stars_var_confusion_metrics_100$TP[,1])/100, length(which(stars_var_confusion_metrics_100$FP[,1] >= 1))/100 )
mean_barplot_vect_hmm <- c( sum(hmm_mean_confusion_metrics$TP)/100, length(which(hmm_mean_confusion_metrics$FP >= 1))/100 )
var_barplot_vect_hmm <- c( sum(hmm_var_confusion_metrics$TP)/100, length(which(hmm_var_confusion_metrics$FP >= 1))/100 )

full_stars_hmm_df <- data.frame(value = c(mean_barplot_vect_stars, 
                                    var_barplot_vect_stars,
                                    mean_barplot_vect_hmm,
                                    var_barplot_vect_hmm), 
                          method = c(rep('STARS', 4),
                                    rep('HMM', 4)),
                          type = c(rep('mean', 2),
                                   rep('variance', 2),
                                   rep('mean', 2),
                                   rep('variance', 2)),
                          value_type = rep(c('Correct','Incorrect'), 4))

ggplot(full_stars_hmm_df, aes(x = factor(method), y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'STARS vs HMM Proportion of Simulations \n that Correctly/Incorrectly Detected Regime Shifts', 
       x = 'Significance Level',
       y = 'Proportion') + 
  theme_bw() +
  facet_grid( ~ type) +
  theme(plot.title = element_text(size=26, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 25, vjust = -4), 
        axis.title.y = element_text(size = 25, vjust = 2), 
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        aspect.ratio = 1.67, 
        panel.grid.major = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.key.width= unit(2, 'cm'), 
        legend.title = element_text(size=15, face = 'bold'),
        legend.text = element_text(size=14),
        strip.background = element_rect(colour="black", fill="light grey", size = 0.55),
        strip.text = element_text(size = 27),
        plot.background = element_rect(size = 2)) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'METRIC') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 6,vjust=-0.5)
