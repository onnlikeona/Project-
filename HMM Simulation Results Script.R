####### HMM RESULTS ########

#### Load Libraries ####

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
library(ggplot2)
library(ggthemes)


### Read in Excel Simulation Datasets ####

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "A1:CV1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
all_datasets <- read_excel_allsheets("Simulations.xlsx")


#### Applying HMM ####

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


#### Performance Metrics ####

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
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  
  for(k in 1:1){
    
    method_detections <- as.data.frame(.data) 
    
    for(j in 1:100){
      
    # Turnning RSIs into zeros and ones (for STARS method)
    # for (i in 1000) {
    #  ifelse(.data[[k]][i,j] == 0, method_detections[i, j] <- 0, method_detections[i, j] <- 1) 
    #}
    
    # Where the method detected shifts
    method_detections_index <- which(method_detections[, j] == 1)
    
    # Where the simulation detected a shift , i.e. at time t = 500
    simulated_shifts <- rep(0, 1000)
    simulated_shifts[500] <- 1
    simulated_shift_index <- which(simulated_shifts == 1)
    
    # Accommodating early/late shifts (within 10 units)
    if(method_detections[simulated_shift_index, j] == 1){ # if the method has the shift at t = 500, leave it
      
      method_detections[simulated_shift_index, j] <- 1
      
    } else if(method_detections[simulated_shift_index, j] != 1 & any(method_detections_index <= (simulated_shift_index + 10) & method_detections_index >= (simulated_shift_index - 10))){ 
      
      index_of_shifts_within_10_units <- method_detections_index[which(method_detections_index <= (simulated_shift_index + 10) & method_detections_index >= (simulated_shift_index - 10) )] 
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
   
    TP[j] <- sum(TP_count, na.rm = T)
    TN[j] <- sum(TN_count, na.rm = T)
    FP[j] <- sum(FP_count, na.rm = T)
    FN[j] <- sum(FN_count, na.rm = T)
    
    TPR[j] <- TP[j]/(TP[j] + FN[j])
    TNR[j] <- TN[j]/(FP[j] + TN[j])
    FPR[j] <- FP[j]/(FP[j] + TN[j]) 
    FNR[j] <- FN[j]/(FN[j] + TP[j])
    
 
 } # for j
    # Total number of detections made at time t, i.e. rowsums
    #otal_detections_at_time_t <- rowSums(method_detections)
    
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
                  Predictions = as.data.frame(.data),
                  Preds = method_detections)
    return(object_list)
    
} 

hmm_mean_confusion_metrics <- confusion_metrics(hmm_apply_mean)
hmm_var_confusion_metrics <- confusion_metrics(hmm_apply_var)

#### Plots ####

#### Plotting Number of Detections ####
#### For Change in Mean ####
plot_loc <- hmm_mean_confusion_metrics$total_detections_at_time_t
ggplot(plot_loc, aes(x = index(plot_loc), y = total_detections_at_time_t)) + 
  geom_bar(fill = "steelblue",stat = "identity") +
  ggtitle("HMM Total Detections Over time for \n Regime Shift in the Mean") +
  theme_bw() + 
  coord_fixed(ratio = 10) +
  scale_x_continuous(n.breaks=10,name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10)  + 
  theme(plot.title = element_text(size=26, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5), 
        axis.title.x = element_text(size = 25, 
                                    vjust = -4), 
        axis.title.y = element_text(size = 25, 
                                    vjust = 2),
        axis.text.x = element_text(size = 20,
                                   vjust = 0),
        axis.text.y = element_text(size = 20),
        aspect.ratio = 1, panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank())

#### For Change in Variance ####
plot_loc <- hmm_var_confusion_metrics$total_detections_at_time_t
ggplot(plot_loc, aes(x = index(plot_loc), y = total_detections_at_time_t)) + 
  geom_bar(fill = "steelblue",stat = "identity") +
  ggtitle("HMM Total Detections Over time for \n Regime Shift in the Variance") +
  theme_bw() + 
  coord_fixed(ratio = 10) +
  scale_x_continuous(n.breaks=10,name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10)  + 
  theme(plot.title = element_text(size=26, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5), 
        axis.title.x = element_text(size = 25, 
                                    vjust = -4), 
        axis.title.y = element_text(size = 25, 
                                    vjust = 2),
        axis.text.x = element_text(size = 20,
                                   vjust = 0),
        axis.text.y = element_text(size = 20),
        aspect.ratio = 1, panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


mean_barplot_vect_hmm <- c( sum(hmm_mean_confusion_metrics$TP)/100, length(which(hmm_mean_confusion_metrics$FP >= 1))/100 )
var_barplot_vect_hmm <- c( sum(hmm_var_confusion_metrics$TP)/100, length(which(hmm_var_confusion_metrics$FP >= 1))/100 )
full_hmm_df <- data.frame(value = c(mean_barplot_vect_hmm, 
                                    var_barplot_vect_hmm), 
                           type = c(rep('mean', 2),
                                   rep('variance', 2)),
                          value_type = rep(c('Correct','Incorrect'), 2))
ggplot(full_hmm_df, aes(x = type, y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'HMM Proportion of Correctly and \n Incorrectly Detected Shifts', 
       x = 'Regime Shift Type',
       y = 'Proportion') + 
  theme_bw()  +
  theme(plot.title = element_text(size=26, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 25, vjust = -4), 
        axis.title.y = element_text(size = 25, vjust = 2), 
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 20),
        aspect.ratio = 1.2, 
        panel.grid.major = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.key.width= unit(2, 'cm'), 
        legend.title = element_text(size=15, face = 'bold'),
        legend.text = element_text(size=14),
        strip.background = element_rect(colour="black", fill="light grey", size = 0.55),
        strip.text = element_text(size = 12),
        plot.background = element_rect(size = 2)) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'METRIC') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 6,vjust=-0.5)
#### Percentage correctly/incorrectly detected ####
# mean_barplot_vect <- c( sum(hmm_mean_confusion_metrics$TP)/100, length(which(hmm_mean_confusion_metrics$FP >= 1))/100 )
# var_barplot_vect <- c( sum(hmm_var_confusion_metrics$TP)/100, length(which(hmm_var_confusion_metrics$FP >= 1))/100 )
# plot_barplot_1 <- matrix(c(mean_barplot_vect, var_barplot_vect), nrow = 2, byrow = F)
# colnames(plot_barplot_1) <- c("Mean", "Variance")
# rownames(plot_barplot_1) <- c("Proportion Correctly Detected", "Proportion Incorrectly Detected")
# barplot(plot_barplot_1, 
#         col=colors()[c(23,12)], 
#         border="white", 
#         font.axis=2, 
#         beside=T, 
#         legend=rownames(plot_barplot_1), 
#         xlab="Regime Shift Type", 
#         font.lab=2, ylim = c(0,1))
# title("HMM Proportion of Correctly and Incorrectly Detected for Simulations")


