####### STARS RESULTS ########

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
library(reshape)
library(ggthemes)


### Read in Excel Simulation Datasets ####

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "B1:CW1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
stars_mean_data <- read_excel_allsheets("STARS Mean Results.xlsm")
stars_variance_data <- read_excel_allsheets("STARS Variance Results.xlsm")




#### Performance Metrics ####

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
        
      } else if(method_detections[simulated_shift_index, j] != 1 & any(method_detections_index <= (simulated_shift_index + 10) & method_detections_index >= (simulated_shift_index - 10))){ 
        
        index_of_shifts_within_10_units <- method_detections_index[which(method_detections_index <= (simulated_shift_index + 10) & method_detections_index >= (simulated_shift_index - 10) )] # choose the observations that are nearest to the simulated shift, within a 10 unit threshold. 
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



################################ PLOTS ###################################

###################### ggplot Bar Graphs ########################

par(mfrow = c(2, 1))
#### Cutoff Length = 100 ####

mean_barplot_vect_1 <- c( sum(stars_mean_confusion_metrics_100$TP[,1])/100, length(which(stars_mean_confusion_metrics_100$FP[,1] >= 1))/100 )
var_barplot_vect_1 <- c( sum(stars_var_confusion_metrics_100$TP[,1])/100, length(which(stars_var_confusion_metrics_100$FP[,1] >= 1))/100 )
mean_barplot_vect_2 <- c( sum(stars_mean_confusion_metrics_100$TP[,2])/100, length(which(stars_mean_confusion_metrics_100$FP[,2] >= 1))/100 )
var_barplot_vect_2 <- c( sum(stars_var_confusion_metrics_100$TP[,2])/100, length(which(stars_var_confusion_metrics_100$FP[,2] >= 1))/100 )
mean_barplot_vect_3 <- c( sum(stars_mean_confusion_metrics_100$TP[,3])/100, length(which(stars_mean_confusion_metrics_100$FP[,3] >= 1))/100 )
var_barplot_vect_3 <- c( sum(stars_var_confusion_metrics_100$TP[,3])/100, length(which(stars_var_confusion_metrics_100$FP[,3] >= 1))/100 )

full_df_100 <- data.frame(value = c(mean_barplot_vect_1, 
                                mean_barplot_vect_2,
                                mean_barplot_vect_3,
                                var_barplot_vect_1,
                                var_barplot_vect_2, 
                                var_barplot_vect_3), 
                      alpha = c(rep('0.01', 2),
                                rep('0.05', 2),
                                rep('0.10', 2),
                                rep('0.01', 2),
                                rep('0.05', 2),
                                rep('0.10', 2)),
                      type = c(rep('mean', 6),
                               rep('variance', 6)),
                      value_type = rep(c('Correct','Incorrect'), 6))
ggplot(full_df_100, aes(x = factor(alpha), y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'STARS Proportion of Correctly and Incorrectly Detected Shifts', 
       x = 'Regime Shift Type',
       y = 'Proportion') + 
  theme_bw() +
  facet_wrap( ~ type ) +
  theme(plot.title = element_text(size=11, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 14, vjust = -4), 
        axis.title.y = element_text(size = 14, vjust = 2), 
        aspect.ratio = 1.50, 
        panel.grid.major = element_blank()) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'Metric') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 4,vjust=-0.5)

#### Cutoff Length = 250 ####

mean_barplot_vect_1 <- c( sum(stars_mean_confusion_metrics_250$TP[,1])/100, length(which(stars_mean_confusion_metrics_250$FP[,1] >= 1))/100 )
var_barplot_vect_1 <- c( sum(stars_var_confusion_metrics_250$TP[,1])/100, length(which(stars_var_confusion_metrics_250$FP[,1] >= 1))/100 )
mean_barplot_vect_2 <- c( sum(stars_mean_confusion_metrics_250$TP[,2])/100, length(which(stars_mean_confusion_metrics_250$FP[,2] >= 1))/100 )
var_barplot_vect_2 <- c( sum(stars_var_confusion_metrics_250$TP[,2])/100, length(which(stars_var_confusion_metrics_250$FP[,2] >= 1))/100 )
mean_barplot_vect_3 <- c( sum(stars_mean_confusion_metrics_250$TP[,3])/100, length(which(stars_mean_confusion_metrics_250$FP[,3] >= 1))/100 )
var_barplot_vect_3 <- c( sum(stars_var_confusion_metrics_250$TP[,3])/100, length(which(stars_var_confusion_metrics_250$FP[,3] >= 1))/100 )

full_df_250 <- data.frame(value = c(mean_barplot_vect_1, 
                                    mean_barplot_vect_2,
                                    mean_barplot_vect_3,
                                    var_barplot_vect_1,
                                    var_barplot_vect_2, 
                                    var_barplot_vect_3), 
                          alpha = c(rep('0.01', 2),
                                    rep('0.05', 2),
                                    rep('0.10', 2),
                                    rep('0.01', 2),
                                    rep('0.05', 2),
                                    rep('0.10', 2)),
                          type = c(rep('mean', 6),
                                   rep('variance', 6)),
                          value_type = rep(c('Correct','Incorrect'), 6))
ggplot(full_df_250, aes(x = factor(alpha), y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'STARS Proportion of Correctly and Incorrectly Detected Shifts', 
       x = 'Regime Shift Type',
       y = 'Proportion') + 
  theme_bw() +
  facet_wrap( ~ type ) +
  theme(plot.title = element_text(size=11, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 14, vjust = -4), 
        axis.title.y = element_text(size = 14, vjust = 2), 
        aspect.ratio = 1.50, 
        panel.grid.major = element_blank()) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'Metric') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 4,vjust=-0.5)

#### Cutoff Length = 500 ####

mean_barplot_vect_1 <- c( sum(stars_mean_confusion_metrics_500$TP[,1])/100, length(which(stars_mean_confusion_metrics_500$FP[,1] >= 1))/100 )
var_barplot_vect_1 <- c( sum(stars_var_confusion_metrics_500$TP[,1])/100, length(which(stars_var_confusion_metrics_500$FP[,1] >= 1))/100 )
mean_barplot_vect_2 <- c( sum(stars_mean_confusion_metrics_500$TP[,2])/100, length(which(stars_mean_confusion_metrics_500$FP[,2] >= 1))/100 )
var_barplot_vect_2 <- c( sum(stars_var_confusion_metrics_500$TP[,2])/100, length(which(stars_var_confusion_metrics_500$FP[,2] >= 1))/100 )
mean_barplot_vect_3 <- c( sum(stars_mean_confusion_metrics_500$TP[,3])/100, length(which(stars_mean_confusion_metrics_500$FP[,3] >= 1))/100 )
var_barplot_vect_3 <- c( sum(stars_var_confusion_metrics_500$TP[,3])/100, length(which(stars_var_confusion_metrics_500$FP[,3] >= 1))/100 )

full_df_500 <- data.frame(value = c(mean_barplot_vect_1, 
                                    mean_barplot_vect_2,
                                    mean_barplot_vect_3,
                                    var_barplot_vect_1,
                                    var_barplot_vect_2, 
                                    var_barplot_vect_3), 
                          alpha = c(rep('0.01', 2),
                                    rep('0.05', 2),
                                    rep('0.10', 2),
                                    rep('0.01', 2),
                                    rep('0.05', 2),
                                    rep('0.10', 2)),
                          type = c(rep('mean', 6),
                                   rep('variance', 6)),
                          value_type = rep(c('Correct','Incorrect'), 6))
ggplot(full_df_500, aes(x = factor(alpha), y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  labs(title = 'STARS Proportion of Correctly and Incorrectly Detected Shifts', 
       x = 'Regime Shift Type',
       y = 'Proportion') + 
  theme_bw() +
  facet_wrap( ~ type ) +
  theme(plot.title = element_text(size=11, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 14, vjust = -4), 
        axis.title.y = element_text(size = 14, vjust = 2), 
        aspect.ratio = 1.50, 
        panel.grid.major = element_blank()) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'Metric') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 4,vjust=-0.5)


complete_df <- rbind(full_df_100, full_df_250, full_df_500)
complete_df$cutoff <- c(rep('l = 100', 12),
                               rep('l = 250', 12), 
                               rep('l = 500', 12))


ggplot(complete_df, aes(x = factor(alpha), y = value, fill = value_type)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  ggtitle('STARS\' Proportion of Simulations that \n Correctly and Incorrectly Detected Regime Shifts') +
  labs(x = 'Significance Level',
       y = 'Proportion') + 
  theme_bw()  +
  facet_grid(type ~ cutoff) +
  theme(plot.title = element_text(size=26, 
                                  face = "bold", 
                                  vjust = 5, 
                                  hjust = 0.5,), 
        axis.title.x = element_text(size = 25, vjust = -4), 
        axis.title.y = element_text(size = 25, vjust = 2),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        aspect.ratio = 1.4, 
        panel.grid.major = element_blank(),
        legend.key.size = unit(1.35, 'cm'),
        legend.key.width= unit(1.35, 'cm'), 
        legend.title = element_text(size=15, face = 'bold'),
        legend.text = element_text(size=15),
        strip.background = element_rect(colour="black", fill="light grey", size = 0.55),
        strip.text = element_text(size = 23),
        plot.background = element_rect(size = 2)) +
  scale_fill_discrete(
    labels = c("Correct", "Incorrect"), name = 'METRIC') + guides(colour = FALSE) +
  geom_text(aes(label = value),
            position=position_dodge(width=0.9), size = 4,vjust=-0.5)


# theme(plot.title = element_text(size=16, 
#                                 face = "bold", 
#                                 vjust = 5, 
#                                 hjust = 0.5,), 
#       axis.title.x = element_text(size = 14, vjust = -4), 
#       axis.title.y = element_text(size = 14, vjust = 2), 
#       aspect.ratio = 1.34, 
#       panel.grid.major = element_blank()) +
#   scale_fill_discrete(
#     labels = c("Correct", "Incorrect"), name = 'Metric') + guides(colour = FALSE) +
#   geom_text(aes(label = value),
#             position=position_dodge(width=0.9), size = 4,vjust=-0.5)
#### STARS DBN ####

stars_gg <- data.frame(Time =  seq(1, 1000, 1),
                          freq =  stars_mean_confusion_metrics_100$total_detections_at_time_t$X0.01)

ggplot(stars_gg, aes(x = Time, y = freq)) +
  geom_bar(fill = "steelblue",stat = "identity") +
  ggtitle("STARS' Total Detections Over time for \n Regime Shift in the Mean (Cutoff Length = 100, Significance Level = 0.01)") +
  theme_bw() + 
  coord_fixed(ratio = 10) +
  scale_x_continuous(n.breaks=10,name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10)  + 
  theme(plot.title = element_text(size=20, 
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





# Full for ggplot

# plot_loc_gg <- data.frame(Time =  rep(seq(1, 1000, 1), 3),
#                           freq =  c(stars_mean_confusion_metrics_100$total_detections_at_time_t[,1],
#                                     stars_mean_confusion_metrics_100$total_detections_at_time_t[,2],
#                                     stars_mean_confusion_metrics_100$total_detections_at_time_t[,3]),
#                           alpha = c(rep('0.01', 1000),
#                                     rep('0.05', 1000),
#                                     rep('0.10', 1000)))
# ggplot(plot_loc_gg, aes(x = Time, y = freq, fill = factor(alpha))) +
#   geom_bar(stat = 'identity', position = 'dodge')


# #### Plotting Number of Detections ####
#   #### For Change in Mean ####
# 
# 
# plot_loc <- data.frame(Time = seq(1, 1000, 1), 
#                        "0.01" = stars_mean_confusion_metrics$total_detections_at_time_t[,1],
#                        "0.05" = stars_mean_confusion_metrics$total_detections_at_time_t[,2],
#                        "0.10" = stars_mean_confusion_metrics$total_detections_at_time_t[,3])
# molten <- melt(plot_loc, id.vars = "Time")
# ggplot(plot_loc, aes(x = Time, y = value, colour = variable, fill = variable)) + 
#   geom_bar(stat = "identity") +
#   ggtitle("STARS Total Detections Over time for Shift in Mean") +
#   theme_bw() + 
#   coord_fixed(ratio = 10) +
#   scale_x_continuous(n.breaks=10,name="Time") +
#   scale_y_continuous(name="Frequency",n.breaks = 10) + 
#   theme(plot.title = element_text(size=11, 
#                                   face = "bold", 
#                                   vjust = 5, 
#                                   hjust = 0.5,), 
#                                   axis.title.x = element_text(size = 13, vjust = -4), 
#                                   axis.title.y = element_text(size = 13, vjust = 2), 
#                                   aspect.ratio = 0.75, 
#                                   panel.grid.major = element_blank(),
#                                   panel.grid.minor = element_blank()) +
#   scale_fill_discrete(
#     labels = c("0.01", "0.05", "0.10"), name = 'Significance Level') + guides(colour = FALSE)
# ggplot(as.data.frame(plot_loc[,2]), aes(x = Time, y = X0.01)) + geom_bar()
# 
#                                   
                                       

#   #### For Change in Variance ####
# plot_loc <- data.frame(Time = seq(1, 1000, 1), 
#                        "0.01" = stars_var_confusion_metrics$total_detections_at_time_t[,1],
#                        "0.05" = stars_var_confusion_metrics$total_detections_at_time_t[,2],
#                        "0.10" = stars_var_confusion_metrics$total_detections_at_time_t[,3])
# molten <- melt(plot_loc, id.vars = "Time")
# ggplot(molten, aes(x = Time, y = value, colour = variable, fill = variable)) + 
#   geom_bar(stat = "identity") +
#   ggtitle("STARS Total Detections Over time for Shift in Variance") +
#   theme_bw() + 
#   coord_fixed(ratio = 10) +
#   scale_x_continuous(n.breaks=10,name="Time") +
#   scale_y_continuous(name="Frequency",n.breaks = 10) + 
#   theme(plot.title = element_text(size=11, 
#                                   face = "bold", 
#                                   vjust = 5, 
#                                   hjust = 0.5,), 
#         axis.title.x = element_text(size = 13, vjust = -4), 
#         axis.title.y = element_text(size = 13, vjust = 2), 
#         aspect.ratio = 0.75, 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_fill_discrete(
#     labels = c("0.01", "0.05", "0.10"), name = 'Significance Level') + guides(colour = FALSE)
# 
# 
#   ##### Percentage correctly/incorrectly detected ####
# #### Significance level of 0.01 ####
# mean_barplot_vect_1 <- c( sum(stars_mean_confusion_metrics$TP[,1])/100, length(which(stars_mean_confusion_metrics$FP[,1] >= 1))/100 )
# var_barplot_vect_1 <- c( sum(stars_var_confusion_metrics$TP[,1])/100, length(which(stars_var_confusion_metrics$FP[,1] >= 1))/100 )
# plot_barplot_1 <- matrix(c(mean_barplot_vect, var_barplot_vect), nrow = 2, byrow = F)
# colnames(plot_barplot_1) <- c("Mean", "Variance")
# rownames(plot_barplot_1) <- c("Proportion Correctly Detected", "Proportion Incorrectly Detected")
# bp_1 <- barplot(plot_barplot_1, 
#         col=colors()[c(23,12)], 
#         border="white", 
#         font.axis=2, 
#         beside=T, 
#         legend=rownames(plot_barplot_1), 
#         xlab="Regime Shift Type", 
#         font.lab=2, 
#         ylim = c(0,1), 
#         args.legend = list(x = 'center', cex = 0.8))
# title("Significance Level of 0.01 ")
# text(bp_1, 0, cex = 0.5, pos = 3)
# 
# #### Significance level of 0.05 ####
# mean_barplot_vect_2 <- c( sum(stars_mean_confusion_metrics$TP[,2])/100, length(which(stars_mean_confusion_metrics$FP[,2] >= 1))/100 )
# var_barplot_vect_2 <- c( sum(stars_var_confusion_metrics$TP[,2])/100, length(which(stars_var_confusion_metrics$FP[,2] >= 1))/100 )
# plot_barplot_2 <- matrix(c(mean_barplot_vect, var_barplot_vect), nrow = 2, byrow = F)
# colnames(plot_barplot_2) <- c("Mean", "Variance")
# rownames(plot_barplot_2) <- c("Proportion Correctly Detected", "Proportion Incorrectly Detected")
# barplot(plot_barplot_2, 
#         col=colors()[c(23,12)], 
#         border="white", 
#         font.axis=2, 
#         beside=T, 
#         legend=rownames(plot_barplot_2), 
#         xlab="Regime Shift Type", 
#         font.lab=2, 
#         ylim = c(0,1), 
#         args.legend = list(x = 'center', cex = 0.8))
# title("Significance Level of 0.05 ")
# 
# #### Significance level of 0.10 ####
# mean_barplot_vect_3 <- c( sum(stars_mean_confusion_metrics$TP[,3])/100, length(which(stars_mean_confusion_metrics$FP[,3] >= 1))/100 )
# var_barplot_vect_3 <- c( sum(stars_var_confusion_metrics$TP[,3])/100, length(which(stars_var_confusion_metrics$FP[,3] >= 1))/100 )
# plot_barplot_3 <- matrix(c(mean_barplot_vect, var_barplot_vect), nrow = 2, byrow = F)
# colnames(plot_barplot_3) <- c("Mean", "Variance")
# rownames(plot_barplot_3) <- c("Proportion Correctly Detected", "Proportion Incorrectly Detected")
# barplot(plot_barplot_3, 
#         col=colors()[c(23,12)], 
#         border="white", 
#         font.axis=2, 
#         beside=T, 
#         legend=rownames(plot_barplot_3), 
#         xlab="Regime Shift Type", 
#         font.lab=2, 
#         ylim = c(0,1), 
#         args.legend = list(x = 'center', cex = 0.8))
# title("Significance Level of 0.10 ")
# 
# 
# #### FOR GGPLOT ####
# 
# full_df <- data.frame(value = c(mean_barplot_vect_1, 
#                                 mean_barplot_vect_2,
#                                 mean_barplot_vect_3,
#                                 var_barplot_vect_1,
#                                 var_barplot_vect_2, 
#                                 var_barplot_vect_3), 
#                       alpha = c(rep('0.01', 2),
#                                 rep('0.05', 2),
#                                 rep('0.10', 2),
#                                 rep('0.01', 2),
#                                 rep('0.05', 2),
#                                 rep('0.10', 2)),
#                       type = c(rep('mean', 6),
#                                rep('variance', 6)),
#                       value_type = rep(c('Correct','Incorrect'), 6))
# ggplot(full_df, aes(x = factor(alpha), y = value, fill = value_type)) + 
#   geom_bar(stat = 'identity', position = 'dodge') +
#   labs(title = 'STARS Proportion of Correctly and Incorrectly Detected Shifts', 
#        x = 'Regime Shift Type',
#        y = 'Proportion') + 
#   theme_bw() +
#   facet_wrap( ~ type ) +
#   theme(plot.title = element_text(size=11, 
#                                   face = "bold", 
#                                   vjust = 5, 
#                                   hjust = 0.5,), 
#         axis.title.x = element_text(size = 14, vjust = -4), 
#         axis.title.y = element_text(size = 14, vjust = 2), 
#         aspect.ratio = 1.50, 
#         panel.grid.major = element_blank()) +
#   scale_fill_discrete(
#     labels = c("Correct", "Incorrect"), name = 'Metric') + guides(colour = FALSE) +
#   geom_text(aes(label = value),
#             position=position_dodge(width=0.9), size = 4,vjust=-0.5)
# 
# 
