#############################################################
### STARS RESULTS: SHIFT IN THE MEAN
#############################################################

##############################
## 0 - Load libraries
##############################
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

##############################
## 1 - Function for importing all sheets in an Excel workbook.
##############################

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "B1:CW1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
STARS_SHIFTS_IN_MEAN <- read_excel_allsheets("STARS_SHIFTS_IN_MEAN_RESULTS.xlsx")

##############################
## 2 - Function for working with STARS Results.
##############################

STARS_DETECTION_COUNTER <- function(.data){
  # .data (list) is the complete data set from the excel file that has all the RSI values and time series data for all parameters
  # N (numerical) is the interval of the simulated shifts, every N = 100 or N = 200 observations
  # dataset_ (character), look at the dataset_'' objects above. eg. 'dataset_ones'.  
  S <- c()
  U <- c()
  D <- matrix(NA,nrow = 100,ncol = 9)
  FP <- matrix(NA,nrow = 100,ncol = 9)
  FPR <- matrix(NA,nrow = 100,ncol = 9)
  DD <-matrix(NA,nrow = 9,ncol = 100)
  .FP <-matrix(NA,nrow = 9,ncol = 100)
  .FPR <-matrix(NA,nrow = 100,ncol = 9)
  LOC <- matrix(NA,nrow = 1000,ncol = 9)
  SFT <- list()
  for (k in 1:9){
    N <- c()
    M <- c()
    shifts <- matrix(NA,nrow = 1000,ncol = 100)
    for(j in 1:100){
      # Turns RSIs into zeros and ones
      for (i in 1:1000){
        if(.data[[k]][i,j] == 0){
          shifts[i,j] <- 0
        } else {
          shifts[i,j] <- 1
        }
      }
      # IDENTIFYING LOCATIONS OF DETECTIONS
      index_shifts <- which(shifts[,j] == 1)
      # DEFINING LOCATION OF KNOWN SHIFT
      know <- rep(0,1000)
      know[500] <- 1
      index_know <- which(know == 1)
      
      # ACCOMODATING LATE DETECTIONS OF KNOWN SHIFT WITHIN 10 UNITS
      if (shifts[index_know,j] == 1) {
        shifts[index_know,j] <- 1
      } else {
        index_near_shifts <- index_shifts[index_shifts >= (index_know-10) & index_shifts<=(index_know+10)]
        # IF THERE ARE MULTIPLE DETECTION NEAR BY CHOOSE THE NEAREST ONE TO REPESENET THE DETECTION OF THE KNOWN SHIFT
        if(length(index_near_shifts)==1){
          shifts[index_know,j] <- 1
          shifts[index_near_shifts,j] <- 0
        } else if(length(index_near_shifts)>1){
          closest <- which(abs(index_near_shifts-index_know)==min(abs(index_near_shifts-index_know)))
          closest_index <- index_near_shifts[closest[1]]
          shifts[closest_index,j] <- 0
          shifts[index_know,j] <- 1
        } else{
          shifts[index_know,j] <- 0
        }
      }
      
      # INDICATES WHETHER THE KNOWN SHIFT WAS DETECTED IN THIS VERSION
      # TO KEEP TRACK OF VERSIONS WHERE THE KNOWN SHIFT WAS DETECTED
      if(shifts[index_know,j]==1){
        N[j] <- 1
      } else {
        N[j] <- 0
      }
      # KEEP TRACK OF THE NUMBER OF DETECTIONS MADE IN EACH VERSION OF EACH TEST INCLUDING THE DETECTION OF THE KNOWN SHIFT
      D[j,k] <- sum(shifts[,j])
      # KEEP TRACK OF THE NUMBER OF UNKNOWN DETECTIONS MADE IN EACH VERSION OF EACH TEST EXCLUDING THE DETECTION OF THE KNOWN SHIFT
      if(D[j,k] == 0){
        FP[j,k] <- 0
        FPR[j,k] <- FP[j,k]/(1000-FP[j,k])
      } else if(D[j,k]>=1 && N[j]==0){
        FP[j,k] <- D[j,k]
        FPR[j,k] <- FP[j,k]/1000
      } else if(D[j,k]>=1 && N[j]==1){
        FP[j,k] <- D[j,k] - 1
        FPR[j,k] <- FP[j,k]/999
      }
      # INDICATES WHETHER ANY UNKNOWN SHIFTS WERE DETECTED IN THIS VERSION
      # TO KEEP TRACK OF VERSIONS WHERE MULTIPLE UNKNOWN DETECTIONS WERE MADE
      if(FP[j,k]==0){
        M[j] <- 0
      } else if(FP[j,k]>=1){
        M[j] <- 1
      }
    }
    # INDICATES THE PROPERTION OF VERSIONS WHERE THE KNOWN SHIFT WAS DETECTED
    # AND INDICATES THE PROPERTION OF VERSIONS WHERE MULTIPLE UNKNOWN SHIFTS WERE DETECTED
    # KEEPING THE LOCATIONS OF THE DETECTIONS
    LOC[,k] <- rowSums(shifts)
    SFT[[k]] <- shifts
    S[k] <- sum(N)/100
    U[k] <- sum(M)/100
    DD[k,] <- round(D[,k],0) 
    .FP[k,] <- round(FP[,k],0) 
    .FPR[,k] <- round(FPR[,k],2)
  }
  # COMBINE INTO ONE OBJECT
  SS <- data.frame(Known=S,row.names = c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1"))
  UU <- data.frame(Unknown=U,row.names = c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1"))
  DDD <- data.frame(DD,check.names = TRUE,row.names = c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1"))
  FP_ <- data.frame(.FP,check.names = TRUE,row.names = c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1"))
  FPR_ <- data.frame(.FPR,check.names = TRUE)
  names(FPR_) <- c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1")
  LOC_ <- data.frame(LOC,row.names = as.character(seq(1,1000,1)))
  names(LOC_) <- c("20-0.01","20-0.05","20-0.1","100-0.01","100-0.05","100-0.1","500-0.01","500-0.05","500-0.1")
  FINAL <- list(KNOWN=SS,UNKNOWN=UU,NO.DECT=DDD,FPRate=FPR_,FalsePos=FP_,SHIFTS=SFT,Locations=LOC_)
  return(FINAL)
}
a <- STARS_DETECTION_COUNTER(STARS_SHIFTS_IN_MEAN)

##############################
## 3 - Plotting Results.
##############################

##############################
# 3.1 - Plotting number of detection made.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
CUTOFF_LENGTH <- c(rep("20" , 300) , rep("100" , 300) , rep("500" , 300))
PVALUE <- rep(c(rep("0.01",100) ,rep("0.05",100) ,rep("0.1",100)) , 3)
NUMBERS <- rbind(transpose(a$NO.DECT[1,1:100]),transpose(a$NO.DECT[2,1:100]),transpose(a$NO.DECT[3,1:100]),
                 transpose(a$NO.DECT[4,1:100]),transpose(a$NO.DECT[5,1:100]),transpose(a$NO.DECT[6,1:100]),
                 transpose(a$NO.DECT[7,1:100]),transpose(a$NO.DECT[8,1:100]),transpose(a$NO.DECT[9,1:100]))
colnames(NUMBERS) <- "DETECTIONS"
PLOT_DATA <- data.frame(CUTOFF_LENGTH=as.factor(CUTOFF_LENGTH),'Probability Level'=as.factor(PVALUE),NUMBERS)

# GROUPED AND STACKED BOXPLOTS
ggplot(PLOT_DATA, aes(x=CUTOFF_LENGTH, y=DETECTIONS, fill=PVALUE)) + 
  geom_boxplot() + 
  ggtitle("Number of Detections Made By STARS Method") +
  theme_bw() +
  facet_wrap(~Probability.Level,labeller = "label_both") + 
  scale_y_continuous(name="Number of Detections Made", limits=c(0, 10),n.breaks = 10) + 
  scale_x_discrete(limits=c("20","100","500"),name="Cutoff Length") +
  scale_fill_discrete(name = "Probability\nLevel",labels = c("0.01", "0.05", "0.10"))

##############################
# 3.2 - Plotting FREQUENCY DISTRIBUTION.
##############################

# FREQUENCY GRAPH FOR VARIOUS CUTOFF LENGTHS
PLOT_DATA <- data.frame(L20P0.05=a$Locations[,2],L100P0.05=a$Locations[,5],L500P0.05=a$Locations[,8],check.names = TRUE)
# FREQUENCY GRAPH
colors <- c("20" = "steelblue", "100" = "salmon", "500" = "khaki")
a$Locations$`100-0.01`
ggplot(data = PLOT_DATA) + 
  geom_bar(aes(x = index(PLOT_DATA), y = L20P0.05,fill ="20"),stat = "identity",width=7.1) +
  geom_bar(aes(x = index(PLOT_DATA), y = L100P0.05,fill="100"),stat = "identity",width=7.1) +
  geom_bar(aes(x = index(PLOT_DATA), y = L500P0.05,fill ="500"),stat = "identity",width=7.1) +
  ggtitle("Frequency Distribution of Regime Shift Detections By STARS") +
  theme_bw() +
  scale_x_continuous(name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10) + 
  labs(fill = "Cut-off Length") +
  scale_fill_manual(values = colors)

# FREQUENCY GRAPH FOR VARIOUS PVALUES

PLOT_DATA <- data.frame(L100P0.01=a$Locations[,4],L100P0.05=a$Locations[,5],L100P0.10=a$Locations[,6],check.names = TRUE)

# FREQUENCY GRAPH
colors <- c("0.01" = "steelblue", "0.05" = "salmon", "0.10" = "khaki")

ggplot(data = PLOT_DATA) + 
  geom_bar(aes(x = index(PLOT_DATA), y = L100P0.01,fill = "0.01"),stat = "identity",width=7.1) +
  geom_bar(aes(x = index(PLOT_DATA), y = L100P0.05,fill="0.05"),stat = "identity",width=7.1) +
  geom_bar(aes(x = index(PLOT_DATA), y = L100P0.10,fill ="0.10"),stat = "identity",width=7.1) +
  ggtitle("Frequency Distribution of Regime Shift Detections By STARS") +
  theme_bw() +
  scale_x_continuous(name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10) + 
  labs(fill = "Probability\nLevel") +
  scale_fill_manual(values = colors)




##############################
# 3.2 - Plotting Proportion of Experiments that Detected the True Regime Shift (SENSITIVITY)
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
CUTOFF_LENGTH <- c(rep("20" , 3) , rep("100" , 3) , rep("500" , 3))
PVALUE <- rep(c("0.01" , "0.05" , "0.1") , 3)
PROPORTIONS <- a$KNOWN$Known 
PLOT_DATA <- data.frame(as.factor(CUTOFF_LENGTH),as.factor(PVALUE),PROPORTIONS)

# STACKED HISTOGRAM 
ggplot(PLOT_DATA, aes(fill=PVALUE, y=PROPORTIONS, x=CUTOFF_LENGTH)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Proportion of Experiments that Detected the True Regime Shift ") +
  theme_bw() +
  scale_x_discrete(limits=c("20","100","500"),name="Cutoff Length") +
  scale_y_continuous(name="Proportion", limits=c(0, 1),n.breaks = 10) +
  scale_fill_discrete(name = "Probability\nLevel",labels = c("0.01", "0.05", "0.10"))



##############################
# 3.4 - Plotting probability of at least one false detection.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
CUTOFF_LENGTH <- c(rep("20" , 3) , rep("100" , 3) , rep("500" , 3))
PVALUE <- rep(c("0.01" , "0.05" , "0.1") , 3)
PROPORTIONS <- a$UNKNOWN$Unknown 
PLOT_DATA <- data.frame(as.factor(CUTOFF_LENGTH),as.factor(PVALUE),PROPORTIONS)

# STACKED HISTOGRAM 
ggplot(PLOT_DATA, aes(fill=PVALUE, y=PROPORTIONS, x=CUTOFF_LENGTH)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Proportion of Experiments with At Least One False Detectioin") +
  theme_bw() +
  scale_x_discrete(limits=c("20","100","500"),name="Cutoff Length") +
  scale_y_continuous(name="Probability", limits=c(0, 1),n.breaks =10) +
  scale_fill_discrete(name = "Probability\nLevel",labels = c("0.01", "0.05", "0.10"))


