#############################################################
### STARS RESULTS: SHIFT IN THE VARIANCE
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
library(ggthemes)

##############################
## 1 - Function for importing all sheets in an Excel workbook.
##############################

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, range = "A1:CV1001",col_names = TRUE))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
DATA <- read_excel_allsheets("SIMULATIONS SPREADSHEETS.xlsx")

##############################
## 2 - Function for Applying HMM.
##############################

detections <- function(posterior.states){ 
  # outputs a vector indicating where shifts occurred
  # 1 indicates no shift, 0 indicates a shift
  
  N <- length(posterior.states)
  pred_shifts <- c()
  pred_shifts[1] <- 0
  pred_shifts[N] <- 0
  for(i in 2:(N-1)){
    ifelse(posterior.states[i] - posterior.states[i-1] == 0, pred_shifts[i] <- 0, pred_shifts[i] <- 1)
  }
  return(pred_shifts)
}

application_hmm <- function(.data,type){
  D <- matrix(data=NA,nrow = 1000,ncol = 100)
  DDD <- list()
  states <- c(2)
  if (type=="mean") {
    df <- data.frame(matrix(unlist(.data[[1]]), ncol=length(.data[[1]]), byrow=F))
    L <- dim(df)[2]
    for (k in 1:1) {
      for (i in 1:L) {
        hmm <- depmix(df[,i] ~ 1, family = gaussian(), nstates = states[k], data=data.frame(returns=df[,i]))
        hmmfit <- fit(hmm, verbose = FALSE)
        post_probs <- posterior(hmmfit)
        D[,i]<- detections(post_probs[,1])
      }
      DD <- data.frame(D,check.names = TRUE)
      DDD[[k]] <- DD
    }
    return(DDD)
  } else if (type=="variance") {
    df <- data.frame(matrix(unlist(.data[[2]]), ncol=length(.data[[2]]), byrow=F))
    L <- dim(df)[2]
    for (k in 1:1) {
      for (i in 1:L) {
        hmm <- depmix(df[,i] ~ 1, family = gaussian(), nstates = states[k], data=data.frame(returns=df[,i]))
        hmmfit <- fit(hmm, verbose = FALSE)
        post_probs <- posterior(hmmfit)
        D[,i]<- detections(post_probs[,1])
      }
      DD <- data.frame(D,check.names = TRUE)
      DDD[[k]] <- DD
    }
    return(DDD)
  } else if (type=="both") {
    df <- data.frame(matrix(unlist(.data[[3]]), ncol=length(.data[[3]]), byrow=F))
    L <- dim(df)[2]
    for (k in 1:1) {
      for (i in 1:L) {
        hmm <- depmix(df[,i] ~ 1, family = gaussian(), nstates = states[k], data=data.frame(returns=df[,i]))
        hmmfit <- fit(hmm, verbose = FALSE)
        post_probs <- posterior(hmmfit)
        D[,i]<- detections(post_probs[,1])
      }
      DD <- data.frame(D,check.names = TRUE)
      DDD[[k]] <- DD
    }
    return(DDD)
  }
}
o <- application_hmm(DATA,"mean")

##############################
## 2 - Function for working with STARS Results.
##############################

hmm_detection_counter <- function(.data){
  # .data (list) is the complete data set from the excel file that has all the RSI values and time series data for all parameters
  # N (numerical) is the interval of the simulated shifts, every N = 100 or N = 200 observations
  # dataset_ (character), look at the dataset_'' objects above. eg. 'dataset_ones'.  
  S <- c()
  U <- c()
  D <- matrix(NA,nrow = 100,ncol = 1)
  FP <- matrix(NA,nrow = 100,ncol = 1)
  FPR <- matrix(NA,nrow = 100,ncol = 1)
  DD <-matrix(NA,nrow = 1,ncol = 100)
  .FP <-matrix(NA,nrow = 1,ncol = 100)
  .FPR <-matrix(NA,nrow = 1,ncol = 100)
  LOC <- matrix(NA,nrow = 1000,ncol = 1)
  shifts <- matrix(NA,nrow = 1000,ncol = 100)
  for (k in 1:1){
    N <- c()
    M <- c()
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
    LOC <- rowSums(shifts)
    S <- N
    U[k] <- sum(M)/100
    DD[k,] <- round(D[,k],0) 
    .FP[k,] <- round(FP[,k],0) 
    .FPR[k,] <- round(FPR[,k],2)
  }
  # COMBINE INTO ONE OBJECT
  SS <- data.frame(Known=S)
  UU <- data.frame(Unknown=U)
  DDD <- data.frame(DD,check.names = TRUE,row.names = c("2"))
  FP_ <- data.frame(.FP,check.names = TRUE,row.names = c("2"))
  FPR_ <- data.frame(.FPR,check.names = TRUE,row.names = c("2"))
  LOC_ <- data.frame(LOC,row.names = as.character(seq(1,1000,1)))
  FINAL <- list(KNOWN=SS,UNKNOWN=UU,NO.DECT=DDD,FPRate=FPR_,FalsePos=FP_,Locations=LOC_)
  return(FINAL)
}

r <- hmm_detection_counter(o)

##############################
## 3 - Plotting Results.
##############################

##############################
## 3.1 - Plotting number of detection made.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
STATES <- c(rep("2" , 100))
NUMBERS <- rbind(transpose(r$NO.DECT[1,1:100]))
colnames(NUMBERS) <- "DETECTIONS"
PLOT_DATA <- data.frame(STATES=as.factor(STATES),NUMBERS)
ggplot(PLOT_DATA, aes(x=DETECTIONS)) + geom_density()
ggplot(PLOT_DATA, aes(x=DETECTIONS)) +
  geom_histogram(binwidth=.5, colour="black", fill="white")
# GROUPED AND STACKED BOXPLOTS
ggplot(PLOT_DATA, aes(x=STATES, y=DETECTIONS,fill = STATES)) + 
  geom_boxplot() + 
  ggtitle("Number of Detections Made by an Experiment") +
  theme_bw() +
  scale_x_discrete(limits=c("2"),name="Number of States") +
  scale_y_continuous(name="Number of Detections",n.breaks = 10,limits = c(0,300))

# FREQUENCY GRAPH
PLOT_DATA <- r$Locations
ggplot(PLOT_DATA, aes(x = index(PLOT_DATA), y = LOC)) + 
  geom_bar(fill = "steelblue",stat = "identity") +
  ggtitle("Frequency of Detections Made by an Experiment at each Time") +
  theme_bw() +
  scale_x_continuous(n.breaks=10,name="Time") +
  scale_y_continuous(name="Frequency",n.breaks = 10)

ggplot(dat, aes(x=rating)) + geom_density()
##############################
# 3.2 - Plotting probability of achieving 100% sensitivity.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
STATES <- c(rep("2" , 1))
PROPORTIONS <- r$KNOWN$Known 
PLOT_DATA <- data.frame(STATES=as.factor(STATES),PROPORTIONS)

# STACKED HISTOGRAM 
ggplot(PLOT_DATA, aes(fill=STATES, y=PROPORTIONS, x=STATES)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Probability of HMM achieving 100% Sensitivity") +
  theme_bw() +
  scale_x_discrete(limits=c("2"),name="Number of States") +
  scale_y_continuous(name="Probabilities", limits=c(0, 1),n.breaks = 10) +
  scale_fill_discrete(name = "Number of\nStates",labels = c("2"))

# HISTOGRAM OF SENSITIVITY
SENSITIVITY <- r$KNOWN$Known 
PLOT_DATA <- data.frame(SENSITIVITY)
ggplot(PLOT_DATA, aes(x=SENSITIVITY)) + 
  geom_histogram(binwidth=100,colour="black", fill="steelblue") +
  ggtitle("Histogram of Sensitivity of HMM Experiments") +
  theme_bw() +
  scale_x_continuous(n.breaks=2,name="Sensitivity") +
  scale_y_continuous(name="Frequency",n.breaks = 15)

##############################
# 3.3 - Plotting number of false positives made.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
STATES <- c(rep("2" , 100))
NUMBERS <- rbind(transpose(r$FalsePos[1,1:100]))
colnames(NUMBERS) <- "DETECTIONS"
PLOT_DATA <- data.frame(STATES=as.factor(STATES),NUMBERS)

# GROUPED AND STACKED BOXPLOTS
ggplot(PLOT_DATA, aes(x=STATES, y=DETECTIONS, fill=STATES)) + 
  geom_boxplot() + 
  ggtitle("Frequency of False Positives across Two-State HMM Experiments") +
  theme_bw() +
  scale_y_continuous(name="Frequency of False Positives",n.breaks = 10) + 
  scale_x_discrete(limits=c("2"),name="Number of States") +
  scale_fill_discrete(name = "Number of\nStates",labels = c("2"))

##############################
# 3.4 - Plotting FPR.
##############################

# CREATING DATASET THAT WILL BE USED TO PLOT
PLOT_DATA <- transpose(r$FPRate)
# Plotting FPR.
ggplot(PLOT_DATA, aes(x=V1)) + 
  geom_histogram(binwidth=.005,colour="black", fill="steelblue") +
  ggtitle("Histogram of the False Positive Rates") +
  theme_bw() +
  scale_x_continuous(n.breaks=10,name="False Positive Rate (FPR)") +
  scale_y_continuous(name="Frequency",n.breaks = 15,limits = c(0,15))
