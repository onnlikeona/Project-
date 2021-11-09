rm(list=ls())
library(tidyverse)
library(openxlsx) 
library(Rcpp) 
library(timeSeries) 
library(zoo)
library(xts)
library(latex2exp)
library(rbenchmark)
library(nloptr)
library(quadprog)
library(gdata)
library(quantmod)
library(depmixS4)
library(hsmm)
library(seqHMM)
library(lubridate)
library(anytime)
library(ggplot2)
library(forecast)
library(dygraphs)
library(sweep)
library(timekit)
library(googleVis)
citation("hsmm")

setwd("/Users/ndu/OneDrive - University of Cape Town/THE GRIND/ACADEMIC GRIND/2020/PROJECT/NEW STUFF/Project--master")

#######################
#LOADING_DATA_AND_FUNCTION
#######################
source("rSTARS.R")

#######################
#PLOTTING_ALSH_DATA
#######################
alsh_original <- read.csv("/Users/ndu/OneDrive - University of Cape Town/THE GRIND/ACADEMIC GRIND/2020/PROJECT/NEW STUFF/Project--master/AllShareHistoricalPrices.csv", header = T)
alsh_original[,1] <- as.Date(alsh_original[,1],format = "%m/%d/%y")
alsh <- xts(alsh_original[, -1], order.by=as.POSIXct(alsh_original$Date))

chartSeries(alsh$Close,theme = chartTheme("white"),
            name = "All Share Index Prices",up.col="black", )
diff_alsh = diff(log( Cl( alsh ) ))
diff_alsh = na.omit(diff_alsh)
chartSeries(diff_alsh,theme = chartTheme("white"),
            name = "First Difference of the Log Prices For All Share Index",up.col="black")
plot(diff_alsh,plot.type = c("single"), 
     format = "auto", 
     at=pretty(diff_alsh),
     ylab = "Log Prices", 
     main = "First Difference of the Log Prices For All Share Index", lwd = 1)
chartSeries(diff_alsh,theme = chartTheme("white"),name = "First Difference of the Log Prices For All Share Index",up.col="black")
plot(alsh$Close,plot.type = c("single"), 
     format = "auto", 
     at=pretty(alsh$Close),
     ylab = "Prices", 
     main = "All Share Index Prices", lwd = 1)


###############################################
#       Hidden_Markov
##############################################
#Hidden_Markov_Using_Logged_Differenced_Data
##############################################
hmm <- depmix(Close ~ 1, family =gaussian(), nstates = 3, data=diff_alsh)
set.seed(1)
hmmfit <- fit(hmm, verbose = FALSE)
hmmvit <- viterbi(hmmfit)
hmmpred <- cbind(alsh$Close[-1,],hmmvit$state)
hmmpostprobs <- posterior(hmmfit, type = 'viterbi') ####
#hmm.post.df <- melt(hmmpostprobs, measure.vars = c("S1","S2","S3"))
#hmm.post.df$Time <- date(alsh)[-1]
mycols <- c("salmon", "khaki","lightgreen")
## PLOT 1
date__ <- alsh_original$Date[-1]
#par(mar=c(6, 5, 4, 7) + 0.1, xpd=T) # used to be 2
par(mar=c(5, 5.5, 4, 10), xpd=TRUE)
matplot(date__,hmmpostprobs[,-1], 
        type='l', 
        lty = 1, 
        lwd = 1.6, 
        main='Regime Viterbi Posterior Probabilities', 
        ylab='Probability', 
        col=mycols, 
        xlab = 'Date',
        cex.axis = 2.2, 
        cex.lab = 2.5, 
        cex.main = 2.5)
legend("right", 
       inset=c(-0.22, 0), 
       legend=c("State 1","State 2", "State 3"), 
       col = c("salmon", "khaki","lightgreen"), 
       lty = c(1, 1, 1),
       lwd = 3,
       title="Regimes",
       cex = 1.5)

plot(alsh$Close,format="auto",main="All Share Index Prices", xlab="", ylab="Prices")


## PLOT 2
hmmpostprobs <- as.data.frame(hmmpostprobs)
hmmpostprobs <- cbind(date__, hmmpostprobs) %>% as.data.frame()
hmmpostprobs_ts <- xts(hmmpostprobs[, -1], order.by=as.POSIXct(hmmpostprobs$date__))
par(mar=c(6, 5, 4, 2) + 0.1)
chartSeries(alsh$Close[-1],theme = chartTheme("white"),name = "All Share Index Prices",up.col="black")
State_1 <- hmmpred[hmmpred[,2] == 1, 1]
State_2 <- hmmpred[hmmpred[,2] == 2, 1]
State_3 <- hmmpred[hmmpred[,2] == 3, 1]
par(mar=c(6, 5, 4, 2) + 0.1)
addTA(State_1,on=1,type="p",col=mycols[1],pch=20)
addTA(State_2,on=1,type="p",col=mycols[3],pch=20)
addTA(State_3,on=1,type="p",col=mycols[2],pch=20)
# Parameter distriution
summary(hmmfit)

##############################################
#Hidden_Markov_Using_Original_Data
##############################################
hmm <- depmix(Close ~ 1, family =gaussian(), nstates = 3, data=alsh)
hmmfit <- fit(hmm, verbose = FALSE)
hmmvit <- viterbi(hmmfit,alsh$Close)
hmmpred <- cbind(alsh$Close,hmmvit$state)
hmmpostprobs <- posterior(hmmfit)
hmm.post.df <- melt(hmmpostprobs, measure.vars = c("S1","S2","S3"))
hmm.post.df$Time <- date(alsh$Close)
mycols <- c("salmon", "khaki","lightgreen")
## PLOT 1
matplot(hmmpostprobs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability',col=mycols)
plot(alsh$Close,format="auto",main="All Share Index Prices", xlab="", ylab="Prices")
## PLOT 2
chartSeries(alsh$Close,theme = chartTheme("white"),name = "All Share Index Prices",up.col="black")
addTA(hmmpred[hmmpred[,2]==1,1],on=1,type="p",col=mycols[1],pch=20)
addTA(hmmpred[hmmpred[,2]==2,1],on=1,type="p",col=mycols[3],pch=20)
addTA(hmmpred[hmmpred[,2]==3,1],on=1,type="p",col=mycols[2],pch=20)

##############################
#STARS
##############################

a <- data.frame(alsh)
Time <-as.character(row.names(a))
x <- data.frame(cbind(Time,Close=a$Close))
x$Close <- as.numeric(x$Close)
r = diff(log( Cl( alsh ) ))
r = na.omit(r)
r <- data.frame(r)
x <- data.frame(cbind(Time[-1],Close=r$Close))
x$Close <- as.numeric(x$Close)

rsi_31 <- rstars(x,l.cutoff = 31,pValue = 0.05,save.path = getwd())
rsi_14 <- rstars(x,l.cutoff = 14,pValue = 0.05,save.path = getwd())
rsi_24 <- rstars(x,l.cutoff = 24,pValue = 0.05,save.path = getwd())


table_plot = cbind(as.data.frame(rsi_31$Close),as.data.frame(rsi_24$Close),as.data.frame(rsi_14$Close),diff(log( Cl( alsh ) ))[-1])
A <- as.xts(table_plot,order.by = as.Date(Time[-1]))
plot(A,col = c("#EB4C43","khaki","lightgreen","#3081B5"),
     main = paste("Regime Shift detection in","Daily Closing ALSH Price Index","using STARS"),
     xlab = "Time", ylab = "Log Prices")
addLegend('topright',legend.names = c("l = 31 days", "l = 24 days" , "l = 14 days"),
          col = c("#EB4C43","khaki","lightgreen"),lwd = c(2,2))


rsi <- cbind(rsi_31$Close.1,rsi_24$Close.1,rsi_14$Close.1)
rsi <- abs(data.frame(rsi, row.names = NULL))
rsi <- xts(rsi, order.by = as.Date(Time[-1]))

plot(rsi, plot.type = c("multiple"),col = c("#EB4C43","khaki","lightgreen"),
     main = paste("Regime Shift Index for", "Daily Closing ALSH Index"),
     xlab = "Time", ylab = "RSI")
addLegend('topright',legend.names = c("l = 31 days", "l = 24 days" , "l = 14 days"),
          col = c("#EB4C43","khaki","lightgreen"),lwd = c(2,2))

rsi_p01 <- rstars(x,l.cutoff = 31,pValue = 0.1,save.path = getwd())
rsi_p005 <- rstars(x,l.cutoff = 31,pValue = 0.05,save.path = getwd())
rsi_p001<- rstars(x,l.cutoff = 31,pValue = 0.01,save.path = getwd())


table_plot = cbind(as.data.frame(rsi_p01$Close),as.data.frame(rsi_p005$Close),as.data.frame(rsi_p001$Close),diff(log( Cl( alsh ) ))[-1])
B <- as.xts(table_plot,order.by = as.Date(Time[-1]))
plot(B,col = c("khaki","lightgreen","#EB4C43","#3081B5"),
     main = paste("Regime Shift detection in","Daily Closing ALSH Price Index","using STARS"),
     xlab = "Time", ylab = "Log Prices")
addLegend('topright',legend.names = c("p = 0.1", "p = 0.05" , "p = 0.01"),
          col = c("khaki","lightgreen","#EB4C43"),lwd = c(2,2))


rsi <- cbind(rsi_p001$Close.1,rsi_p005$Close.1,rsi_p01$Close.1)
rsi <- abs(data.frame(rsi, row.names = NULL))
rsi <- xts(rsi, order.by = as.Date(Time[-1]))

plot(rsi, plot.type = c("multiple"),col = c("khaki","lightgreen","#EB4C43"),
     main = paste("Regime Shift Index for", "Daily Closing ALSH Index"),
     xlab = "Time", ylab = "RSI")
addLegend('topright',legend.names = c("p = 0.1", "p = 0.05" , "p = 0.01"),
          col = c("khaki","lightgreen","#EB4C43"),lwd = c(2,2))

