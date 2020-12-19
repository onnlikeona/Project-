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
setwd("~/Documents/2020/Semester 2/Project/Application of Methods/REALDATA")
#######################
LOADING_DATA_AND_FUNCTION
#######################

source("rSTARS.R")
load(file = "All_Share_Daily_1995_2020")

#######################
PLOTTING_ALSH_DATA
#######################

chartSeries(alsh$Close,theme = chartTheme("white"),
            name = "All Share Index Prices",up.col="black")
chartSeries(diff_alsh,theme = chartTheme("white"),
            name = "First Difference of the Log Prices For All Share Index",up.col="black")

###############################################
Hidden_Markov
##############################################
Hidden_Markov_Using_Logged_Differenced_Data
##############################################
hmm <- depmix(Close ~ 1, family =gaussian(), nstates = 3, data=diff_alsh)
hmmfit <- fit(hmm, verbose = FALSE)
hmmvit <- viterbi(hmmfit,r)
hmmpred <- cbind(alsh$Close[-1,],hmmvit$state)
hmmpostprobs <- posterior(hmmfit)
hmm.post.df <- melt(hmmpostprobs, measure.vars = c("S1","S2","S3"))
hmm.post.df$Time <- date(alsh)[-1]
mycols <- c("salmon", "khaki","lightgreen")
## PLOT 1
matplot(hmmpostprobs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability',col=mycols)
plot(alsh$Close,format="auto",main="All Share Index Prices", xlab="", ylab="Prices")
## PLOT 2
chartSeries(alsh$Close,theme = chartTheme("white"),name = "All Share Index Prices",up.col="black")
addTA(hmmpred[hmmpred[,2]==1,1],on=1,type="p",col=mycols[1],pch=20)
addTA(hmmpred[hmmpred[,2]==2,1],on=1,type="p",col=mycols[3],pch=20)
addTA(hmmpred[hmmpred[,2]==3,1],on=1,type="p",col=mycols[2],pch=20)

##############################################
Hidden_Markov_Using_Original_Data
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
STARS
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
plot(xts(table_plot,order.by = Time[-1]),col = c("#EB4C43","khaki","lightgreen","#3081B5"),
     main = paste("Regime Shift detection in","Daily Closing ALSH Price Index","using STARS"))
addLegend('topright',legend.names = c("l = 31 days", "l = 24 days" , "l = 14 days"),
          col = c("#EB4C43","khaki","lightgreen"),lwd = c(2,2))


rsi <- cbind(rsi_31$Close.1,rsi_24$Close.1,rsi_14$Close.1)
rsi <- abs(data.frame(rsi, row.names = NULL))

rsi <- xts(rsi, order.by = Time[-1])

plot(rsi, plot.type = c("multiple"),col = c("#EB4C43","khaki","lightgreen"),
     main = paste("Regime Shift Index for", "Daily Closing ALSH Index"))
addLegend('topright',legend.names = c("l = 31 days", "l = 24 days" , "l = 14 days"),
          col = c("#EB4C43","khaki","lightgreen"),lwd = c(2,2))

rsi_p01 <- rstars(x,l.cutoff = 31,pValue = 0.1,save.path = getwd())
rsi_p005 <- rstars(x,l.cutoff = 31,pValue = 0.05,save.path = getwd())
rsi_p001<- rstars(x,l.cutoff = 31,pValue = 0.01,save.path = getwd())


table_plot = cbind(as.data.frame(rsi_p01$Close),as.data.frame(rsi_p005$Close),as.data.frame(rsi_p001$Close),diff(log( Cl( alsh ) ))[-1])
plot(xts(table_plot,order.by = Time[-1]),col = c("khaki","lightgreen","#EB4C43","#3081B5"),
     main = paste("Regime Shift detection in","Daily Closing ALSH Price Index","using STARS"))
addLegend('topright',legend.names = c("p = 0.1", "p = 0.05" , "p = 0.01"),
          col = c("khaki","lightgreen","#EB4C43"),lwd = c(2,2))


rsi <- cbind(rsi_p001$Close.1,rsi_p005$Close.1,rsi_p01$Close.1)
rsi <- abs(data.frame(rsi, row.names = NULL))

rsi <- xts(rsi, order.by = Time[-1])

plot(rsi, plot.type = c("multiple"),col = c("khaki","lightgreen","#EB4C43"),
     main = paste("Regime Shift Index for", "Daily Closing ALSH Index"))
addLegend('topright',legend.names = c("p = 0.1", "p = 0.05" , "p = 0.01"),
          col = c("khaki","lightgreen","#EB4C43"),lwd = c(2,2))

