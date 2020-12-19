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
IMPORTING_ALSH_DATA
#######################
alsh <- read.csv("AllShareHistoricalPrices.csv")
alsh[,1] <- as.Date(alsh[,1],format = "%m/%d/%y")
alsh <- xts(alsh[, -1], order.by=as.POSIXct(alsh$Date))


#######################
PLOTTING_ALSH_DATA
#######################

chartSeries(alsh$Close,theme = chartTheme("white"),name = "All Share Index Prices",up.col="black")

#######################
PLOTTING_DIFFERENCED_DATA
#######################

diff_alsh = diff(log( Cl( alsh ) ))
diff_alsh = na.omit(diff_alsh)
plot(diff_alsh,plot.type = c("single"), 
     format = "auto", 
     at=pretty(diff_alsh),
     ylab = "Log Prices", 
     main = "First Difference of the Log Prices For All Share Index")
chartSeries(diff_alsh,theme = chartTheme("white"),name = "First Difference of the Log Prices For All Share Index",up.col="black")

keep(alsh,diff_alsh,sure = TRUE)
save(alsh,diff_alsh,file = "All_Share_Daily_1995_2020")
unlink("AllShare.RData")
unlink(".RData")