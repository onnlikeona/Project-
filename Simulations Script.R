#### LOAD LIBRARIES ####

rm(list=ls())
library(tidyverse)
library(stats)
library(astsa)
library(gdata)
library(openxlsx)

#### Data Set 1 ####

dataset_1 <- list()
for(i in 1:100){
  set.seed(i)
  mu1 <- 0
  mu2 <- 1
  var <- 1
  x1 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu1,var))
  
  x2 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu2,var))
  
  data <- ts(c(x1[101:600],x2[101:600]))
  dataset_1[[i]] <- data
}
plot.ts(dataset_1[[9]], main="Simulated Daily Returns Over Time for Mean Regime Shift", xlab="Days", ylab="Daily Return")
abline(v=c(500),lty=2,col="red")
legend(1, 5, legend=c("Daily Return", "Regime Shift in the Mean"),
       col=c("black", "red"), lty=1:2, cex=0.8,
       title="Line Types", text.font=4, bg="white")


#### Data Set 2 ####

dataset_2 <- list()
for(k in 1:100){
  set.seed(k)
  var1 <- 0.1
  var2 <- 0.3
  mu <- 0
  x1 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu,var1))
  
  x2 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu,var2))
  
  data <- ts(c(x1[101:600],x2[101:600]))
  dataset_2[[k]] <- data
}
plot.ts(dataset_2[[6]], main="Simulated Daily Returns Over Time for Variance Regime Shift", xlab="Days", ylab="Daily Return")
abline(v=c(500),lty=2,col="red")
legend(1, 1, legend=c("Daily Return", "Regime Shift in the Variance"),
       col=c("black", "red"), lty=1:2, cex=0.8,
       title="Line Types", text.font=4, bg="white")

#### Data Set 3 ####

# dataset_3 <- list()
# for(j in 1:100){
#   set.seed(j)
#   var1 <- 0.3
#   var2 <- 0.6
#   mu1 <- 0.5
#   mu2 <- -0.1
#   x1 <- arima.sim(list(order=c(2,0,0),
#                        ar=c(0.6,-0.5)),
#                   n=700,
#                   innov=rnorm(700,mu1,var1))
#   
#   x2 <- arima.sim(list(order=c(2,0,0),
#                        ar=c(0.6,-0.5)),
#                   n=700,
#                   innov=rnorm(700,mu2,var2))
#   
#   data <- ts(c(x1[101:600],x2[101:600]))
#   dataset_3[[j]] <- data
# }
# plot.ts(dataset_3[[9]], main="Simulated Daily Returns Over Time", xlab="Days", ylab="Daily Return")
# abline(v=c(500),lty=2,col="red")
# legend(1, -0.5, legend=c("Daily Return", "Regime Shift"),
#        col=c("black", "red"), lty=1:2, cex=0.8,
#        title="Line Types", text.font=4, bg="white")

#### EXPORT TO SIMULATIONS TO EXCEL WORKBOOK ####
Simulated_data_workbook <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(Simulated_data_workbook, "Data Set 1")
addWorksheet(Simulated_data_workbook, "Data Set 2")
#addWorksheet(Simulated_data_workbook, "Data Set 3")


# Write the data to the sheets
writeData(Simulated_data_workbook, sheet = "Data Set 1", x = dataset_1)
writeData(Simulated_data_workbook, sheet = "Data Set 2", x = dataset_2)
#writeData(Simulated_data_workbook, sheet = "Data Set 3", x = dataset_3)

# Export the file
saveWorkbook(Simulated_data_workbook, "Simulations.xlsx")

