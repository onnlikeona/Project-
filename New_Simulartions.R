library(tidyverse)
library(stats)
library(astsa)
library(gdata)
library(openxlsx)
#Dataset1
dataset1 <- list()
for(j in 1:100){
  set.seed(j)
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
  dataset1[[j]] <- data
}
plot.ts(dataset1[[9]], main="Simulated Daily Returns Over Time", xlab="Days", ylab="Daily Return")
abline(v=c(500),lty=2,col="red")
legend(1, 5, legend=c("Daily Return", "Regime Shift in the Mean"),
       col=c("black", "red"), lty=1:2, cex=0.8,
       title="Line Types", text.font=4, bg="white")


#Dataset2
dataset2 <- list()
for(j in 1:100){
  set.seed(j)
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
  dataset2[[j]] <- data
}
plot.ts(dataset2[[6]], main="Simulated Daily Returns Over Time", xlab="Days", ylab="Daily Return")
abline(v=c(500),lty=2,col="red")
legend(1, 1, legend=c("Daily Return", "Regime Shift in the Variance"),
       col=c("black", "red"), lty=1:2, cex=0.8,
       title="Line Types", text.font=4, bg="white")

#Dataset3
dataset3 <- list()
for(j in 1:100){
  set.seed(j)
  var1 <- 0.3
  var2 <- 0.6
  mu1 <- 0.5
  mu2 <- -0.1
  x1 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu1,var1))
  
  x2 <- arima.sim(list(order=c(2,0,0),
                       ar=c(0.6,-0.5)),
                  n=700,
                  innov=rnorm(700,mu2,var2))
  
  data <- ts(c(x1[101:600],x2[101:600]))
  dataset3[[j]] <- data
}
plot.ts(dataset3[[9]], main="Simulated Daily Returns Over Time", xlab="Days", ylab="Daily Return")
abline(v=c(500),lty=2,col="red")
legend(1, -0.5, legend=c("Daily Return", "Regime Shift"),
       col=c("black", "red"), lty=1:2, cex=0.8,
       title="Line Types", text.font=4, bg="white")

Simulated_data_workbook <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(Simulated_data_workbook, "Data Set 1")
addWorksheet(Simulated_data_workbook, "Data Set 2")
addWorksheet(Simulated_data_workbook, "Data Set 3")


# Write the data to the sheets
writeData(Simulated_data_workbook, sheet = "Data Set 1", x = dataset1)
writeData(Simulated_data_workbook, sheet = "Data Set 2", x = dataset2)
writeData(Simulated_data_workbook, sheet = "Data Set 3", x = dataset3)

# Export the file
saveWorkbook(Simulated_data_workbook, "SIMULATIONS SPREADSHEETS.xlsx")

