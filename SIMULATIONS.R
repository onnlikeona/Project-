library(tidyverse)
library(astsa)
library(gdata)
library(openxlsx)
set.seed(696)

# Data Set 1
f <- function(x){
  if((x %% 2) == 0) {
    return(1)
  } else {
    return(0)
  }
}

dataset_1 <- list()

for (r in 1:100) {
  set.seed(r)
  mu = 0
  mu1 = -0.1
  mu2 = 1.1
  var = 0.1
  data.set <- c()
  
for (i in 0:9) {
  mu1 = mu + 0.2*i
  mu2 = mu + 0.1*i
  mu = f(i)*mu1+(1-f(i))*mu2
  x <- arima.sim(list(order=c(2,0,0),ar=c(0.6,-0.5)),n=100,innov=rnorm(100,mu,1),n.start=100)
  data.set[((i*100)+1):(100*(i+1))] <- ts(x)
}
  
  dataset_1[[r]] <- data.set
  
}


# Data Set 2
dataset_2 <- list()
for(r in 1:100){
  set.seed(696)
  mu = 0
  mu1 = -0.1
  mu2 = 1.1
  var = 0.1
  
for (i in 0:9) {
  x <- arima.sim(list(order=c(2,0,0), ar=c(0.6,-0.5)), n = 100, innov = rnorm(100,mu,var),n.start = 100)
  data.set[((i*100)+1):(100*(i+1))] <- ts(x)
  var = var + 0.1
}

  dataset_2[[r]] <- data.set

}

# Data Set 3
dataset_3 <- list()
for(r in 1:100){
  
set.seed(r)
mu1 <- 0.1
mu2 <- -0.01
var1 <- 0.1
var2 <- 0.2
x1 <- rnorm( 100,mu1,var1) 
x2 <- rnorm( 100,mu2,var2) 
x3 <- rnorm( 100,mu1,var1) 
x4 <- rnorm( 100,mu2,var2) 
x5 <- rnorm( 100,mu1,var1)
x6 <- rnorm( 100,mu1,var2)
x7 <- rnorm( 100,mu1,var1)
x8 <- rnorm( 100,mu1,var2)
x9 <- rnorm( 100,mu1,var1)
x10 <- rnorm( 100,mu1,var2)
data.set <- ts(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
dataset_3[[r]] <- data.set

}


a <- cbind(dataset_1, dataset_2, dataset_3)


# write.csv(a,"SIMULATIONS AFTER COMMENTS.csv")

# Create a blank workbook
data_wb <- createWorkbook()

# Add some sheets to the workbook
addWorksheet(data_wb, "Data Set 1")
addWorksheet(data_wb, "Data Set 2")
addWorksheet(data_wb, "Data Set 3")


# Write the data to the sheets
writeData(data_wb, sheet = "Data Set 1", x = dataset_1)
writeData(data_wb, sheet = "Data Set 2", x = dataset_2)
writeData(data_wb, sheet = "Data Set 3", x = dataset_3)

# Export the file
saveWorkbook(data_wb, "SIMULATION SHEETS.xlsx")


