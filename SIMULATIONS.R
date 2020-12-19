library(tidyverse)
library(astsa)
library(gdata)
set.seed(696)
setwd("~/Documents/2020/Semester 2/Project/Application of Methods/SIMULATIONS")
# Data Set 1
f <- function(x){
  if((x %% 2) == 0) {
    return(1)
  } else {
    return(0)
  }
}

set.seed(696)
mu = 0
mu1 = -0.1
mu2 = 1.1
var = 0.1
data.set1 <- c()
for (i in 0:9) {
  mu1 = mu + 0.2*i
  mu2 = mu + 0.1*i
  mu = f(i)*mu1+(1-f(i))*mu2
  x <- arima.sim(list(order=c(2,0,0),ar=c(0.6,-0.5)),n=100,innov=rnorm(100,mu,1),n.start=100)
  data.set1[((i*100)+1):(100*(i+1))] <- ts(x)
}

# Data Set 2
set.seed(696)
mu = 0
mu1 = -0.1
mu2 = 1.1
data.set2 <- c()
var = 0.1
for (i in 0:9) {
  x <- arima.sim(list(order=c(2,0,0), ar=c(0.6,-0.5)), n = 100, innov = rnorm(100,mu,var),n.start = 100)
  data.set2[((i*100)+1):(100*(i+1))] <- ts(x)
  var = var + 0.1
}
data.set1 <- as.ts(data.set1)
data.set2 <- as.ts(data.set2)


# Data Set 3
set.seed(696)
mu1 <- 0.1
mu2 <- -0.01
var1 <- 0.1
var2 <- 0.2
x1 <- rnorm( 200,mu1,var1) 
x2 <- rnorm( 200,mu2,var2) 
x3 <- rnorm( 200,mu1,var1) 
x4 <- rnorm( 200,mu2,var2) 
x5 <- rnorm( 200,mu1,var1)
data.set3 <- ts(c(x1,x2,x3,x4,x5))


a <- cbind(data.set1,data.set2,data.set3)
colnames(a) <- c("AR(2)","AR(4)","WiteNoise")
write.csv(a,"SIMULATIONS.csv")

keep(data.set1,
     data.set2,
     data.set3,
     sure = TRUE)
