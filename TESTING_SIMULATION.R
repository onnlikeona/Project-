library(tidyverse)
library(astsa)
library(timeSeries)
library(forecast)
library(depmixS4)
source("SIMULATIONS.R")
source("rSTARS.R")

# Data Set 1 Shift in Mean
hmm <- depmix(data.set1 ~ 1, family = gaussian(), nstates = 5, data=data.frame(returns=data.set1))
hmmfit <- fit(hmm, verbose = FALSE)
post_probs <- posterior(hmmfit)
a <- AIC(hmmfit,k=1:120)
b <- BIC(hmmfit)

layout(1:2)
plot.ts(data.set1, main="Simulated Data with Shifts in Mean", xlab="", ylab="")
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')

Time <- as.character(seq(from = as.Date("2000-01-01"),by = "day",length.out = 1000))
x <- data.frame(cbind(Time,Returns=data.set1[1:1000]))
x$Returns <- as.numeric(x$Returns)
rstars(x,l.cutoff = 10,pValue = 0.05,save.path = getwd(),preWhitening = TRUE)

# Data Set 2 Shift in Variance
hmm <- depmix(data.set2 ~ 1, family = gaussian(), nstates = 10, data=data.frame(returns=data.set2))
hmmfit <- fit(hmm, verbose = FALSE)
post_probs <- posterior(hmmfit)

layout(1:2)
plot.ts(data.set2, main="Simulated Data with Shifts in Variance", xlab="", ylab="")
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')

Time <- as.character(seq(from = as.Date("2000-01-01"),by = "day",length.out = 1000))
x <- data.frame(cbind(Time,Returns=data.set2))
x$Returns <- as.numeric(x$Returns)
rstars(x,l.cutoff = 10,pValue = 0.01,save.path = getwd(),preWhitening = TRUE)

# Data Set 3 Shift in Mean and Variance
hmm <- depmix(data.set3 ~ 1, family = gaussian(), nstates = 2, data=data.frame(returns=data.set3))
hmmfit <- fit(hmm, verbose = FALSE)
post_probs <- posterior(hmmfit)

layout(1:2)
plot.ts(data.set3, main="Simulated Data with Shifts in Variance", xlab="", ylab="")
abline(v=c(100,200,300,400,500,600,700,800,900),lty=2,col="steelblue")
matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')

Time <- as.character(seq(from = as.Date("2000-01-01"),by = "day",length.out = 1000))
x <- data.frame(cbind(Time,Returns=data.set3))
x$Returns <- as.numeric(x$Returns)
rstars(x,l.cutoff = 10,pValue = 0.01,save.path = getwd(),preWhitening = TRUE)
