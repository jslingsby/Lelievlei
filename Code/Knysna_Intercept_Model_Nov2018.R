##############################################################################
######## Script to run JAGS "intercept" model for mortality at Lilyvlei
##############################################################################
######## Compiled by Jasper Slingsby 2018
######## Last edited: 20 Nov 2018
######## Data: /Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/Lelievlei PSPs g01a&b 201208.xlsx
##############################################################################
######## 
###Steps:
###1) 
##############################################################################


#/usr/local/bin/jags

library(MCMCpack)
library(rjags)
library(R2jags)
library(ggplot2)
library(GGally)
library(reshape2)
library(dplyr)


#Set working directory
if(Sys.getenv("USER")=='jasper') setwd("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest")

load("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_21Dec2016.Rdata")

###############################################################################
###Simple logistic regression model to explain mortality
###############################################################################

cat(
  "model {
  
  #Model
  for( i in 1 : N ) {
  Mortality[i] ~ dbern(mu[i])
  mu[i]<-1/(1+exp(-alpha))
  }
  
  
  #Priors
  #intercept(s)
  alpha ~ dnorm( 0 , 1.0E-3 )
  
  }", file="intlogistic.jag"
)

###Set up data

#Create data list
datInt<-list(Mortality = MortDat$Mortality, 
           N = nrow(MortDat))

###Set initial values based on a GLM
estInits<-with(MortDat, glm(Mortality~1, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5),
  list(alpha = estInits$coef[1] - .5)       
)

###Run model
parameters<-c("alpha")

logmod<- jags(data = datInt,
              inits = inits,
              parameters.to.save = parameters,
              model.file = "intlogistic.jag",
              n.chains = 2,
              n.iter = 5000,
              n.burnin = 2000,
              n.thin = 1)


###Explore output
logmod
plot(logmod)
plot(as.mcmc(logmod), ask=T)

lapply(logmod$BUGSoutput$mean, FUN=exp) #exponentiated point estimates
exp(logmod$BUGSoutput$summary[,c(1,3,7)]) #exponentiated point estimates with their 95% credible intervals
