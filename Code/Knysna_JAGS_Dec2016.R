##############################################################################
######## Script to run JAGS model for mortality at Lilyvlei
##############################################################################
######## Compiled by Jasper Slingsby 2016
######## Last edited: 21 Dec 2016
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
        mu[i]<-1/(1+exp(-(alpha + betaD*Diameter[i] + betaY[Year[i]] + betaC[Condition[i]] + betaP[Plot[i]] + betaS[Species[i]])))
      }


  #Priors
      #intercept(s)
      alpha ~ dnorm( 0 , 1.0E-3 )
      betaY[1] <- 0
      betaC[1] <- 0
      betaP[1] <- 0
      betaS[1] <- 0

      #Diameter
      betaD ~ dnorm( 0 , 1.0E-3 )

      #Year
      betaY[2] ~ dnorm( 0 , 1.0E-3 )

      #Condition
      for( j in 2 : 5 ) {
      betaC[j] ~ dnorm( 0 , 1.0E-3 )
      }

      #Plot
      for( k in 2 : 108 ) {
      betaP[k] ~ dnorm( 0 , 1.0E-3 )
      }

      #Species
      for( l in 2 : 26 ) {
      betaS[l] ~ dnorm( 0 , 1.0E-3 )
      }

            }", file="logistic.jag"
)

###Set up data

#Create data list
dat1<-list(Mortality = MortDat$Mortality, 
           Diameter = (MortDat$Diameter - mean(MortDat$Diameter)), 
           Year = MortDat$Year, 
           Condition = MortDat$Condition, 
           Plot = MortDat$Plot, 
           Species = MortDat$Species, 
           N = nrow(MortDat))

###Set initial values based on a GLM
estInits<-with(MortDat, glm(Mortality~Diameter+Year+Condition+Plot+Species, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5,
       betaD = estInits$coef[2] + .1,
       betaY = c(NA, estInits$coef[3] + .1),
       betaC = c(NA, estInits$coef[4:7] + .1),
       betaP = c(NA, estInits$coef[8:114] + .1),
       betaS = c(NA, estInits$coef[115:139] + .1)),
  list(alpha = estInits$coef[1] - .5,
       betaD = estInits$coef[2] - .1,
       betaY = c(NA, estInits$coef[3] - .1),
       betaC = c(NA, estInits$coef[4:7] - .1),
       betaP = c(NA, estInits$coef[8:114] - .1),
       betaS = c(NA, estInits$coef[115:139] - .1))       
)

###Run model
parameters<-c("alpha", "betaD", "betaY", "betaC", "betaP", "betaS")

logmod<- jags(data = dat1,
              inits = inits,
              parameters.to.save = parameters,
              model.file = "logistic.jag",
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
