##############################################################################
######## Script to run JAGS model for mortality at Lilyvlei
##############################################################################
######## Compiled by Jasper Slingsby 2016
######## Last edited: 14 Nov 2018
######## Data: /Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/Lelievlei PSPs g01a&b 201208.xlsx
##############################################################################
###Steps:
###1) 
##############################################################################

#/usr/local/bin/jags

library(MCMCpack)
library(rjags)
library(R2jags)
library(coda)
library(ggplot2)
library(GGally)
library(reshape2)
library(dplyr)


#Set working directory
if(Sys.getenv("USER")=='jasper') setwd("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest")

load("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_21Dec2016.Rdata")

###############################################################################
###Hierarchical logistic regression model to explain mortality across both intervals
###############################################################################

cat(
  "model {
  
  ##Model
  for( i in 1 : N ) {
  Mortality[i] ~ dbern(mu[i])
  mu[i]<-1/(1+exp(-(alpha + betaD*Diameter[i] + betaY[Year[i]] + betaC[Condition[i]] + betaP[Plot[i]] + betaS[Species[i]])))
  }
  
  ##Priors
  #intercept(s)
  alpha ~ dnorm( 0 , 1.0E-3 )
  betaY[1] <- 0
  betaC[1] <- 0
  betaP[1] <- 0
  betaS[1] <- 0
  #plot.mu[1] <-0
  #spp.mu[1] <-0

  #Diameter
  betaD ~ dnorm( 0 , 1.0E-3 )
  
  #Year
  betaY[2] ~ dnorm( 0 , 1.0E-3 )
  
  #Condition
  for( j in 2 : 5 ) {
  betaC[j] ~ dnorm( 0 , 1.0E-3 )
  }
  
  #Plot - hierarchical - specify precision at the same time
  for( k in 2 : 108 ) {
  betaP[k] ~ dnorm( plot.mu[k] , plot.tau )
  }
  
  #Species - hierarchical
  for( l in 2 : 19 ) {
  betaS[l] ~ dnorm( spp.mu[l] , spp.tau )
  }

  ##Hierarchical regression
  plot.mu <- plotdat %*% plot.beta
  spp.mu <- sppdat %*% spp.beta

  #Prior on precisions
  plot.tau ~ dgamma( 0.01, 0.01 )
  spp.tau ~ dgamma( 0.01, 0.01 )

  ##Hierarchical beta priors
  for (m in 1:nplotBeta) {
  plot.beta[m] ~ dnorm( 0, 1.0E-3 )
  }

  for (n in 1:nsppBeta) {
  spp.beta[n] ~ dnorm( 0, 1.0E-3 )
  }
  
  }", file="hlogistic.jag"
)

###Set up data

#Create a trimmed dataset with only the abundant species
MortDat_trim <- MortDat[-which(MortDat$Species%in%names(which(summary(as.factor(datAll$Species))<15))),] #only 39 stems!!!
MortDat_trim <- droplevels(MortDat_trim)
tmeans_trim <- tmeans[-which(tmeans$Species%in%names(which(summary(as.factor(datAll$Species))<15))),]
tmeans_trim <- droplevels(tmeans_trim)

#Calculate plot means - how to code plot*year? - or just make plot a random effect?
StemD <- rowMeans(cbind(datSum$Stem_Density[1:108], datSum$Stem_Density[109:216]))
BA <- rowMeans(cbind(datSum$Basal_Area[1:108], datSum$Basal_Area[109:216]))

#Create data list
hdat1<-list(Mortality = MortDat_trim$Mortality, 
           Diameter =  (MortDat_trim$Diameter - mean(MortDat_trim$Diameter)), 
           Year = MortDat_trim$Year, 
           Condition = MortDat_trim$Condition, 
           Plot = MortDat_trim$Plot, 
           Species = MortDat_trim$Species, 
           N = nrow(MortDat_trim),
           plotdat = data.frame(
             Stem_Density = (StemD-mean(StemD))/sd(StemD),
             Basal_Area = (BA-mean(BA))/sd(BA)),
           sppdat = data.frame(
             Leaf_Area = (tmeans_trim$Leaf_Area - mean(tmeans_trim$Leaf_Area))/sd(tmeans_trim$Leaf_Area),
             SLA = (tmeans_trim$SLA - mean(tmeans_trim$SLA))/sd(tmeans_trim$SLA),
             Max_Height = (tmeans_trim$Maximum_Height - mean(tmeans_trim$Maximum_Height))/sd(tmeans_trim$Maximum_Height)),
#             MAT = (tmeans_trim$bio1 - mean(tmeans_trim$bio1))/sd(tmeans_trim$bio1)),
           nplotBeta = 2,
           nsppBeta = 3)

#Wood_Density = (tmeans_trim$Wood_Density - mean(tmeans_trim$Wood_Density))/sd(tmeans_trim$Wood_Density),

###Set initial values based on a GLM (allow hierarchical initials to be estimated by JAGS?)
estInits<-with(MortDat_trim, glm(Mortality~Diameter+Year+Condition+Plot+Species, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5,
       betaD = estInits$coef[2] + .1,
       betaY = c(NA, estInits$coef[3] + .1),
       betaC = c(NA, estInits$coef[4:7] + .1),
       betaP = c(NA, estInits$coef[8:114] + .1),
       betaS = c(NA, estInits$coef[115:132] + .1)), #to 139 for all species
  list(alpha = estInits$coef[1] - .5,
       betaD = estInits$coef[2] - .1,
       betaY = c(NA, estInits$coef[3] - .1),
       betaC = c(NA, estInits$coef[4:7] - .1),
       betaP = c(NA, estInits$coef[8:114] - .1),
       betaS = c(NA, estInits$coef[115:132] - .1)) #to 139 for all species      
)

###Run model
parameters<-c("alpha", "betaD", "betaY", "betaC", "betaP", "betaS", "plot.beta", "spp.beta") #, "plot.mu", "plot.tau", "spp.mu", "spp.tau")

hlogmod<- jags(data = hdat1,
              inits = inits,
              parameters.to.save = parameters,
              model.file = "hlogistic.jag",
              n.chains = 2,
              n.iter = 5000,
              n.burnin = 2000,
              n.thin = 1)


###Explore output
hlogmod
plot(hlogmod)
#plot(as.mcmc(hlogmod), ask=T)
#gelman.plot(as.mcmc(hlogmod))

lapply(hlogmod$BUGSoutput$mean, FUN=exp) #exponentiated point estimates
exp(hlogmod$BUGSoutput$summary[,c(1,3,7)]) #point estimates with their 95% credible intervals

###############################################################################
###HLM to explain mortality for each interval separately
###############################################################################

cat(
  "model {
  
  ##Model
  for( i in 1 : N ) {
  Mortality[i] ~ dbern(mu[i])
  mu[i]<-1/(1+exp(-(alpha + betaD*Diameter[i] + betaC[Condition[i]] + betaP[Plot[i]] + betaS[Species[i]])))
  }
  
  ##Priors
  #intercept(s)
  alpha ~ dnorm( 0 , 1.0E-3 )
  betaC[1] <- 0
  betaP[1] <- 0
  betaS[1] <- 0
  
  #Diameter
  betaD ~ dnorm( 0 , 1.0E-3 )

  #Condition
  for( j in 2 : 5 ) {
  betaC[j] ~ dnorm( 0 , 1.0E-3 )
  }
  
  #Plot - hierarchical - specify precision at the same time
  for( k in 2 : 108 ) {
  betaP[k] ~ dnorm( plot.mu[k] , plot.tau )
  }
  
  #Species - hierarchical
  for( l in 2 : 19 ) {
  betaS[l] ~ dnorm( spp.mu[l] , spp.tau )
  }
  
  ##Hierarchical regression
  plot.mu <- plotdat %*% plot.beta
  spp.mu <- sppdat %*% spp.beta
  
  #Prior on precisions
  plot.tau ~ dgamma( 0.01, 0.01 )
  spp.tau ~ dgamma( 0.01, 0.01 )
  
  ##Hierarchical beta priors
  for (m in 1:nplotBeta) {
  plot.beta[m] ~ dnorm( 0, 1.0E-3 )
  }
  
  for (n in 1:nsppBeta) {
  spp.beta[n] ~ dnorm( 0, 1.0E-3 )
  }
  
  }", file="hylogistic.jag"
)

###Set up and run 2001 data
MortDat_trim2001 <- MortDat_trim[which(MortDat_trim$Year==2001),]

#Create data list
hdat2001<-list(Mortality = MortDat_trim2001$Mortality, 
            Diameter =  (MortDat_trim2001$Diameter - mean(MortDat_trim2001$Diameter)), 
            Condition = MortDat_trim2001$Condition, 
            Plot = MortDat_trim2001$Plot, 
            Species = MortDat_trim2001$Species, 
            N = nrow(MortDat_trim2001),
            plotdat = data.frame(
              Stem_Density = (datSum$Stem_Density[1:108]-mean(datSum$Stem_Density[1:108]))/sd(datSum$Stem_Density[1:108]),
              Basal_Area = (datSum$Basal_Area[1:108]-mean(datSum$Basal_Area[1:108]))/sd(datSum$Basal_Area[1:108])),
            sppdat = data.frame(
              Leaf_Area = (tmeans_trim$Leaf_Area - mean(tmeans_trim$Leaf_Area))/sd(tmeans_trim$Leaf_Area),
              Wood_Density = (tmeans_trim$Wood_Density - mean(tmeans_trim$Wood_Density))/sd(tmeans_trim$Wood_Density),
              SLA = (tmeans_trim$SLA - mean(tmeans_trim$SLA))/sd(tmeans_trim$SLA)),
            nplotBeta = 2,
            nsppBeta = 3)

###Set initial values based on a GLM (allow hierarchical initials to be estimated by JAGS?)
estInits<-with(MortDat_trim2001, glm(Mortality~Diameter+Condition+Plot+Species, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5,
       betaD = estInits$coef[2] + .1,
       betaC = c(NA, estInits$coef[3:6] + .1),
       betaP = c(NA, estInits$coef[7:113] + .1),
       betaS = c(NA, estInits$coef[114:131] + .1)), #to 139 for all species
  list(alpha = estInits$coef[1] - .5,
       betaD = estInits$coef[2] - .1,
       betaC = c(NA, estInits$coef[3:6] - .1),
       betaP = c(NA, estInits$coef[7:113] - .1),
       betaS = c(NA, estInits$coef[114:131] - .1)) #to 139 for all species      
)

###Run model
parametersy<-c("alpha", "betaD", "betaC", "betaP", "betaS", "plot.beta", "spp.beta") #, "plot.mu", "plot.tau", "spp.mu", "spp.tau")

hlogmod2001<- jags(data = hdat2001,
               inits = inits,
               parameters.to.save = parametersy,
               model.file = "hylogistic.jag",
               n.chains = 2,
               n.iter = 5000,
               n.burnin = 2000,
               n.thin = 1)


###Explore output
hlogmod2001
plot(hlogmod2001)
#plot(as.mcmc(hlogmod), ask=T)

#lapply(hlogmod$BUGSoutput$mean, FUN=exp) #exponentiated point estimates
#exp(hlogmod$BUGSoutput$summary[,c(1,3,7)]) #point estimates with their 95% credible intervals

###Set up and run 2011 data
MortDat_trim2011 <- MortDat_trim[which(MortDat_trim$Year==2011),]

#Create data list
hdat2011<-list(Mortality = MortDat_trim2011$Mortality, 
               Diameter =  (MortDat_trim2011$Diameter - mean(MortDat_trim2011$Diameter)), 
               Condition = MortDat_trim2011$Condition, 
               Plot = MortDat_trim2011$Plot, 
               Species = MortDat_trim2011$Species, 
               N = nrow(MortDat_trim2011),
               plotdat = data.frame(
                 Stem_Density = (datSum$Stem_Density[109:216]-mean(datSum$Stem_Density[109:216]))/sd(datSum$Stem_Density[109:216]),
                 Basal_Area = (datSum$Basal_Area[109:216]-mean(datSum$Basal_Area[109:216]))/sd(datSum$Basal_Area[109:216])),
               sppdat = data.frame(
                 Leaf_Area = (tmeans_trim$Leaf_Area - mean(tmeans_trim$Leaf_Area))/sd(tmeans_trim$Leaf_Area),
                 Wood_Density = (tmeans_trim$Wood_Density - mean(tmeans_trim$Wood_Density))/sd(tmeans_trim$Wood_Density),
                 SLA = (tmeans_trim$SLA - mean(tmeans_trim$SLA))/sd(tmeans_trim$SLA)),
               nplotBeta = 2,
               nsppBeta = 3)

###Set initial values based on a GLM (allow hierarchical initials to be estimated by JAGS?)
estInits<-with(MortDat_trim2011, glm(Mortality~Diameter+Condition+Plot+Species, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5,
       betaD = estInits$coef[2] + .1,
       betaC = c(NA, estInits$coef[3:6] + .1),
       betaP = c(NA, estInits$coef[7:113] + .1),
       betaS = c(NA, estInits$coef[114:131] + .1)), #to 139 for all species
  list(alpha = estInits$coef[1] - .5,
       betaD = estInits$coef[2] - .1,
       betaC = c(NA, estInits$coef[3:6] - .1),
       betaP = c(NA, estInits$coef[7:113] - .1),
       betaS = c(NA, estInits$coef[114:131] - .1)) #to 139 for all species      
)

###Run model
parametersy<-c("alpha", "betaD", "betaC", "betaP", "betaS", "plot.beta", "spp.beta") #, "plot.mu", "plot.tau", "spp.mu", "spp.tau")

hlogmod2011<- jags(data = hdat2011,
                   inits = inits,
                   parameters.to.save = parametersy,
                   model.file = "hylogistic.jag",
                   n.chains = 2,
                   n.iter = 5000,
                   n.burnin = 2000,
                   n.thin = 1)


###Explore output
hlogmod2011
plot(hlogmod2011)
#plot(as.mcmc(hlogmod), ask=T)

###############################################################################

###############################################################################


###Explore some of the ideas
par(mfrow=c(3,1))
hist(MortDat$Diameter, xlim = c(0,200), ylim=c(0,3000), main = "All", xlab = "")
hist(MortDat$Diameter[which(MortDat$Mortality==0)], xlim = c(0,200), ylim=c(0,3000), main = "Survived", xlab = "")
hist(MortDat$Diameter[which(MortDat$Mortality==1)], xlim = c(0,200), ylim=c(0,300), main = "Died", xlab = "")

pd <- droplevels(Plot_dynamics[which(Plot_dynamics$Period %in% c("1991-2001", "2001-2011")),])
boxplot(Mortality ~ Period, data = pd)
summary(aov(Mortality ~ Period, data = pd))
