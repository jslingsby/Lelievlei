##############################################################################
######## Script to run JAGS model for mortality at Lilyvlei and predict for 2021
##############################################################################
######## Compiled by Jasper Slingsby 2018
######## Last edited: 19 Nov 2018
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
library(ggpubr)


#Set working directory
if(Sys.getenv("USER")=='jasper') setwd("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest")

#load("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_21Dec2016.Rdata")

load("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_19Nov2018.Rdata")


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
  for( o in 2 : 3 ) {
  betaY[o] ~ dnorm( 0 , 1.0E-3 )
  }
  
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
MortDatAll_trim <- MortDatAll[-which(MortDatAll$Species%in%names(which(summary(as.factor(datAll$Species))<15))),] #only 39 stems!!!
MortDatAll_trim <- droplevels(MortDatAll_trim)
tmeans_trim <- tmeans[-which(tmeans$Species%in%names(which(summary(as.factor(datAll$Species))<15))),]
tmeans_trim <- droplevels(tmeans_trim)

#Censor 30% of 2011 data for validation
samp <- sample(which(MortDatAll_trim$Year == 2011), round(length(which(MortDatAll_trim$Year == 2011))*.3, digits = 0), replace = F)

sampMort <- MortDatAll_trim$Mortality[samp]

MortDatAll_trim$Mortality[samp] <- NA

#Calculate plot means - how to code plot*year? - or just make plot a random effect?
StemD <- rowMeans(cbind(datSum$Stem_Density[1:108], datSum$Stem_Density[109:216]))
BA <- rowMeans(cbind(datSum$Basal_Area[1:108], datSum$Basal_Area[109:216]))

#Create data list
hdat1<-list(Mortality = MortDatAll_trim$Mortality, 
            Diameter =  (MortDatAll_trim$Diameter - mean(MortDatAll_trim$Diameter)), 
            Year = MortDatAll_trim$Year, 
            Condition = MortDatAll_trim$Condition, 
            Plot = MortDatAll_trim$Plot, 
            Species = MortDatAll_trim$Species, 
            N = nrow(MortDatAll_trim),
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
estInits<-with(MortDatAll_trim, glm(Mortality~Diameter+Year+Condition+Plot+Species, family=binomial(logit)))
#estInits

inits<-list(
  list(alpha = estInits$coef[1] + .5,
       betaD = estInits$coef[2] + .1,
       betaY = c(NA, estInits$coef[3:4] + .1),
       betaC = c(NA, estInits$coef[5:8] + .1),
       betaP = c(NA, estInits$coef[9:115] + .1),
       betaS = c(NA, estInits$coef[116:133] + .1)), #to 139 for all species
  list(alpha = estInits$coef[1] - .5,
       betaD = estInits$coef[2] - .1,
       betaY = c(NA, estInits$coef[3:4] - .1),
       betaC = c(NA, estInits$coef[5:8] - .1),
       betaP = c(NA, estInits$coef[9:115] - .1),
       betaS = c(NA, estInits$coef[116:133] - .1)) #to 139 for all species      
)

###Run model
#parameters<-c("alpha", "betaD", "betaY", "betaC", "betaP", "betaS", "plot.beta", "spp.beta") #, "plot.mu", "plot.tau", "spp.mu", "spp.tau", "mu") #add/remove "mu"?

parameters<-c("alpha", "betaD", "betaY", "betaC", "betaP", "betaS", "plot.beta", "spp.beta", "Mortality", "mu")

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
hldat <- as.data.frame(hlogmod$BUGSoutput$summary[,c(1,3,7)]) #point estimates with their 95% credible intervals

##########################################################
###Plotting

###Validation against test data

# boxplot(colMeans(hlogmod$BUGSoutput$sims.list$Mortality[,samp]) ~ sampMort)
# accuracy <- table(hlogmod$BUGSoutput$sims.list$Mortality[1,samp], sampMort)
# sum(diag(accuracy))/sum(accuracy)

auchb <- function(x, samp, sampMort){auc(roc(x[samp], sampMort))}
#auchb(hlogmod$BUGSoutput$sims.list$Mortality[1,], samp, sampMort)

mortAUC <- apply(hlogmod$BUGSoutput$sims.list$Mortality, MARGIN = 1, FUN = "auchb", samp = samp, sampMort = sampMort)

a <- ggplot(data.frame(AUC = mortAUC), aes(AUC)) 
a + geom_histogram()

###Predicted mortality by 2021

MT <- function(nm, t) {(1-((length(nm) - sum(nm, na.rm = T))/length(nm))^(1/t))*100} #Where n0 = length(nm) is the initial number of trees and nm is the vector indicating whether the tree died in this interval

mort2021 <- apply(hlogmod$BUGSoutput$sims.list$Mortality[,which(MortDatAll_trim$Year == 2021)], MARGIN = 1, FUN = "MT", t = 10)

m <- ggplot(data.frame(Mortality = mort2021), aes(Mortality)) 
m + geom_histogram()

# ##[Mortality]
# ## select betas of interest
# plot.dat <- hldat[grep("Mortality", rownames(hldat)),]
# #Obs <- MortDatAll_trim %>% group_by(Condition) %>% summarise(mort = mean(Mortality))
# #plot.dat <- log(plot.dat)
# plot.dat <- cbind(plot.dat, MortDatAll_trim[,c("Year","Mortality","Species","Condition")])
# 
# pd2021 <- plot.dat[plot.dat$Year == "2021",]
# 
# boxplot(mean ~ Species, data = pd2021)
# 
# ## order the data by the factor scores, for better visualization
# plot.dat <- plot.dat[order(plot.dat$mean), ]

# ## order the observation IDs as well (seems redundant, but is necessary)
# plot.dat$Condition <- factor(x = as.character(plot.dat$Condition), levels = as.character(plot.dat$Condition))
# 
# ## make plots
# sck <- ggplot(plot.dat, aes(x = Condition, y = mean, ymin = `2.5%`, ymax = `97.5%`))
# (sck <- sck + geom_pointrange() + ylab("Estimate"))
# 
# scl <- ggplot(plot.dat, aes(x = mort, y = mean))
# (scl <- scl + geom_point() + geom_segment(aes(x = mort, xend = mort, y = `2.5%`, yend = `97.5%`)) + ylab("Estimate") + xlab("Average mortality by condition"))
# 
# ggarrange(sck,scl, ncol = 2)