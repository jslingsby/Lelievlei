##############################################################################
######## Script to analyse Lilyvlei data
##############################################################################
######## Compiled by Jasper Slingsby 2016
######## Last edited: 21 Dec 2016
######## Data: /Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/Lelievlei PSPs g01a&b 201208.xlsx
##############################################################################
######## 
###Steps:
###1) 
##############################################################################

library(gdata)
library(picante)
library(FD)
library(gplots)
library(MCMCpack)
library(MCMCglmm)
library(ggplot2)
library(GGally)
library(reshape2)
library(dplyr)


#Set working directory
if(Sys.getenv("USER")=='jasper') setwd("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest")

load("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_21Dec2016.Rdata")

###############################################################################
###Explore stand level stats
###############################################################################

t(AnnualDynamics)

StemCondition

###############################################################################
###Run analyses of plot summaries
###############################################################################

#Plot Dynamics
ggpairs(Plot_dynamics, columns = c(3:5,2), ggplot2::aes(colour=Period))

#Species Dynamics - THESE NEED CHECKING...
ggpairs(Species_dynamics, columns = c(3:5,2), ggplot2::aes(colour=Period))

#Plot Stem density, BA, AGB
ggpairs(datSum, columns = c(3:6,2), ggplot2::aes(colour=Year))

#Analyze stem density, BA, AGB with Bayesian repeated measures
fitStem_Density <- MCMCregress(Stem_Density ~ Year, random = ~Plot, data = datSum, mcmc=10000, burnin=3000, thin = 10)
summary(fitStem_Density)
plot(fitStem_Density)

fitBA <- MCMCregress(Basal_Area ~ Year, random = ~Plot, data = datSum, mcmc=10000, burnin=3000, thin = 10)
summary(fitBA)
plot(fitBA)

fitAGB_Constant <- MCMCregress(AGB_Constant ~ Year, random = ~Plot, data = datSum, mcmc=10000, burnin=3000)
summary(fitAGB_Constant)
plot(fitAGB_Constant)

fitAGB_Specific <- MCMCregress(AGB_Specific ~ Year, random = ~Plot, data = datSum, mcmc=10000, burnin=3000)
summary(fitAGB_Specific)
plot(fitAGB_Specific)


###############################################################################
###Run mortality analyses (do for recruitment too???)
###############################################################################

###Explore evidence for change in mortality through time - i.e. a survival analysis or similar
#THIS SHOULD BE SOMETHING LIKE A NEGATIVE BINOMIAL - SEE van Mantgem et al. 2009 or Peng et al. 2011
MortDat <- merge(MortDat, datSum, all.x = T)

#THIS IS WRONG...#
fit0 <- MCMCglmm(Mortality ~ Year, random = ~Plot, data = MortDat, family = "categorical", prior=prior, verbose=TRUE, pr=TRUE, nitt=13000, burnin=3000)
summary(fit0)
plot(fit0)


###Explore influences on mortality (by individual, species and site)

#Create a trimmed dataset with only the abundant species
MortDat_trim <- MortDat[-which(MortDat$Species%in%names(which(summary(as.factor(datAll$Species))<15))),] #only 39 stems!!!
MortDat_trim <- droplevels(MortDat_trim)

fit1 <- MCMClogit(Mortality ~ Condition + Diameter + Year + SLA + Wood_Density, random = ~Plot, data = MortDat, mcmc=13000, burnin=3000) #, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
summary(fit1)
plot(fit1)

fit2 <- MCMClogit(Mortality ~ Condition + Diameter + Year + SLA + Wood_Density + Basal_Area + Stem_Density, random = ~Plot, data = MortDat, mcmc=13000, burnin=3000, thin = 10) #, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
summary(fit2)
plot(fit2)

fit3 <- MCMClogit(Mortality ~ Condition + Diameter + SLA + Wood_Density + Basal_Area + Stem_Density, random = ~Plot + Year, data = MortDat, mcmc=13000, burnin=3000, thin = 10) #, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
summary(fit3)
plot(fit3)

###Attempt at Jimmying a hierarchical model for species and traits (using data trimmed to more abundant species)

#Fit model with site and individual effects and species as a factor
fitSpecies <- MCMClogit(Mortality ~ Condition + Diameter + Species + Basal_Area + Stem_Density + Year, random = ~Plot, data = MortDat_trim, mcmc=55000, burnin=30000, thin = 55) #, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
summary(fitSpecies)
plot(fitSpecies)

#Extract coefficients for species
sppeff <- as.data.frame(fitSpecies[,7:24])
spf <- melt(sppeff, value.name = "Effect")
names(spf)[1] <- "Species"
spf$Species <- gsub("Species", "", spf$Species)
spf <- merge(spf, tmeans, all.x=T)

#Check boxplot of species effects
par(mar = c(5, 12, 4, 2) + 0.1)
boxplot(Effect ~ Species, data = spf, horizontal = T, las=1, cex.axis = .75)
abline(v=0, col="grey50")
par(mar = c(5, 4, 4, 2) + 0.1)

#Check ggpairs of species means ~ traits
x <- aggregate(spf$Effect, by = list(spf$Species), FUN = "mean")
names(x)[1] <- "Species"
x <- merge(x, tmeans, all.x = T)
ggpairs(x, columns = 2:13)

#Fit model with traits as predictors of the posterior distribution of species coefficients on mortality
fitspf <- MCMCregress(Effect ~ SLA + Wood_Density + Leaf_Area + LDMC, data = spf, mcmc=13000, burnin=3000, thin = 10) #, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
summary(fitspf)
plot(fitspf)

###############################################################################
###################################THE END#####################################
###############################################################################

#WARNING!!! This takes ~3 hours...
#fit2 <- MCMChlogit(Mortality ~ Condition + Diameter + Year + SLA + Wood_Density, random = ~Plot, group = "Species", data = MortDat, mcmc=10000, burnin=7500, r=108, R=diag(c(1,rep(0.1,107)))) #Species removed ( )
#summary(fit2$mcmc)

###############################################################################
###Draw out community data matrices
###############################################################################

###Draw out 1991 community data matrix
dat<-datAll[!is.na(datAll$DBH2),]

#Samp by abundance
samp<-data.frame(dat$PLID, rep(1,dim(dat)[1]), dat$SPECIES) #select species occurrence and plot names
samp91<-sample2matrix(samp)

#Samp by basal area
samp<-data.frame(dat$PLID, pi*(dat$DBH2/2)^2, dat$SPECIES) #select species occurrence and plot names
samp91ba<-sample2matrix(samp)

#Samp by biomass
samp<-data.frame(dat$PLID, apply(as.matrix(dat$DBH2),2,constantB), dat$SPECIES) #select species occurrence and plot names
samp91bm<-sample2matrix(samp)


###Draw out 2001 community data matrix
#Samp by abundance
dat<-datAll[!is.na(datAll$DBH3),]
samp<-data.frame(dat$PLID, rep(1,dim(dat)[1]), dat$SPECIES) #select species occurrence and plot names
samp01<-sample2matrix(samp)

#Samp by basal area
samp<-data.frame(dat$PLID, pi*(dat$DBH3/2)^2, dat$SPECIES) #select species occurrence and plot names
samp01ba<-sample2matrix(samp)

#Samp by biomass
samp<-data.frame(dat$PLID, apply(as.matrix(dat$DBH3),2,constantB), dat$SPECIES) #select species occurrence and plot names
samp01bm<-sample2matrix(samp)

###Draw out 2011 community data matrix
dat<-datAll[!is.na(datAll$DBH4),]
samp<-data.frame(dat$PLID, rep(1,dim(dat)[1]), dat$SPECIES) #select species occurrence and plot names
samp11<-sample2matrix(samp)

#Samp by basal area
samp<-data.frame(dat$PLID, pi*(dat$DBH4/2)^2, dat$SPECIES) #select species occurrence and plot names
samp11ba<-sample2matrix(samp)

#Samp by biomass
samp<-data.frame(dat$PLID, apply(as.matrix(dat$DBH4),2,constantB), dat$SPECIES) #select species occurrence and plot names
samp11bm<-sample2matrix(samp)

#Calculate combined abundances/basal area/biomass by survey
sampAll<-rbind(colSums(samp91), colSums(samp01), colSums(samp11))
sampAllba<-rbind(colSums(samp91ba), colSums(samp01ba), colSums(samp11ba))
sampAllbm<-rbind(colSums(samp91bm), colSums(samp01bm), colSums(samp11bm))

###############################################################################
###Get trait data
###############################################################################
lf1<-read.xls("Data/datacomplete/leaf data_Knysna forest_Thomas Morris data.xlsx", sheet=1, stringsAsFactors =FALSE)
for(i in 1:dim(lf1)[1]) {lf1$X[i]<-paste(lf1[i,1:2], sep="", collapse=" ")} #Combine collection labels to allow merging

lf2<-read.xls("Data/datacomplete/leaf data_Knysna forest_Thomas Morris data.xlsx", sheet=2, stringsAsFactors =FALSE)
for(i in 1:dim(lf2)[1]) {lf2$X[i]<-paste(lf2[i,1:2], sep="", collapse=" ")} #Combine collection labels to allow merging

trt<-read.xls("Data/datacomplete/stem_data.xlsx", sheet=1, stringsAsFactors =FALSE)
for(i in 1:dim(trt)[1]) {trt$X[i]<-paste(trt[i,1:2], sep="", collapse=" ")}

lf1means<-aggregate(lf1[,5:9],by=list(lf1$X),"mean") #calculate species means for all quantitative data
colnames(lf1means)[1]<-"X"

#Merge 2 leaf spreadsheets and wood density
leaf<-merge(lf2[,3:7], lf1means)
trait<-merge(leaf, trt[,c(14,16)])

#Combine genus and species names
for(i in 1:dim(trait)[1]) {trait$X[i]<-paste(trait[i,2:3], sep="", collapse=" ")}

#Fix names to match stem data with plot data
trait$X<-gsub("Psydrax obovata","Psydrax obovata subsp. obovata", trait$X)
trait$X<-gsub("Psydrax ovata","Psydrax obovata subsp. obovata", trait$X)
trait$X<-gsub("Olea capensis subsp. Capensis", "Olea capensis subsp. capensis", trait$X) 
trait$X<-gsub("Olea capensis subsp. Macrocarpa", "Olea capensis subsp. macrocarpa", trait$X)
trait$X<-gsub("Apodytes dimidiata", "Apodytes dimidiata subsp. dimidiata", trait$X)
trait$X<-gsub("Cassine eucliformis", "Robsonodendron eucleiforme", trait$X)
trait$X<-gsub("Eleaodendron croceum", "Elaeodendron croceum", trait$X)
trait$X<-gsub("Diospyrus whyteana", "Diospyros whyteana", trait$X)
trait$X<-gsub("Occotea bullata", "Ocotea bullata", trait$X)

#remove alien Guava and juvenile Elaeodendron data
crap<-sort(unique(trait$X[-which(trait$X%in%spp[,3])])) #unmatched names - to remove
trait<-trait[-which(trait$X%in%crap),] #remove alien Guava and juvenile Elaeodendron data

#calculate additional traits
trait$SLA<-trait$Area..cm./trait$leaf_dry_wt
trait$LDMC<-trait$leaf_dry_wt/trait$leaf_wet_wt..g.

#calculate means and rename traits
tmeans<-aggregate(trait[,c(6:8,10:13)],by=list(trait$X),"mean") #calculate species means for all quantitative data
rownames(tmeans)<-tmeans[,1]
tmeans<-tmeans[,-1]
colnames(tmeans)<-c("Area_cm", "Length_cm", "Ave_Lf_Width_cm", "Leaf_Thickness_mm", "Wood_Density_kg/m3", "SLA", "LDMC")

###############################################################################
###Calculate FD metrics and plot results
###############################################################################

#source("Code/FDplots.R")
#source("Code/FDplots_ba.R")
#source("Code/FDplots_bm.R")

###############################################################################
###Species climate data analyses
###############################################################################

#source("Code/clim.R")

