##############################################################################
######## Script to prep Lilyvlei data
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

###############################################################################
###Get species names
###############################################################################

spp<-read.xls("Data/Lelievlei PSPs g01a&b 201208.xlsx", sheet=3, pattern="Number Code", strip.white=TRUE, stringsAsFactors=F)
spp<-spp[1:102,] #get rid of crap below data
spp[,3]<-trim(as.character(spp[,3])) #remove white spaces from end of names

###############################################################################
###Get abundance data
###############################################################################
datAll<-read.xls("Data/Lelievlei PSPs g01a&b 201208.xlsx", sheet=1, stringsAsFactors=F)

#Change tree codes for species names
datAll$SPECIESCODE <- datAll$SPECIES
datAll$SPECIES<-spp$Scientific.Name[match(as.numeric(datAll$SPECIES), as.numeric(spp$Number.Code))]
datAll<-datAll[-which(datAll$SPECIES=="Chionanthus foveolata subsp. foveolata"),] #Remove singleton

###############################################################################
###Calculate mortality and recruitment
###############################################################################
datRM <- data.frame(PLID = datAll$PLID, SPECIES = datAll$SPECIES, p1991 = !is.na(datAll$DBH2), p2001 = !is.na(datAll$DBH3), p2011 = !is.na(datAll$DBH4))
datRM$State <- apply(datRM[, 3:5], MARGIN = 1, FUN = paste0, collapse = ", ")

datAll$Mortality2001 <- NA
datAll$Mortality2001[which(datRM$State %in% c("TRUE, TRUE, TRUE", "TRUE, TRUE, FALSE"))] <- 0
datAll$Mortality2001[which(datRM$State  == "TRUE, FALSE, FALSE")] <- 1

datAll$Mortality2011 <- NA
datAll$Mortality2011[which(datRM$State %in% c("TRUE, TRUE, TRUE", "FALSE, TRUE, TRUE"))] <- 0
datAll$Mortality2011[which(datRM$State  == "TRUE, TRUE, FALSE")] <- 1

#datAll$Recruitment2001 <- NA
#datAll$Recruitment2001[which(datRM$State  %in% c("FALSE, TRUE, TRUE", "FALSE, TRUE, FALSE"))] <- 1
#datAll$Recruitment2011 <- NA
#datAll$Recruitment2011[which(datRM$State  == "FALSE, FALSE, TRUE")] <- 1

datAll$recruit2001 <- as.numeric(datRM$State == "FALSE, TRUE, TRUE")
datAll$recruit2011 <- as.numeric(datRM$State == "FALSE, FALSE, TRUE")
datAll$recruit1991_2011 <- rowSums(cbind(datAll$recruit2001, datAll$recruit2011))

datAll$died2001 <- as.numeric(datRM$State == "TRUE, FALSE, FALSE")
datAll$died2011 <- as.numeric(datRM$State == "TRUE, TRUE, FALSE")
datAll$died1991_2011 <- rowSums(cbind(datAll$died2001, datAll$died2011))

###############################################################################
###Bin stem conditions - Dying > Unhealthy > Damaged > Leaning > Healthy
###############################################################################

###1991
datAll$Condition1991[which(rowSums(cbind(datAll$COND2A %in% c(0, 60:69), datAll$COND2B %in% c(0, 60:69)))==2)] <- "Healthy"
datAll$Condition1991[which(rowSums(cbind(datAll$COND2A %in% c(40:49), datAll$COND2B %in% c(40:49)))>0)] <- "Leaning"
datAll$Condition1991[which(rowSums(cbind(datAll$COND2A %in% c(31:33, 36, 38), datAll$COND2B %in% c(31:33, 36, 38)))>0)] <- "Damaged"
datAll$Condition1991[which(rowSums(cbind(datAll$COND2A %in% c(10:19, 34, 35, 37, 50:59, 70:90), datAll$COND2B %in% c(10:19, 34, 35, 37, 50:59, 70:90)))>0)] <- "Unhealthy"
datAll$Condition1991[which(rowSums(cbind(datAll$COND2A %in% c(20:29), datAll$COND2B %in% c(20:29)))>0)] <- "Dying"

###2001
datAll$Condition2001[which(rowSums(cbind(datAll$COND3A %in% c(0, 60:69), datAll$COND3B %in% c(0, 60:69)))==2)] <- "Healthy"
datAll$Condition2001[which(rowSums(cbind(datAll$COND3A %in% c(40:49), datAll$COND3B %in% c(40:49)))>0)] <- "Leaning"
datAll$Condition2001[which(rowSums(cbind(datAll$COND3A %in% c(31:33, 36, 38), datAll$COND3B %in% c(31:33, 36, 38)))>0)] <- "Damaged"
datAll$Condition2001[which(rowSums(cbind(datAll$COND3A %in% c(10:19, 34, 35, 37, 50:59, 70:90), datAll$COND3B %in% c(10:19, 34, 35, 37, 50:59, 70:90)))>0)] <- "Unhealthy"
datAll$Condition2001[which(rowSums(cbind(datAll$COND3A %in% c(20:29), datAll$COND3B %in% c(20:29)))>0)] <- "Dying"

###2011
datAll$Condition2011[which(rowSums(cbind(datAll$COND4A %in% c(0, 60:69), datAll$COND4B %in% c(0, 60:69)))==2)] <- "Healthy"
datAll$Condition2011[which(rowSums(cbind(datAll$COND4A %in% c(40:49), datAll$COND4B %in% c(40:49)))>0)] <- "Leaning"
datAll$Condition2011[which(rowSums(cbind(datAll$COND4A %in% c(31:33, 36, 38), datAll$COND4B %in% c(31:33, 36, 38)))>0)] <- "Damaged"
datAll$Condition2011[which(rowSums(cbind(datAll$COND4A %in% c(10:19, 34, 35, 37, 50:59, 70:90), datAll$COND4B %in% c(10:19, 34, 35, 37, 50:59, 70:90)))>0)] <- "Unhealthy"
datAll$Condition2011[which(rowSums(cbind(datAll$COND4A %in% c(20:29), datAll$COND4B %in% c(20:29)))>0)] <- "Dying"


###############################################################################
###Get trait data
###############################################################################

lf1<-read.xls("Data/datacomplete/leaf data_Knysna forest_Thomas Morris data.xlsx", sheet=1, stringsAsFactors =FALSE)
lf1$X<-apply(trim(lf1[,1:2]), MARGIN=1, FUN="paste0", collapse="_")

lf2<-read.xls("Data/datacomplete/leaf data_Knysna forest_Thomas Morris data.xlsx", sheet=2, stringsAsFactors =FALSE)
lf2$X<-apply(lf2[,1:2], MARGIN=1, FUN="paste0", collapse="_")
lf2$leaf_wet_wt..g. <- lf2$leaf_wet_wt..g./3
lf2$leaf_dry_wt <- lf2$leaf_dry_wt/3

trt<-read.xls("Data/datacomplete/stem_data.xlsx", sheet=1, stringsAsFactors =FALSE)
trt$X<-apply(trt[,1:2], MARGIN=1, FUN="paste0", collapse="_")
trt$StemWaterContent <- (trt$stem__wet_wt..g. - trt$stem_dry_wt..g.)/trt$stem__wet_wt..g.

lf1means<-aggregate(lf1[,5:9],by=list(lf1$X),"mean") #calculate species means for all quantitative data
colnames(lf1means)[1]<-"X"

#Merge 2 leaf spreadsheets and wood density
leaf<-merge(lf2[,3:7], lf1means)
trait<-merge(leaf, trt[,c(14,16,17)])

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
trait$LMA<-trait$leaf_dry_wt/trait$Area..cm.
trait$LDMC<-trait$leaf_dry_wt/trait$leaf_wet_wt..g.
trait$LeafWaterContent<-(trait$leaf_wet_wt..g. - trait$leaf_dry_wt)/trait$leaf_wet_wt..g.
trait$Density..Immertion..kg.m..3. <- trait$Density..Immertion..kg.m..3./1000 #convert wood density to g/cm^3

#calculate means and rename traits
tmeans<-aggregate(trait[,6:ncol(trait)],by=list(trait$X),"mean") #calculate species means for all quantitative data
#rownames(tmeans)<-tmeans[,1]
#tmeans<-tmeans[,-1]
colnames(tmeans)<-c("Species", "Leaf_Area", "Leaf_Length", "Ave_Leaf_Width", "Max_Leaf_Width", "Leaf_Thickness", "Wood_Density", "Stem_Water_Content", "SLA", "LMA", "LDMC", "Leaf_Water_Content")
#Units cm, cm, cm, cm, mm, kg/m^3, %, cm/g, g/cm, %

#add height data from POSA
height <- read.csv("heights.csv")
tmeans <- merge(tmeans, height, all.x=T)

#ggpairs(tmeans[,2:ncol(tmeans)])

###############################################################################
###Make final dataframes for analysis
###############################################################################

###Mortality
MortDat <- data.frame(Plot = datAll$PLID, Species	= factor(datAll$SPECIES), Year = factor(c(rep(2001, nrow(datAll)), rep(2011, nrow(datAll)))), Mortality = c(datAll$Mortality2001, datAll$Mortality2011), Condition = factor(c(datAll$Condition1991, datAll$Condition2001), levels = c("Healthy", "Leaning", "Damaged", "Unhealthy", "Dying")),	Diameter = c(datAll$DBH2, datAll$DBH3))
#Recruitment = c(datAll$Recruitment2001, datAll$Recruitment2011),

MortDat <- na.omit(MortDat)
MortDat <- merge(MortDat, tmeans)

###Mortality with 2011 data for 2021 prediction
MortDatAll <- data.frame(Plot = datAll$PLID, Species	= factor(datAll$SPECIES), Year = factor(c(rep(2001, nrow(datAll)), rep(2011, nrow(datAll)), rep(2021, nrow(datAll)))), Mortality = c(datAll$Mortality2001, datAll$Mortality2011, rep(NA, nrow(datAll))), Condition = factor(c(datAll$Condition1991, datAll$Condition2001,datAll$Condition2011), levels = c("Healthy", "Leaning", "Damaged", "Unhealthy", "Dying")),	Diameter = c(datAll$DBH2, datAll$DBH3, datAll$DBH4))
#Recruitment = c(datAll$Recruitment2001, datAll$Recruitment2011),

MortDatAll <- MortDatAll[-which(rowSums(is.na(MortDatAll)) > 1), ]
MortDatAll <- merge(MortDatAll, tmeans)

###AllDat
names(datAll)[which(names(datAll)=="SPECIES")] <- "Species"
datAll <- merge(datAll, tmeans)


###############################################################################
###Stem condition summary
###############################################################################

StemCondition <- data.frame(
  Condition_1991 = summary(factor(datAll$Condition1991, levels = c("Healthy", "Leaning", "Damaged", "Unhealthy", "Dying"))),
  Condition_2001 = summary(factor(datAll$Condition2001, levels = c("Healthy", "Leaning", "Damaged", "Unhealthy", "Dying"))),
  Condition_2011 = summary(factor(datAll$Condition2011, levels = c("Healthy", "Leaning", "Damaged", "Unhealthy", "Dying"))))

StemCondition

#Ave # of stems by condition by plot?

###############################################################################
###Calculate Dynamics
###############################################################################

###Set functions - Needs fixing? Or use other method?
RC <- function(n0, nr, t) {((((sum(decostand(n0, "pa", na.rm=T), na.rm=T) + sum(nr, na.rm = T))/sum(decostand(n0, "pa", na.rm=T), na.rm=T))^(1/t))-1)*100} #Where n0 are the initial dbh values for each stem (including NAs) and nr is the vector indicating whether the stem recruited in this interval
MT <- function(n0, nm, t) {(1-((sum(decostand(n0, "pa", na.rm=T), na.rm=T) - sum(nm, na.rm = T))/sum(decostand(n0, "pa", na.rm=T), na.rm=T))^(1/t))*100} #Where n0 are the initial dbh values for each stem (including NAs) and nm is the vector indicating whether the stem died in this interval

##Overall
AnnualDynamics <- data.frame(Recruitment = c(RC(n0 = datAll$DBH2, nr = datAll$recruit2001, t=10), RC(n0 = datAll$DBH3, nr = datAll$recruit2011, t=10), RC(n0 = datAll$DBH2, nr = datAll$recruit1991_2011, t=20)), Mortality = c(MT(n0 = datAll$DBH2, nm = datAll$died2001, t = 10), MT(n0 = datAll$DBH3, nm = datAll$died2011, t = 10), MT(n0 = datAll$DBH2, nm = datAll$died1991_2011, t = 20)))

AnnualDynamics$Dynamism <- (AnnualDynamics$Recruitment + AnnualDynamics$Mortality)/2
rownames(AnnualDynamics) <- c("1991 - 2001", "2001 - 2011", "1991 - 2011")

t(AnnualDynamics)

##By plot
Recruitment2001 <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH2, nr = datAll$recruit2001, t=10)}, INDICES = list(datAll$PLID))
Recruitment2011 <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH3, nr = datAll$recruit2011, t=10)}, INDICES = list(datAll$PLID))
Recruitment1991_2011  <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH2, nr = datAll$recruit1991_2011, t=20)}, INDICES = list(datAll$PLID))

Mortality2001 <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH2, nm = datAll$died2001, t=10)}, INDICES = list(datAll$PLID))
Mortality2011 <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH3, nm = datAll$died2011, t=10)}, INDICES = list(datAll$PLID))
Mortality1991_2011  <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH2, nm = datAll$died1991_2011, t=20)}, INDICES = list(datAll$PLID))

Plot_dynamics <- data.frame(Plot = as.factor(rep(levels(as.factor(datAll$PLID)),3)), Period = as.factor(c(rep("1991-2001", length(Recruitment2001)), rep("2001-2011", length(Recruitment2001)), rep("1991-2011", length(Recruitment2001)))), Recruitment = c(Recruitment2001, Recruitment2011, Recruitment1991_2011), Mortality = c(Mortality2001, Mortality2011, Mortality1991_2011))

Plot_dynamics$Dynamism <- (Plot_dynamics$Recruitment + Plot_dynamics$Mortality)/2

##By species
Recruitment2001 <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH2, nr = datAll$recruit2001, t=10)}, INDICES = list(datAll$Species))
Recruitment2011 <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH3, nr = datAll$recruit2011, t=10)}, INDICES = list(datAll$Species))
Recruitment1991_2011  <- by(datAll, FUN = function(datAll) {RC(n0 = datAll$DBH2, nr = datAll$recruit1991_2011, t=20)}, INDICES = list(datAll$Species))

Mortality2001 <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH2, nm = datAll$died2001, t=10)}, INDICES = list(datAll$Species))
Mortality2011 <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH3, nm = datAll$died2011, t=10)}, INDICES = list(datAll$Species))
Mortality1991_2011  <- by(datAll, FUN = function(datAll) {MT(n0 = datAll$DBH2, nm = datAll$died1991_2011, t=20)}, INDICES = list(datAll$Species))

Species_dynamics <- data.frame(Plot = as.factor(rep(levels(as.factor(datAll$Species)),3)), Period = as.factor(c(rep("1991-2001", length(Recruitment2001)), rep("2001-2011", length(Recruitment2001)), rep("1991-2011", length(Recruitment2001)))), Recruitment = c(Recruitment2001, Recruitment2011, Recruitment1991_2011), Mortality = c(Mortality2001, Mortality2011, Mortality1991_2011))

Species_dynamics$Dynamism <- (Species_dynamics$Recruitment + Species_dynamics$Mortality)/2

###############################################################################
###Calculate Stem density, BA and AGB
###############################################################################

### Basal Area
BA <- function(x) {pi*(x/(2*100))^2} #Where x = diameter (Note the *100 to convert cm to m)

datAll$BA1991 <- BA(datAll$DBH2)
datAll$BA2001 <- BA(datAll$DBH3)
datAll$BA2011 <- BA(datAll$DBH4)

###AGB according to Baker et al. 2004, using wood density from Midgley and Seydack 2006
#constantAGB <- function(x){.933/.58*exp(2.42*log(x)-2)} #.933 = wood density, x=diameter #Function to calculate biomass per plot by year (eq 2.3 nd 2.4 in Baker et al. 2004, using wood density from Midgley and Seydack 2006)
constantAGB <- function(x){exp(.33*log(x) + .933*log(x)^2 - .122*log(x)^3 - .37)}

datAll$BakerAGB1991 <- constantAGB(datAll$DBH2)
datAll$BakerAGB2001 <- constantAGB(datAll$DBH3)
datAll$BakerAGB2011 <- constantAGB(datAll$DBH4)

###AGB with species-specific wood density following the method of Chave et al. 2005 - NEED TO FIX!!!
chaveAGB <- function(x, d) {d * exp(1.78*log(x) + .207*log(x)^2 - .0281*log(x)^3 - .677)} #Where d = species specific wood density and x = diameter

datAll$ChaveAGB1991 <- chaveAGB(x = datAll$DBH2, d = datAll$Wood_Density)
datAll$ChaveAGB2001 <- chaveAGB(x = datAll$DBH3, d = datAll$Wood_Density)
datAll$ChaveAGB2011 <- chaveAGB(x = datAll$DBH4, d = datAll$Wood_Density)

###Summarize and plot results
#Stand level summary
colSums(datAll[, c("BA1991", "BA2001", "BA2011", "BakerAGB1991", "BakerAGB2001", "BakerAGB2011", "ChaveAGB1991", "ChaveAGB2001", "ChaveAGB2011")], na.rm=T)/(0.04*108) #/(0.04*108) converts from /400m^2 to Ha

#Plot level summary
datStemDens <- aggregate(decostand(datAll[, c("DBH2","DBH3","DBH4")], "pa", na.rm=T), by = list(datAll$PLID), FUN = sum, na.rm=T)
datSum <- aggregate(datAll[, c("BA1991", "BA2001", "BA2011", "BakerAGB1991", "BakerAGB2001", "BakerAGB2011", "ChaveAGB1991", "ChaveAGB2001", "ChaveAGB2011")], by = list(datAll$PLID), FUN = "sum", na.rm=T)

datSum <- data.frame(Plot = as.factor(rep(datSum$Group.1, 3)), Year = as.factor(c(rep(1991, nrow(datSum)), rep(2001, nrow(datSum)), rep(2011, nrow(datSum)))), Stem_Density = c(datStemDens$DBH2, datStemDens$DBH3, datStemDens$DBH4), Basal_Area = c(datSum$BA1991, datSum$BA2001, datSum$BA2011), AGB_Constant = c(datSum$BakerAGB1991, datSum$BakerAGB2001, datSum$BakerAGB2011), AGB_Specific = c(datSum$ChaveAGB1991, datSum$ChaveAGB2001, datSum$ChaveAGB2011))

###############################################################################
###Get distribution data and fix accuracy
###############################################################################
acc <- read.csv("/Users/jasper/Documents/Databases/PRECIS and SAPIA_Feb2013/precis.csv", stringsAsFactors=F)#[,1:2]

###Match records
x <- which(acc$TAXNAME %in% tmeans$Species)
###Match records with another name
y <- which(acc$TAXNAME %in% c("Afrocanthium mundianum", "Ilex mitis var. mitis","Ochna arborea var. arborea"))

###Retrieve records
acc<-acc[c(x,y),]

###Fix names to match samp and trait dataframes
acc <- acc[-which(acc$TAXNAME=="Ochna arborea"),]
acc$TAXNAME[which(acc$TAXNAME=="Ochna arborea var. arborea")]<-"Ochna arborea"
acc$TAXNAME[which(acc$TAXNAME=="Ilex mitis var. mitis")]<-"Ilex mitis"
acc$TAXNAME[which(acc$TAXNAME=="Afrocanthium mundianum")]<-"Canthium mundianum"

#sort(unique(acc$TAXNAME)) #to check names

###############################################################################
###Read in climate data and draw out for species localities
###############################################################################
#source("/Users/jasper/Documents/Dell/GIS/WorldClimCurrent/clim2QDS.R")

load("/Users/jasper/Documents/GIS/WorldClimCurrent/Bio/BioSA/bioclimQDS.RData")

###Trim climdat and acc to common QDS cells
climdat <- climdat[which(names(climdat)%in%unique(acc$GRIDREF))] 
acc<-acc[which(acc$GRIDREF%in%names(climdat)),]

###Vector of species names
sp <- tmeans$Species

###Draw out climate data into a list of tables for each species
spclim<-list()
for(i in 1:length(sp))
{
  z<-unique(acc[which(acc$TAXNAME==sp[i]),]$GRIDREF) #get QDS names for spp i
  p<-which(names(climdat)%in%z) #match QDS names in z with positions in climdat
  
  spclim[[i]]<-unlist(climdat[[p[1]]]) #place data from first QDS in output
  for(j in 2:length(p)) #loop tp draw out QDS data 1 by 1 and bind into single table 
  {
    spclim[[i]]<-rbind(spclim[[i]],unlist(climdat[[p[j]]]))
  }
}


###Summarize tables
climmeans <- as.data.frame(t(sapply(spclim,"colMeans", na.rm=TRUE)))
climmeans$Species <- sp

tmeans <- merge(tmeans, climmeans, all.x = T)

###############################################################################
###Clean up and save
###############################################################################

#save(list = c("AnnualDynamics", "StemCondition", "Plot_dynamics", "Species_dynamics", "MortDat", "datAll", "datSum", "tmeans"), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_3Oct2017.Rdata")

save(list = c("AnnualDynamics", "StemCondition", "Plot_dynamics", "Species_dynamics", "MortDat", "MortDatAll", "datAll", "datSum", "tmeans"), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/KnysnaData_19Nov2018.Rdata")
