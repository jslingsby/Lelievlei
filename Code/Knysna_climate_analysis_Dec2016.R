##############################################################################
######## Script to analyse CRU data, monthly summaries from Adam Wilson's interpolated daily weather data for the CFR, and rainfall record for Diepwalle
##############################################################################
######## Compiled by Jasper Slingsby 2016
######## Last edited: 30 Nov 2016
######## Data: weather_fine.nc - summary of Adam's data with 6 vars: 
######## tmax_mean, tmax_sd, tmin_mean, tmin_sd, ppt_mean, ppt_sd for the period 1990 to 2009
##############################################################################
######## 
###Steps:
###1) setup CDO commands for system call to MacOS
###2) read in and process outputs
##############################################################################

library(raster); library(gdata); library(SPEI); library(dplyr); library(ggplot2); library(GGally);library(MCMCglmm)

###Make speiplot function
speiplot <- function (x, tit, ncl, pcl, ...) {
  ser <- ts(as.matrix(x$fitted[-c(1:x$scale), ]), end = end(x$fitted), 
            frequency = frequency(x$fitted))
  ser[is.nan(ser - ser)] <- 0
  se <- ifelse(ser == 0, ser, NA)
  tit <- tit
  if (start(ser)[2] == 1) {
    ns <- c(start(ser)[1] - 1, 12)
  } else {
    ns <- c(start(ser)[1], start(ser)[2] - 1)
  }
  if (end(ser)[2] == 12) {
    ne <- c(end(ser)[1] + 1, 1)
  } else {
    ne <- c(end(ser)[1], end(ser)[2] + 1)
  }
  n <- ncol(ser)
  if (is.null(n)) 
    n <- 1
  par(mar = c(4, 4, 2, 1) + 0.1)
  if (n > 1 & n < 5) 
    par(mfrow = c(n, 1))
  if (n > 1 & n >= 5) 
    par(mfrow = c({
      n + 1
    }%/%2, 2))
  for (i in 1:n) {
    datt <- ts(c(0, ser[, i], 0), frequency = frequency(ser), 
               start = ns, end = ne)
    datt.pos <- ifelse(datt > 0, datt, 0)
    datt.neg <- ifelse(datt <= 0, datt, 0)
    plot(datt, type = "n", main = tit, ...)
    if (!is.null(x$ref.period)) {
      k <- ts(5, start = x$ref.period[1, ], end = x$ref.period[2, 
                                                               ], frequency = 12)
      k[1] <- k[length(k)] <- -5
      polygon(k, col = "light grey", border = NA, density = 20)
      abline(v = x$ref.period[1, 1] + (x$ref.period[1, 
                                                    2] - 1)/12, col = "grey")
      abline(v = x$ref.period[2, 1] + (x$ref.period[2, 
                                                    2] - 1)/12, col = "grey")
    }
    #grid(col = "black")
    polygon(datt.pos, col = pcl, border = NA)
    polygon(datt.neg, col = ncl, border = NA)
    lines(datt, col = "dark grey")
    abline(h = 0)
    points(se, pch = 21, col = "white", bg = "black")
  }
}

##############################################################################
###1) setup CDO command
##############################################################################

#cdo monmean /Users/jasper/Documents/PostDoc/Connecticut/Adam/weatherdat/weather_fine.nc /Users/jasper/Documents/PostDoc/Connecticut/Adam/weatherdat/monmean_17June2015.nc


##############################################################################
###2) Read in and process Wilson data
##############################################################################

tmaxmean=stack("/Users/jasper/Documents/PostDoc/Connecticut/Adam/weatherdat/monmean_17June2015.nc", varname="tmax_mean")
tminmean=stack("/Users/jasper/Documents/PostDoc/Connecticut/Adam/weatherdat/monmean_17June2015.nc", varname="tmin_mean")
pptmean=stack("/Users/jasper/Documents/PostDoc/Connecticut/Adam/weatherdat/monmean_17June2015.nc", varname="ppt_mean")

###Get site data and limit to populations we want displayed
sites = data.frame(x = -33.922, y = 23.043)
coordinates(sites) = ~ y + x
proj4string(sites) <- CRS(proj4string(tmaxmean))

###Extract data to sites and bind site and climate data
tmaxW=extract(tmaxmean,sites)
tminW=extract(tminmean,sites)
tmeanW= data.frame(Year = as.numeric(substr(colnames(tminW),2,5)), Month = month.name[as.numeric(substr(colnames(tminW),7,8))], MeanTW = t((tmaxW+tminW)/2)) # Date = gsub("\\.", "-", substr(colnames(tminW),2,11)), 
pptW=extract(pptmean,sites)
pptW = data.frame(Year = as.numeric(substr(colnames(pptW),2,5)), Month = month.name[as.numeric(substr(colnames(pptW),7,8))], PptW = t(pptW*as.numeric(substr(colnames(pptW), 10, 11)))) # Date = gsub("\\.", "-", substr(colnames(pptW),2,11)), 

###Calculate SPEI with 1,3,6 and 12 month base periods
#PET <- thornthwaite(tmeanW$MeanTW,coordinates(sites)[,2])
#spei12W <- spei(pptW$PptW-PET,12)

#wdat <- data.frame(Date = gsub("\\.", "-", substr(colnames(tmaxW),2,11)), Year = as.numeric(substr(colnames(tmaxW),2,5)), Month = month.name[as.numeric(substr(colnames(tmaxW),7,8))], SPEIW = spei12W$fitted, stringsAsFactors = F, MeanTW = tmeanW, PptW = pptW)
#wdat$Month <- factor(wdat$Month, levels = month.name[1:12])

#wldat <- data.frame(Date = gsub("\\.", "-", substr(colnames(tmaxW),2,11)), SPEI = spei12W$fitted, stringsAsFactors = F, MeanTW = tmeanW, PptW = pptW)

PETW <- thornthwaite(tmeanW$MeanTW,coordinates(sites)[,2])
spei12W <- spei(pptW$PptW-PETW,12)

##############################################################################
###3) Read in CRU data - http://www.cru.uea.ac.uk/data accessed 30 Nov 2016
##############################################################################

tmeanC <- read.table("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/weather_data/CRU_Knysna_1901_2015_Temp.txt", header = F, skip = 5)
colnames(tmeanC) <- c("Year", "Month", "MeanTC")
#tmeanC$Date = month.name[tmeanC$Month]
tmeanC$Month = month.name[tmeanC$Month]

pptC <- read.table("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/weather_data/CRU_Knysna_1901_2015_Precip.txt", header = F, skip = 5)
colnames(pptC) <- c("Year", "Month", "PptC")
pptC$Month = month.name[pptC$Month]

PETC <- thornthwaite(tmeanC$MeanTC,coordinates(sites)[,2])
spei12C <- spei(pptC$PptC-PETC,12)

##############################################################################
###5) Get Diepwalle data
##############################################################################

pptD <- read.csv("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/weather_data/Diepwalle.csv", header = T, stringsAsFactors = F, row.names=1)
temp <- data.frame(expand.grid(dimnames(pptD))[1:2], as.vector(as.matrix(pptD)))
pptD <- data.frame(Year = temp[, 1], Month = month.name[as.numeric(temp[, 2])], PptD = as.numeric(as.character(temp[, 3])))

##############################################################################
###6) Compare datasets
##############################################################################

###Merge T datasets
tmean <- merge(tmeanC, tmeanW, all.x=T)

###Merge Ppt datasets
ppt <- merge(pptC, pptD, all.x=T)
ppt <- merge(ppt, pptW, all.x=T)

pdf("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_precipTS.pdf", width = 9, height = 7)
par(mfrow = c(3,1), mar = c(2,4,1,0) + 0.1)
plot(ppt[,3], xaxt = "n", bty = "n", ylab = "", xlab = "", type = "l", ylim = c(0,600), cex.axis = 1.25)
abline(mean(ppt[,3], na.rm = T), 0, col="grey")
abline(quantile(ppt[,3], probs = .95, na.rm = T), 0, col="grey", lty = "dashed")
abline(quantile(ppt[,3], probs = .05, na.rm = T), 0, col="grey", lty = "dashed")
abline(v=c(1092, 1212, 1332), col="grey")
legend("topleft", "CRU", bty = "n", cex = 1.25)
plot(ppt[,4], xaxt = "n", bty = "n", ylab = "Precipitation", xlab = "", type = "l", ylim = c(0,600), cex.lab=1.45, cex.axis = 1.25)
abline(mean(ppt[,4], na.rm = T), 0, col="grey")
abline(quantile(ppt[,4], probs = .95, na.rm = T), 0, col="grey", lty = "dashed")
abline(quantile(ppt[,4], probs = .05, na.rm = T), 0, col="grey", lty = "dashed")
abline(v=c(1092, 1212, 1332), col="grey")
legend("topleft", "Diepwalle rain gauge", bty = "n", cex = 1.25)
plot(ppt[,5], xaxt = "n", bty = "n", ylab = "", xlab = "", type = "l", ylim = c(0,600), cex.axis = 1.25)
abline(mean(ppt[,5], na.rm = T), 0, col="grey")
abline(quantile(ppt[,5], probs = .95, na.rm = T), 0, col="grey", lty = "dashed")
abline(quantile(ppt[,5], probs = .05, na.rm = T), 0, col="grey", lty = "dashed")
abline(v=c(1092, 1212, 1332), col="grey")
legend("topleft", "Wilson et al. 2014", bty = "n", cex = 1.25)
axis(1, at=seq(0,1380,120), labels=seq(1900,2010,10), cex.lab=1.5, cex.axis = 1.25)
dev.off()

#par(mfrow = c(1,1), mar = mardef)

pdf("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_precip_pairs.pdf", width = 7, height = 7)
pairs(ppt[,3:5], labels = c("CRU", "Diepvalle", "Wilson et al. 2014"))
dev.off()

capture.output(summary(lm(ppt[,4] ~ ppt[,3])), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_precip_pairs_DC.txt")
capture.output(summary(lm(ppt[,5] ~ ppt[,3])), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_precip_pairs_WC.txt")
capture.output(summary(lm(ppt[,5] ~ ppt[,4])), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_precip_pairs_WD.txt")

###Plot temp figures
pdf("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_tempTS.pdf", width = 9, height = 9)
par(mfrow = c(2,1), mar = c(2,4,1,0) + 0.1)
plot(tmean[,3], xaxt = "n", bty = "n", ylab = "", xlab = "", type = "l", ylim = c(8,22))
abline(mean(tmean[,3], na.rm = T), 0, col="grey")
abline(quantile(tmean[,3], probs = .95, na.rm = T), 0, col="grey", lty = "dashed")
abline(quantile(tmean[,3], probs = .05, na.rm = T), 0, col="grey", lty = "dashed")
abline(v=c(1092, 1212, 1332), col="grey")
legend("topleft", "CRU", bty = "n")
plot(tmean[,4], xaxt = "n", bty = "n", ylab = "Mean Temperature", xlab = "Year", type = "l", ylim = c(8,22))
abline(mean(tmean[,4], na.rm = T), 0, col="grey")
abline(quantile(tmean[,4], probs = .95, na.rm = T), 0, col="grey", lty = "dashed")
abline(quantile(tmean[,4], probs = .05, na.rm = T), 0, col="grey", lty = "dashed")
abline(v=c(1092, 1212, 1332), col="grey")
legend("topleft", "Wilson et al. 2014", bty = "n", cex = 1.25)
axis(1, at=seq(0,1380,120), labels=seq(1900,2010,10))
dev.off()

#par(mfrow = c(1,1), mar = mardef)

pdf("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_temp_pairs.pdf", width = 5, height = 5)
plot(tmean[,4] ~ tmean[,3], xlab = "CRU", ylab = "Wilson et al. 2014")
#abline(1,1, col = "grey")
#pairs(ppt[,3:5], labels = c("CRU", "Diepvalle", "Wilson et al. 2014"))
dev.off()

capture.output(summary(lm(tmean[,4] ~ tmean[,3])), file = "/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SM_temp_pairs.txt")

###Plot SPEI figure
pdf("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Paper/Figures/SPEI.pdf", width = 9, height = 7)
par(mfrow = c(2,1), mar = c(2,4,1,0) + 0.1)
speiplot(x = spei12C, ylab="SPEI", tit="", ncl = "#D6604D", pcl = "#4393C3", xaxt = "n", xlab = "", bty = "n")
axis(1, at=seq(0,115,5), labels=seq(1900,2015,5))
abline(v=c(91,101,111), lty = "dashed", lwd=1.5) #
legend("topleft", "CRU", bty = "n")
speiplot(x = spei12W, ylab="SPEI", tit="", ncl = "#D6604D", pcl = "#4393C3", xaxt = "n", xlab = "Year", bty = "n", xlim = c(0,22))
axis(1, at=seq(1,21,5), labels=c("", seq(1995,2010,5)))
abline(v=c(2,12,22), lty = "dashed", lwd=1.5)
legend("topleft", "Wilson et al. 2014", bty = "n")
dev.off()

##############################################################################
###7) Explore NDVI-climate relationship
##############################################################################

# #Get NDVI
# ndvinames <- read.delim("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/NDVI/LT5_L1T_TOA.txt", sep = ":", header = F)
# ndvi <- stack("/Users/jasper/Dropbox/SAEON/Projects/Knysna forest/Data/NDVI/20161212_v1_LT5_L1T_TOA_Lilyvlei_daily__1984-06-20-2011-02-07.tif")
# #NAvalue(ndvi) <- 0
# offs(ndvi) <- -2
# gain(ndvi) <- .001
# names(ndvi) <- as.character(ndvinames[,2])
# ndvi <- setZ(ndvi,as.Date(as.character(ndvinames[,2])))
# 
# ndvimeans <- cellStats(ndvi, stat='mean')
# ndvimins <- cellStats(ndvi, stat='min')
# ndvimaxs <- cellStats(ndvi, stat='max')
# 
# ndvit <- data.frame(Date = as.Date(as.character(ndvinames[,2])), MeanNDVI = ndvimeans, MinNDVI = ndvimins, MaxNDVI = ndvimaxs)
# 
# ndvit <- ndvit[-which(ndvit$MinNDVI<0),]
# 
# ndvit$Month <- months(ndvit$Date)
# ndvit$Year <- substr(ndvit$Date,1,4)
# 
# ###Merge all
# dat <- merge(ndvit, ppt, all.x = T)
# dat <- merge(dat, tmean, all.x = T)
# 
# #Add SPEI
# speiW <- data.frame(Year = tmeanW$Year, Month = tmeanW$Month, SPEIW = as.vector(spei12W$fitted))
# speiC <- data.frame(Year = tmeanC$Year, Month = tmeanC$Month, SPEIC = as.vector(spei12C$fitted))
# 
# dat <- merge(dat, speiC, all.x=T)
# dat <- merge(dat, speiW, all.x=T)
# dat$Month <- factor(dat$Month, levels = month.name)
# dat$YearMonth <- paste(dat$Year, dat$Month, sep="")
# 
# 
# boxplot(SPEIW ~ Month, data = dat)
# boxplot(SPEIC ~ Month, data = dat)
# boxplot(MeanNDVI ~ Month, data = dat)
# boxplot(MeanNDVI ~ Month + Year, data = dat)
# 
# 
# density <- diwish(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.3,.3,1),2,2))
# draw <- riwish(3, matrix(c(1,.3,.3,1),2,2))
# 
# datn <- na.omit(dat)
# 
# fit1 <- MCMCglmm(MeanNDVI ~ SPEIW + PptW + MeanTW, data = datn, burnin=1000, nitt=10000, thin=100)
# 
# summary(fit1)
# 
# #fit1 <- MCMChregress(MeanNDVI ~ SPEIC + PptC + MeanTC, random = ~Month, group = "YearMonth", data = dat, burnin=1000, mcmc=10000, thin=100, r=12, R=diag(c(1,rep(0.1,11))))#, R=diag(c(1,0.1,0.1)))
# 
# density(fit1$mcmc[,2])
# 
# ggpairs(dat[,c(4,7:13)])
# ggpairs(dat[which(dat$Month == "July"),c(4,7:13)])

