#This is a series of function for dendrochronology using dplR, treeclim, and TRADER.
#By Stockton Maxwell
#NADEF2017
#betweentwotrees

#load libraries
library(dplR)
library(treeclim)
library(graphics)
library(utils)
library(TRADER)
library(burnr)
library(ggplot2)

###############################################
#Cofecha type stuff in dplR

#read in raw ring widths, change file name to run your own data
grow.rwl <- read.rwl(fname = "wa091.rwl", format = "auto") #just change fname to run stuff below
rwl.stats(grow.rwl) #summary and stats of raw ring width file
seg.plot(grow.rwl) #plot of series time spans
skel.plot(grow.rwl[1]) #skeleton plot of an individual series

colSums(grow.rwl, na.rm = TRUE, dims = 1) # get radii length in mm
wa091.common <- common.interval(grow.rwl, type=c("years"), make.plot=TRUE) #find common time interval

#Calculate mean sensitivity
sens1(grow.rwl)

#crossdating
corr.rwl.seg(rwl = grow.rwl, seg.length = 50, bin.floor = 100, n = NULL, prewhiten = TRUE, pcrit = 0.05, 
             biweight = TRUE, method = c("spearman"), make.plot = TRUE, label.cex = 1, floor.plus1 = FALSE,
             master = NULL) #cofecha essentially
series.rwl.plot(grow.rwl, series = "TCT202A", series.yrs = as.numeric(names(series)), #look at an individual series
                seg.length = 50, bin.floor = 100, n = NULL,
                prewhiten = TRUE, biweight = TRUE, floor.plus1 = FALSE)

#spaghetti plot of raw ring widths
spag.plot(rwl = grow.rwl, zfac = 1, useRaster = FALSE, res = 300)

#Calculate marker rings from raw ring width series
markers <- pointer(grow.rwl)

###########
#Arstan stuff in dplR

#power transformation of raw ring widths to stablize variance
#grow.pow <- powt(grow.rwl) #optional, skip if problems

#interactive detrending
grow.rwi.int <- i.detrend(rwl = grow.rwl, nyrs = NULL, f = 0.5,pos.slope = FALSE) #allows you to see a variety of fits
spag.plot(rwl = grow.rwi.int, zfac = 1, useRaster = FALSE, res = 300)

#detrend all series at once
grow.rwi <- detrend(rwl = grow.rwl, method = c("Spline"), nyrs = NULL, f = 0.5, pos.slope = FALSE) 
rwi.stats(grow.rwi) #stats for entire crn
rwi.stats.running(grow.rwi) #running stats - time periods can be adjusted, see help

#building crn without AR model, this produces a residual crn
grow.crn <- chron(x = grow.rwi, prefix = "TAH", biweight = TRUE, prewhiten = TRUE)

#building crn without AR model, this produces a standardized crn
grow.crn <- chron(x = grow.rwi, prefix = "TAH", biweight = TRUE, prewhiten = FALSE)

#plot crn
crn.plot(crn = grow.crn, add.spline = TRUE, nyrs = NULL, f = 0.5, crn.line.col='grey50',
         spline.line.col='red', samp.depth.col='grey90', samp.depth.border.col='grey80',
         crn.lwd=1, spline.lwd=2.0, abline.pos=1, abline.col='black', abline.lty=1,abline.lwd=1,
         xlab="Time", ylab="RWI")

#wavelet transform
Years <- as.numeric(rownames(grow.crn))
PANstd <- grow.crn[, 1]
out.wave <- morlet(y1 = PANstd, x1 = Years, p2 = 9, dj = 0.1,
                   siglvl = 0.99)
wavelet.plot(out.wave, useRaster = NA)


##########################################################################
###Bones for Treeclim

#bring in data frame from dplR, run a summary on it
summary(grow.crn)

#bring in climate data, change file name to run your own data
climate <- read.csv("WA_div5_climate_cascades_west.csv", header = TRUE)
names(climate)
summary(climate)
ym <- climate[,1:2] #pull year and month columns
var1 <- climate[3] #pull climate variables, 1 or 2 at a time to avoid N problems
clim <- data.frame(c(ym, var1)) #build climate data frame 
summary(clim)

#Response function analysis in treeclim - modeled after Dendroclim2002. 
#Now can take dynamic = "static", "moving", "evolving"
resp <- dcc(chrono = grow.crn, climate = clim, selection = -5:10, 
            method = "response", dynamic = "evolving", win_size = 35, win_offset = 1, start_last = TRUE,
            timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE)
coef <- coef(resp)
plot(resp)
resp #show model results
traceplot(resp, variables = NULL, facet = FALSE)
#save the output
write.csv(coef, file = ("pcp_coef.csv")) #output coeff
tiff("resp_coef.tiff", width = 8, height = 4, units = 'in', res = 300) #write plot to file
plot.new()
plot(resp)
title(main = "Climate", xlab = "Month")
dev.off()

#evaluate recon skill with split calibration, requires 2 climate variables
recon <- dcc(chrono = grow.crn, climate = clim, selection = 1:2, #use a selection with recon variable of interest
             method = "response", dynamic = "static", win_size = 35, win_offset = 1, start_last = TRUE,
             timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE)
recon_coef <- coef(recon)
plot(recon)
recon #show model results
#this does the evaluation
skillz <- skills(object = recon, target = .mean(1:2), model = "ols", calibration = "50%", timespan = NULL)
plot(skillz)
skillz #show model results

#evaluate recon skill with split calibration with single variable
recon_dlm <- dlm(chrono = grow.crn, climate = clim, selection = 2, timespan = NULL, var_names = NULL,
                 param_names = NULL, intercept = TRUE, scale = FALSE)
recon_coef <- coef(recon_dlm)
plot(recon_dlm)

recon_dlm #show model results
#this does the evaluation
skillz <- skills(object = recon_dlm, target = 2, model = "ols", calibration = "50%", timespan = NULL)
plot(skillz)
skillz #show model results

#seasonal correlation - try it
climate <- read.csv("WA_div5_climate_cascades_west.csv", header = TRUE)
names(climate)
summary(climate)
ym <- climate[,1:2] #pull year and month columns
var1 <- climate[3:4] #pull climate variables, 1 or 2 at a time to avoid N problems
clim <- data.frame(c(ym, var1)) #build climate data frame 
summary(clim)
seas <- seascorr(grow.crn, climate, var_names = NULL, timespan = NULL, complete = 9,
                 season_lengths = c(1, 3, 6), primary = 1, secondary = 2, ci = 0.05)
plot(seas)
seas

############################################################################################
#TRADER Code for Growth Release Detection

#read in crns, change file name to run your own data
thedata <- read.rwl('wa091.rwl') #Use your file in your WD

#Abrams and Nowacki technique creates a large number of files in your wd
growthAveragingALL(thedata, releases = NULL, m1 = 10, m2 = 10,buffer = 10, drawing = TRUE, criteria = -0.25, criteria2 = -0.50,prefix = "ga", gfun = mean, length = 5, storedev = jpeg)

#plot raw data
spag.plot(thedata, zfac = 1, useRaster = FALSE, res = 300)
thedata.raw.crn <- chron(thedata, prefix = "CAM", prewhiten=FALSE)
plot(thedata.raw.crn,abline.pos=NULL,ylab='mm',xlab='Year')

##############################################################################
#Basal Area Increment calculation in dplR

grow.rwl <- read.rwl(fname = "wa091.rwl", format = "auto") #just change fname to run stuff below
basal <- bai.out(grow.rwl, diam = NULL)
basal_p <- print(basal)
spag.plot(basal[1], zfac = 1, useRaster = FALSE, res = 300)
write.csv(basal_p, "basal_bai.csv")

###############################################################################
#From Chris Gentry, NADEF 2017

Zion <- read_fhx('Zion.fhx')
Sites <- read.csv('ZionSiteIDs.csv')
facetplot <- plot_demograph(Zion, facet_group = Sites$SiteID, facet_id = Sites$series, plot_legend = TRUE)
print(facetplot)
rugplot <- plot_demograph(Zion, composite_rug = TRUE, plot_legend = TRUE)
compositerug <- rugplot + 
  annotate('rect', xmin = 1721, xmax = 1723, ymin = 0, ymax = 21, alpha = 0.4) + 
  annotate('rect', xmin = 1734, xmax = 1736, ymin = 0, ymax = 21, alpha = 0.4) + 
  annotate('rect', xmin = 1748, xmax = 1750, ymin = 0, ymax = 21, alpha = 0.4) + 
  annotate('rect', xmin = 1777, xmax = 1779, ymin = 0, ymax = 21, alpha = 0.4) + 
  annotate('rect', xmin = 1793, xmax = 1795, ymin = 0, ymax = 21, alpha = 0.4) + 
  scale_x_continuous(limits=c(1450, 2005), breaks = seq(1450,2005,25)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.2), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(compositerug)
