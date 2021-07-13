#This is a series of functions for dendrochronological analysis
#By Stockton Maxwell, NADEF Co-organizer
#North American Dendroecological Fieldweek
#dendrotilyoudie

##############################################################################

#GETTING STARTED

#Set working directory or use Session menu
setwd("C:/Users/rmaxwell2/Google Drive/Research Projects/NADEF_SM/NADEF 2021 Dendroclimatology/Arrr_Demo_NADEF_2021")#make sure those slashes face the correct way

#install libraries
install.packages("dplR")
install.packages("treeclim")
install.packages("TRADER")
install.packages("DendroSync")
install.packages("dendroTools")
install.packages("burnr")
install.packages("graphics")
install.packages("utils")

#load libraries
library(dplR)
library(treeclim)
library(TRADER)
library(DendroSync)
library(dendroTools)
library(burnr)
library(graphics)
library(utils)


########################################################################################
#Dendrochronology Program Library in R 
#Cofecha type stuff in dplR - this section of the package helps the user crossdate tree ring series

#read in raw ring widths, change file name to run your own data
grow.rwl <- read.tucson(fname = "Yellowstone_PSME_format.txt") #just change fname to run stuff below or use the appropriate "read" function - see dplR help
rwl.stats(grow.rwl) #summary and stats of raw ring width file
seg.plot(grow.rwl) #plot of series time spans
skel.plot(grow.rwl[1]) #skeleton plot of an individual series
colSums(grow.rwl, na.rm = TRUE, dims = 1) # get radii length in mm
crn.common <- common.interval(grow.rwl, type=c("years"), make.plot=TRUE) #find common time interval
sens1(grow.rwl) #Calculate mean sensitivity
rwl.report(grow.rwl)

#crossdating- you can use this but I highly recommend Andy Bunn's xdater app
corr.rwl.seg(rwl = grow.rwl, seg.length = 50, bin.floor = 100, n = NULL, prewhiten = TRUE, pcrit = 0.05, 
             biweight = TRUE, method = c("spearman"), make.plot = TRUE, label.cex = 1, floor.plus1 = FALSE,
             master = NULL) #cofecha essentially
series.rwl.plot(grow.rwl, series = "BTK05A", series.yrs = as.numeric(names(series)), #look at an individual series
                seg.length = 50, bin.floor = 100, n = NULL,
                prewhiten = TRUE, biweight = TRUE, floor.plus1 = FALSE)
interseries.cor(grow.rwl, n = NULL, prewhiten = TRUE, biweight = TRUE, method = "spearman")#calculate interseries correlations for each series


#spaghetti plot of raw ring widths
spag.plot(rwl = grow.rwl, zfac = 1, useRaster = FALSE, res = 300)

#Calculate marker rings from raw ring width series
markers <- pointer(grow.rwl)


###################################################################################
#Arstan stuff in dplR - this section of the package detrends or standardizes series into a site chronology

#interactive detrending - this allows you to explore curve fits for each tree ring series
grow.rwi.int <- i.detrend(rwl = grow.rwl, nyrs = NULL, f = 0.5,pos.slope = FALSE) #allows you to see a variety of fits
spag.plot(rwl = grow.rwi.int, zfac = 1, useRaster = FALSE, res = 300) #again but with the detrended series

#detrend all series at once - after you know which option is best for your data. Just adjust the method.
grow.rwi <- detrend(rwl = grow.rwl, method = c("Spline"), nyrs = NULL, f = 0.5, pos.slope = FALSE) 
rwi.stats(grow.rwi) #stats for entire crn
rwi.stats.running(grow.rwi) #running stats - time periods can be adjusted, see help

#building crn with AR model, this produces a residual crn
grow.crn <- chron(x = grow.rwi, prefix = "BTK", biweight = TRUE, prewhiten = TRUE)
#plot crn
crn.plot(crn = grow.crn, add.spline = TRUE, nyrs = NULL, f = 0.5, crn.line.col='grey50',
         spline.line.col='red', samp.depth.col='grey90', samp.depth.border.col='grey80',
         crn.lwd=1, spline.lwd=2.0, abline.pos=1, abline.col='black', abline.lty=1,abline.lwd=1,
         xlab="Time", ylab="RWI")

#building crn without AR model, this produces a standardized crn
grow.crn <- chron(x = grow.rwi, prefix = "BTK", biweight = TRUE, prewhiten = FALSE)
#plot crn
crn.plot(crn = grow.crn, add.spline = TRUE, nyrs = NULL, f = 0.5, crn.line.col='grey50',
         spline.line.col='red', samp.depth.col='grey90', samp.depth.border.col='grey80',
         crn.lwd=1, spline.lwd=2.0, abline.pos=1, abline.col='black', abline.lty=1,abline.lwd=1,
         xlab="Time", ylab="RWI")

#wavelet transform - this allows you to look at frequencies or temporal patterns in your crn. It's good for paleoclimatology.
Years <- as.numeric(rownames(grow.crn))
rings <- grow.crn[, 1]
tubular <- morlet(y1 = rings, x1 = Years, p2 = 9, dj = 0.1,
                   siglvl = 0.99)
wavelet.plot(tubular, useRaster = NA)

#subset and save crn as .csv for later analysis in DendroTools package
btk_std <- grow.crn[1] #subset only year (already as.numeric) and index columns
write.csv(btk_std, file = "BTK_std.csv")


#####################################################################################
#Treeclim - this package allow for the assessment of growth-climate relationships

#bring in data frame from dplR, run a summary on it
summary(grow.crn)

#bring in climate data, change file name to run your own data
climate <- read.csv("WY_yellowstone_climate.csv", header = TRUE)
names(climate)
summary(climate)
ym <- climate[,1:2] #pull year and month columns
var1 <- climate[3] #pull climate variables, 1 or 2 at a time to avoid N problems
clim <- data.frame(c(ym, var1)) #build climate data frame 
summary(clim)

#Response function analysis in treeclim - modeled after Dendroclim2002. 
#Can take dynamic = "static", "moving", "evolving"
resp <- dcc(chrono = grow.crn, climate = clim, selection = -5:10, 
                method = "correlation", dynamic = "moving", win_size = 35, win_offset = 1, start_last = TRUE,
                timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE) #this is the main function in treeclim
coef <- coef(resp) #model coefficients 
plot(resp) #plot the model coefficients
resp #show model results if you'd like
traceplot(resp, variables = c("PCP.prev.aug", "PCP.curr.jun", "PCP.curr.jul"), facet = FALSE) #shows correlations over time if moving or evolutionary selected
#save the output
write.csv(coef, file = ("pcp_coef.csv")) #output coeff
#write plot to file
tiff("resp_coef.tiff", width = 8, height = 4, units = 'in', res = 300) 
plot.new()
plot(resp)
title(main = "Climate", xlab = "Month")
dev.off()

#evaluate recon skill with split calibration, requires 2 climate variables/months
recon <- dcc(chrono = grow.crn, climate = clim, selection = 6:7, #use a selection with recon variable of interest - modifiers like .mean or .sum can be used to average across months
                method = "response", dynamic = "static", win_size = 35, win_offset = 1, start_last = TRUE,
                timespan = NULL, var_names = NULL, ci = 0.05, boot = "std", sb = FALSE)
recon_coef <- coef(recon)
plot(recon)
recon #show model results
#this does the evaluation
skillz <- skills(object = recon, target = .mean(1:2), model = "ols", calibration = "50%", timespan = NULL)
plot(skillz)
skillz #show model results

#evaluate recon skill with split calibration with single variable - is your model time stable and does it verify?
recon_dlm <- dlm(chrono = grow.crn, climate = clim, selection = 6, timespan = NULL, var_names = NULL,
    param_names = NULL, intercept = TRUE, scale = FALSE)
recon_coef <- coef(recon_dlm)
plot(recon_dlm)
recon_dlm #show model results
#this does the evaluation
skillz <- skills(object = recon_dlm, target = 6, model = "ols", calibration = "50%", timespan = NULL)
plot(skillz)
skillz #show model results

#seasonal correlation - developed from Dave Meko orginally in Matlab
climate <- read.csv("WY_yellowstone_climate.csv", header = TRUE)
names(climate)
summary(climate)
ym <- climate[,1:2] #pull year and month columns
var1 <- climate[3:4] #pull climate variables, 1 or 2 at a time to avoid N problems
clim <- data.frame(c(ym, var1)) #build climate data frame 
summary(clim)
seas <- seascorr(grow.crn, climate, var_names = NULL, timespan = NULL, complete = 9,
         season_lengths = c(1, 3, 6), primary = 1, secondary = 2, ci = 0.05)#this is the main function
plot(seas)
seas


############################################################################################
#TRADER Code for Growth Release Detection

#load libraries needed
library(dplR)
library(TRADER)

#set working directory
setwd("C:/Users/rmaxwell2/Google Drive/Research Projects/NADEF_SM/NADEF 2021 Dendroclimatology/Arrr_Demo_NADEF_2021")

#read in crns, change file name to run your own data
thedata <- read.tucson('Yellowstone_PSME_format.txt') #Use your file in your WD

#Abrams and Nowacki technique - this is just one of many techniques
growthAveragingALL(thedata, releases = NULL, m1 = 15, m2 = 15,buffer = 10, drawing = TRUE, criteria = 0.25, criteria2 = 0.50,prefix = "ga", gfun = mean, length = 5, storedev = jpeg)

#plot raw data
spag.plot(thedata, zfac = 1, useRaster = FALSE, res = 300)
thedata.raw.crn <- chron(thedata, prefix = "BTK", prewhiten=FALSE)
plot(thedata.raw.crn,abline.pos=NULL,ylab='mm',xlab='Year')

#This function calculates the synchronous growth changes (sgc), semi synchronous growth changes (ssgc) and the length of the compared overlap for a given set of tree-ring records.
changes <- sgc(thedata,overlap = 50, prob = TRUE)
mean(changes$sgc_mat, na.rm = TRUE)
mean(changes$ssgc_mat, na.rm = TRUE)


##############################################################################
#Basal Area Increment calculation in dplR

grow.rwl <- read.tucson(fname = "Yellowstone_PSME_format.txt") #just change fname to run stuff below
basal <- bai.out(grow.rwl, diam = NULL)
basal_p <- print(basal)
spag.plot(basal[5], zfac = 1, useRaster = FALSE, res = 300)
write.csv(basal_p, "basal_bai.csv")


###############################################################################
#BURNR
#How about a little fire history graphing - there's some superposed epoch analysis in there too if you fancy
#From Chris Gentry, NADEF Co-organizer

setwd("C:/Users/rmaxwell2/Google Drive/Research Projects/NADEF_SM/NADEF 2021 Dendroclimatology/Arrr_Demo_NADEF_2021")
library(burnr)
library(dplR)
library(ggplot2)
Zion <- read_fhx('Zion.fhx')
Sites <- read.csv('ZionSiteIDs.csv')
facetplot <- plot_demograph(Zion, facet_group = Sites$SiteID, facet_id = Sites$series, plot_legend = TRUE)
print(facetplot)
rugplot <- plot_demograph(Zion, composite_rug = TRUE, plot_legend = TRUE)
compositerug <- rugplot + annotate('rect', xmin = 1721, xmax = 1723, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1734, xmax = 1736, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1748, xmax = 1750, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1777, xmax = 1779, ymin = 0, ymax = 21, alpha = 0.4) + annotate('rect', xmin = 1793, xmax = 1795, ymin = 0, ymax = 21, alpha = 0.4) + scale_x_continuous(limits=c(1450, 2005), breaks = seq(1450,2005,25)) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(compositerug)


#############################################################################
#DENDROTOOLS
#Alright ya'll, let's step up the game
#Let's try some climate response with daily climate data 'cause trees don't know what a month is
#You can get daily climate data for the US from the PRISM Climate Group

setwd("C:/Users/rmaxwell2/Google Drive/Research Projects/NADEF_SM/NADEF 2021 Dendroclimatology/Arrr_Demo_NADEF_2021")
library(dendroTools)
#download daily data from a single point for a single variable on the PRISM website
cdata <- read.table(file = "BTK_daily_pcp_original.csv", skip = 10, header = TRUE, sep = ",") #skips reading the header
head(cdata)
cdata <- data_transform(input = cdata, format = "daily",
  monthly_aggregate_function = "auto", date_format = "ymd")#transform data

#have a look at the daily climate data
glimpse_daily_data(env_data = cdata, tidy_env_data = FALSE, na.color = "white")

#analyze growth vs climate with fixed window width
fixed_width <- daily_response(response = crn_d, env_data = cdata,
                                      method = "cor", fixed_width = 60,
                                      row_names_subset = TRUE, remove_insignificant = TRUE,
                                      alpha = 0.05)
fixed_width$plot_extreme #creates a plot showing best correlated period

#Compare the response across two periods of analysis to assess time stability, or you can leave the subset to all the years
btk_past <- daily_response(response = crn_d, env_data = cdata,
                                   method = "cor", lower_limit = 50, upper_limit = 70,
                                   row_names_subset = TRUE, previous_year = TRUE,
                                   remove_insignificant = TRUE, alpha = 0.05, 
                                   plot_specific_window = 60, subset_years = c(1982, 1998))
btk_present <- daily_response(response = crn_d, env_data = cdata,
                                      method = "cor", lower_limit = 50, upper_limit = 70,
                                      row_names_subset = TRUE, previous_year = TRUE,
                                      remove_insignificant = TRUE, alpha = 0.05, 
                                      plot_specific_window = 60, subset_years = c(1999, 2016))
#plot the results
btk_past$plot_heatmap
btk_present$plot_heatmap
btk_past$plot_specific #choose a specific window length to plot if you set this above
btk_present$plot_specific


##############################################################################
#Principal Components Analysis then climate response on the PCs
# Load data
data(example_proxies_individual)
data(LJ_daily_temperatures)
# Example PCA - just create a data frame with multiple crns
example_PCA <- daily_response(response = example_proxies_individual, 
                              env_data = LJ_daily_temperatures, method = "lm", 
                              lower_limit = 60, upper_limit = 70,
                              row_names_subset = TRUE, remove_insignificant = TRUE,
                              alpha = 0.01, PCA_transformation = TRUE,
                              components_selection = "manual", N_components = 2)
summary(example_PCA$PCA_output)
example_PCA$plot_heatmap

#oh, you want a quick common period reconstruction?
data(data_TRW)
data(KRE_daily_temperatures)
example_reconstruction_lin <- daily_response(response = data_TRW, 
                                             env_data = KRE_daily_temperatures, 
                                             method = "lm", metric = "r.squared", 
                                             lower_limit = 30, upper_limit = 40,
                                             row_names_subset = TRUE, 
                                             temporal_stability_check = "progressive",
                                             cross_validation_type = "randomized", k = 3)
example_reconstruction_lin$plot_extreme
example_reconstruction_lin$temporal_stability
example_reconstruction_lin$cross_validation
example_reconstruction_lin$transfer_function
linear_model <- lm(Optimized_return ~ TRW, data = example_reconstruction_lin$optimized_return)
reconstruction <- data.frame(predictions = predict(linear_model, newdata = data_TRW))
linear_model <- lm(Optimized_return ~ TRW, data = example_reconstruction_lin$optimized_return)
reconstruction <- data.frame(predictions = predict(linear_model, newdata = data_TRW))
plot(row.names(data_TRW), reconstruction$predictions, type = "l", xlab = "Year", ylab = "Mean temperature May 15 - Jun 27 [ºC]")


############################################################################################
#MONTHLY ANALYSIS
#load data
flow_d <- read.csv("BighornXavier_r.csv", header = TRUE)
row.names(flow_d) <- as.numeric(flow_d$year)
flow_d <- flow_d[,2:13]#subset months, not year column
#run analysis with split period
flow_past <- monthly_response(response = crn_d, env_data = flow_d,
                                     method = "cor", row_names_subset = TRUE, previous_year = TRUE,
                                     remove_insignificant = TRUE, alpha = 0.05,
                                     subset_years = c(1936, 1976), aggregate_function = 'mean')

flow_present <- monthly_response(response = crn_d, env_data = flow_d,
                                        method = "cor", row_names_subset = TRUE, alpha = 0.05,
                                        previous_year = TRUE, remove_insignificant = TRUE,
                                        subset_years = c(1977, 2016), aggregate_function = 'mean')

flow_past$plot_heatmap
flow_present$plot_heatmap
flow_past$plot_extreme
flow_present$plot_extreme


########### NOT READY YET #######################
#Monthly analysis using a PCA of chronologies
#bring in crns
###############Still need to get crns in column format
flow_PCA <- monthly_response(response = example_proxies_individual,
                                env_data = flow_d, method = "lm",
                                row_names_subset = TRUE, remove_insignificant = TRUE,
                                alpha = 0.01, PCA_transformation = TRUE, previous_year = TRUE,
                                components_selection = "manual", N_components = 2)

summary(flow_PCA$PCA_output)
flow_PCA$plot_heatmap
flow_PCA$plot_extreme


#################################################################################
#DendroSync - Provides functions for the calculation and plotting of synchrony in 
#tree growth from tree-ring width chronologies (TRW index).
library(DendroSync)
data(conifersIP) #note the format of the data if you choose to use this analysis
head(conifersIP)

## Calculate synchrony for null.model (broad evaluation, mBE) and homoscedastic variant
# of unstructured model (or full, mUN) for conifersIP data, 
# and heteroscedastic variant for 1970-1999 period.

##Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN) 
#using taxonomic grouping criteria (i.e. Species)
ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
                       data = conifersIP, homoscedastic = TRUE)

summary(ModHm)# Class and length of list elements

#Synchrony for mBE and mUN models
sync(ModHm, modname = "mBE")
sync(ModHm, modname = "mUN")

##Chop the data from 1970 to 1999.
conif.30 <- conifersIP[conifersIP$Year>1969 & conifersIP$Year<2000,]
summary(conif.30$Year)

#Fit the heteroscedastic set of variance covariance mixed models (mBE, mHeNE, mHeCS, mHeUN)
# using taxonomic grouping criteria (ie. Species)
ModHt30 <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
                         data = conif.30, homoscedastic = FALSE)
sync(ModHt30, modname = "mBE")
sync(ModHt30, modname = "mHeUN")
