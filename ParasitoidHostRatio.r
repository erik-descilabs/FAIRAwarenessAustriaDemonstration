####################STATISTICAL ANALYSIS################################

# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4d - Effect of meadows (500 m radius) on parasitoid-host ratio####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 4a (based on model predictions, not the raw data)

# https://datadryad.org/stash/dataset/doi:10.5061/dryad.484tt
Landscape.affects.pest.and.crop.yield_2023<-read.csv("https://datadryad.org/stash/downloads/file_stream/28280", na.strings=c("NA", "", "N/A"))

#shorten dataset name
LandscapeData <- Landscape.affects.pest.and.crop.yield_2023
#remove the long name dataset, so we accidentally don't end up 
#working with the wrong data
rm(Landscape.affects.pest.and.crop.yield_2023)
#check dataset: should yield 43 obs. of 20 variables
head(LandscapeData)
#Dataset looks fine! Let's check the summary
summary(LandscapeData)
#year is read as a numeric variable. same with Plot_ID
#converting Year and Plot_ID from numeric to character
LandscapeData$Year <- as.character(LandscapeData$Year)
LandscapeData$Plot_ID <- as.character(LandscapeData$Plot_ID)
#check summary again
summary(LandscapeData)
#looks good now!

############ LOAD PACKAGES ##################### 
library(nlme)#nlme_3.1-152  
library(MuMIn)#MuMIn_1.43.17  
library (ggplot2)#ggplot2_3.3.6  
library (visreg)#visreg_2.7.0 
library(effects)#effects_4.2-0  
library(orkg)
library (ggpubr)#vesion 0.4.0
library(dplyr)#version 1.08

######## REMOVING MISSING VALUES ########
#Removing the rows with missing data from each variable individually (this way we don't delete more data than we need to)

names(LandscapeData)
FleaBeetlesIncidence <- LandscapeData[!is.na(LandscapeData$FleaBeetles_incidence), ]
FleaBeetlesAbundance <- LandscapeData
LepidopteranIncidence <- LandscapeData[!is.na(LandscapeData$Lepidoptera_incidence), ]
LepidopteranAbundance <-LandscapeData[!is.na(LandscapeData$Lepidoptera_abundance), ]
AphidsIncidence <- LandscapeData[!is.na(LandscapeData$Aphid_incidence), ]
ParasitoidHostRatio <-LandscapeData[!is.na(LandscapeData$ParasitoidHostRatio), ]
PlantDamageIndex <-LandscapeData[!is.na(LandscapeData$Plant_damage), ]


#######CREATING PLOTS BASED ON LME MODEL PREDICTIONS ########
#making model predictions and calculating confidence intervals


####Parasitoid Host-ratio - figure 4d#########
fitlme.Ph <- lme(ParasitoidHostRatio~  mead_500, data = ParasitoidHostRatio, random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme.Ph$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'


#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 500 meter radius of the field plot) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme.Ph,type='marginal'))
anovaOutput


#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme.Ph)$tTable, check.names=FALSE)
sum1 

newdat.lme.Ph  = data.frame(mead_500 = ParasitoidHostRatio$mead_500,
                            ParasitoidHostRatio=ParasitoidHostRatio$ParasitoidHostRatio
)
head(newdat.lme.Ph )
newdat.lme.Ph$predlme.Ph = predict(fitlme.Ph, newdata = newdat.lme.Ph, level = 0)
#ggplot(ParasitoidHostRatio, aes(x = mead_500, y = ParasitoidHostRatio) ) +
# geom_rug(sides = "b", size = 1) +
#geom_line(data = newdat.lme.Ph, aes(y = predlme.Ph), size = 1)
des.Ph = model.matrix(formula(fitlme.Ph)[-2], newdat.lme.Ph)
predvar.Ph = diag( des.Ph %*% vcov(fitlme.Ph) %*% t(des.Ph) )
newdat.lme.Ph$lower = with(newdat.lme.Ph, predlme.Ph - 2*sqrt(predvar.Ph) )
newdat.lme.Ph$upper = with(newdat.lme.Ph, predlme.Ph + 2*sqrt(predvar.Ph) )

ggplot(ParasitoidHostRatio, aes(x = mead_500, y = ParasitoidHostRatio) ) +
  geom_point()+
  geom_rug(sides = "b", size = 1) +
  geom_ribbon(data = newdat.lme.Ph, aes(y = NULL, ymin = lower, ymax = upper,
                                        color = NULL),
              alpha = .15) +
  geom_line(data = newdat.lme.Ph, aes(y = predlme.Ph), size = .75)+
  theme_classic()+xlab('Proportion of meadows at 500-m')+ylab('Parasitoid Host-ratio')




############Parasitoid-host ratios- FIGURE 4d #############################

#Figure 4d--Keeping the default settings in ggplot
ParasitoidHostRatio.default<- ggplot(LandscapeData, aes(x=mead_500, y=ParasitoidHostRatio)) +   geom_point(shape=1) + geom_smooth (aes (x=mead_500, y=ParasitoidHostRatio), method=lm, se=TRUE, fullrange=TRUE) + theme(panel.background = element_rect(fill='white', colour='black'))+xlab('Proportion of meadows at 500-m')+ylab('ParasitoidHostRatio')
ParasitoidHostRatio.default#we can ignore the warning. It is just telling us that there is missing data
#Figure4c--Reproducing the plot theme selected in the paper (the way the plot looks like in the paper)
ParasitoidHostRatio.papertheme<- ggplot(LandscapeData, aes(x=mead_500, y=ParasitoidHostRatio)) +  geom_point(size = 2) + geom_smooth (aes (x=mead_500, y=ParasitoidHostRatio), method=lm, se=TRUE, fullrange=TRUE,color = "black") +
  scale_color_grey (start=0.7, end=0.2)+theme(panel.background = element_rect(fill='white', colour='black'))+
  xlab('Proportion of meadows at 500-m')+ ylab('ParasitoidHostRatio')
ParasitoidHostRatio.papertheme
#add equation to the plot---

ParasitoidHostRatio.papertheme+ stat_regline_equation()

#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models
ParasitoidHostRatio.model <-lm (ParasitoidHostRatio~mead_500, data=LandscapeData,na.action=na.omit )
sumLR <- data.frame(summary(ParasitoidHostRatio.model)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4d based on the model predictions.
PredictedValuesParasitoidHostRatio <- newdat.lme.Ph
PredictedValuesParasitoidHostRatio <- subset(PredictedValuesParasitoidHostRatio, select = -c(ParasitoidHostRatio))
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled ParasitoidHostRatio as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_500: same as described above
#predlme.Ph: predicted (based on model) parasitoid-host ratio of a field
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI



###End of script

#Input data set
inputDF <- ParasitoidHostRatio[, c("Year", "Farm_ID", "Plot_ID", "ParasitoidHostRatio", "mead_500")]
