####################STATISTICAL ANALYSIS################################
# Contrasting effects of landscape composition on crop yield mediated by specialist herbivores
# https://doi.org/10.1002/eap.1695
# Authors: Ricardo Perez-Alvarez, Brian Nault, and Katja Poveda
# Analysis conducted with R version 4.0.4 (2021-02-15)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)

########### LOAD DATA  ##############
#####This script is for Fig. 4a - Effect of meadows (250 m radius) on aphid incidence####
##There are three key sets of data to extract:
#1 - Output from ANOVA (Type III sum of squares) 
#2 - Output of summary function on lme (fixed effects)
#3 - Predicted model values to recreate figure 4a (based on model predictions, not the raw data)

# https://datadryad.org/stash/dataset/doi:10.5061/dryad.484tt
Landscape.affects.pest.and.crop.yield <- read.csv("https://datadryad.org/stash/downloads/file_stream/28280", na.strings=c("NA", "", "N/A"))

#shorten dataset name
LandscapeData <- Landscape.affects.pest.and.crop.yield
#remove the long name dataset, so we accidentally don't end up 
#working with the wrong data
rm(Landscape.affects.pest.and.crop.yield)
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
library(ggplot2)#ggplot2_3.3.6  
library(visreg)#visreg_2.7.0 
library(effects)#effects_4.2-0  
library(ggpubr)#vesion 0.4.0
library(dplyr)#version 1.08

######## REMOVING MISSING VALUES ########
#Removing the rows with missing data from each variable individually (this way we don't delete more data than we need to)

names(LandscapeData)
FleaBeetlesIncidence <- LandscapeData[!is.na(LandscapeData$FleaBeetles_incidence), ]
FleaBeetlesAbundance <- LandscapeData
LepidopteranIncidence <- LandscapeData[!is.na(LandscapeData$Lepidoptera_incidence), ]
LepidopteranAbundance <- LandscapeData[!is.na(LandscapeData$Lepidoptera_abundance), ]
AphidsIncidence <- LandscapeData[!is.na(LandscapeData$Aphid_incidence), ]
ParasitoidHostRatio <- LandscapeData[!is.na(LandscapeData$ParasitoidHostRatio), ]
PlantDamageIndex <-LandscapeData[!is.na(LandscapeData$Plant_damage), ]


#######CREATING PLOTS BASED ON LME MODEL PREDICTIONS ########
#making model predictions and calculating confidence intervals

##########Aphid Incidence -Figure 4a##############

#The first step is to fit a simplified mixed-effect model (nlme)
fitlme <- lme(sqrt(Aphid_incidence)~  mead_250+Year, data = AphidsIncidence, random=~1|Farm_ID/Plot_ID)
LMMOutput <- data.frame(fitlme$coefficients$fixed)
colnames(LMMOutput)[1] <- 'Value'

#First set of data to EXTRACT from anova. Here, we run an ANOVA (Type III sum of squares) and would like to capture: numDF (degrees of freedom of the numerator), denDF (degrees of freedom of the denominator), the F-value (test statistic from the F test), and the associated p-value for all three rows - Intercept, Mead_250 (proportion of meadows within a 250 meter radius of the field plot) and Year (sampling year)
anovaOutput <- data.frame(anova(fitlme,type='marginal'))
anovaOutput

#Second set of data to EXTRACT from summary output. Here we would like to extract the information associated with the fixed effects: Value (slope estimates), Std.Error (approximate standard error of the slope estimates), DF (denominator degrees of freedom), t- value (ratios between slope estimates and their standard errors), p-value (associated p-value from a t-distribution)
sum1 <- data.frame(summary(fitlme)$tTable, check.names=FALSE)
sum1 
sampleSize <- fitlme$dims$N

#Next, create a new dataset from which we will make predictions
newdat.lme <- data.frame(Year = AphidsIncidence$Year,
                         mead_250 = AphidsIncidence$mead_250,
                         Aphid_incidence=AphidsIncidence$Aphid_incidence
)
head(newdat.lme)

#I can add the predicted values to the dataset.
newdat.lme$predlme <- predict(fitlme, newdata = newdat.lme, level = 0)
#And then use these in geom_line() to add fitted lines based on the new predlm variable. (figure code just provided for reference)
#ggplot(AphidsIncidence, aes(x = mead_250, y = Aphid_incidence, color = Year) ) +
# geom_rug(sides = "b", size = 1) +
#geom_line(data = newdat.lme, aes(y = predlme), size = 1)

#create confidence intervals for lme objects
des <- model.matrix(formula(fitlme)[-2], newdat.lme)
#Then we use matrix multiplication on the model matrix and variance-covariance matrix
#extracted from the model with vcov(). We pull out the values on the diagonal, which are the
#variances of the predicted values.
predvar <- diag( des %*% vcov(fitlme) %*% t(des) )
#I add the confidence interval limits to the dataset for plotting.
newdat.lme$lower <- with(newdat.lme, predlme - 2*sqrt(predvar) )
newdat.lme$upper <- with(newdat.lme, predlme + 2*sqrt(predvar) )

#Here's the code to generate the figure, with a confidence envelope for each line. (Again, figure just provided for reference, but not integral to code to generate data for extraction)
p1 <- ggplot(AphidsIncidence, aes(x = mead_250, y = Aphid_incidence, color = Year) ) +
  geom_point()+
  geom_rug(sides = "b", size = 1) +
  geom_ribbon(data = newdat.lme, aes(y = NULL, ymin = lower, ymax = upper,
                                     color = NULL, fill = Year),
              alpha = .15) +
  geom_line(data = newdat.lme, aes(y = predlme), size = .75)+
  theme_classic()+xlab('Proportion of meadows at 250-m')+ylab('Proportion of plants infested by aphids')
p1

# Save ggplot figure as png
# ggsave("Fig.4a.png", plot = p1, scale=0.5)

# Remove line 17 to replicate Fig. 4a
LandscapeData <- LandscapeData[-17,]
###########APHID incidence--FIGURE 4a########################
#Figure 4a--Keeping the default settings in ggplot
PlotAphidsIncidence.default <- ggplot(LandscapeData, aes(x=mead_250, y=Aphid_incidence, colour = Year)) +   geom_point(shape=1) + geom_smooth (aes (x=mead_250, y=Aphid_incidence, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE) + theme(panel.background = element_rect(fill='white', colour='black'))+scale_y_continuous(breaks=c(0,0.3, 0.6, 0.9))+xlab('Proportion of meadows at 250-m')+ylab('Proportion of plants infested by aphids')
PlotAphidsIncidence.default#we can ignore the warning. It is just telling us that there is missing data
#Figure4a--Reproducing the plot theme selected in the paper (the way the plot looks like in the paper)
PlotAphidsIncidence.papertheme <- ggplot(LandscapeData, aes(x=mead_250, y=Aphid_incidence, colour = Year, linetype = Year)) +  geom_point(aes(shape = factor(Year)), size = 2) + geom_smooth (aes (x=mead_250, y=Aphid_incidence, colour=factor(Year)), method=lm, se=TRUE, fullrange=TRUE,show.legend = NA) + 
  scale_color_grey (start=0.7, end=0.2)+ scale_linetype_manual(values = c("dashed", "solid"))+ scale_shape_manual(name = "Year",values = c(16, 17))+theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(legend.position = "top")+xlab('Proportion of meadows at 250-m')+ ylab('Proportion of plants infested by aphids')
PlotAphidsIncidence.papertheme
#add equation to the plot---legend location need to be adjusted 
PlotAphidsIncidence.papertheme + stat_regline_equation()

# Save ggplot figure as png
# ggsave("Fig.4a.2.png", plot = PlotAphidsIncidence.papertheme, scale=0.5)

#not sure how to extract equation from the ggplot figure. Alternatively, we could run the individual models (lm) 
#First, we need to split the dataset per year
Data2014 <- filter (LandscapeData, Year=="2014")
Data2015 <- filter(LandscapeData, Year=="2015")
#now we can create individual models (lm)

AphidsIncidence2014 <- lm (Aphid_incidence~mead_250, data=Data2014 )
sum2014 <- data.frame(summary(AphidsIncidence2014)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope 
rownames(sum2014) <- c("(Intercept)(2014)","mead_250 (2014)")
AphidsIncidence2015 <- lm (Aphid_incidence~mead_250, data=Data2015 )
sum2015 <- data.frame(summary(AphidsIncidence2015)$coefficients, check.names=FALSE) #the estimate of mead_250 represents the slope
rownames(sum2015) <- c("(Intercept)(2015)","mead_250 (2015)")
sumLR <- rbind(sum2014,sum2015)
#keep in mind that the model significance could be different from the significance
#obtained from the linear mixed-effect models (lme)

#Third set of data to EXTRACT. This would allow someone to recreate FIG. 4a based on the model predictions.
PredictedValuesAphid_incidence <- newdat.lme

#Drop Aphid_Incidence
PredictedValuesAphid_incidence <- subset(PredictedValuesAphid_incidence, select = -c(Aphid_incidence))
PredictedValuesAphid_incidence 
# We would like to extract the data from all rows and columns, EXCEPT for the column labeled Aphid_Incidence as this is just the raw data. Description of data in each column: 
#Year: sampling year (either 2014 or 2015)
#mead_250: same as described above
#predlme-: predicted (based on model) percentage of plants in a field that would have 10 or more aphids on it 
#lower: lower limit of the 95% CI
#upper: Upper limit of the 95% CI

###End of script
