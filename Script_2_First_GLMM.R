##############################################################
# ESA 2024 GLMM Deep Dive Workshop
#
# Exercise 2 :Fit a GLMM
#
# Objective: Specify random and fixed effects structure appropriately
##############################################################

#rm(list=ls())

##############################
# Load packages

library(glmmTMB)
library(DHARMa)
library(tidyr)
library(ggeffects)

##############################
# Load Data
path <- r'[C:\Analysis\glmm\Workshop_Aug_2024_ESA]'
dataPath <- file.path(path, "Data_workshop")

Shrub <- readRDS(file.path(dataPath, "Shrub.rds"))

###############################
# Specify and Fit the Full Model
###############################

# Reminder:
# What we got from data exploration:
# Families to try: Binomial & Beta-Binomial
# We want to see if there is a trend over time
# We want to account for Variation between years and Sites, although we are not interested in estimates


# We are starting with a full model of random effects based on on Piepho & Ogutu (2002)
# Fixed:
# we are interested in a trend so we will include year as fix effect to get an estimate

# Random:
# random site intercept effect to account for the random sample of sites and correlation among all 
# (1|Site) -years visited in the same site

# random year intercept effect to account for correlation among all sites visited in the same year
# (1|Year) - random intercept for each year 

# we didn't see a lot of variation between the trends for each site
# however to be sure we can check random site slope effect to account for variation among site-level trend lines
# (1 + WYear|Site) - this assumes correlation between intercept and slope


# When doing a binomial model with success rates ( in this case presence of shrub is success and absence is failure)
# We layout the model specifying each of them 

# success: Shrub_cover
# failure: 100-Shrub_cover

# (success, failure)

# lets fit a binomial first:
fit_Bin <- glmmTMB(cbind(Shrub$Shrub_cover, 100-Shrub$Shrub_cover) ~ WYear + 
                     (1|Site) + 
                     (1|Year), 
                   data = Shrub,
                   REML = T,                        # this is T when looking at residuals
                   family=binomial(link='logit'))  # start with binomial

summary(fit_Bin)

# Lets look at the residuals using the DHARMa package

testResiduals(fit_Bin, plot = T) 
testDispersion(fit_Bin,plot = T)
testQuantiles(fit_Bin, plot = T)
testUniformity(fit_Bin, plot = T)

# these all look like there is overdispersion
# lets fit a beta binomial instead

fit_betaBin <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + (1|Site) + (1|Year), 
                       data = Shrub,
                       REML = T, 
                       family=betabinomial(link='logit'))  

summary(fit_betaBin)


# Lets look at the residuals 

# When interpreting DHARMA residuals, the residuals are expected to follow a uniform distribution instead of the 
# normal distribution, and are standardized to values between 0 and 1


testResiduals(fit_betaBin, plot = T) 
testDispersion(fit_betaBin,plot = T)
testQuantiles(fit_betaBin, plot = T)
testUniformity(fit_betaBin, plot = T)

# these all look good!
# Beta-Binomial is probably a good choice!

# What if we wanted to check if there is random site slope effect to account for variation among site-level trend lines
# (1 + WYear|Site) - this assumes correlation between intercept and slope

fit_betaBin_v1 <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + 
                            (1|Year)+
                            (1+WYear|Site), # this means random intercept for site and slope time correlated
                          data = Shrub,
                          REML = T, 
                          family=betabinomial(link='logit'))  


#did not converge
# check what's going on:
summary(fit_betaBin_v1)
# if we look at the estimates correlation between site intercept and site slope = 1
# Overparametrization: Model is over fitted when slope and intercept are correlated, correlation = 1

# we can try to fit one without the random intercept and slope correlated

fit_betaBin_v2 <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + 
                            (1|Year)+
                            (1+WYear||Site), # double bar means random intercept for site and slope time not correlated
                          data = Shrub,
                          REML = T, 
                          family=betabinomial(link='logit'))  

# Let's look at the estimates:
summary(fit_betaBin_v2)

# Both these models fit_betaBin_v2 and fit_betaBin: are beta-binomial, and have random intercepts for site and year as well as trend over time
# the difference is that fit_betaBin_v2 accounts for variation in trend for each site, but if we look at the coefficient this value is really small
# we could probably remove the slope for Site because its so small however, we can use AIC to compare models:

# To compare models we need to run them with ML
fit_betaBin <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + (1|Site) + (1|Year), 
                       data = Shrub,
                       REML = F, 
                       family=betabinomial(link='logit'))  

fit_betaBin_v2 <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + 
                            (1|Year)+
                            (1+WYear||Site), # double bar means random intercept for site and slope time not correlated
                          data = Shrub,
                          REML = F, 
                          family=betabinomial(link='logit'))  



AIC(fit_betaBin, fit_betaBin_v2)

# We have very similar AICs, however second model is more complex. We would continue with the most parsimonious model.

# ok! lets continue with our beta-binomial with random intercept model for site and Year


##################################################################
##################################################################
# Check Model Assumptions 
##################################################################
##################################################################

# We did much of this above.  Let's examine a few more assumptions:

#################################################
# Additional residual diagnostics for final model
# Based on Pinheiro and Bates
################################################

# Residual diagnostics
# with REML=T for more accurate estimates
fit_betaBin <- glmmTMB(cbind(Shrub_cover,100-Shrub_cover) ~ WYear + (1|Site) + (1|Year), 
                       data = Shrub,
                       REML = T, 
                       family=betabinomial(link='logit')) 

fit = fit_betaBin 

# we want to see if there are any patterns in the residuals which might indicate we have missing covariates

# Resids v Time
modelDF = Shrub
modelDF$pResid = resid(fit,type="pearson")
modelDF$fittedMS = fitted(fit)  # ,type = "mean_subject")
modelDF$Shrub_cover = modelDF$Shrub_cover

# Plot residuals by year to look for patterns
residByYear <- data.frame(WYear=modelDF$WYear,resid=modelDF$pResid)
residByYear = residByYear[order(residByYear$WYear),]
plot(residByYear$WYear,residByYear$resid,xlab="Year",ylab="Pearson residual")
# looks pretty unform and centered around 0


# Variance uniformity across sites? 
boxplot(modelDF$pResid~modelDF$Site)  

# Variance uniformity across years? 
boxplot(modelDF$pResid~modelDF$WYear)
# looks good


# examine temporal autocorrelation in residuals
range(table(modelDF$Site))  # all Sites visited 24 times
Sites <- unique(modelDF$Site)
length(Sites) # 100 sites

# Plot autocorrelation plots for each site.
# We expect our random effects structure to account for this!
# Autocorrelation in the residuals suggests that there’s a relationship or dependency between 
# current and past errors in the time series. 

par(mfrow=c(3,3))
for(i in 1:200) {
  print(pacf(ts(residuals(fit)[modelDF$Site==Sites[i]]),main=Sites[i]))
}
dev.off()
#...

# no autocorrelation in the residuals!

#################################################
# Assess that random effects are distributed normally
#################################################
ranef(fit)$cond$Year
ranef(fit)$cond$Site


qqnorm(ranef(fit)$cond$Year[,1], main = "Year Intercept Random Effects")
qqline(ranef(fit)$cond$Year[,1])

qqnorm(ranef(fit)$cond$Site[,1], main = "Site Intercept Random Effects")
qqline(ranef(fit)$cond$Site[,1])

dev.off()


#################################################
# Assess that residuals are independent from 
# levels of random effects 
#################################################


random_effects <- ranef(fit)
SiteRE = data.frame(Site = rownames(random_effects$cond$Site),
                    SiteRE = random_effects$cond$Site)
names(SiteRE)[2] = 'SiteRE_Intercept'

rownames(SiteRE)=NULL


YearRE = data.frame(Year=rownames(ranef(fit)$cond$Year),YearRE=ranef(fit)$cond$Year)
names(YearRE)[2] = 'YearRE'
rownames(YearRE)=NULL
YearRE <- tidyr::separate(data = YearRE, col = Year, into = c("Year"), sep = ":")

modelDF = merge(modelDF, SiteRE)
modelDF = merge(modelDF, YearRE)

# Is the residual value correlated with the year random effect?
cor.test_Year = cor.test(modelDF$pResid, modelDF$YearRE, method='kendall')
print(cor.test_Year)#Some correlation between residuals and slope but very small
# cor = 0.03837626
# p-value = 0.00575

# How about site?
cor.test_Site = cor.test(modelDF$pResid, modelDF$SiteRE_Intercept, method='kendall')
print(cor.test_Site) #Some correlation between residuals and slope but very small
# cor = 0.02375734
# p-value = 0.08252

# Some significant correlation test, but the correlation is very small. 
# We could consider other model structure, but in this simulated data set 
# we know that these variables are independent. 

#################################################################
##################################################################
# Biological Conclusions
##################################################################
##################################################################

# This uses z-tests:
# Examine fixed effects
summary(fit)
#positive significant trend

#               Estimate  Std. Error  z value     Pr(>|z|)
# (Intercept) 0.78371581 0.098943188 7.920867 2.358605e-15
# WYear       0.02314495 0.005188291 4.460997 8.157938e-06

# Examine whether confidence intervals contain zero for FE:
confint(fit) # Do CIs include 0?

#                           2.5 %       97.5 %   Estimate
# (Intercept)              0.58979073 0.97764089 0.78371581
# WYear                    0.01297609 0.03331382 0.02314495
# Std.Dev.(Intercept)|Site 0.59961275 0.82458339 0.70315767
# Std.Dev.(Intercept)|Year 0.07773768 0.22368498 0.13186641

# Binomial and beta-binomial GLMMs model the cover proportion with a logit link function
# trend is generated as linear on the scale of the odds of the cover probability
# The odds is p/(1-p)



# WYear
# obtain estimate and 95%-CI of trend on the odds scale
(exp(confint(fit)[2,])-1)[c(3,1,2)]

# 0.0234% trend over time


# Marginal effects plots # type="fe" accounting for fixed effect variation
# Predicted Annual Status
margEffects_WYear <- ggpredict(fit, terms = c("WYear"),ci.lvl = 0.95, type="fe")
margEffects_WYear
plot(margEffects_WYear,rawdata=FALSE,jitter=0)




