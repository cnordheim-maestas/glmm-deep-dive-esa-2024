#######################################################################
# ESA 2024 GLMM Deep Dive Workshop
#
# Exercise 3 :Fit a GLMM with a design variable
#
# Objective: Specify random and fixed effects structure appropriately
#   and identify sources of heterogeneity
#######################################################################

#rm(list=ls())

#################################
# Load packages

library(glmmTMB)
library(DHARMa)
library(tidyr)
library(ggeffects)


#################################
#set paths
#################################
# set your working directory to ".../Workshop_Aug_2024_ESA"
#path <- r'[C:\Analysis\glmm\Workshop_Aug_2024_ESA]'
#codePath <- file.path(path, 'Code_workshop')
dataPath <- file.path(path, 'Data_workshop')

######################################
# Define functions
######################################
logit <- function(x) log(x/(1-x))
expit <- plogis


########################################################
# Example 3 data set
########################################################

# read data
Vegetation <- readRDS(file.path(dataPath,'Vegetation.rds'))


###############################
# Specify and Fit the Full Model
###############################

# Reminder:
# What we got from data exploration:
# Families to try: NB
# Trend over time
# Variation between years and Sites
# Heterogeneity between Islands

# We are starting with a full model of random effects based on on Piepho & Ogutu (2002)
# random site intercept effect to account for the random sample of sites and correlation among all 
#    years visited in the same site
# random site slope effect to account for variation among site-level trend lines
# (1 + WYear|Site) - this assumes correlation between intercept and slope

# random year intercept effect to account for correlation among all sites visited in the same year
#  (1|Year) - random intercept for each year 

# We will fit a few candidate distributions and evaluate the fit using simulation via 
# the DHARMa package.
# we have NBI and NB2 as options we would check both of these

# We fit the final model with REML to get accurate variance estimates. 
NB1_fit_a <- glmmTMB(Y ~ WYear + (1+WYear|Site)+(1|Year), 
                data = Vegetation,
                REML = T, 
                family=nbinom1(link='log')) 
summary(NB1_fit_a)


NB2_fit_a <- glmmTMB(Y ~ WYear + (1+WYear|Site)+(1|Year), 
                      data = Vegetation,
                      REML = T, 
                      family=nbinom2(link='log')) 
summary(NB2_fit_a)


# simulate residuals and save them
NB1_fit_a_simResids <- simulateResiduals(NB1_fit_a, plot=T)
NB2_fit_a_simResids <- simulateResiduals(NB2_fit_a, plot=T)
# nbinom2 has better fit

# When interpreting DHARMA residuals, the residuals are expected to follow a uniform distribution instead of the 
# normal distribution, and are standardized to values between 0 and 1

testResiduals(NB1_fit_a_simResids, plot = T)   #  
testResiduals(NB2_fit_a_simResids, plot = T)   #  
#  

testZeroInflation(NB1_fit_a_simResids)        
testZeroInflation(NB2_fit_a_simResids)      
# 

testQuantiles(NB1_fit_a_simResids, plot = T)  #  
testQuantiles(NB2_fit_a_simResids, plot = T)  #  

# Proceed with NB2

# check heterogeneity by island 
dat <- data.frame(Vegetation,resids=NB2_fit_a_simResids$scaledResiduals)
ggplot(dat, aes( WYear, resids, colour = Island)) +
  geom_smooth(method = "lm")
ggplot(dat, aes(x=Island, y=resids)) +
  geom_boxplot()

##############################################################################
# Now we are going to add:   separate trends by Island
##############################################################################
NB2_fit_b <- glmmTMB(Y ~ WYear*Island + (1+WYear|Site)+(1|Year), 
                        data = Vegetation,
                        REML = T, 
                        family=nbinom2(link='log')) 
summary(NB2_fit_b)

# simulate residuals and save them
NB2_fit_b_simResids <- simulateResiduals(NB2_fit_b, plot=T)

# lets compare both models with the same function.
testResiduals(NB2_fit_b_simResids, plot = T)  
# 

testZeroInflation(NB2_fit_b_simResids)        
# 

testQuantiles(NB2_fit_b_simResids, plot = T) 


# check heterogeneity by island 
dat <- data.frame(Vegetation,resids=NB2_fit_b_simResids$scaledResiduals)
ggplot(dat, aes( WYear, resids, colour = Island)) +
  geom_smooth(method = "lm")
ggplot(dat, aes(x=Island, y=resids)) +
  geom_boxplot()


###################################################################################
# Model c:  Trend by Island
# model random site effect independently by Island  
# random variation in intercept among Years within Islands (nested)
###################################################################################

NB2_fit_c <- glmmTMB(Y ~ WYear*Island + (-1+Island|Site)+ (-1+WYear|Site)+(1|Year), 
                        data = Vegetation,
                        REML = T, 
                        family=nbinom2(link='log')) 

summary(NB2_fit_c)  #  
# did not converge

# simulate residuals and save them
NB2_fit_c_simResids <- simulateResiduals(NB2_fit_c, plot=T)

testResiduals(NB2_fit_c_simResids, plot = T)  #  


# check heterogeneity by island 
dat <- data.frame(Vegetation,resids=NB2_fit_c_simResids$scaledResiduals)
ggplot(dat, aes( WYear, resids, colour = Island)) +
  geom_smooth(method = "lm")
ggplot(dat, aes(x=Island, y=resids)) +
  geom_boxplot()


###################################################################################
# Model d: Trend by Island
# different random intercept for year by island  
# Site intercept and slope modeled across Islands
###################################################################################

NB2_fit_d <- glmmTMB(Y ~ WYear*Island + (1+WYear|Site)+(-1+Island|Year),
                        data = Vegetation,
                        REML = T, 
                        family=nbinom2(link='log')) 

summary(NB2_fit_d)  #  
# Year variation by Island is very different

# simulate residuals and save them
NB2_fit_d_simResids <- simulateResiduals(NB2_fit_d, plot=T)

testResiduals(NB2_fit_d_simResids, plot = T)  # sign. outlier test now

# check heterogeneity by island 
dat <- data.frame(Vegetation,resids=NB2_fit_d_simResids$scaledResiduals)
ggplot(dat, aes( WYear, resids, colour = Island)) +
  geom_smooth(method = "lm")
ggplot(dat, aes(x=Island, y=resids)) +
  geom_boxplot()

##############################################################
# Several models have good residuals how do we choose?
# fit with REML=FALSE to check AIC
bbmle::AICtab(NB2_fit_a,NB2_fit_b,NB2_fit_c,NB2_fit_d)

#          dAIC  df
#NB2_fit_d   0.0 11 <-- true model
#NB2_fit_b  42.8 9 
#NB2_fit_a 186.1 7 
#NB2_fit_c    NA 10

# What makes more biological sense, whats the question we are asking?
summary(NB2_fit_d) # best candidate, Nbinom2 better than Nbinom1

# Note: AIC doesn't tell us how well the model fits
# We still need to check assumptions of the mixed model

##################################################################
##################################################################
# 4. Check Model Assumptions 
##################################################################
##################################################################

# We did much of this above.  Let's examine a few more assumptions:

#################################################
# Additional residual diagnostics for final model
# Based on Pinheiro and Bates
#################################################
# Residual diagnostics

fit = NB2_fit_d
summary(fit)

# Resids v Time 
modelDF = Vegetation
modelDF$pResid = resid(fit,type="pearson")
modelDF$fittedMS = fitted(fit)  # ,type = "mean_subject")
modelDF$y = modelDF$Y

# Plot residuals by year to examine any patterns and correlation
residByYear <- data.frame(WYear=modelDF$WYear,resid=modelDF$pResid)
residByYear <- data.frame(WYear=modelDF$WYear,resid=modelDF$pResid,resid_response=modelDF$pResid_response,resid_working=modelDF$pResid_working,resid_deviance=modelDF$pResid_deviance)
residByYear <- residByYear[order(residByYear$WYear),]
plot(residByYear$WYear,residByYear$resid,xlab="Year",ylab="Pearson residual")

# Variance uniformity across sites? 
boxplot(modelDF$pResid~modelDF$Site)  

# Variance uniformity across years? 
boxplot(modelDF$pResid~modelDF$WYear)
# looks good

# Variance uniformity across islands? 
boxplot(modelDF$pResid~modelDF$Island)
# looks good

# examine temporal autocorrelation in residuals
range(table(modelDF$Site))  # all Sites visited 24 times
Sites <- unique(modelDF$Site)
length(Sites)

# Plot autocorrelation plots for each site.
# We expect our random effects structure to account for this!
# We just want to make sure there is no pattern here for each site,
# though we probably could not pick it up given how few years of data we have.
par(mfrow=c(3,3))
for(i in 1:200) {
  print(pacf(ts(residuals(fit)[modelDF$Site==Sites[i]]),main=Sites[i]))
}
dev.off()
#...

# no autocorrelation in the residuals 

#################################################
# Assess that random effects are distributed normally
#################################################
ranef(fit)$cond$Year
ranef(fit)$cond$Site

par(mfrow=c(2,2))

qqnorm(ranef(fit)$cond$Year[,1], main = "Random Year Intercept Effects")
qqline(ranef(fit)$cond$Year[,1])

qqnorm(ranef(fit)$cond$Site[,1], main = "Random Site Intercept Effects")
qqline(ranef(fit)$cond$Site[,1])

qqnorm(ranef(fit)$cond$Site[,2], main = "Random Site Slope Effects")
qqline(ranef(fit)$cond$Site[,2])

dev.off()

#################################################
# Assess that residuals are independent from 
# levels of random effects 
#################################################

random_effects <- ranef(fit)
SiteRE = data.frame(Site = rownames(random_effects$cond$Site),
                    SiteRE = random_effects$cond$Site)
names(SiteRE)[2] = 'SiteRE_Intercept'
names(SiteRE)[3] = 'SiteRE_Slope'
rownames(SiteRE)=NULL

YearRE <- rbind(data.frame(Year=rownames(ranef(fit)$cond$`Year`),Island=1, YearRE=ranef(fit)$cond$`Year`[,1]),
                data.frame(Year=rownames(ranef(fit)$cond$`Year`),Island=2, YearRE=ranef(fit)$cond$`Year`[,2]))

modelDF = merge(modelDF, SiteRE)
modelDF = merge(modelDF, YearRE)

# Is the residual value correlated with the Year random effect?
cor.test_Year = cor.test(modelDF$pResid, modelDF$YearRE, method='kendall')
print(cor.test_Year) # p-value = 0.03095
# cor = 0.03140531

# Is the residual value correlated with the Site random effect?
cor.test_Site = cor.test(modelDF$pResid, modelDF$SiteRE_Intercept, method='kendall')
print(cor.test_Site) #  p-value = 0.01747
# cor = 0.0229327

# Is the residual value correlated with the Site Slope random effect?
cor.test_SiteSlope = cor.test(modelDF$pResid, modelDF$SiteRE_Slope, method='kendall')
print(cor.test_SiteSlope) # p-value = 0.3688
# cor = 0.008671914


#################################################################
##################################################################
# Biological Conclusions
##################################################################
##################################################################


# Interpret results
summary(fit)
# Obtain confidence intervals for FE & RE:
confint(fit)

# Obtain estimates of trend by Island with CIs
# define contrasts
c_Island1 <- matrix(c(0,1,0,0),4,1)
c_Island2 <- matrix(c(0,1,0,1),4,1)

# obtain FE coefficients and variance of FE coefficients
beta <- matrix(fixef(fit)$cond,4,1)
VarBeta <- vcov(fit)$cond

trendEst <- c(t(c_Island1) %*% beta, t(c_Island2) %*% beta)
trendEst_Var <- c(t(c_Island1) %*% VarBeta %*% c_Island1, 
                 t(c_Island2) %*% VarBeta %*% c_Island2)
trendEst_SE <- sqrt(trendEst_Var)
trendEst
trendEst_SE 

# Conduct two-sided Z-test for the trend on each Island
2*pnorm(abs(trendEst/trendEst_SE), lower.tail=F)


# Marginal effects trend plots with confidence intervals
# FE variance only
margEffects_WYear <- ggpredict(fit, terms = c("WYear", "Island"),ci.lvl = 0.95, type="fe")  
margEffects_WYear
plot(margEffects_WYear,rawdata=TRUE,jitter=0)

# FE & RE variance 
margEffects_WYear_RE <- ggpredict(fit, terms = c("WYear", "Island"),ci.lvl = 0.95, type="re")  
margEffects_WYear_RE
plot(margEffects_WYear_RE,rawdata=TRUE,jitter=0)

