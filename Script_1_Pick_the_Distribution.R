######################################################################################################
# ESA 2024 GLMM Deep Dive
#
# Exercise 1: Identify Candidate Conditional Distributions of Response Variables and Link Functions
#
# Objective: Use graphical tools to determine candidate conditional distributions & potential structure for random & fixed effects
###################################################################################################


rm(list=ls())

#################################
# Load packages
require(ggplot2)
require(dplyr)


## set your path ##
getwd()
path <- r'[C:\Analysis\glmm\Workshop_Aug_2024_ESA]'# Set your own path
dataPath <- file.path(path, "/Data_workshop")
  
#################################
# Study 1 - Marine Debris - 
#################################
# Background on the data: 
# Marine debris was weighed over time throughout a coastal area

Marine_debris <- readRDS(file.path(dataPath, "Marine_debris.rds"))

######################
# Examine data
######################

# What type of response variable is debris?

range(Marine_debris$debris) # 0.000 304.663
table(Marine_debris$debris)

# Type of data?
# Continuous, Count, presence/absence, percentage....? 
# continuous
# We can use the flow chart in the presentation to help us decide what distributions we can consider.

# Look at how distribution of data looks like:
hist(Marine_debris$debris)
# not normal

# What type of distribution represents continuous data?
# ________ & _________


# zeros
table(Marine_debris$debris==0)
# FALSE  TRUE 
# 28754    46 

# This makes our life easy! gamma distribution represents data >0 so in this case tweedie would be our best choice

# Lets explore the data a bit more to confirm and look at potential structures for random and fixed effects

# boxplot - year
plot(Marine_debris$Year, Marine_debris$debris)
# variability between years looks inconsistent: meaning including random intercept for year might be a good to account to differences in dispersion by year

# boxplot - site
plot(Marine_debris$Site, Marine_debris$debris)
# Here is clear that variation is not consistent between sites: including random site intercept might be good to account for differences in dispersion by site


# look at the means by site
# If we were interested in estimating effect for each site we could include it as fixed effect
meanY <- aggregate(list(meanY=Marine_debris$debris), 
                   list(Site=Marine_debris$Site),
                   mean)

plot(meanY[,1:2],
     pch = 19,
     ylab = 'mean(debris)',
     xlab = 'Site')


# look at the means by year
meanYr <- aggregate(list(meanYr=Marine_debris$debris), 
                    list(WYear=Marine_debris$WYear),
                    mean)

plot(meanYr[,1:2],
     pch = 19,
     ylab = 'mean(debris)',
     xlab = 'WYear')
lines(meanYr)
# nice trend

# tweedie uses log link to linearize the mean so lets check that:
plot(log(meanYr[,1:2]),
     pch = 19,
     ylab = 'mean(debris)',
     xlab = 'WYear')
lines(log(meanYr))


# Examine variance
# For each distribution there is an assumed relationship between the mean and the variance.
# We can use that relationship to help select our distribution.
# Tweedie variance = mean^P  (P = power parameter between 1 and 2 )
# gamma variance = mean*scale
# we can use our data to do a rough scale calculation for gamma  & tweedie

varY <- aggregate(list(varY=Marine_debris$debris), 
                  list(Site=Marine_debris$Site),
                  var)
meanVar <- merge(meanY,varY)

scaleY <- var(meanVar$varY)/mean(meanVar$meanY)

ggplot(meanVar, aes(meanY, varY)) +
  geom_smooth(method='lm', formula = y ~I(x*scaleY), se=F, aes(color='blue')) + 
  geom_smooth(method='lm', formula = y ~ I(x^1.6) + offset(x) -1, se=F, aes(color='red')) + # you can play with different values for p (between 1 and 2) 
  geom_point() +  
  scale_color_identity(name = "Model fit",
                       breaks = c("blue","red"),
                       labels = c("gamma", "tweedie"),
                       guide = "legend")


# This further supports that Tweedie is likely the most suitable distribution

                                          # THE END #



#################################
# Study 2 - Shrub Cover - 
#################################
# Background on the data: 
# These data are from a sample of point-line transects. 
# Each transect has 100 points at which a species of interest has a hit (1) or a miss (0)
# The response variable Shrub is the number of hits for the transect at site i (shrub coverage)

######################
# load data
######################

Shrub <- readRDS(file.path(dataPath, "Shrub.rds"))


######################
# Examine data
######################

# What type of response variable is Shrub?

range(Shrub$Shrub_cover) # 0 - 100.
table(Shrub$Shrub_cover) # All integers.

# Continuous, Count, presence/absence, percentage....?

# Since we have presence/absence on the same site 100 times our response variable comes in percent

# Look at distribution
hist(Shrub$Shrub_cover)

# What are distributions that represent percentages? Try to guess before moving forward:
# __________ or _________________


# Binomial data can either be modeled at the individual (binary response) or group (proportion) level.
# The binomial distribution arises when n independent Bernoulli trials are conducted, 
# each with the same probability of success p. 
# The number of successes is a binomial(n, p) random variable.

# The beta-binomial(n, α, β) distribution is  generated by choosing the probability p for a binomial(n, p) 
#  distribution from a beta(α, β) distribution.

# One major difference between a binomial and beta-binomial is that in a binomial distribution,
# p is fixed for a set number of trials; in a beta-binomial, p is not fixed and changes from trial to trial.
# Usually beta-binomial have higher dispersion than a binomial

# Lets keep exploring the data:

# A large number of zeroes may indicate that we need to use a zero-inflated model
table(Shrub$Shrub_cover==0)
# FALSE  TRUE 
# 2396     4 

# Now lets look at our data vs potential predictors, this will help us build the model later.
# Looking at responses over predictors can guide us on what to input as random or fixed effects.


# plot response variable as a function of Year
ggplot(Shrub, aes(WYear, Shrub_cover)) +
  geom_point(alpha = 0.05, position = position_jitter(height = 0.15)) +
  ylab("Percent Cover") +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2)
# Looks like Percent Cover increases over time, so if we are interested in looking at trend, time could be included as fixed effect to get an estimate.


# box plot response variable as a function of Site
ggplot(Shrub, aes(Site, Shrub_cover)) +
  geom_boxplot() +
  ylab("Percent Cover") 
# There is high variability between sites. 
# This suggests that Site might be a good variable to incorporate as random variable to account for heterogeneity between Sites.

# rough trend for each site over the years.
# just plotting a few so we can visualize
  ggplot(Shrub[Shrub$Site%in% c(1:50),], aes(WYear, y=Shrub_cover/100, colour=Site )) +
    geom_point() +
    geom_smooth(
      method="glm",
      method.args=list(family=binomial(link='logit')),
                        se = FALSE)

# Just plotting a few we can see that each site varies in intercept but slopes are pretty similar
# We would want to include random intercept for site and maybe check random slope for site as well when we are fitting the model
  
  
# box plot response variable as a function of Year
ggplot(Shrub, aes(Year, Shrub_cover)) +
  geom_boxplot() +
  ylab("Percent Cover") 
# There is high variability between Years
# Incorporating a random intercept for Year might be a good idea to account for heterogeneity between years
# However this high variability could also indicate that this data could be fitted with a beta binomial which accounts for higher dispersion


# lets look at the means #

# means by site
mean_cover <- aggregate(list(Shrub$Shrub_cover), 
                        list(Site=Shrub$Site),
                        mean)

plot(mean_cover[,1:2],
     pch = 19,
     ylab = 'mean(Shrub_cover)',
     xlab = 'Site')
# we could try to include site as fixed effect however in this study we are interested in getting a coefficient per site, we want to look at an overall trend


# means by Year
mean_coverYr <- aggregate(list(mean_coverYr = Shrub$Shrub_cover), 
                        list(WYear = Shrub$WYear),
                        mean)

plot(mean_coverYr[,1:2],
     pch = 19,
     ylab = 'mean(Shrub_cover)',
     xlab = 'WYear')

lines(mean_coverYr)
# looks like positive trend over time again

# Examine logit which is the log link for binomial and beta-binomial
logit <- function(x) log(x/(1-x))

plot(mean_coverYr$WYear, logit(mean_coverYr$mean_coverYr/100),
     pch=19,
     ylab='logit[mean(Shrub_cover)]',
     xlab='WYear')  
lines(mean_coverYr$WYear, logit(mean_coverYr$mean_coverYr/100))



# What did we gather?
# Percent cover data so could be binomial or beta-binomial
# High variability between sites and years - potential for random intercepts for site and year
# sites seem to have similar slope but different intercepts - potential for random slope for site
# positive trend over time - could incorporate time as a fixed effect since we are interested in trend

# Sometimes just by using data exploration we can't figure out the distribution and we need to fit the model 
# and look at the residuals and coefficients
# We will continue with this example during model selection in script 2

# TO BE CONTINUED ......
#########################################


# NEXT EXERCISE #


########################################
# Study 3 - Live native shrub density 
########################################
# Background on the data: 
# These data represent number of live native shrubs counted in 30 m^2 plots in two of the Channel Islands 
# off the California Coast. We are interested in determining trends in live native shrub density over time
# in the two Islands

######################
# load data
######################

Vegetation <- readRDS(file.path(dataPath, "Vegetation.rds"))


######################
# Examine data
######################

# what type of response variable do we have? 
Vegetation$Y
range(Vegetation$Y) #0-26



# What kind of data do we have? What are candidate distributions?
# Continuous, Count, presence/absence, percentage....?
# tweedie, nbinom or poisson

# Response variable by site by year
Sites <- sort(unique(Vegetation$Site))
plot(jitter(Vegetation$WYear),Vegetation$Y,col=Vegetation$Site)
for (i in 1:length(Sites)) lines(Vegetation$WYear[Vegetation$Site==Sites[i]],
                                 Vegetation$Y[Vegetation$Site==Sites[i]],
                                 col=i)
# Are there any relationships?

# Does it look like any distribution?                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
ggplot(Vegetation, aes(x=Y, fill=Island)) +
  geom_histogram(alpha=0.5, position="identity") +
  facet_wrap(~Island)
# Could be ____________, ________ or __________


# zero inflated?
table(Vegetation$Y==0)
#could be ZI, there are lots of zeros

# We could consider a zero inflated model, however zero inflated models should not be the first model to fit:
# First the overall mean is low covariates might help account for 0s
# Also certain distributions can account for a large number of zeros

#########################################################
# Look at the means of the response variable by SITE
meanY2 <- Vegetation %>% group_by(Island,Site) %>% summarise(meanY2 =mean(Y))
meanY2 

ggplot(meanY2, aes(Site, meanY2, colour = Island)) + geom_point()
# Looks like there are higher values for Island 2
# if we want to compare islands, incorporating Island as fixed effect could be a good idea


# We can get a feeling for what distribution might work best by plotting the mean and variance for each group level:
# For Negative Binomial 1 Variance =  μ+σμ so for NB2 variance is quadratic function of the mean
# For Negative Binomial 2 Variance =  μ+μ2/θ= μ+σμ2

 
varY2 <- Vegetation %>% group_by(Island,Site) %>% summarise(varY2 = var(Y))
meanVarY2 <- merge(meanY2,varY2)

ggplot(meanVarY2, aes(meanY2, varY2)) + 
  geom_smooth(method='lm', formula = y ~x -1, se=F, aes(color='blue')) + 
  geom_smooth(method='lm', formula = y ~ I(x^2) + offset(x) -1, se=F, aes(color='red')) +
  geom_abline(aes(intercept =0, slope=1, color='black')) +
  geom_point() +  
  facet_wrap(~Island) +
  scale_color_identity(name = "Model fit",
                       breaks = c("blue","red",'black'),
                       labels = c("NB1", "NB2", "Poisson"),
                       guide = "legend")

# In the above plot, the black line represents the Poisson (one-to-one), 
# the blue line represents the NB1 negative binomial parameterization,
# and the red line represents the NB2 negative binomial parameterization.

# could be NBINOM1 or NBINOM2, although looks a lot like nbinom2
# looks like poisson is a bad fit

#########################################################
# look at the means by YEAR #
meanY2Yr <- Vegetation %>% group_by(Island,WYear) %>% summarise(meanY2 =mean(Y))
meanY2Yr

ggplot(meanY2Yr, aes(WYear, meanY2, colour = Island)) + geom_point()+ geom_smooth(method = "lm",se = FALSE)


# Examine log link
ggplot(meanY2Yr, aes(WYear, log(meanY2), colour = Island)) + geom_point() + geom_smooth(method = "lm", se = FALSE)
# we can see a trend over time, and a separate trend per Island

# Examine variance again
varY2Yr <- Vegetation %>% group_by(Island,WYear) %>% summarise(varY2 = var(Y))
meanVarY2Yr <- merge(meanY2Yr,varY2Yr)

ggplot(meanVarY2Yr, aes(meanY2, varY2)) + geom_smooth(method='lm', formula = y ~x -1, se=F, aes(color='blue')) + 
  geom_smooth(method='lm', formula = y ~ I(x^2) + offset(x) -1, se=F, aes(color='red')) +
  geom_abline(aes(intercept =0, slope=1, color='black')) +
  geom_point(aes(shape=Island)) +  
  scale_color_identity(name = "Model fit",
                       breaks = c("blue","red",'black'),
                       labels = c("NB1", "NB2", "Poisson"),
                       guide = "legend") +
  facet_wrap(~Island)

# NBINOM2 is looking pretty good

#########################################################
# We can fit a quick linear regression to see response by Year
# plot Y+1 so that 0's can be visualized
ggplot(Vegetation, aes( WYear, log(Y+1), colour = Island)) +
  geom_point(position = position_jitter(height = 0.4, width = 0)) +
  geom_smooth(method = "lm")
# what do we see? Variation in slope and intercept 

# remove points to check trend
ggplot(Vegetation, aes( WYear, log(Y+1), colour = Island)) +
  geom_smooth(method = "lm")
# (Y+1) to account for zeros 

# we can fit a quick linear regression to see response by Site
ggplot(Vegetation[Vegetation$Island==2,], aes( WYear, log(Y+1), colour = Site)) +
  geom_point(position = position_jitter(height = 0.4, width = 0)) +
  geom_smooth(method = "lm", se = FALSE)+
  theme(legend.position = "none")
# There is variation in intercept and slope over time for each site
# There is variation in the trends by site, so a random slope for site

# check variation by Island #

# Variance uniformity across Islands? 
ggplot(meanVarY2Yr, aes(x=Island, y=varY2)) +
  geom_boxplot()
# Difference in variance by Island?

# Variance uniformity across Years? 
ggplot(meanVarY2Yr, aes(x=as.factor(WYear), y=varY2)) +
  geom_boxplot() 
# Difference in the variance by year?

# The information gathered from data exploration that we will take to model fitting:

# Candidate distributions are: Negative Binomial 2 (most likely), although we might want to check Negative Binomial 1
# We will use log link for this

# Plotting Year and Site with the response variable we saw that:
# We want to account for: 

### Variation between years and Sites ###
# Heterogeneity between Islands
# Potential for random effects to account for these if we are not interested in estimating them #

### Trend over time ###
# Since we are interested in estimating trend, we might want to include it as fixed effect

# To be continued............



