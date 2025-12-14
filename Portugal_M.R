##########
# AS3204 #
#####################################
# mortality forecasting with StMoMo #
# Pietro Millossovich #
#######################

#################
# Preliminaries #
#################

# what is StMoMo?
# St = Stochastic, Mo = Mortality, Mo = Modelling

# a package for fitting, analysing, forecasting, simulating, bootstrapping stochastic mortality models 

# by Villegas, Millossovich, Kaishev (2016)

# the package has many existing facilities for implementing standard models, but give freedom to users to build up new models 

# see the vignette at
# https://cran.r-project.org/web/packages/StMoMo/vignettes/StMoMoVignette.pdf
# and the full paper at
# https://www.jstatsoft.org/article/view/v084i03

# in this session we will focus on the Lee-Carter and Cairns-Blake-Dowd (cbd) models

# set your directory here: in RStudio -> Session -> Set Working Directory -> To Source File Location

# load the package
install.packages("StMoMo")
library(StMoMo)
library(help = "StMoMo")

# load the data set for Spain (females)
load(file = "PORTUGAL-males-mortality.RData")

# "PORTUGALStMoMom", "PORTUGALStMoMom": StMoModata objects
# these are lists with several components
print(PORTUGALStMoMom)

# matrices with death counts and expsoures, ages on rows, years on columns
PORTUGALStMoMom$Dxt
PORTUGALStMoMom$Dxt[1:5, 1:5]
PORTUGALStMoMom$Ext[1:5, 1:5]

# ages, years, type, series, label
PORTUGALStMoMom$ages # 0 - 110
PORTUGALStMoMom$years # 1940 - 2022
PORTUGALStMoMom$type # central
PORTUGALStMoMom$series # male
PORTUGALStMoMom$label # portugal

#######################################
# Define a stochastic mortality model #
#######################################

# Lee-Carter (LC)

LC <- lc() # sama je as below!! bc below is the default setting! below is to just show u the ropes

print(LC)

# or, more in detail
LC <- lc(link = "log", const = "sum")
# link = "log" -> bc log(mxt)!! -> Poisson modelling of death counts
# * ??? recall: model is based in death count (Dxt), and Dxt follows a Poisson distrubtion with param = mxt * Ext like below
# Dxt ~ Poisson(mxt * Ext) with log(mxt) = ax+bx*kt 
# * mxt -> central death rate // Ext -> central exposure 
# * // by LC; log(mxt) = ax+bx*kt -> idea behind LC is separating effect of age and effect of time -> two separate effects!
# const = "sum" -> standard constraints on parameters to make model identifiable: \sum b1[x] = 1, \sum k1[t] = 0 -> bc we don't want eg two diff sets of values of params leading to same value of log(mxt)
# * other models will have a similar structure!! just diff expressions for log(mxt)!!! (see CBD later)

# Cairns-Blake-Dowd (cbd, or M5)

CBD <- cbd(link = "log")

# ??? link = "log" -> Poisson modelling of death counts -> this is now necessary as by default the cbd model works with a Binomial distribution and logit function of one year intitial rates instead of the central mortality

# ??? Dxt ~ Poisson(mxt * Ext) with log(mxt) = k1t+(x-xbar) * k2t

print(CBD)

# no constraint needed on parameters; ??? here f2[x] = x - \bar{x}, where \bar{x} is the average age in the data set


####################################
# Fit a stochastic mortality model #
####################################

# set the ages and years to use in fitting (we want to skip periods of time like world war etc)
# ages 55 to 89 (35 ages), years 1960 to 2012 (53 years)
# must be a subset of the available ages and years
ages.fit <- 55:89
years.fit <- 1963:2022

# * try fitting diff sets of years eg 1950 to 2012 etc to see robustness of model
# unlike Lee-Carter, the cbd model can only be fitted on adult mortality data -> cbd has a linear age function (x- xbar), so it doesn't adress stuff like infant mortality and accidental hump

# fit Lee-Carter with ML
# this uses an StMoMo object, ??? and GNM (generalised non-linear model) to fit the data!
LCfit <- fit(object = LC, # model
             data = PORTUGALStMoMom, # a StMoModata object
             ages.fit = ages.fit, # ages in fit
             years.fit = years.fit # years in fit
)

# alternative: specify each element separately -> for if u don't have an StMoMo object (data not from mortality website or wtv) 
LCfit1 <- fit(LC, # model -> uses standard structure; ages on rows and calendar years on columns
              Dxt = PORTUGALStMoMom$Dxt, # death counts
              Ext = PORTUGALStMoMom$Ext, # exposures
              ages = PORTUGALStMoMom$ages, # total ages (that corresponds the the Dxt and Ext matrix)
              years = PORTUGALStMoMom$years, # total years (that corresponds the the Dxt and Ext matrix)
              ages.fit = ages.fit, # ages in fit (fitting process - only a subset of ages corresponding to whole matrix above!)
              years.fit = years.fit, # calendar years in fit (fitting process - only a subset of calendar yrs corresponding to whole matrix above!)
)

# plot actual vs fitted rates
# actual central death rates = death counts / exposure
mxt <- PORTUGALStMoMom$Dxt / PORTUGALStMoMom$Ext
mxt[56:60, 24:29]

# fitted rates
mxt.hat.lc <- fitted(LCfit, type = "rates")
mxt.hat.lc[1:5, 1:5]

# alternatives are 
# type = "link" for log death rates (log(mxt))
# type = "deaths" for death counts

# plot actual rates
plot(years.fit,
     mxt["55", as.character(years.fit)],
     xlab = "year",
     ylab = "central death rate",
     main = "Fitted vs. observed rates at age 55 - Males"
)

# plot fitted rates
lines(years.fit, mxt.hat.lc["55", ], col = "red")

# fit Cairns-Blake-Dowd with ML

CBDfit <- fit(object = CBD,
              data = PORTUGALStMoMom,
              ages.fit = ages.fit,
              years.fit = years.fit
)

# to avoid issues with non integer death counts (bc poisson distribution is discrete!), specify one by one
CBDfit <- fit(object = CBD,
              Dxt = round(PORTUGALStMoMom$Dxt),
              Ext = PORTUGALStMoMom$Ext,
              ages = PORTUGALStMoMom$ages,
              years = PORTUGALStMoMom$years,
              ages.fit = ages.fit,
              years.fit = years.fit
)

# plot fitted rates
mxt.hat.cbd <- fitted(CBDfit, type = "rates")
lines(years.fit, mxt.hat.cbd["55", ], col = "blue")

legend("topright", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

###############################
# Analyse the estimated model #
###############################

# convergence:
LCfit$conv # algo behind $conv is just MLE, and returns if convergence testing via (black box) numerical optimisation methods was successful!
CBDfit$conv
# LCfit and CBDfit contain all the data! - ori data, fitted parameters, info abt fitting process, etc

# output: for instance, for LC
# extract estimated parameters using coef()
coef(LCfit)
coef(LCfit)$ax # a vector
coef(LCfit)$bx # a matrix
coef(LCfit)$kt # a matrix

LCfit$loglik # log-likelihood

# sample size = length(ages.fit) * length(years.fit)
LCfit$nobs # number of observations (aka sample size)

# dof = length(ages.fit) + length(ages.fit) + length(years.fit) - 2
# number of ax (one for each age), bx (one for each age), kt (one for each calendar year) - number of constraints (= 2 in LC model!)
LCfit$npar 

# plot estimated parameters

# plot LC model
plot(LCfit)
plot(LCfit, nCol = 3)
plot(LCfit, type = "p", pch = "*")

# * theory recap!: 
# ax -> time avg of log mortality rates at age x over set of calendar years we r considering -> shape str8 line bc on log scale, central death rates are approximately linear!
# bt -> modulates the effect kt term on different ages -> highest at ages 70-75! and lower for younger ages and older ages
# kt -> decreasing! -> dictates time behaviour of our central death rates. the time index of morality level! (goes down over time of mortality improves!) -> pushes whole central mortality down! good!

# plot cbd model
plot(CBDfit, parametricbx = FALSE)
# ??? bx? what bx?? the option parametricbx = FALSE excludes bx that are parametric functions: check
# excludes plot of bx bx bc bx is just constant
# k1t -> main driver!
# ??? k2t -> multiplied by linear function (x - xbar) -> modulates, allows to have different effect by age!!
# advantage of having two time indices k1t and k2t that change w mortality rates in cbd -> not perfectly correlated in cbd, good! but perfectlt correlated in LC, bad!
plot(CBDfit)

##################################
# Analyse residuals of the models#
##################################

# to analyse goodness of fit of model! how good is model, does it represent data? 
# is there any effect that is not captured by the model chosen?

# we will see that LC does better job than CBD which is too simple
# will also see that none of them is actually perfect
# will see how to forecast

# Compute deviance residuals (residuals = actual - fitted)
# we don't wanna see a pattern in these residuals! we only want random errors!
# ??? residuals are calculated according to underlying dist - Poisson dist!
LCres <- residuals(LCfit) 
CBDres <- residuals(CBDfit)

# three type of plots:
# 1) colour maps - set type = "colourmap" -> strength of colour reflects extent of the error! we don't want to see a pattern in these colours appearing!
# 1') sign plot, just white and grey - set type "signplot"
# 2) scatter plots by age, year and cohort (= year - age) - set type = "scatter"

par(mfrow = c(1, 1))

# colour map of residuals for LC
plot(LCres,
     type = "colourmap",
     reslim = c(-3.5, 3.5), # reslim helps painting a range of colours by setting the limit of residuals
     # ??? how to know what limits to choose??
     main = "Lee-Carter")
# some patterns visible:
# - the dark blue diagonal line appearing at age 60 in 1960 then going up diagonally
# - some patches of red/blue clustered tgt -> those areas always negative/positive
# so not entirely random! - fit is not too good!

# colour map of residuals for CBD
plot(CBDres,
     type = "colourmap",
     reslim = c(-3.5, 3.5),
     main = "Cairns-Blake-Dowd")
# strong pattern visible - way worse
# pattern kinda follows diagonal approach - the red area getting stronger and wider left to right, and going up
# -> suggests that it follows a given combination of ages and calendar years which are increasing together
# -> WHICH MEANS it's specific to some cohorts!!!!
# -> some group of individuals that share same period of time of birth, have sth in common, in particular, their mortality improvements have been faster than previous generations, and to the following generations... weird...

# diagonal patterns: combinations of ages and years -> cohort effect -> not explicitly captured in the Lee-Carter and the cbd models
# eg: this golden generation in uk is the one belonging to 1920s - 1940s (20-25 yrs)! they have lower mortality rates/higher mortality improvement rates so their mortality has been getting better than ppl who were born before 1920 and after 1945
# no clear explanation of why this effect is present! but insurers who are selling annuities hv to be aware of this and account for that!
# LC and CBD don't account for the cohort effect bc they separate age and calendar year, there is no cohort effect allowed! but there are other models that can be used, and the simplest way to do it is by adding a term that depends on (t - x), but forecasting is more complicated, we wil not be discussing that here

# scatter map of residuals for LC
plot(LCres,
     type = "scatter",
     cex = 0.3,
     reslim = c(-3.5, 3.5),
     main = "Lee-Carter")
# clear pattern by cohort (negative, then positive, then negative-ish again, then positive)

# scatter map of residuals for CBD
plot(CBDres,
     type = "scatter",
     cex = 0.3,
     reslim = c(-3.5, 3.5),
     main = "Cairns-Blake-Dowd")
# wayyyy worse than LC. strong pattern by age, year and cohort
# for calendar year it would've been okay if it weren't for the explosion at the end ard 2010

# AIC (Akaike Information Criterion) and BIC (Bayes Information Criterion) 
# metrics that account for goodness of fit (log-likelihood -> the higher the better) and model complexity (n. of parameters -> the lower the better u don't want it to be too complicated) -> but overall, there's a goodness of fit and model complexity trade-off
# -> the lower value of index, the better the model

# LC
AIC(LCfit) # the surface definition of AIC n BIC is just the log likelihood minus some function of the number of parameters
BIC(LCfit)

# CBD 
AIC(CBDfit)
BIC(CBDfit)

# clear winner is Lee-Carter. CBD index values twice as much for LC. LC more complicated, but it pays off!


###############
# Forecasting #
###############

# Forecast the models 20 years ahead: h = 20 (default option)
# we are forecasting kt
# period indices kt are forecasted using a MRW with drift (with correlated shocks)

# LC (only one period index kt)
LCfor <- forecast(LCfit)
# same as
LCfor <- forecast(LCfit, h = 20)
plot(LCfor, only.kt = TRUE) # can see that confidence bounds r getting wider n wider!
# if u want to forecast mortality rates next, just sub the forecasted kt in the original (ln ?) mxt eqn!

# confidence bounds are controlled by the argument 
# level = c(80, 95)

# access results of the forecast
LCfor$kt.f # kappa t index! period index (central estimate and lower / upper bounds)
LCfor$kt.f$mean # central estimate -> the black line straight line in the plot!
LCfor$kt.f$lower # ??? array corresponding to 20th and 5th percentiles
LCfor$kt.f$upper # array corresponding to 80th and 95th percentiles

# forecast central death rates mxt for 2023-2042 (20 years) and the ages 55-89
# -> can obtain death probs qxt and survival probs pxt
LCfor$rates
LCfor$rates[1:5, 1:5]

LCfor_qxt <- 1 - exp(- LCfor$rates) # ??? getting qxt from mxt like this, bc of the assumption that force of mortality is constant
LCfor_qxt[, "2042"] # life table (in terms of qx) for t=2042

# as a data frame -> life table!
data.frame(x = ages.fit,
           qx = as.numeric(LCfor_qxt[, "2042"]))

# cbd (two period indices k1t, k2t)
CBDfor <- forecast(CBDfit)
plot(CBDfor, only.kt = TRUE)

CBDfor$kt.f$mean # central estimate (a matrix for the two period indices)
CBDfor$kt.f$lower # array
CBDfor$kt.f$upper # array

# plot fitted and forecast rates at age 55
m55t <- mxt["55", as.character(years.fit)] # actual rates
m55t_LC <- LCfor$rates["55", ] # forecasted rates past 2022

m55t_CBD <- CBDfor$rates["55", ] # forecasted rates past 2022

plot(years.fit, m55t,
     xlim = range(years.fit, LCfor$years), # ??? 1960 - 2060
     # ylim = range(m55t, m55t_LC), # ???
     ylim = range(m55t, m55t_CBD),
     type = "p", xlab = "year", ylab = "rate",
     main = "Forecast of mortality rates at age 55 - Males")

lines(years.fit, mxt.hat.lc["55", ], col = "red") # fitted rates
lines(LCfor$years, m55t_LC, col = "red") # forecasted rates

lines(years.fit, mxt.hat.cbd["55", ], col = "blue") # fitted rates
lines(CBDfor$years, m55t_CBD, col = "blue") # forecasted rates

legend("topright", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

#########################################
# Temporary curtate expectation of life #
#########################################

# We want, for each calendar year t = 2023,...,2042:
#   e_{55,t:35} = sum_{k=1}^{35} {}_k p_{55,t},
# where {}_k p_{55,t} is the k-year survival probability for a life aged 55
# in year t, under a *period* life table using the forecast forces mu_{x,t}.
#
# Under constant force between integer ages:
#   p_{x,t} = exp(-mu_{x,t}),
# and {}_k p_{55,t} = product_{j=0}^{k-1} p_{55+j,t}.

# Forecast years we care about
years.forecast <- 2023:2042
age.range      <- 55:89  # 35 ages

# 1) Extract the forecast forces mu_{x,t} for ages 55–89, years 2023–2042

mu_LC  <- LCfor$rates[as.character(age.range), as.character(years.forecast)]
mu_CBD <- CBDfor$rates[as.character(age.range), as.character(years.forecast)]

# 2) Helper to compute e_{55,t:35} from a mu_{x,t} matrix
compute_e55_t35 <- function(mu_mat, x_start = 55, term = 35) {
  # mu_mat: rows = ages ("55",...,"89"), cols = forecast years ("2023",...,"2042")
  
  # 1-year survival probs p_{x,t}
  p_xt <- exp(-mu_mat)
  
  years <- colnames(p_xt)
  e_vec <- numeric(length(years))
  names(e_vec) <- years
  
  # index of age 55 in the rownames
  ages_num <- as.integer(rownames(p_xt))
  idx0 <- which(ages_num == x_start)
  
  for (j in seq_along(years)) {
    # p_55,t, p_56,t, ..., up to p_(55+term-1),t
    p_age <- p_xt[, j]
    
    # cumulative products => {}_k p_{55,t}, k = 1,...,term
    surv_k <- cumprod(p_age[idx0:(idx0 + term - 1)])
    
    # temporary curtate expectation = sum_k {}_k p_{55,t}
    e_vec[j] <- sum(surv_k)
  }
  
  e_vec
}

# 3) Compute e_{55,t:35} for LC and CBD (females, in this script)

e55_t35_LC  <- compute_e55_t35(mu_LC)
e55_t35_CBD <- compute_e55_t35(mu_CBD)

# Have a look
e55_t35_LC
e55_t35_CBD

# 4) Optional: put into a data frame for tables / plots
e55_t35_male <- data.frame(
  year        = as.numeric(names(e55_t35_LC)),
  e55_t35_LC  = as.numeric(e55_t35_LC),
  e55_t35_CBD = as.numeric(e55_t35_CBD)
)

print(e55_t35_male)

# Example plot – LC vs CBD temporary curtate life expectancy at age 55

plot(
  e55_t35_male$year, e55_t35_male$e55_t35_LC,
  col = "red", type = "l", lwd = 2,
  xlab = "Calendar year t",
  ylab = expression(e[55,t:35]),
  main = "Temporary curtate e[55,t:35] – Males"
)

lines(
  e55_t35_male$year, e55_t35_male$e55_t35_CBD,
  col = "blue", lwd = 2, lty = 2
)

legend("topleft", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
