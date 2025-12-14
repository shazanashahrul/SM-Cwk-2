#######################################
# AS3204 Survival Models Coursework 2 #
#############################################################
# Mortality forecasting FEMALES Portugal mortality data with StMoMo #
# Group 5 #
###########



#################
# Preliminaries #
#################

# set directory: Session -> Set Working Directory -> To Source File Location

# load the package
install.packages("StMoMo")
library(StMoMo)

# load the data set for Portugal (females)
load(file = "PORTUGAL-females-mortality.RData")

# "PORTUGALStMoMof": StMoModata object
# these are lists with several components
print(PORTUGALStMoMof)

# matrices with death counts and expsoures, ages on rows, years on columns
PORTUGALStMoMof$Dxt[1:5, 1:5]
PORTUGALStMoMof$Ext[1:5, 1:5]

# ages, years, type, series, label
PORTUGALStMoMof$ages # 0 - 110
PORTUGALStMoMof$years # 1940 - 2023
PORTUGALStMoMof$type # central
PORTUGALStMoMof$series # female
PORTUGALStMoMof$label # portugal



#######################################
# Define a stochastic mortality model #
#######################################

# Lee-Carter (LC)

LC <- lc() # default settings same as below

print(LC)

# or, more in detail
LC <- lc(link = "log", const = "sum")
# link = "log" -> log(mxt) -> Poisson modelling of death counts
# recall: model is based in death count (Dxt), and Dxt follows a Poisson distrubtion with param = mxt * Ext like below
# Dxt ~ Poisson(mxt * Ext) with log(mxt) = ax+bx*kt 
# mxt -> central death rate, Ext -> central exposure 
# by LC; log(mxt) = ax+bx*kt -> idea behind LC is separating effect of age and effect of time -> two separate effects
# const = "sum" -> standard constraints on parameters to make model identifiable: sum b1[x] = 1, sum k1[t] = 0 -> we don't want two different sets of values of parameters leading to same value of log(mxt)
# other models will have a similar structure, just different expressions for log(mxt) (see CBD later)


# Cairns-Blake-Dowd (cbd, or M5)

CBD <- cbd(link = "log")

# link = "log" -> Poisson modelling of death counts -> this is now necessary as by default the cbd model works with a Binomial distribution and logit function of one year intitial rates instead of the central mortality

# Dxt ~ Poisson(mxt * Ext) with log(mxt) = k1t + (x-xbar) * k2t

print(CBD)

# no constraint needed on parameters; here f2[x] = x - xbar, where xbar is the average age in the data set



####################################
# Fit a stochastic mortality model #
####################################

# set the ages and years to use in fitting (we want to skip periods of time like world war and plagues etc)
# ages 55 to 89 (35 ages), years 1963 to 2022 (60 years)
# must be a subset of the available ages and years
ages.fit <- 55:89
years.fit <- 1963:2022

# unlike Lee-Carter, the CBD model can only be fitted on adult mortality data -> CBD has a linear age function (x- xbar), so it doesn't address instances like infant mortality and accidental hump

# fit Lee-Carter with ML
# fit() uses a GNM (generalised non-linear model) to fit the data
LCfit <- fit(object = LC, # model
             data = PORTUGALStMoMof, # a StMoModata object
             ages.fit = ages.fit, # ages in fit
             years.fit = years.fit # years in fit
)

# plot actual vs fitted rates
# actual central death rates = death counts / exposure
mxt <- PORTUGALStMoMof$Dxt / PORTUGALStMoMof$Ext

# LC fitted rates
mxt.hat.lc <- fitted(LCfit, type = "rates")

# fit Cairns-Blake-Dowd with ML

# CBDfit <- fit(object = CBD,
#               data = PORTUGALStMoMof,
#               ages.fit = ages.fit,
#               years.fit = years.fit
# )

# to avoid issues with non integer death counts (since poisson distribution is discrete), specify one by one
CBDfit <- fit(object = CBD,
              Dxt = round(PORTUGALStMoMof$Dxt),
              Ext = PORTUGALStMoMof$Ext,
              ages = PORTUGALStMoMof$ages,
              years = PORTUGALStMoMof$years,
              ages.fit = ages.fit,
              years.fit = years.fit
)


# CBD fitted rates
mxt.hat.cbd <- fitted(CBDfit, type = "rates")

# plot actual rates
plot(years.fit,
     mxt["55", as.character(years.fit)],
     ylim = range(mxt["55", as.character(years.fit)], mxt.hat.lc["55", ], mxt.hat.cbd["55", ]), # makes sure everything is in frame
     xlab = "year",
     ylab = "central death rate",
     main = "Fitted vs. observed rates at age 55 - Females"
)

# plot LC fitted rates
lines(years.fit, mxt.hat.lc["55", ], col = "red")

# overlay CBD fitted rates
lines(years.fit, mxt.hat.cbd["55", ], col = "blue")

# add legend
legend("topright", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:1, cex=0.8)



###############################
# Analyse the estimated model #
###############################

# convergence:
LCfit$conv 
CBDfit$conv
# algorithm behind $conv is just MLE, and returns if convergence testing via numerical optimisation methods was successful
# LCfit and CBDfit contain all the data - original data, fitted parameters, information about fitting process, etc

# output: for instance, for LC
# extract estimated parameters using coef()
coef(LCfit)
coef(LCfit)$ax # vector
coef(LCfit)$bx # matrix
coef(LCfit)$kt # matrix

LCfit$loglik # log-likelihood

# sample size = length(ages.fit) * length(years.fit)
LCfit$nobs # number of observations (aka sample size)

# dof = length(ages.fit) + length(ages.fit) + length(years.fit) - 2
# number of ax (one for each age), bx (one for each age), kt (one for each calendar year) - number of constraints (= 2 in LC model)
LCfit$npar 

# plot estimated parameters

# plot LC model
plot(LCfit, nCol = 3)

# theory recap: 
# ax -> time average of log mortality rates at age x over set of calendar years we are considering -> shape is straight line because on log scale, central death rates are approximately linear
# bt -> modulates the effect of kt term on different ages -> highest at late 60s/early 70s, and lower for younger ages and older ages
# kt -> dictates time behaviour of our central death rates. This is the time index of morality level (goes down over time if mortality improves) -> here it is decreasing -> pushes whole log of central mortality down, which is good

# plot cbd model
plot(CBDfit, parametricbx = FALSE)
# the option parametricbx = FALSE excludes bx that are parametric functions. Excludes plot of bx because bx is just constant
# k1t -> main driver
# k2t -> multiplied by linear function (x - xbar) -> allows to have different effect by age
# advantage is having two time indices k1t and k2t that change with mortality rates -> not perfectly correlated in CBD, which is good, but perfectly correlated in LC, which is bad



###################################
# Analyse residuals of the models #
###################################

# to analyse goodness of fit of model. How good is the model, does it represent data? 
# is there any effect that is not captured by the model chosen?

# we will see that LC does a better job than CBD which is too simple

# Compute deviance residuals (residuals = actual - fitted)
# we do not want to see a pattern in these residuals. We only want random errors
# residuals are calculated according to underlying distribution (Poisson distribution)
LCres <- residuals(LCfit) 
CBDres <- residuals(CBDfit)

# we will be using two types of plots:
# 1) colour maps -> set type = "colourmap" -> strength of colour reflects extent of the error. We do not want to see a pattern in these colours appearing
# 2) scatter plots by age, year and cohort (= year - age) -> set type = "scatter"

par(mfrow = c(1, 1))

# colour map of residuals for LC
plot(LCres,
     type = "colourmap",
     reslim = c(-3.5, 3.5), # reslim helps painting a range of colours by setting the limit of residuals
     main = "LC Residuals Heatmap - Females")
# some patterns visible:
# - some diagonal lines appearing at age 60 in 1960 then going up diagonally
# - some patches of red/blue clustered together
# so not entirely random - fit is not too good

# colour map of residuals for CBD
plot(CBDres,
     type = "colourmap",
     reslim = c(-3.5, 3.5),
     main = "CBD Residuals Heatmap - Females")
# strong pattern visible - much worse
# pattern follows a slightly diagonal approach - the red area is getting stronger and wider left to right, and going up
# -> suggests that it follows a given combination of ages and calendar years which are increasing together
# -> which means it is specific to some cohorts
# -> some group of individuals that share same period of time of birth, have something in common, in particular, their mortality improvements have been faster than previous generations, as well as the following generations

# diagonal patterns: combinations of ages and years -> cohort effect -> is not explicitly captured in the Lee-Carter and the CBD models
# LC and CBD do not account for the cohort effect since they separate age and calendar year, there is no cohort effect allowed. But there are other models that can be used, and the simplest way to do it is by adding a term that depends on (t - x), but forecasting is more complicated, so we will not be carrying this procedure out,

# scatter map of residuals for LC
plot(LCres,
     type = "scatter",
     cex = 0.3,
     reslim = c(-3.5, 3.5),
     main = "Lee-Carter - F")
# fairly random scatterplot

# scatter map of residuals for CBD
plot(CBDres,
     type = "scatter",
     cex = 0.3,
     reslim = c(-3.5, 3.5),
     main = "Cairns-Blake-Dowd - F")
# much worse than LC. strong pattern by age, year and cohort

# AIC (Akaike Information Criterion) and BIC (Bayes Information Criterion) 
# metrics that account for goodness of fit (log-likelihood) and model complexity (n. of parameters) -> there is a goodness of fit and model complexity trade-off
# -> the lower value of index, the better the model

# LC
AIC(LCfit)
BIC(LCfit)

# CBD 
AIC(CBDfit)
BIC(CBDfit)



###############
# Forecasting #
###############

# Forecast the models 20 years ahead: h = 20
# Firstly, forecasting of kt
# period indices kt are forecasted using a MRW with drift (with correlated shocks)

# LC (only one period index kt)
LCfor <- forecast(LCfit, h = 20)
plot(LCfor, only.kt = TRUE) 
# we can see that confidence bounds are getting wider n wider

# confidence bounds are controlled by the argument 
# level = c(80, 95)

# access results of the forecast
LCfor$kt.f # kappa t index
LCfor$kt.f$mean # central estimate -> the black straight line in the plot
LCfor$kt.f$lower # array corresponding to 20th and 5th percentiles
LCfor$kt.f$upper # array corresponding to 80th and 95th percentiles

# forecast central death rates mxt for 2023-2042 (20 years) and the ages 55-89
# -> can obtain death probability qxt and survival probability pxt
LCfor$rates[1:5, 1:5]

LCfor_qxt <- 1 - exp(- LCfor$rates) # derive qxt from mxt this way, because of the assumption that force of mortality is constant
LCfor_qxt[, "2042"] # life table (in terms of qx) for t = 2042

# as a data frame -> qxt displayed in life table form
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
m55t_LC <- LCfor$rates["55", ] # forecasted LC rates past 2022
m55t_CBD <- CBDfor$rates["55", ] # forecasted CBD rates past 2022

plot(years.fit, m55t,
     xlim = range(years.fit, LCfor$years), 
     ylim = range(m55t, m55t_CBD),
     type = "p", xlab = "year", ylab = "rate",
     main = "Forecast of mortality rates at age 55 - Females")

lines(years.fit, mxt.hat.lc["55", ], col = "red") # fitted LC rates
lines(LCfor$years, m55t_LC, col = "red") # forecasted LC rates

lines(years.fit, mxt.hat.cbd["55", ], col = "blue") # fitted CBD rates
lines(CBDfor$years, m55t_CBD, col = "blue") # forecasted CBD rates

legend("topright", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:1, cex=0.8)



#########################################
# Temporary curtate expectation of life #
#########################################

# (formulas used in comments below are in LaTeX format)

# we want, for each calendar year t = 2023,...,2042:
#   e_{55,t:35} = sum_{k=1}^{35} {}_k p_{55,t},
# where {}_k p_{55,t} is the k-year survival probability for a life aged 55
# in year t, under a period life table using the forecast forces mu_{x,t}.

# under constant force between integer ages:
#   p_{x,t} = exp(-mu_{x,t}),
# and {}_k p_{55,t} = product_{j=0}^{k-1} p_{55+j,t}.

# forecast years we want
years.forecast <- 2023:2042
age.range <- 55:89  # 35 ages

# 1) Extract the forecast forces mu_{x,t} for ages 55–89, years 2023–2042
mu_LC <- LCfor$rates[as.character(age.range), as.character(years.forecast)]
mu_CBD <- CBDfor$rates[as.character(age.range), as.character(years.forecast)]

# 2) Helper to compute e_{55,t:35} from a mu_{x,t} matrix
compute_e55_t35 <- function(mu_mat, x_start = 55, term = 35) {
  # mu_mat:
  #   - rows correspond to integer ages (55, 56, ..., 89)
  #   - columns correspond to calendar years (2023, ..., 2042)
  #
  # x_start:
  #   - starting age (55 in coursework)
  #
  # term:
  #   - number of whole years we consider (35 in coursework)
  #
  # OUTPUT:
  #   - named numeric vector:
  #       names = calendar years (e.g. "2023", ..., "2042")
  #       values = e_{55,t:35} for each year t
  
  # convert mu to 1-year survival probabilities p_{x,t} = exp(-mu_{x,t})
  p_xt <- exp(-mu_mat)
  
  # get the calendar years from the column names
  years <- colnames(p_xt)
  
  # prepare a vector to store e_{55,t:35} for each year
  e_vec <- numeric(length(years))
  names(e_vec) <- years
  
  # find the row index corresponding to age 55
  # (we assume rownames are "55", "56", ..., "89")
  ages_num <- as.integer(rownames(p_xt)) # turn "55" -> 55, etc.
  idx0 <- which(ages_num == x_start)
  
  # loop over each forecast year t
  for (j in seq_along(years)) {
    # extract the column of survival probabilities for that year:
    # p_55,t, p_56,t, ..., p_89,t
    p_age <- p_xt[, j]
    
    # we only need the first 'term' ages starting at 55:
    # (55, 56, ..., 55 + term - 1)
    # here term = 35, so this is ages 55,...,89.
    p_55_to_89 <- p_age[idx0:(idx0 + term - 1)]
    
    # cumulative product:
    #   surv_k[1] = p_55,t = {}_1 p_{55,t}
    #   surv_k[2] = p_55,t * p_56,t = {}_2 p_{55,t}
    #   ...
    #   surv_k[k] = product_{j=0..k-1} p_{55+j,t} = {}_k p_{55,t}
    surv_k <- cumprod(p_55_to_89)
    
    # temporary curtate expectation:
    #   e_{55,t:35} = sum_{k=1}^{35} {}_k p_{55,t}
    e_vec[j] <- sum(surv_k)
  }
  
  # return the vector of e_{55,t:35} across years
  e_vec
}

# 3) compute e_{55,t:35} for LC and CBD (females, in this script)
e55_t35_LC  <- compute_e55_t35(mu_LC)
e55_t35_CBD <- compute_e55_t35(mu_CBD)

# have a look
e55_t35_LC
e55_t35_CBD

# 4) put into a data frame for tables / plots
e55_t35_female <- data.frame(
  year        = as.numeric(names(e55_t35_LC)),
  e55_t35_LC  = as.numeric(e55_t35_LC),
  e55_t35_CBD = as.numeric(e55_t35_CBD)
)

print(e55_t35_female)

# example plot – LC vs CBD temporary curtate life expectancy at age 55

plot(
  e55_t35_female$year, e55_t35_female$e55_t35_LC,
  col = "red", type = "l", lwd = 2,
  xlab = "Calendar year t",
  ylab = expression(e[55,t:35]),
  main = "Temporary curtate e[55,t:35] – Females"
)

lines(
  e55_t35_female$year, e55_t35_female$e55_t35_CBD,
  col = "blue", lwd = 2, lty = 2
)

legend("topleft", legend=c("LC", "CBD"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
