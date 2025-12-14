##########################################
# AS3204 Coursework 2 - StMoMo Template  #
# Country: PORTUGAL (example)           #
##########################################

# Install packages if needed:
# install.packages("StMoMo")

library(StMoMo)

################################
# 1. Load and inspect the data #
################################

# Set working directory to where this script and RData files live:
# Session -> Set Working Directory -> To Source File Location

# Load male and female mortality data (StMoMoData objects)
load("PORTUGAL-males-mortality.RData")
load("PORTUGAL-females-mortality.RData")

# After loading, check what objects you got:
ls()

# *** IMPORTANT ***
# Adjust these names to match what you see from ls().
# I assume below they’re called PORTUGAL.m and PORTUGAL.f,
# similar to SPAIN.m / SPAIN.f in the lecturer’s script.

PORT.m <- PORTUGAL.m   # rename if needed
PORT.f <- PORTUGAL.f   # rename if needed

# Quick sanity check:
print(PORT.m)
print(PORT.f)

# Check ages and years ranges:
PORT.m$ages
PORT.m$years

#########################################
# 2. Define ages/years and StMoMo models #
#########################################

# Coursework spec: ages 55–89, years 1963–2022
ages.fit  <- 55:89
years.fit <- 1963:2022

# Lee–Carter with Poisson/log, standard constraints
LC  <- lc(link = "log", const = "sum")

# CBD with Poisson/log (override default logit/Binomial)
CBD <- cbd(link = "log")

LC
CBD

##############################################
# 3. Fit LC and CBD – males and females     #
##############################################

## 3.1 LC fits

LCfit.m <- fit(
  object    = LC,
  data      = PORT.m,
  ages.fit  = ages.fit,
  years.fit = years.fit
)

LCfit.f <- fit(
  object    = LC,
  data      = PORT.f,
  ages.fit  = ages.fit,
  years.fit = years.fit
)

## 3.2 CBD fits (round Dxt to avoid non-integer issues)

CBDfit.m <- fit(
  object    = CBD,
  Dxt       = round(PORT.m$Dxt),
  Ext       = PORT.m$Ext,
  ages      = PORT.m$ages,
  years     = PORT.m$years,
  ages.fit  = ages.fit,
  years.fit = years.fit
)

CBDfit.f <- fit(
  object    = CBD,
  Dxt       = round(PORT.f$Dxt),
  Ext       = PORT.f$Ext,
  ages      = PORT.f$ages,
  years     = PORT.f$years,
  ages.fit  = ages.fit,
  years.fit = years.fit
)

# Check convergence flags
LCfit.m$conv; LCfit.f$conv
CBDfit.m$conv; CBDfit.f$conv

############################################
# 4. Basic parameter plots (optional)      #
############################################

# These plots are handy for the report
par(mfrow = c(2, 2))

plot(LCfit.m, main = "LC – Males")
plot(LCfit.f, main = "LC – Females")
plot(CBDfit.m, parametricbx = FALSE, main = "CBD – Males")
plot(CBDfit.f, parametricbx = FALSE, main = "CBD – Females")

par(mfrow = c(1, 1))

############################################
# 5. Residuals and heatmaps (part iii a)   #
############################################

LCres.m  <- residuals(LCfit.m)
LCres.f  <- residuals(LCfit.f)
CBDres.m <- residuals(CBDfit.m)
CBDres.f <- residuals(CBDfit.f)

# Heatmaps of deviance residuals
par(mfrow = c(2, 2))

plot(LCres.m,
     type   = "colourmap",
     reslim = c(-3.5, 3.5),
     main   = "LC residuals – Males")

plot(LCres.f,
     type   = "colourmap",
     reslim = c(-3.5, 3.5),
     main   = "LC residuals – Females")

plot(CBDres.m,
     type   = "colourmap",
     reslim = c(-3.5, 3.5),
     main   = "CBD residuals – Males")

plot(CBDres.f,
     type   = "colourmap",
     reslim = c(-3.5, 3.5),
     main   = "CBD residuals – Females")

par(mfrow = c(1, 1))

# You can also use type = "scatter" for cohort stripes if you want.

###########################################
# 6. AIC and BIC (used in report)        #
###########################################

AIC(LCfit.m); BIC(LCfit.m)
AIC(LCfit.f); BIC(LCfit.f)

AIC(CBDfit.m); BIC(CBDfit.m)
AIC(CBDfit.f); BIC(CBDfit.f)

###########################################
# 7. Forecast kappas and mortality (iii b)#
###########################################

# Forecast horizon: 2023–2042 from last observed year 2022
h <- 2042 - max(years.fit)  # should be 20

LCfor.m  <- forecast(LCfit.m,  h = h)
LCfor.f  <- forecast(LCfit.f,  h = h)
CBDfor.m <- forecast(CBDfit.m, h = h)
CBDfor.f <- forecast(CBDfit.f, h = h)

# Plot kappa paths only
par(mfrow = c(2, 2))
plot(LCfor.m,  only.kt = TRUE, main = "LC κ_t – Males")
plot(LCfor.f,  only.kt = TRUE, main = "LC κ_t – Females")
plot(CBDfor.m, only.kt = TRUE, main = "CBD κ_t – Males")
plot(CBDfor.f, only.kt = TRUE, main = "CBD κ_t – Females")
par(mfrow = c(1, 1))

# Access forecast rates (μ or m_x,t) for ages 55–89, years 2023–2042
LC_rates.m  <- LCfor.m$rates
LC_rates.f  <- LCfor.f$rates
CBD_rates.m <- CBDfor.m$rates
CBD_rates.f <- CBDfor.f$rates

# Sanity check dimensions and names
dim(LC_rates.m)
rownames(LC_rates.m)
colnames(LC_rates.m)  # should be "2023", ..., "2042"

#############################################
# 8. Helper functions for e^{:35}_{55,t}    #
#############################################

# We treat LCfor$rates / CBDfor$rates as (approx) forces μ_{x,t}
# under constant force of mortality between ages x and x+1.
# Then 1-year survival at age x in year t is:
#   p_{x,t} = exp(-μ_{x,t}).  (Unit 6 notes, CFM => μ ≈ m_x,t). 
#
# For each calendar year t, we build a life table for ages 55–89
# using {μ_{55,t}, ..., μ_{89,t}} and compute:
#   e^{:35}_{55,t} = Σ_{k=1}^{35} {}_k p_{55,t}
# where {}_k p_{55,t} is the product of 1-year survival probs
# from age 55 up to age 55 + k - 1 (STAYING in calendar year t
# for the purposes of this life table).

compute_e55_35 <- function(rate_matrix) {
  # rate_matrix: matrix with rows named "55"..."89",
  # columns = forecast years "2023"... "2042"
  
  # 1-year survival probabilities
  pxt <- exp(-rate_matrix)
  
  years <- colnames(rate_matrix)
  e_vec <- numeric(length(years))
  names(e_vec) <- years
  
  for (j in seq_along(years)) {
    # p_55,t ... p_89,t for that calendar year
    p_age <- pxt[, j]
    # Cumulative products => {}_k p_{55,t}
    surv_k <- cumprod(p_age)  # length = 35
    e_vec[j] <- sum(surv_k)
  }
  
  return(e_vec)
}

# Compute e^{:35}_{55,t} for each combination (iii c)

e55_35_LC_m  <- compute_e55_35(LC_rates.m)
e55_35_LC_f  <- compute_e55_35(LC_rates.f)
e55_35_CBD_m <- compute_e55_35(CBD_rates.m)
e55_35_CBD_f <- compute_e55_35(CBD_rates.f)

# Quick plots for intuition
plot(as.numeric(names(e55_35_LC_m)),  e55_35_LC_m,  type = "l",
     xlab = "Calendar year t", ylab = "e^{:35}_{55,t}",
     main = "Temporary curtate e^{:35}_{55,t} – LC males")

lines(as.numeric(names(e55_35_LC_f)),  e55_35_LC_f)
# add legend, etc., as you like

# You can repeat with CBD, or overlay LC vs CBD, and males vs females.

#############################################
# 9. Save results for report (tables/plots) #
#############################################

# For example, put life expectancies in a data.frame:
e55_35_df <- data.frame(
  year         = as.numeric(names(e55_35_LC_m)),
  LC_males     = e55_35_LC_m,
  LC_females   = e55_35_LC_f,
  CBD_males    = e55_35_CBD_m,
  CBD_females  = e55_35_CBD_f
)

print(e55_35_df)

# You can write this to CSV for use in plots / LaTeX / Word:
# write.csv(e55_35_df, "e55_35_results.csv", row.names = FALSE)
