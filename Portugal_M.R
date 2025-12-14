load("C:/Users/maanu/Downloads/PORTUGAL-males-mortality.RData")
load("C:/Users/maanu/Downloads/PORTUGAL-females-mortality.RData")
ls()
str(PORTUGALStMoMof)
str(PORTUGALStMoMom)
library(StMoMo)
install.packages("StMoMo")
install.packages("StMoMo")
library(StMoMo)

# Define Lee–Carter model: log-Poisson (as required)
LC <- lc(link = "log")   # or just LC <- lc() since "log" is the default
ages_fit  <- 55:89
years_fit <- 1963:2022

# Females
LC_females <- fit(
  LC,
  data      = PORTUGALStMoMof,
  ages.fit  = ages_fit,
  years.fit = years_fit
)

# Males
LC_males <- fit(
  LC,
  data      = PORTUGALStMoMom,
  ages.fit  = ages_fit,
  years.fit = years_fit
)

summary(LC_females)
summary(LC_males)
plot(LC_females, param = "residuals")
plot(LC_males,   param = "residuals")


