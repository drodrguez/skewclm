# Cargar librerías necesarias
library(MASS)   # Para comparar con el modelo ordinal estándar
library(ordinal)
data(housing)
# Cargar nuestras funciones implementadas
source("R/distributions.R")
source("R/model.R")
source("R/estimation.R")
source("R/methods.R")

model <- clm_stn(Sat ~ Infl + Type, data = housing)
fit <- fit_clm_stn(model)
print(fit)
fit$gradient


model <- clm_stn(Sat ~ Infl + Cont, data = housing)
fit <- fit_clm_stn(model)
print(fit)
fit$gradient


model <- clm_stn(Sat ~ Type + Cont, data = housing)
fit <- fit_clm_stn(model)
print(fit)
fit$gradient
