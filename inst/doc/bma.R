## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StanMoMo)

## ----eval=FALSE---------------------------------------------------------------
#  #We extract deaths and exposures
#  ages.fit<-50:90
#  years.fit<-1970:2017
#  deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
#  exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
#  #We fit the three mortality models
#  fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 10, validation=10,family = "poisson",cores=4)
#  fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 10, validation=10,family = "poisson",cores=4)
#  fitAPC=apc_stan(death = deathFR,exposure=exposureFR, forecast = 10, validation=10,family = "poisson",cores=4)
#  #We compute the model weights
#  model_weights<-mortality_weights(list(fitLC,fitRH,fitAPC))

