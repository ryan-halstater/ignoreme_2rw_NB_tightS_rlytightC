## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev = "png",dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center")

## ----setup--------------------------------------------------------------------
library(StanMoMo)

## -----------------------------------------------------------------------------
ages.fit<-60:90
years.fit<-1980:2010
deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]

## -----------------------------------------------------------------------------
fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",chains=1,iter=1000,cores=1)

## ----warning=FALSE------------------------------------------------------------
#Extract model parameters from the fitted model
params<-rstan::extract(fitLC)
print(names(params))

## ----eval=FALSE---------------------------------------------------------------
#  library(shinystan)
#  launch_shinystan(fitLC)

## ----echo=TRUE,warning=FALSE--------------------------------------------------
boxplot_post_dist(fitLC, "a", ages.fit, years.fit)
boxplot_post_dist(fitLC, "b", ages.fit, years.fit)
boxplot_post_dist(fitLC, "k", ages.fit, years.fit)

## ----eval=FALSE---------------------------------------------------------------
#  fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",cores=4)
#  fitAPC=apc_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",cores=4)
#  fitCBD=cbd_stan(death = deathFR,exposure=exposureFR, age=ages.fit, forecast=10,family = "poisson",cores=4)
#  fitM6=m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit,forecast = 10, family = "poisson",cores=4)

