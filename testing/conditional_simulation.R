library(tidyverse)
library(crawl)
library(sf)

ak_base = nPacMaps::alaska(resolution = "f", simplify = FALSE)

data("harborSeal")

levels(harborSeal$Argos_loc_class) = c("3","2","1","0","A","B")

harborSeal_sf <- harborSeal %>% 
  filter(!is.na(latitude)) %>%
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  sf::st_transform(3338)

#######################################################
# fit crawl model
######################################################

initial = list(a = c(sf::st_coordinates(harborSeal_sf)[[1,1]], 0,
                     sf::st_coordinates(harborSeal_sf)[[1,2]], 0),
  P=diag(c(10000^2,5400^2,10000^2,5400^2))
)

fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)

constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal_sf, Time.name="Time", 
  initial.state=initial, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)#,
  #initialSANN=list(maxit=200, trace=1, REPORT=1)
)

predTimes <- seq(min(harborSeal_sf$Time), max(harborSeal_sf$Time), by = 0.5)


hs_pred <- crawl::crwPredict(fit1, predTime=predTimes)
hs_sim <- crawl::crwSimulator(fit1, predTime = predTimes)
hs_crwIS <- crwPostIS(hs_sim)
