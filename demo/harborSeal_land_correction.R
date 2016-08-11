
library(dplyr)
library(sp)
library(crawl)
library(raster)
library(gdistance)

data("harborSeal")

harborSeal %>% dplyr::filter(!is.na(latitude)) %>% 
  as.data.frame() -> harborSeal

sp::coordinates(harborSeal) <- ~longitude + latitude
sp::proj4string(harborSeal) <- CRS("+init=epsg:4326")

harborSeal <- sp::spTransform(harborSeal, CRS("+init=epsg:3338"))
harborSeal <- harborSeal[order(harborSeal$Time),]

# fit crawl model

initial = list(
  a=c(coordinates(harborSeal)[1,1],0,
      coordinates(harborSeal)[1,2],0),
  P=diag(c(10000^2,5400^2,10000^2,5400^2))
)

fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)

constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("x","y"), Time.name="Time", 
  initial.state=initial, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)#,
  #initialSANN=list(maxit=200, trace=1, REPORT=1)
)

predTimes <- seq(min(harborSeal$Time), max(harborSeal$Time), by = 0.5)

# Created predicted path
hs_pred <- crawl::crwPredict(fit1, predTime=predTimes)
# coordinates(hs_pred) <- ~mu.x+mu.y
# proj4string(hs_pred) <- CRS("+init=epsg:3338")

# Create simulated posterior track
# simObj <- crwSimulator(fit1, predTimes, parIS=0)
# samp = crwPostIS(simObj, fullPost = FALSE)

# Get land raster
land = raster(system.file("raster/land_res.grd", package="crawl"))
plot(land)
lines(hs_pred$mu.x, hs_pred$mu.y)

# Convert to water areas
water = asFactor(1-land)
#Create transition matrix
trans = transition(water, "areas", directions = 16)[[1]]
# Move path based on shortest distance around land
new_hs_pred = fix_path(hs_pred, land, trans)
plot(land)
lines(new_hs_pred, col="red")

# new_samp = fix_path(samp, land, trans)
# plot(land)
# lines(new_samp, col="red")
