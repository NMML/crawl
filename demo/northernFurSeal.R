data(northernFurSeal)

northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class, 
                                          levels=c("3", "2", "1", "0", "A"))
## Project data ##
library(rgdal)
coordinates(northernFurSeal) = ~longitude+latitude
proj4string(northernFurSeal) <- CRS("+proj=longlat")
northernFurSeal <- spTransform(northernFurSeal, CRS("+proj=aea +lat_1=30 +lat_2=70 +lat_0=52 +lon_0=-170 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

initial = list(
  a=c(coordinates(northernFurSeal)[1,1],0,coordinates(northernFurSeal)[1,2],0),
  P=diag(c(10000^2,54000^2,10000^2,5400^2))
)

##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)
fixPar = c(log(250), log(500), log(1500), rep(NA,3), NA)
displayPar( mov.model=~1, err.model=list(x=~Argos_loc_class-1),
            data=northernFurSeal,fixPar=fixPar)
constr=list(
  lower=c(rep(log(1500),2), rep(-Inf,2)),
  upper=rep(Inf,4)
)

ln.prior = function(theta){-abs(theta[4])/0.5}

#set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), 
  data=northernFurSeal, Time.name="Time", 
  initial.state=initial, fixPar=fixPar, constr=constr, prior=ln.prior,
  control=list(maxit=30, trace=1, REPORT=1),
  initialSANN=list(maxit=200, trace=1, REPORT=1)
)
fit1
##Make hourly location predictions
predTime <- seq(ceiling(min(northernFurSeal$Time)), floor(max(northernFurSeal$Time)), 1)
predObj <- crwPredict(object.crwFit=fit1, predTime, speedEst=TRUE, flat=TRUE)
head(predObj)
crwPredictPlot(predObj, "map")

##Create simulation object with 100 parameter draws
set.seed(123)
simObj <- crwSimulator(fit1, predTime, parIS=100, df=Inf, scale=18/20)

## Examine IS weight distribution
w <- simObj$thetaSampList[[1]][,1]
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')

##Approximate number of independent samples
round(100/(1+(sd(w)/mean(w))^2))

dev.new(bg=gray(0.75))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
crwPredictPlot(predObj, 'map')

## Sample 20 tracks from posterior predictive distribution
iter <- 20
cols <- jet.colors(iter)
for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim[,'mu.x'], samp$alpha.sim[,'mu.y'],col=cols[i])
}

