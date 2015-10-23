## ----message=FALSE-------------------------------------------------------
library(crawl)

## ------------------------------------------------------------------------
data("northernFurSeal")
head(northernFurSeal)

## ------------------------------------------------------------------------
northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,levels=c("3", "2", "1","0","A"))

## ---- message=FALSE------------------------------------------------------
library(rgdal)
coordinates(northernFurSeal) = ~longitude+latitude

## ------------------------------------------------------------------------
proj4string(northernFurSeal) <- CRS("+proj=longlat")

## ----message=FALSE-------------------------------------------------------
northernFurSeal <- spTransform(northernFurSeal, CRS("+proj=aea +lat_1=30 +lat_2=70 +lat_0=52 +lon_0=-170 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

## ----message=FALSE-------------------------------------------------------
initial = list(a=c(coordinates(northernFurSeal)[1,1],0,coordinates(northernFurSeal)[1,2],0),P=diag(c(10000^2,54000^2,10000^2,5400^2)))

## ----message=FALSE-------------------------------------------------------
fixPar = c(log(250), log(500), log(1500), rep(NA,3), NA)

## ----message=FALSE-------------------------------------------------------
displayPar(mov.model=~1,err.model=list(x=~Argos_loc_class-1),data=northernFurSeal,fixPar=fixPar)

## ---- message=FALSE------------------------------------------------------
constr=list(lower=c(rep(log(1500),2), rep(-Inf,2)),upper=rep(Inf,4))

## ---- message=FALSE------------------------------------------------------
ln.prior = function(theta){-abs(theta[4]-3)/0.5}

## ----message=FALSE-------------------------------------------------------
set.seed(123)
fit1 <- crwMLE(mov.model=~1, 
               err.model=list(x=~Argos_loc_class-1),
               data=northernFurSeal, 
               Time.name="Time",
               initial.state=initial,fixPar=fixPar, constr=constr, prior=ln.prior,
               control=list(maxit=30, trace=0,REPORT=1),
               initialSANN=list(maxit=200, trace=0, REPORT=1))

## ------------------------------------------------------------------------
fit1

## ----message=FALSE-------------------------------------------------------
predTime <- seq(ceiling(min(northernFurSeal$Time)), floor(max(northernFurSeal$Time)), 1)

## ----message=FALSE-------------------------------------------------------
predObj <- crwPredict(object.crwFit=fit1, predTime, speedEst=TRUE, flat=TRUE)

## ----message=FALSE-------------------------------------------------------
crwPredictPlot(predObj, "map")

## ----message=FALSE-------------------------------------------------------
set.seed(123)
simObj <- crwSimulator(fit1, predTime, method="IS", parIS=100, df=5, scale=18/20)

## ----message=FALSE-------------------------------------------------------
w <- simObj$thetaSampList[[1]][,1]
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')

## ----message=FALSE-------------------------------------------------------
round(100/(1+(sd(w)/mean(w))^2))

## ------------------------------------------------------------------------
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## ------------------------------------------------------------------------
iter <- 20
cols <- jet.colors(iter)

## ----message=FALSE-------------------------------------------------------
crwPredictPlot(predObj, 'map')
for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim[,'mu.x'], samp$alpha.sim[,'mu.y'],col=cols[i]) 
  }

## ----message=FALSE-------------------------------------------------------
# library(kotzeb0912)
library(sp)
library(rgdal)
library(trip)
library(crawl)
library(argosfilter)
library(parallel)
library(dplyr)
library(lubridate)

## ----message=FALSE-------------------------------------------------------
# locs <- dplyr::filter(kotzeb0912_locs,instr=="Mk10") %>% 
#   dplyr::arrange(deployid,unique_posix)
# # speedfilter using the paralell package for multi-core speed
# cfilter<-mclapply(split(locs,locs$deployid),function(x) sdafilter(
#   lat=x$latitude, lon=x$longitude, dtime=x$unique_posix,
#   lc=x$quality, ang=-1,vmax=3.5),mc.preschedule=F,mc.cores=3)
# cfilter<-do.call("c",cfilter)
# cfilter<-as.vector(cfilter)
# locs$filtered <- cfilter
# 
# data <- filter(locs,filtered=="not", !is.na(error_semimajor_axis)) %>% arrange(deployid,unique_posix) %>% as.data.frame(.)


## ----message=FALSE-------------------------------------------------------
# coordinates(data) = ~longitude+latitude
# proj4string(data) = CRS("+proj=longlat +datum=WGS84")
# 
# data=spTransform(data, CRS("+init=epsg:3571"))

## ----message=FALSE-------------------------------------------------------
# ids = unique(data@data$deployid)              #define seal IDs
# model.fits = vector("list",length(ids))       #make an list of the seal IDs
# names(model.fits) = ids                       #define the names of each model as the seal ID
# 
# for(i in 1:length(ids)){                      
#   idData = subset(data,deployid==ids[i])      #subset the data using the first PTT (seal ID)
#   diagData = model.matrix(~error_semimajor_axis+error_semiminor_axis+error_ellipse_orientation, idData@data)[,-1] 
#   idData@data = cbind(idData@data, argosDiag2Cov(diagData[,1], diagData[,2], diagData[,3]))
#   init = list(
#     a=c(coordinates(idData)[1,1],0,coordinates(idData)[1,2],0),
#     P=diag(c(5000^2,10*3600^2, 5000^2, 10*3600^2))
#   )
#   model.fits[[i]] = crwMLE(mov.model=~1, 
#                            err.model=list(x=~ln.sd.x-1, y=~ln.sd.y-1, rho=~error.corr), 
#                            data=idData, 
#                            Time.name="unique_posix", 
#                            initial.state=init,
#                            fixPar = c(1,1,NA,NA), 
#                            theta=c(log(10), 3), 
#                            initialSANN=list(maxit=2500),
#                            #prior=function(par){-(abs(par[2]-3)/10)^2/2},
#                            control=list(REPORT=10, trace=1))
#   print(model.fits[[i]])
# }


## ----message=FALSE-------------------------------------------------------
# predData=NULL
# for(i in 1:length(model.fits)){
#   if(any(is.nan(model.fits[[i]]$se))) {
#     next
#   }
#   model.fits[[i]]$data$unique_posix <- lubridate::with_tz(
#     model.fits[[i]]$data$unique_posix,"GMT")
#   predTimes <- seq(
#     lubridate::ceiling_date(min(model.fits[[i]]$data$unique_posix),"hour"),
#     lubridate::floor_date(max(model.fits[[i]]$data$unique_posix),"hour"),
#     "1 hour")
#   tmp = crwPredict(model.fits[[i]], predTime=predTimes)
#   crwPredictPlot(tmp,plotType="map")
#   predData = rbind(predData, tmp)
# }
# 
# predData$predTimes <- intToPOSIX(predData$TimeNum)
# predData_sp <- predData
# coordinates(predData_sp) <- ~mu.x+mu.y
# proj4string(predData_sp) <- CRS("+init=epsg:3571")

## ----plot-1--------------------------------------------------------------
# p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + geom_path(aes(colour=deployid)) + 
#   fte_theme()
# p1

