## ----message=FALSE-------------------------------------------------------
library(crawl)

## ------------------------------------------------------------------------
data("northernFurSeal")
head(northernFurSeal)

## ------------------------------------------------------------------------
northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,
                                          levels=c("3", "2", "1","0","A"))

## ---- message=FALSE------------------------------------------------------
library(sp)
library(rgdal)
coordinates(northernFurSeal) = ~longitude+latitude

## ------------------------------------------------------------------------
proj4string(northernFurSeal) <- CRS("+proj=longlat")

## ----message=FALSE-------------------------------------------------------
northernFurSeal <- spTransform(northernFurSeal, 
                               CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                         "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                         "+ellps=GRS80 +datum=NAD83",
                                         "+units=m +no_defs"))
)

## ----message=FALSE-------------------------------------------------------
initial = list(a=c(coordinates(northernFurSeal)[1,1],0,
                   coordinates(northernFurSeal)[1,2],0),
               P=diag(c(10000^2,54000^2,10000^2,5400^2)))

## ----message=FALSE-------------------------------------------------------
fixPar = c(log(250), log(500), log(1500), rep(NA,3), NA)

## ----message=FALSE-------------------------------------------------------
displayPar(mov.model=~1,
           err.model=list(x=~Argos_loc_class-1),
           data=northernFurSeal,
           fixPar=fixPar)

## ---- message=FALSE------------------------------------------------------
constr=list(lower=c(rep(log(1500),2), rep(-Inf,2)),
            upper=rep(Inf,4))

## ---- message=FALSE------------------------------------------------------
ln.prior = function(theta){-abs(theta[4]-3)/0.5}

## ----message=FALSE-------------------------------------------------------
set.seed(123)
fit1 <- crwMLE(mov.model=~1, 
               err.model=list(x=~Argos_loc_class-1),
               data=northernFurSeal, 
               Time.name="Time",
               initial.state=initial,
               fixPar=fixPar, 
               constr=constr, 
               prior=ln.prior,
               control=list(maxit=30, trace=0,REPORT=1),
               initialSANN=list(maxit=200, trace=0, REPORT=1))

## ------------------------------------------------------------------------
fit1

## ----message=FALSE-------------------------------------------------------
predTime <- seq(ceiling(min(northernFurSeal$Time)), 
                floor(max(northernFurSeal$Time)), 1)

## ----message=FALSE-------------------------------------------------------
predObj <- crwPredict(object.crwFit=fit1, 
                      predTime, 
                      speedEst=TRUE, 
                      flat=TRUE)

## ----message=FALSE-------------------------------------------------------
crwPredictPlot(predObj, "map")

## ----message=FALSE-------------------------------------------------------
set.seed(123)
simObj <- crwSimulator(fit1, 
                       predTime, 
                       method="IS", 
                       parIS=100, 
                       df=5, 
                       scale=18/20)

## ----message=FALSE-------------------------------------------------------
w <- simObj$thetaSampList[[1]][,1]
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')

## ----message=FALSE-------------------------------------------------------
round(100/(1+(sd(w)/mean(w))^2))

## ------------------------------------------------------------------------
my.colors <-colorRampPalette(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))

## ------------------------------------------------------------------------
iter <- 20
cols <- my.colors(iter)

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
library(argosfilter)
library(parallel)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)

## ----message=FALSE-------------------------------------------------------
data("beardedSeals")
beardedSeals

## ----message=FALSE-------------------------------------------------------
beardedSeals %>% 
  group_by(deployid,date_time) %>% 
  filter(n()>1)

## ----message=FALSE-------------------------------------------------------
library(xts)
date_unique <-beardedSeals %>% 
  group_by(deployid) %>%
  do(unique_date = xts::make.time.unique(.$date_time,eps=1)) %>%
  tidyr::unnest(unique_date) %>%
  mutate(unique_posix = as.POSIXct(.$unique_date,origin='1970-01-01 00:00:00',tz='UTC')) %>%
  dplyr::arrange(deployid,unique_posix) %>% 
  dplyr::select(unique_posix)

beardedSeals <- beardedSeals %>% arrange(deployid,date_time) %>%
  bind_cols(date_unique)

## ----message=FALSE-------------------------------------------------------
beardedSeals %>% 
  group_by(deployid,unique_posix) %>% 
  filter(n()>1)

## ----message=FALSE-------------------------------------------------------
beardedSeals <- beardedSeals %>%  
  dplyr::arrange(deployid,unique_posix)

# speedfilter using the paralell package for multi-core speed
cfilter<-mclapply(split(beardedSeals, beardedSeals$deployid),
                  function(x) sdafilter(
                    lat=x$latitude, 
                    lon=x$longitude, 
                    dtime=x$unique_posix,
                    lc=x$quality, 
                    ang=-1,
                    vmax=5),
                  mc.preschedule=F,mc.cores=3)

cfilter<-do.call("c",cfilter)
cfilter<-as.vector(cfilter)
beardedSeals$filtered <- cfilter

beardedSeals <- beardedSeals %>% 
  dplyr::filter(filtered=="not", !is.na(error_semimajor_axis)) %>%
  arrange(deployid,unique_posix)

## ----message=FALSE-------------------------------------------------------
beardedSeals <- as.data.frame(beardedSeals)
coordinates(beardedSeals) = ~longitude+latitude
proj4string(beardedSeals) = CRS("+proj=longlat +datum=WGS84")

beardedSeals <- spTransform(beardedSeals, CRS("+init=epsg:3571"))

## ----message=FALSE-------------------------------------------------------
ids = unique(beardedSeals@data$deployid)      #define seal IDs

library(doParallel)
n.cores <- detectCores()
registerDoParallel(cores=3)

model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(beardedSeals,deployid == ids[i])
    
    diag_data = model.matrix(
      ~ error_semimajor_axis + error_semiminor_axis + error_ellipse_orientation,
      id_data@data
    )[,-1]
    
    id_data@data = cbind(id_data@data, 
                         argosDiag2Cov(
                           diag_data[,1], 
                           diag_data[,2], 
                           diag_data[,3]))
    
    init = list(a = c(coordinates(id_data)[1,1],0,
                      coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                           5000 ^ 2, 10 * 3600 ^ 2)))
    
    fit <- crwMLE(
      mov.model =  ~ 1,
      err.model = list(
        x =  ~ ln.sd.x - 1, 
        y =  ~ ln.sd.y - 1, 
        rho =  ~ error.corr
      ),
      data = id_data,
      Time.name = "unique_posix",
      initial.state = init,
      fixPar = c(1,1,NA,NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }

names(model_fits) <- ids

print(model_fits)

## ----message=FALSE-------------------------------------------------------
predData <- foreach(i = 1:length(model_fits)) %dopar% {

  model_fits[[i]]$data$unique_posix <- lubridate::with_tz(
    model_fits[[i]]$data$unique_posix,"GMT")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix),"hour"),
    lubridate::floor_date(max(model_fits[[i]]$data$unique_posix),"hour"),
    "1 hour")
  tmp = crwPredict(model_fits[[i]], predTime=predTimes)
}

predData <- dplyr::bind_rows(predData) %>% as.data.frame(.)

predData$predTimes <- intToPOSIX(predData$TimeNum)

## ----plot-1--------------------------------------------------------------
theme_map = function(base_size=9, base_family="")
{
    require(grid)
    theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.title.x=element_text(vjust=0),
          axis.title.y=element_text(angle=90, vjust=1.25),
          axis.text.y=element_text(angle=90),
          axis.ticks=element_line(colour="black", size=0.25),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text=element_text(),
          legend.title=element_text(face="bold", hjust=0),
          panel.border=element_rect(fill=NA, colour="black"),
          panel.grid.major=element_line(colour="grey92", size=0.3, linetype=1),
          panel.grid.minor=element_blank(),
          plot.title=element_text(vjust=1),
          strip.background=element_rect(fill="grey90", colour="black", size=0.3),
          strip.text=element_text()
          )
}

p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + 
  geom_path(aes(colour=deployid)) + xlab("easting (meters)") +
  ylab("northing (meters)") + theme_map()
p1

