require(crawl)
data(harborSeal)
head(harborSeal)
## Calculate Log multipliers for Argos error
argosClasses <- c("3", "2", "1", "0", "A", "B")
ArgosMultFactors <- data.frame(Argos_loc_class=argosClasses, 
                               errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)), 
                               errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
hsNew <- merge(harborSeal, ArgosMultFactors, by=c("Argos_loc_class"), all=TRUE)
hsNew <- hsNew[order(hsNew$Time), ]
head(hsNew)


## Project data ##
library(sp)
library(rgdal)
# get EPSG code for AK projection
epsg = make_EPSG()
epsg[grep("Alaska", epsg$note),1:2]

toProj = hsNew[!is.na(hsNew$latitude),]
coordinates(toProj) = ~longitude+latitude
proj4string(toProj) <- CRS("+proj=longlat")
toProj <- spTransform(toProj, CRS("+init=epsg:3338"))
toProj = as.data.frame(toProj)
#names(toProj)[6:7] = c("x","y")
hsNew = merge(toProj, hsNew, all=TRUE)
hsNew = hsNew[order(hsNew$Time),]

initial.cpp = list(
  a=c(hsNew$x[1],0,hsNew$y[1],0),
  P=diag(c(100000,1,100000,1))
)

##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)

set.seed(123)
fit1 <- crwMLE_cpp(
  mov.model=~1, err.model=list(x=~errX, y=~errY), activity=~I(1-DryTime),
  data=hsNew, coord=c("x","y"), Time.name="Time", 
  initial.state=initial.cpp, fixPar=c(NA, 1, NA, 1, NA, NA, NA), 
  theta=c(5,5,0,0,0),
  control=list(maxit=2000, trace=1, REPORT=1),
  initialSANN=list(maxit=200, trace=1, REPORT=1)
)

print(fit1)

pred1 = crwPredict_cpp(fit1, predTime=NULL, speedEst=FALSE, flat=TRUE, getUseAvail=FALSE)

require(ggplot2)
p1=ggplot(aes(x=mu.x, y=mu.y), data=pred1) + geom_path(col="red", asp=TRUE) + geom_point(aes(x=x, y=y), col="blue") + coord_fixed()
p2=ggplot(aes(x=Time, y=mu.x), data=pred1) + geom_ribbon(aes(ymin=mu.x-2*se.mu.x, ymax=mu.x+2*se.mu.x), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=x), col="blue", size=1)
p3=ggplot(aes(x=Time, y=mu.y), data=pred1) + geom_ribbon(aes(ymin=mu.y-2*se.mu.y, ymax=mu.y+2*se.mu.y), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=y), col="blue", size=1)
print(p1)
print(p2)
print(p2)
# ggsave("map.pdf", p1)
# ggsave("xaxis.pdf", p2, width=10, height=2)
# ggsave("yaxis.pdf", p3, width=10, height=2)

# 
# ##See simulated annealing start values
# fit2$init$par
