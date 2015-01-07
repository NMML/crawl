# library(crawl)
library(ggplot2)
data(harborSeal)
head(harborSeal)
harborSeal$Argos_loc_class = factor(harborSeal$Argos_loc_class, levels=c("3","2","1","0","A","B"))

## Project data ##
library(rgdal)

toProj = harborSeal[!is.na(harborSeal$latitude),]
coordinates(toProj) = ~longitude+latitude
proj4string(toProj) <- CRS("+proj=longlat")
toProj <- spTransform(toProj, CRS("+init=epsg:3338"))
toProj = as.data.frame(toProj)
harborSeal = merge(toProj, harborSeal, all=TRUE)
harborSeal = harborSeal[order(harborSeal$Time),]

initial.cpp = list(
  a=c(harborSeal$x[1],0,harborSeal$y[1],0),
  P=diag(c(10000^2,54000^2,10000^2,5400^2))
)

##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)
fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)
displayPar( mov.model=~1, err.model=list(x=~Argos_loc_class-1),data=harborSeal, 
                activity=~I(1-DryTime),fixPar=fixPar)
constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("x","y"), Time.name="Time", 
  initial.state=initial.cpp, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 3),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1),
  initialSANN=list(maxit=200, trace=1, temp=20, tmax=20, REPORT=1)
)

print(fit1)
pred1 = crwPredict(fit1, predTime=NULL, flat=TRUE)
p1=ggplot(aes(x=mu.x, y=mu.y), data=pred1) + geom_path(col="red", asp=TRUE) + geom_point(aes(x=x, y=y), col="blue") + coord_fixed()
p2=ggplot(aes(x=Time, y=mu.x), data=pred1) + geom_ribbon(aes(ymin=mu.x-2*se.mu.x, ymax=mu.x+2*se.mu.x), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=x), col="blue", size=1)
p3=ggplot(aes(x=Time, y=mu.y), data=pred1) + geom_ribbon(aes(ymin=mu.y-2*se.mu.y, ymax=mu.y+2*se.mu.y), fill="green", alpha=0.5)  + 
  geom_path(, col="red") + geom_point(aes(x=Time, y=y), col="blue", size=1)
print(p1)
print(p2)
print(p3)
# ggsave("map.pdf", p1)
# ggsave("xaxis.pdf", p2, width=10, height=2)
# ggsave("yaxis.pdf", p3, width=10, height=2)


# 
# ##See simulated annealing start values
# fit2$init$par
