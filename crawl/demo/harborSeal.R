
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

## Initial state values
initial.dry <- list(
  a1.x=c(harborSeal$longitude[1],0),
  a1.y=c(harborSeal$latitude[1],0),
  P1.x=diag(c(1,1)),
  P1.y=diag(c(1,1))
)

##Fit model as given in Johnson et al. (2008) Ecology 89:1208-1215
## Start values for theta come from the estimates in Johnson et al. (2008)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
  data=hsNew, coord=c("longitude","latitude"), polar.coord=TRUE, Time.name="Time", 
  initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA), theta=c(-6,-7,-4,-0.5,-1),
  control=list(maxit=2000, trace=1, REPORT=1),
)

fit1
str(fit1)

##Use simulated annealing to obtain start values and place constraints on the parameters

set.seed(123)
fit2 <- crwMLE(
  mov.model=~1, err.model=list(x=~errX, y=~errY), stop.model=~DryTime,
  data=hsNew, coord=c("longitude","latitude"), Time.name="Time", polar.coord=TRUE,
  initial.state=initial.dry, fixPar=c(NA, 1, NA, 1, NA, NA, NA),
  control=list(maxit=2000, trace=1, REPORT=1),
  constr=list(lower=c(-6, -6, -Inf, -Inf, -Inf), upper=Inf),
  initialSANN=list(maxit=100, temp=5, tmax=5, trace=1, REPORT=2)
)

fit2

##See simulated annealing start values
fit2$init$par
