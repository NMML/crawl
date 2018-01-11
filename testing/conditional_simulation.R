library(tidyverse)
library(crawl)
library(sf)

ak_base = nPacMaps::alaska(resolution = "f")

data("harborSeal")

harborSeal_sf <- harborSeal %>% 
  filter(!is.na(latitude)) %>%
  as.data.frame() %>% 
  sf::st_as_sf(coords = c("longitude","latitude")) %>% 
  sf::st_set_crs(4326) %>% 
  sf::st_transform(3338)

harborSeal_sf %>% st_coordinates() %>% as.data.frame() %>% 
  bind_cols(harborSeal_sf %>% st_set_geometry(.,NULL)) %>% 
  right_join(harborSeal) -> harborSeal

levels(harborSeal$Argos_loc_class) = c("3","2","1","0","A","B")

#######################################################
# fit crawl model
######################################################

initial = list(
  a=c(harborSeal$X[1],0,
      harborSeal$Y[1],0),
  P=diag(c(10000^2,5400^2,10000^2,5400^2))
)

fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)

constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal, coord=c("X","Y"), Time.name="Time", 
  initial.state=initial, fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)#,
  #initialSANN=list(maxit=200, trace=1, REPORT=1)
)

predTimes <- seq(min(harborSeal$Time), max(harborSeal$Time), by = 0.5)

hs_pred <- crawl::crwPredict(fit1, predTime=predTimes, flat=FALSE)
state = hs_pred$alpha.hat %>% bind_cols(hs_pred$originalData) %>% filter(locType=='p')
path_sf = state %>% select(starts_with("mu")) %>% 
  st_as_sf(coords=c('mu.x','mu.y')) %>% sf::st_set_crs(3338)

state %>% mutate(
  on_land=sf::st_intersects(path_sf, ak_base) %>% purrr::map_lgl(~ length(.x) > 0)
  ) -> state

b <- sf::st_bbox(path_sf)
# ggplot() +
#   geom_sf(data = ak_base,
#           fill = "grey60", size = 0.2) +
#   geom_sf(data = path_sf[state$on_land,],
#           alpha = 0.1, color = 'blue') +
#   coord_sf(xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))

test_state = state[1239:1274,]
par = tail(fit1$estPar, 2)

### Necessary objects

t0 = test_state[1,"Time"]
alpha0 = test_state[1,1:4] %>% as.numeric()
t2 = tail(test_state,1)[,"Time"]
alpha2 = tail(test_state,1)[,1:4] %>% as.numeric()
t1 = test_state[2,"Time"]
active=1

#### Function

# cond_sim = function(n=500, t0, alpha0, t2, alpha2, t1, par, active=1, inf_fac=1){
#   beta=exp(par[2])
#   sigma2=exp(2*par[1])
#   delta=diff(c(t0, t1, t2))
#   T0 = crawl:::makeT(b=beta, delta=delta[1], active=active)
#   T1 = crawl:::makeT(b=beta, delta=delta[2], active=active)
#   Q0 = crawl:::makeQ(b=beta, sig2=sigma2, delta=delta[1], active=active)
#   Q1 = crawl:::makeQ(b=beta, sig2=sigma2, delta=delta[2], active=active)
#   V_inv = solve(Q0) + t(T1)%*%solve(Q1,T1)
#   v = solve(Q0, T0)%*%alpha0 + t(T1)%*%solve(Q1, alpha2)
#   mu_cond = solve(V_inv, v)
#   V_cond = solve(V_inv)
#   smp = mvtnorm::rmvnorm(n, mu_cond, inf_fac*V_cond)
# }

cond_sim = function(n=500, t0, alpha0, t2, alpha2, t1, par, active=1, inf_fac=1, bm=0){
  if (inherits(t0,"POSIXct")) {
    t0 <- as.numeric(t0)
  }
  if (inherits(t1,"POSIXct")) {
    t1 <- as.numeric(t1)
  }
  if (inherits(t2, "POSIXct")) {
    t2 <- as.numeric(t2)
  }
  if(bm){
    beta=exp(4)
  } else{
    beta=exp(par[2])
  }
  sigma2=exp(2*par[1])
  delta=diff(c(t0, t1, t2))
  T0 = crawl:::makeT(b=beta, delta=delta[1], active=active)
  T1 = crawl:::makeT(b=beta, delta=delta[2], active=active)
  Q0 = crawl:::makeQ(b=beta, sig2=sigma2, delta=delta[1], active=active)
  Q1 = crawl:::makeQ(b=beta, sig2=sigma2, delta=delta[2], active=active)
  V_inv = solve(Q0) + t(T1)%*%solve(Q1,T1)
  v = solve(Q0, T0)%*%alpha0 + t(T1)%*%solve(Q1, alpha2)
  mu_cond = solve(V_inv, v)
  V_cond = solve(V_inv)
  smp = mvtnorm::rmvnorm(n, mu_cond, inf_fac*V_cond)
}


plot(c(alpha0[1], alpha2[1]), c(alpha0[3], alpha2[3]), col=c('red','blue'), asp=1)
points(smp[,c(1,3)])




