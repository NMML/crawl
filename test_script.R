#library(crawl)
library(tidyverse)
library(sf)
library(nPacMaps)
library(ggspatial)

data("harborSeal_sf")
harborSeal_sf$Argos_loc_class = factor(harborSeal_sf$Argos_loc_class, levels=c("3","2","1","0","A","B"))
harborSeal_sf <- harborSeal_sf %>% 
  sf::st_transform(3338) %>% 
  dplyr::slice(550:1000)
fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)
displayPar( mov.model=~1, err.model=list(x=~Argos_loc_class-1),data=harborSeal_sf, 
            activity=~I(1-DryTime),fixPar=fixPar)
constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal_sf, Time.name="Time", fixPar=fixPar, 
  theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr,
  control=list(maxit=2000, trace=1, REPORT=1)#,
  #initialSANN=list(maxit=200, trace=1, REPORT=1)
)

pred1 = crwPredict(fit1, predTime="1 hour")
pred1_sf <- pred1 %>% crawl::crw_as_sf("LINESTRING")
ak_base <- nPacMaps::alaska()

esri_ocean <- 'https://services.arcgisonline.com/arcgis/rest/services/Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpg'

ggplot() +
annotation_map_tile(type = esri_ocean, zoomin = 1) +
  annotation_spatial(data = ak_sf) +
  layer_spatial(data = pred1_sf, lwd = 0.5) +
  ggtitle("Predicted Track of Harbor Seal in Alaska")

sim <- crawl::crwSimulator(fit1, predTime = "1 hour")
sim_track <- crawl::crwPostIS(sim, fullPost = FALSE)

sim_fix <- crawl::fix_path(sim_track,ak_base,fit1)

pred1_fix <- crawl::fix_path(pred1,ak_base,fit1)

pred1_fix_sf <- pred1_fix %>% crawl::crw_as_sf("LINESTRING")

pred1_pts <- pred1 %>% crw_as_sf("POINT","p")
pred1_range <- c(range(st_coordinates(pred1_pts)[,1]),
                 range(st_coordinates(pred1_pts)[,2]))
ak_base <- nPacMaps::alaska()
l_sfc %>% sf::st_sf() %>% 
  ggplot() + 
  geom_sf(data = ak_base, fill = '#c7e9c0', lwd=0.1) +
  geom_sf() +
  coord_sf(xlim = pred1_range[c(1:2)], ylim = pred1_range[c(3:4)]) +
  ggtitle("Predicted Track of Harbor Seal in Alaska")
