#' Coerce to sf/sfc object
#' 
#' Provides reliable conversion of \code{"crwIS"} and \code{"crwPredict"} objects
#' into simple features objects supported in the \code{"sf"} package. Both 
#' \code{"sf"} objects with "POINT" geometry and \code{"sfc_LINESTRING"} objects
#' are created. Coersion of \code{"crwPredict"} objects to \code{"sfc_LINESTRING"}
#' has an option \code{"group"} argument when the \code{"crwPredict"} object
#' includes predictions from multiple deployments. The grouping column will be 
#' used and a tibble of multiple \code{"sfc_LINESTRING"} objects will be returned
#' 
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param epsg integer epsg code specifying the coordinate projection
#' @param ftype character of either "POINT" or "LINESTRING" specifying the feature type
#' @param group (optional) character specifying the column to group by for mulitple LINESTRING features
#' @export

crw_as_sf <- function(crw_object,epsg,ftype,group) {
  UseMethod("crw_as_sf",crw_object)
}

#' @describeIn crw_as_sf coerce crwIS object to sf (POINT geometry) or sfc_LINESTRING
#' @export 
crw_as_sf.crwIS <- function(crw_object,epsg,ftype,group) {
  if(!missing(group)) {
    warning("group argument not applicable to crwIS objects. ignorning")
  }
  stopifnot(!missing(epsg))
  stopifnot(!missing(ftype))
  stopifnot(ftype %in% c("POINT","LINESTRING"))
  if(ftype == "POINT") {
  crw_object <- crw_as_tibble(crw_object) %>% 
    sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
    sf::st_set_crs(epsg) 
  }
  if(ftype == "LINESTRING") {
    crw_object <- crw_as_tibble(crw_object) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
      sf::st_set_crs(epsg) %>% 
      sf::st_coordinates() %>% 
      sf::st_linestring() %>% 
      sf::st_sfc()
  }
  return(crw_object)
}

#' @describeIn crw_as_sf coerce crwPredict object to sf (POINT geometry) or sfc_LINESTRING
#' @export
crw_as_sf.crwPredict <- function(crw_object,epsg,ftype,group) {
  stopifnot(!missing(epsg))
  stopifnot(!missing(ftype))
  stopifnot(ftype %in% c("POINT","LINESTRING"))
  if(ftype == "POINT" && missing(group)) {
  crw_object <- crw_as_tibble(crw_object) %>% 
    sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
    sf::st_set_crs(epsg) 
  }
  if(ftype == "POINT" && !missing(group)) {
   warning("group argument not applicable for 'POINT' type. ignoring")
  }
  if(ftype == "LINESTRING" && missing(group)) {
    crw_object <- crw_as_tibble(crw_object) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
      sf::st_set_crs(epsg) %>% 
      sf::st_coordinates() %>% 
      sf::st_linestring() %>% 
      sf::st_sfc()
  }
  if(ftype == "LINESTRING" && !missing(group)) {
    crw_object <- crw_as_tibble(crw_object) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
      sf::st_set_crs(epsg) %>% 
      dplyr::group_by(group) %>% 
      dplyr::summarize() %>% 
      sf::st_cast("LINESTRING") %>% 
      sf::st_sfc()
  }
  return(crw_object)
}
