#' Coerce to sf/sfc object
#' 
#' Provides reliable conversion of \code{"crwIS"} and \code{"crwPredict"} objects
#' into simple features objects supported in the \code{"sf"} package. Both 
#' \code{"sf"} objects with "POINT" geometry and \code{"sfc_LINESTRING"} objects
#' are created. Coersion of \code{"crwPredict"} objects to \code{"sfc_LINESTRING"}
#' has an option \code{"group"} argument when the \code{"crwPredict"} object
#' includes predictions from multiple deployments. The grouping column will be 
#' used and a tibble of multiple \code{"sf_LINESTRING"} objects will be returned
#' 
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param ftype character of either "POINT" or "LINESTRING" specifying the feature type
#' @param locType character vector of location points to include ("p","o")
#' @param group (optional) character specifying the column to group by for mulitple LINESTRING features
#' @param ... Additional arguments that are ignored
#' @export

crw_as_sf <- function(crw_object,ftype,locType,group) {
  UseMethod("crw_as_sf",crw_object)
}

#' @describeIn crw_as_sf coerce crwIS object to sf (POINT or 
#' LINESTRING geometry)
#' @export 
crw_as_sf.crwIS <- function(crw_object,
                            ftype,
                            locType = c("p", "o"),
                            group = NULL, ...) {
  if (!is.null(group)) {
    warning("group argument not applicable to crwIS objects. ignorning")
  }
  
  stopifnot(!missing(ftype), ftype %in% c("POINT", "LINESTRING"))
  
  crw_crs <- attr(crw_object, "epsg")
  if(is.null(crw_crs) | is.na(crw_crs)) crw_crs <- attr(crw_object, "proj4")
  
  if (ftype == "POINT") {
    crw_object <- crw_as_tibble(crw_object) %>%
      dplyr::filter(.data$locType %in% locType,
                    !is.na(.data$mu.x),
                    !is.na(.data$mu.y)) %>%
      sf::st_as_sf(coords = c("mu.x", "mu.y"))
    crw_object = crw_object %>% sf::st_set_crs(crw_crs)
  }
  if (ftype == "LINESTRING") {
    crw_object <- crw_as_tibble(crw_object) %>%
      dplyr::filter(.data$locType %in% locType,
                    !is.na(.data$mu.x),
                    !is.na(.data$mu.y)) %>%
      sf::st_as_sf(coords = c("mu.x", "mu.y"))
    crw_object = crw_object %>% sf::st_set_crs(crw_crs)
    crw_object = crw_object %>%
      summarise(id = 1, do_union = FALSE) %>%
      sf::st_cast("LINESTRING")
  }
  return(crw_object)
}

#' @describeIn crw_as_sf coerce crwPredict object to sf (POINT or 
#' LINESTRING geometry) 
#' @export
crw_as_sf.crwPredict <- function(crw_object,ftype,
                                 locType = c("p","o"), 
                                 group = NULL, ...) {
  stopifnot(!missing(ftype), ftype %in% c("POINT","LINESTRING"))
  crw_crs <- attr(crw_object, "epsg")
  if(is.null(crw_crs) | is.na(crw_crs)) crw_crs <- attr(crw_object, "proj4")
  if(ftype == "POINT" && is.null(group)) {
    crw_object <- crw_as_tibble(crw_object) %>% 
      dplyr::filter(.data$locType %in% locType) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) 
    crw_object = crw_object %>% sf::st_set_crs(crw_crs)
  }
  if(ftype == "POINT" && !is.null(group)) {
    warning("group argument not applicable for 'POINT' type. ignoring")
  }
  if(ftype == "LINESTRING" && is.null(group)) {
    crw_object <- crw_as_tibble(crw_object) %>% 
      dplyr::filter(.data$locType %in% locType) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y"))
    crw_object = crw_object %>% sf::st_set_crs(crw_crs)
    crw_object = crw_object %>% summarise(id=1,do_union = FALSE) %>% sf::st_cast("LINESTRING")
  }
  if(ftype == "LINESTRING" && !is.null(group)) {
    crw_object <- crw_as_tibble(crw_object) %>% 
      dplyr::filter(.data$locType %in% locType) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) 
    crw_object = crw_object %>% sf::st_set_crs(crw_crs)
    crw_object = crw_object %>% dplyr::group_by(group) %>% 
      dplyr::summarise(do_union = FALSE) %>% 
      sf::st_cast("LINESTRING")
  }
  return(crw_object)
}

#' @export

crw_as_sf.sf <- function(crw_object,ftype,
                         locType = c("p","o"), 
                         group = NULL, ...) {
  message("No conversion between ftypes yet :-(")
  crw_object 
}
