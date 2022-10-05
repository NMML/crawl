#' Coerce to sf/sfc object
#' 
#' Provides reliable conversion of \code{"crwIS"} and \code{"crwPredict"} objects
#' into simple features objects supported in the \code{"sf"} package. Both 
#' \code{"sf"} objects with "POINT" geometry and \code{"sfc_LINESTRING"} objects
#' are created. Coercion of \code{"crwPredict"} objects to \code{"sfc_LINESTRING"}
#' has an option \code{"group"} argument when the \code{"crwPredict"} object
#' includes predictions from multiple deployments. The grouping column will be 
#' used and a tibble of multiple \code{"sf_LINESTRING"} objects will be returned
#' 
#' @param data an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param ftype character of either "POINT" or "LINESTRING" specifying the feature type
#' @param locType character vector of location points to include ("p","o")
#' @param group (optional) character specifying the column to group by for multiple LINESTRING features
#' @param ... Additional arguments that are ignored
#' @export

crw_as_sf <- function(data,ftype,locType,group) {
  UseMethod("crw_as_sf",data)
}

#' @describeIn crw_as_sf coerce crwIS object to sf (POINT or 
#' LINESTRING geometry)
#' @export 
crw_as_sf.crwIS <- function(data,
                            ftype,
                            locType = c("p", "o", "f"),
                            group = NULL, ...) {
  if (!is.null(group)) {
    warning("group argument not applicable to crwIS objects. ignorning")
  }
  stopifnot(!missing(ftype), ftype %in% c("POINT", "LINESTRING"))
  
  crw_crs <- attr(data, "epsg")
  if(is.null(crw_crs) || is.na(crw_crs)) crw_crs <- attr(data, "proj4")  
  if (ftype == "POINT") {
    data <- crw_as_tibble(data) %>%
      dplyr::filter(.data$locType %in% {{ locType }},
                    !is.na(.data$mu.x),
                    !is.na(.data$mu.y)) %>%
      sf::st_as_sf(coords = c("mu.x", "mu.y"))
    data = data %>% sf::st_set_crs(crw_crs)
  }
  if (ftype == "LINESTRING") {
    data <- crw_as_tibble(data) %>%
      dplyr::filter(.data$locType %in% {{ locType }},
                    !is.na(.data$mu.x),
                    !is.na(.data$mu.y)) %>%
      sf::st_as_sf(coords = c("mu.x", "mu.y"))
    data = data %>% sf::st_set_crs(crw_crs)
    data = data %>%
      summarise(id = 1, do_union = FALSE) %>%
      sf::st_cast("LINESTRING")
  }
  return(data)
}

#' @describeIn crw_as_sf coerce crwPredict object to sf (POINT or 
#' LINESTRING geometry) 
#' @export
crw_as_sf.crwPredict <- function(data,ftype,
                                 locType = c("p","o","f"),
                                 group = NULL, ...) {
  stopifnot(!missing(ftype), ftype %in% c("POINT","LINESTRING"))

  crw_crs <- attr(data, "epsg")
  if(is.null(crw_crs) || is.na(crw_crs)) crw_crs <- attr(data, "proj4")
  
  if(ftype == "POINT" && is.null(group)) {
    data <- crw_as_tibble(data) %>% 
      dplyr::filter(.data$locType %in% {{ locType }} ) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) 
    data = data %>% sf::st_set_crs(crw_crs)
  }
  if(ftype == "POINT" && !is.null(group)) {
    warning("group argument not applicable for 'POINT' type. ignoring")
  }
  if(ftype == "LINESTRING" && is.null(group)) {

    data <- crw_as_tibble(data) %>% 
      dplyr::filter(.data$locType %in% {{ locType }} ) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y"))
    data = data %>% sf::st_set_crs(crw_crs)
    data = data %>% summarise(id=1,do_union = FALSE) %>% sf::st_cast("LINESTRING")
  }
  if(ftype == "LINESTRING" && !is.null(group)) {

    data <- crw_as_tibble(data) %>% 
      dplyr::filter(.data$locType %in% {{ locType }} ) %>% 
      dplyr::arrange(.data$TimeNum) %>% 
      sf::st_as_sf(coords = c("mu.x","mu.y")) 
    data = data %>% sf::st_set_crs(crw_crs)
    data = data %>% dplyr::group_by(group) %>% 
      dplyr::summarise(do_union = FALSE) %>% 
      sf::st_cast("LINESTRING")
  }
  return(data)
}

#' @describeIn crw_as_sf coerce list of crwIS objects to sf (LINESTRING or 
#' MULTILINESTRING geometry) 
#' @export
crw_as_sf.list <- function(data,ftype,
                           locType = c("p","o","f"), ...) {
  
  is_list_of_crwis <- data %>% 
    purrr::modify_depth(1, ~purrr::map_lgl(.,inherits,"crwIS")) %>% 
    purrr::map_lgl(all) %>% 
    all()
  
  stopifnot(is_list_of_crwis)
  
  data <- data %>% 
    purrr::modify_depth(1, ~ purrr::map(., crawl::crw_as_sf,
                                        ftype = "LINESTRING", 
                                        locType = {{ locType }} ))
  
  if (ftype == "MULTILINESTRING") {
    make_mls <- function(ll) {
      do.call(rbind,ll) %>% 
        dplyr::group_by(id) %>% 
        dplyr::summarise(do_union = FALSE)
    }
    sf_list <- data %>% purrr::map(make_mls)
  }
  if (ftype == "LINESTRING") {
    make_mls <- function(ll) {
      do.call(rbind,ll) 
    }
    sf_list <- data %>% purrr::map(make_mls)
  }
  
  return(sf_list)
}
#' @export
crw_as_sf.sf <- function(data,ftype,
                         locType = c("p","o","f"),
                         group = NULL, ...) {
  message("No conversion between ftypes yet :-(")
  data 
}
