#' Coerce crawl objects (crwIS and crwPredict) to various objects
#'
#' @author Josh M. London
#' @name coerce-crw-methods
#' @rdname coerce-crw-methods
NULL

#' Coerce to tibble
#' 
#' @name coerce-crw
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"},
#' @rdname coerce-crw-methods
#' @export

crw_as_tibble <- function(crw_object) {
  UseMethod("crw_as_tibble",crw_object)
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_tibble crwIS
#' @export
crw_as_tibble.crwIS <- function(crw_object) {
  tibble::tibble(mu.x = crw_object$alpha.sim[,'mu.x'],
                 mu.y = crw_object$alpha.sim[,'mu.y'],
                 pred_dt = lubridate::as_datetime(crw_object$Time),
                 locType = crw_object$locType
  ) %>% 
    dplyr::arrange(pred_dt)
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_tibble crwPredict
#' @export
crw_as_tibble.crwPredict <- function(crw_object) {
  tibble::as_tibble(crw_object) %>% 
    dplyr::mutate(num_time = lubridate::as_datetime(Time)) %>% 
    dplyr::rename(pred_dt = num_time) %>% 
    dplyr::arrange(pred_dt)
}

#' Coerce to sf point object
#' 
#' @name coerce-crw
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param epsg integer epsg code specifying the coordinate projection
#' @rdname coerce-crw-methods
#' @export

crw_as_sf_point <- function(crw_object,epsg) {
  UseMethod("crw_as_sf_point",crw_object)
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_sf_point crwIS
#' @export 
crw_as_sf_point.crwIS <- function(crw_object,epsg) {
  crw_as_tibble(crw_object) %>% 
    sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
    sf::st_set_crs(epsg)
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_sf_point crwPredict
#' @export
crw_as_sf_point.crwPredict <- function(crw_object,epsg) {
  crw_as_tibble(crw_object) %>% 
    sf::st_as_sf(coords = c("mu.x","mu.y")) %>% 
    sf::st_set_crs(epsg)
}

#' Coerce to sf line object
#' 
#' @name coerce-crw
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"}
#' @param epsg integer epsg code specifying the coordinate projection
#' @param group optional column identifying groups. if provided, separate linestrings will be returned for each group member 
#' @rdname coerce-crw-methods
#' @export
#' 
crw_as_sf_line <- function(crw_object, epsg, group) {
  UseMethod("crw_as_sf_line",crw_object)
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_sf_line crwIS
#' @export
crw_as_sf_line.crwIS <- function(crw_object,epsg,group=NULL) {
  stopifnot(!missing(epsg))
  if(!is.null(group)) {
    warning("group argument not applicable to crwIS objects. ignorning")
  }
    crw_as_sf_point(crw_object,epsg) %>% 
      sf::st_coordinates() %>% 
      sf::st_linestring() %>% 
      sf::st_sfc()
}

#' @return \code{NULL}
#'
#' @rdname coerce-crw-methods
#' @method crw_as_sf_line crwPredict
#' @export
crw_as_sf_line.crwPredict <- function(crw_object,epsg,group) {
  stopifnot(!missing(epsg))
  if(missing(group)) {
    crw_as_sf_point(crw_object,epsg) %>% 
      sf::st_coordinates() %>% 
      sf::st_linestring() %>% 
      sf::st_sfc()
  } else if(!missing(group)) {
    crw_as_sf_point(crw_object, epsg) %>% 
      dplyr::group_by(group) %>% 
      dplyr::summarize() %>% 
      sf::st_cast("LINESTRING") %>% 
      sf::st_sfc()
  }
}
