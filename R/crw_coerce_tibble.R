#' Coerce crawl objects (crwIS and crwPredict) to tibbles
#'
#' @author Josh M. London
#' @param crw_object an object of class \code{"crwIS"} or \code{"crwPredict"},
#' @export

crw_as_tibble <- function(crw_object, ...) {
  UseMethod("crw_as_tibble",crw_object)
}

#' @describeIn crw_as_tibble coerce crwIS object to tibble
#' @export
crw_as_tibble.crwIS <- function(crw_object, ...) {
  tibble::tibble(mu.x = crw_object$alpha.sim[,'mu.x'],
                 nu.x = crw_object$alpha.sim[,'nu.x'],
                 mu.y = crw_object$alpha.sim[,'mu.y'],
                 nu.y = crw_object$alpha.sim[,'nu.y'],
                 pred_dt = lubridate::as_datetime(crw_object$Time),
                 locType = crw_object$locType
  ) %>% 
    dplyr::arrange(pred_dt)
}

#' @describeIn crw_as_tibble coerce crwPredict object to tibble
#' @export
crw_as_tibble.crwPredict <- function(crw_object, ...) {
  tibble::as_tibble(crw_object) %>% 
    dplyr::mutate(num_time = lubridate::as_datetime(crw_object$Time)) %>% 
    dplyr::rename(pred_dt = num_time) %>% 
    dplyr::arrange(pred_dt)
}

#' @describeIn crw_as_tibble 
#' @export
crw_as_tibble.tibble <- function(crw_object, ...) {
  crw_object
}
