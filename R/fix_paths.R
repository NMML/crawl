#' @title Find the sections of a path that pass thorugh a restricted area
#' 
#' @description This function is used to identify sections of a path that pass through 
#' restricted areas. the CTCRW model is crawl cannot actively steer paths away 
#' from restricted areas as it knows nothing of spatial information. So, this function
#' will identify areas that for which the unrestrained path passes through these areas.
#' The user can then use this information to adjust the path as desired. 
#' @param xy A \code{SpatialPoints} object from the \code{sp} package or a 
#' 2-column matrix of x and y locations
#' @param res_raster A \code{raster} object from the raster package that indicates 
#' restricted areas with a 1, else 0 for unrestricted areas.  
#' @return A data.frame with each row associated with each section of the path
#' that crosses a restricted area. The columns provide the start and end row indices of \code{xy} where
#' the section occurs and the previous and post locations that are in unrestricted space.
#' @author Josh M. London (josh.london@noaa.gov)
#' @importFrom dplyr lag lead mutate %>%
#' @importFrom raster extract 
#' @export

get_restricted_segments = function(xy, res_raster){
  restricted <- raster::extract(res_raster, xy)
  in.segment <- (restricted > 0)
  start_idx <- which(c(FALSE, in.segment) == TRUE &
                       dplyr::lag(c(FALSE, in.segment) ==FALSE)) - 1
  end_idx <- which(c(in.segment, FALSE) == TRUE & 
                     dplyr::lead(c(in.segment, FALSE) == FALSE))
  restricted_segments <- data.frame(start_idx, end_idx) %>% 
    rowwise() %>% 
    dplyr::mutate(start_x = xy$mu.x[start_idx-1],
           start_y = xy$mu.y[start_idx-1],
           end_x = xy$mu.x[end_idx+1],
           end_y = xy$mu.y[end_idx+1])
  return(restricted_segments)
}
