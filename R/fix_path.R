#' @title Find the sections of a path that pass thorugh a restricted area
#' @description This function is used to identify sections of a path that pass through 
#' restricted areas. the CTCRW model in crawl cannot actively steer paths away 
#' from restricted areas as it knows nothing of spatial information. So, this function
#' will identify areas that for which the unrestrained path passes through these areas.
#' If the path/points end within the restricted area, those records will be removed.
#' The user can then use this information to adjust the path as desired. 
#' @param xy A \code{SpatialPoints} object from the \code{sp} package or a 
#' 2-column matrix of x and y locations
#' @param res_raster A \code{raster} object from the raster package that indicates 
#' restricted areas with a 1, else 0 for unrestricted areas.  
#' @return A data.frame with each row associated with each section of the path
#' that crosses a restricted area. The columns provide the start and end row indices of \code{xy} where
#' the section occurs and the previous and post locations that are in unrestricted space.
#' @author Josh M. London (josh.london@noaa.gov)
#' @importFrom raster extract 
#' @export

get_restricted_segments = function(xy, res_raster){
  restricted <- raster::extract(res_raster, xy)
  if (any(is.na(restricted))) {
    stop("points in xy fall outside the extent of res_raster")
  }
  if (sum(restricted,na.rm=TRUE) == 0) {return(NULL)}
  
  head_start <- 1
  tail_end <- length(restricted)
  if (min(which(restricted == 1)) == 1) {
    warning(paste0("Path starts in restricted area, first ",
                   min(which(restricted==0)) - 1,
                   " observations removed"))
    head_start = min(which(restricted==0))
  }
  if(max(which(restricted==0)) < length(restricted)){
    warning(paste("Path ends in restricted area, last ", 
                  length(restricted)-max(which(restricted==0)),
                  " observations removed"))
    tail_end <- max(which(restricted==0))
  }
  xy <- xy[head_start:tail_end,]
  restricted <- restricted[head_start:tail_end]
  
  in.segment <- (restricted > 0)
  start_idx <- which(c(FALSE, in.segment) == TRUE &
                       dplyr::lag(c(FALSE, in.segment) == FALSE)) - 1
  end_idx <- which(c(in.segment, FALSE) == TRUE & 
                     dplyr::lead(c(in.segment, FALSE) == FALSE))
  restricted_segments <- data.frame(start_idx, end_idx) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(start_x = xy[start_idx-1,1],
                  start_y = xy[start_idx-1,2],
                  end_x = xy[end_idx+1,1],
                  end_y = xy[end_idx+1,2])
  restricted_segments <- list(
    restricted_segments = restricted_segments,
    fixed_range = c(head_start,tail_end)
  )
  return(restricted_segments)
}


#' @title Project path away from restricted areas
#' @description Corrects a path so that it does not travel through a restricted area.
#' @param xy Coordinate locations for the path. Can be one of the following classes: 
#' (1) a two column matrix, 
#' (2) 'SpatialPoints' or 'SpatialPointsDataFrame' object from the sp package,
#' (3) 'crwPredict' object from the \code{crwPredict} function
#' (4) 'crwIS' object from the \code{crwPostIS} function
#' @param time A vector of times associated with xy locations
#' @param res_raster An indicator raster object with cells = 1 if it is 'off-limits'
#' and 0 elsewise.
#' @param trans A transition matrix object from the gdistance package.
#' @return Either matrix or 'SpatialPoints' object with path projected around
#' restricted areas
#' @importFrom gdistance transition shortestPath
#' @importFrom sp coordinates SpatialPointsDataFrame
#' @importFrom raster cellFromXY 
#' @importFrom stats approx
#' @export
#' 
fix_path = function(xy, time, res_raster, trans){
  if(inherits(xy, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc_data = sp::coordinates(xy)
  } else if(inherits(xy, "matrix")){
    if(ncol(xy)!=2) stop("xy matrix does not have 2 columns")
    loc_data = xy
  } else if(inherits(xy, "crwPredict")){
    loc_data = as.matrix(xy[,c("mu.x","mu.y")]) 
  } else if(inherits(xy, "crwIS")){
    loc_data = xy$alpha.sim[,c("mu.x","mu.y")]
  } else stop("Unrecognized 'xy' format")
  
  if (!missing(time) && inherits(xy,c("crwPredict","crwIS"))) {
    warning("time vector provided for crwPredict or crwIS object. time vector ignored")
  }
  
  rs = get_restricted_segments(loc_data, res_raster)
  if (is.null(rs$restricted_segments) || 
     nrow(rs$restricted_segments) == 0) { 
    return(xy)
    }
  
  seg = rs$restricted_segments
  loc_data = loc_data[rs$fixed_range[1]:rs$fixed_range[2],]

  if (!missing(time) && !inherits(xy,c("crwPredict","crwIS"))) {
    time = time[rs$fixed_range[1]:rs$fixed_range[2]]
  }

  idx = as.matrix(seg[,1:2])
  start_xy = as.matrix(seg[,3:4])
  start_cell = cellFromXY(res_raster, start_xy)
  end_xy = as.matrix(seg[,5:6])
  end_cell = cellFromXY(res_raster, end_xy)
  for(i in 1:nrow(seg)){
    if(start_cell[i] == end_cell[i]){
      path_pts = do.call("cbind",
                         stats::approx(x=c(start_xy[i,1],end_xy[i,1]), 
                                       y=c(start_xy[i,2],end_xy[i,2]),
                                       n=as.integer(seg[i,2]-seg[i,1]+1)
                         ))
    } else{
      path = gdistance::shortestPath(trans, start_xy[i,], end_xy[i,],"SpatialLines")
      path_pts = sp::coordinates(
        sp::spsample(path, n=as.integer(seg[i,2]-seg[i,1]+1), "regular")
      )
    }
    loc_data[idx[i,1]:idx[i,2],] = as.matrix(path_pts)
  }

  if (inherits(xy, "SpatialPoints")){
    loc_data = as.data.frame(loc_data)
    if (inherits(xy,"SpatialPointsDataFrame")) {
      loc_data <- cbind(loc_data,time,
                        xy@data[rs$fixed_range[1]:rs$fixed_range[2],])
    }
    sp::coordinates(loc_data) = c(1,2)
    sp::proj4string(loc_data) = sp::proj4string(xy)
    
  }
  if (inherits(xy,"crwIS")) {
    loc_data <- as.data.frame(loc_data)
    loc_data <- cbind(loc_data,
                      Time = xy$Time[rs$fixed_range[1]:rs$fixed_range[2]],
                      locType = xy$locType[rs$fixed_range[1]:rs$fixed_range[2]])
  }
  if (inherits(xy,"matrix")) {
    loc_data <- cbind(loc_data,time)
  }
  return(loc_data)
  
  if(!missing(time)) {
    loc_data = SpatialPointsDataFrame(
      loc_data,data=data.frame(time))

    return(loc_data)
  } else{
    return(loc_data)
  }
}
