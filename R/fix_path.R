#' @title Identify segments of a path that cross through a restricted area
#' @description This function is used to identify sections of a path that pass through 
#' a restricted area (e.g. for marine mammals or fish, a land mask). the CTCRW model in 
#' crawl cannot actively steer paths away 
#' from land. So, this function
#' will identify path segments from the unrestrained path that pass through these areas.
#' If the path/points end within the land area, those records will be removed.
#' The user can then use this information to adjust the path as desired. 
#' @param crw_object A \code{crwIS} or \code{crwPredict} object from the \code{crawl} package 
#' @param vector_mask A \code{sf} object from sf package that indicates 
#' restricted areas as a polygon feature.  
#' @return A data.frame with each row associated with each section of the path
#' that crosses a restricted area. The columns provide the start and end row indices of \code{xy} where
#' the section occurs and the previous and post locations that are in unrestricted space.
#' @author Josh M. London (josh.london@noaa.gov)
#' @export
#' 

get_mask_segments = function(crw_sf, vector_mask, alpha) {
  # intersect crw_sf with vector mask
  on_mask <- sf::st_intersects(crw_sf, vector_mask) %>% 
    purrr::map_lgl(~ length(.x) > 0)
  if (any(is.na(on_mask))) {
    stop("points in crw_sf fall outside the extent of vector_mask")
  }
  # return NULL if no points within the vector mask
  if (sum(on_mask,na.rm = TRUE) == 0) {return(NULL)}
  
  head_start <- 1
  tail_end <- length(on_mask)
  
  if (min(which(on_mask == TRUE)) == 1) {
    warning(paste0("Path starts within vector mask, first ",
                   min(which(on_mask == 0)) - 1,
                   " observations removed"))
    head_start = min(which(on_mask == FALSE))
  }
  if (max(which(on_mask == 0)) < length(on_mask)) {
    warning(paste("Path ends within vector mask, last ", 
                  length(on_mask) - max(which(on_mask == 0)),
                  " observations removed"))
    tail_end <- max(which(on_mask == 0))
  }
  crw_sf <- crw_sf[head_start:tail_end,]
  alpha$alpha <- alpha$alpha[head_start:tail_end,]
  alpha$times <- alpha$times[head_start:tail_end]
  
  on_mask <- sf::st_intersects(crw_sf, vector_mask) %>% 
    purrr::map_lgl(~ length(.x) > 0)
  
  in.segment <- (on_mask == TRUE)
  
  start_idx <- which(c(FALSE, in.segment) == TRUE &
                       dplyr::lag(c(FALSE, in.segment) == FALSE)) - 2
  end_idx <- which(c(in.segment, FALSE) == TRUE & 
                     dplyr::lead(c(in.segment, FALSE) == FALSE)) + 1
  on_mask_segments <- data.frame(start_idx, end_idx) %>% 
    rowwise() %>% 
    dplyr::mutate(start_alpha = list(alpha$alpha[start_idx, ]),
                  end_alpha = list(alpha$alpha[end_idx,]),
                  times = list(alpha$times[start_idx:end_idx]))
  on_mask_segments <- list(
    on_mask_segments = on_mask_segments,
    fixed_range = c(head_start,tail_end)
  )
  
  return(on_mask_segments)
}


path_simulation_fix <- function(coast_points, sigma, beta, draw_size=500, vector_mask, fix_line=NULL, coast_line=NULL, seg_data=NULL, ...){
  ret_data = coast_points %>% 
    mutate(mu.x=NA, mu.y=NA, nu.x=NA, nu.y=NA)
  ret_data$mu.x[1] = coast_points$mu.x[1]
  ret_data$mu.y[1] = coast_points$mu.y[1]
  coast_points = coast_points %>% 
    tidyr::fill(mu.x, .direction = "up") %>% 
    tidyr::fill(mu.y, .direction = "up") %>% 
    dplyr::mutate(
      delta = c(diff(times),NA),
      Px = c(mu.x[2:n()],NA),
      Py = c(mu.y[2:n()],NA)
    )
  delta = coast_points$delta
  Px = coast_points$Px
  Py = coast_points$Py
  sig = sigma*beta*sqrt(delta)
  for(i in 1:(nrow(coast_points)-2)){
    if(delta[i]!=0){
      #if(i==1){
        nu.x.smp = (Px[i]-ret_data$mu.x[i])/delta[i]+rnorm(1000, 0, sig[i])
        nu.y.smp = (Py[i]-ret_data$mu.y[i])/delta[i]+rnorm(1000, 0, sig[i])
      #} else{
        # nu.x.smp = (1-beta*delta[i])*ret_data$nu.x[i-1] + beta*delta[i]*(Px[i]-ret_data$mu.x[i])/delta[i]+rnorm(1000, 0, sig[i])
        # nu.y.smp = (1-beta*delta[i])*ret_data$nu.y[i-1] + beta*delta[i]*(Py[i]-ret_data$mu.y[i])/delta[i]+rnorm(1000, 0, sig[i])
      #}
      mu.smp = data.frame(
        mu.x=ret_data$mu.x[i] + nu.x.smp*delta[i], 
        mu.y=ret_data$mu.y[i] + nu.y.smp*delta[i]) %>% 
        st_as_sf(coords=c("mu.x","mu.y")) %>% st_set_crs(st_crs(vector_mask))
      off_barrier = which(sapply(st_intersects(mu.smp,vector_mask), 
                                 function(z) if (length(z)==0) TRUE else FALSE))
      if(length(off_barrier)==0){
        NULL
      }
      mu.smp = mu.smp %>% st_coordinates
      ret_data$mu.x[i+1] = mu.smp[off_barrier[1],1]
      ret_data$mu.y[i+1] = mu.smp[off_barrier[1],2]
      ret_data$nu.x[i] = nu.x.smp[off_barrier[1]]
      ret_data$nu.y[i] = nu.y.smp[off_barrier[1]]
    } else{
      ret_data$mu.x[i+1] = ret_data$mu.x[i]
      ret_data$mu.y[i+1] = ret_data$mu.y[i]
      ret_data$nu.x[i] = ret_data$nu.x[i-1]
      ret_data$nu.y[i] = ret_data$nu.y[i-1]
    }
  }
  ret_data$nu.x = ret_data$nu.x
  ret_data$nu.x = ret_data$nu.x
  ret_data <- ret_data %>% 
    dplyr::filter(type %in% c("data","start")) %>% 
    dplyr::select(-type, -times)
  return(ret_data)
} 

path_prediction_fix = function(coast_points, sigma, beta, draw_size=500, 
                               paths=20, ...){
  ret_data = coast_points %>% 
    tidyr::fill(mu.x, .direction = "up") %>% 
    tidyr::fill(mu.y, .direction = "up") %>% 
    dplyr::filter(type %in% c("data","start","end")) %>% 
    mutate(
      nu.x = c(diff(mu.x)/diff(times), coast_points$nu.x[n()]),
      nu.y = c(diff(mu.y)/diff(times), coast_points$nu.y[n()])
    ) %>% 
    dplyr::filter(type %in% c("data","start")) %>%
    dplyr::select(-type, -times)
  return(ret_data)
}


#' @title Identify and re-route segments of a path that cross a restricted area
#' @description This function takes an `sf` object derived from a \code{crwIS} or 
#' \code{crwPredict} and an `sf` polygon object that defines the restricted
#' area. Each segment of the path that crosses the restricted area is identified
#' using the \code{get_mask_segments} function. Each segment begins and ends with
#' a coordinate that is outside the restricted area. The path is re-routed using
#' a mix of brownian movement and correlated parameters from the fitted movement
#' model (\code{crwFit}). A regular distribution of waypoints are established along
#' the restricted area boundary to guide the movement path around the area. 
#' @param crw_sf an `sf` object derived from a \code{crwIS} or \code{crwPredict} 
#' @param vector_mask an 'sf' polygon object that defines the restricted area
#' @param crwFit crwFit object that was used to generate the \code{crwIS} or \code{crwPredict}
#' @return a tibble with each record identifying the segments and pertinant values
#' @export

fix_segments <- function(crw_sf, vector_mask, barrier_buffer=50, crwFit, alpha, crwIS = FALSE) {
  # get par from crwFit
  par <- tail(crwFit$estPar, 2)
  ts <- attr(crwFit,"time.scale")
  
  # identify segments that are within the vector mask
  segments <- get_mask_segments(crw_sf, vector_mask, alpha)
  segments <- segments$on_mask_segments
  # add an empty 'fixed_seg' column to segments
  segments[,"fixed_seg"] <- list(list(NA))
  
  # begin loop through each segment that needs to be fixed
  for (i in 1:nrow(segments)) {
    start_idx <- segments[i,"start_idx"][[1]]
    end_idx <- segments[i,"end_idx"][[1]]
    message(paste('segment',i,'starts at', start_idx,'ends at', end_idx))
    # if (i == 131) {
    #   NULL
    # }
    # length of the segment is determined from the times column
    data_times <- list(times = segments[i, ]$times[[1]],
                       type = "data") %>% 
      tibble::as_tibble() 
    wypt_times <- list(times = c(seq(data_times$times[1], 
                                     data_times$times[nrow(data_times)], 
                                     by= 0.25))) %>% 
      tibble::as_tibble()
    l <- nrow(wypt_times)
    data_times <- data_times %>% 
      dplyr::slice(2:(n()-1))
    
    start_pt = sf::st_point(segments[i,]$start_alpha[[1]][c("mu.x", "mu.y")])
    end_pt = sf::st_point(segments[i,]$end_alpha[[1]][c("mu.x", "mu.y")])
    fix_pts <- sf::st_sfc(start_pt, end_pt) %>% sf::st_set_crs(sf::st_crs(crw_sf))
    fix_line <- fix_pts %>% sf::st_union() %>% sf::st_cast("LINESTRING")
    
    n_polys <- vector_mask %>% 
      sf::st_cast("POLYGON", warn = FALSE) %>% 
      sf::st_intersects(fix_line,sparse = FALSE) %>% 
      sum()
    
    if (n_polys == 0) {
      coast_line <- fix_line
    }
    
    if (n_polys != 0) {
      coast_hull <- vector_mask %>%  
        dplyr::filter(lengths(sf::st_intersects(., fix_line)) > 0) %>% 
        sf::st_difference(sf::st_buffer(sf::st_intersection(fix_line,.),dist = 1)) %>%
        sf::st_collection_extract(type = "POLYGON", warn = FALSE) %>% 
        sf::st_cast("POLYGON", warn = FALSE) %>% 
        sf::st_difference() %>% 
        dplyr::slice(which(!sf::st_is_empty(.))) %>% sf::st_sf() 
      
      coast_hull_over <- coast_hull %>% sf::st_intersects() %>% purrr::map_int(length)
      if(any(coast_hull_over>1)){
        coast_hull <- coast_hull %>% dplyr::slice(-which(coast_hull_over==max(coast_hull_over))) 
      }
      
      coast_hull <- coast_hull %>% # for crwIS segment 51 has 2 very large polys
        dplyr::slice(-which.max(sf::st_area(.))) %>% 
        sf::st_buffer(barrier_buffer) %>% 
        sf::st_union() %>% 
        sf::st_union(fix_line) %>% 
        sf::st_cast("POINT") %>% 
        sf::st_union() %>% 
        sf::st_convex_hull() %>% 
        sf::st_buffer(barrier_buffer) %>% 
        sf::st_sf()
      
      coast_hull_buffer <- coast_hull %>% 
        sf::st_buffer(dist = 1) %>% 
        sf::st_cast("LINESTRING") %>% 
        sf::st_buffer(dist = 1e-5)  
      
      fix_line <- fix_pts %>% 
        sf::st_nearest_points(coast_hull_buffer) %>% 
        sf::st_cast("POINT") %>% sf::st_sf() %>% 
        dplyr::slice(c(2,4)) %>% 
        sf::st_union() %>% 
        sf::st_cast("LINESTRING")
      
      coast_line <- coast_hull %>% sf::st_cast("LINESTRING") %>% 
        sf::st_difference(sf::st_buffer(sf::st_intersection(sf::st_buffer(fix_line,dist=50),.),dist = 1)) %>%
        sf::st_cast("LINESTRING")
      
      new_coast_line <- function(coast_line, vector_mask) {
        no_cross <- sf::st_crosses(coast_line, vector_mask) %>% 
          purrr::map_lgl(~ length(.x) == 0)
        if(sum(no_cross) > 0) {
          coast_line <- coast_line[no_cross,] %>% 
            dplyr::slice(which.max(sf::st_length(.)))
          return(coast_line)
        }
        select_line <- coast_line %>% 
          sf::st_crosses(vector_mask) %>% 
          purrr::map_dbl(~ st_area(vector_mask[.,]) %>% sum()) %>% 
          which.min()
        coast_line[select_line,]
      }
      
      tmp_mask <- vector_mask %>%  
        dplyr::filter(lengths(sf::st_intersects(., fix_line)) > 0)
      coast_line <- coast_line %>% new_coast_line(tmp_mask)
    }
    # 
    # if (n_polys > 1) {
    #   stop(paste("an on-land segment crosses more than one polygon.",
    #              "this scenario is currently not supported. the best",
    #              "solution is to increase the frequency of prediction",
    #              "times or reduce complexity of the vector mask."))
    # }
    
    start_pt <- start_pt %>% 
      sf::st_sfc() %>%
      sf::st_set_crs(sf::st_crs(crw_sf))
    
    end_pt <- end_pt %>% 
      sf::st_sfc() %>% 
      sf::st_set_crs(sf::st_crs(crw_sf))
    
    coast_points <- coast_line %>% 
      sf::st_sample(l-2,type="regular",offset=0.5)
    if (is.null(coast_points)) {
      NULL
    }
    coast_points <- coast_points %>% sf::st_cast("POINT") 
    
    sample_start <- which.min(sf::st_distance(start_pt, coast_points))
    
    if (sample_start != 1) {
      coast_points <- rev(coast_points)
    }
    
    coast_points <- coast_points %>% 
      sf::st_coordinates() %>%
      rbind(sf::st_coordinates(start_pt),.,sf::st_coordinates(end_pt)) %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(mu.x = X, mu.y = Y) %>% 
      dplyr::bind_cols(wypt_times) %>% 
      dplyr::mutate(type = "wypt") %>% 
      dplyr::mutate(type = case_when(
        row_number() == 1L ~ "start",
        row_number() == l ~ "end",
        TRUE ~ type
      )) %>% 
      dplyr::mutate(nu.x = case_when(
        type == "start" ~ crw_sf$nu.x[start_idx],
        type == "end" ~ crw_sf$nu.x[end_idx],
        TRUE ~ NA_real_
      ),
      nu.y = case_when(
        type == "start" ~ crw_sf$nu.y[start_idx],
        type == "end" ~ crw_sf$nu.y[end_idx],
        TRUE ~ NA_real_
      )) %>% 
      dplyr::full_join(data_times, by = c("times", "type")) %>% 
      dplyr::arrange(times) %>% 
      dplyr::select(type,mu.x,nu.x,mu.y,nu.y,times)
    
    if (!crwIS) {
      segments[[i,"fixed_seg"]] <- path_prediction_fix(coast_points, 
                                                       sigma = exp(par[1]), 
                                                       beta = exp(par[2])) 
    }
    if (crwIS) {
      segments[[i,"fixed_seg"]] <- path_simulation_fix(coast_points, 
                                                       sigma = exp(par[1]), 
                                                       beta = exp(par[2]),
                                                       vector_mask = vector_mask, 
                                                       fix_line=fix_line,
                                                       coast_line=coast_line,
                                                       seg_info = segments[i,]
                                                       #seg_data = crw_sf %>% slice(segments$start_idx[i]:segments$end_idx[i])
      )
    }
  }
  return(segments)
}

#' @title Extract alpha values from \code{crwPredict} or \code{crwIS} objects
#' @description extracts the needed alpha values for \code{fix_path}
#' @param crw_object Can be one of the following classes: 
#' (1) 'crwIS' object from the \code{crwPostIS} function
#' (2) 'crwPredict' object from the \code{crwPredict} function
#' @return a list of 2 (matrix of extracted alpha values, times)
#' @export
#'

crw_alpha <- function(crw_object) {
  ts <- attr(crw_object,"time.scale")
  if (inherits(crw_object,"crwIS")) {
    alpha <- crw_object$alpha.sim
    times <- crw_object$Time
  }
  
  if (inherits(crw_object,"crwPredict")) {
    alpha <- data.matrix(crw_object[,c("mu.x","nu.x","mu.y","nu.y")])
    times <- crw_object[,attr(crw_object,"Time.name")]
  }
  out <- list(alpha = alpha,
              times = as.numeric(times)/ts)
  class(out) <- "crwAlpha"
  return(out)
}

#' @title Project path away from restricted areas
#' @description Corrects a path so that it does not travel through a restricted area.
#' @param crw_object Coordinate locations for the path. Can be one of the following classes: 
#' (1) 'crwIS' object from the \code{crwPostIS} function
#' @param vector_mask an 'sf' polygon object that defines the restricted area
#' @param crwFit crwFit object that was used to generate the crw_object
#' @return a new crw_object w/ locType set to "f" for fixed points
#' @export
#'

fix_path <- function(crw_object, vector_mask, crwFit) {
  # check if crwFit used the drift model and stop
  if (inherits(crwFit,"crwFit_drft")) {
    stop("model fits with drift = TRUE are currently not supported within fix_path.")
  }
  alpha <- crw_alpha(crw_object)
  
  if (inherits(crw_object,"crwIS")) {
    crwIS <- TRUE
  } else {
    crwIS <- FALSE
  }
  # convert crw_object to a POINT sf object
  crw_sf <- crawl::crw_as_sf(crw_object,"POINT")
  
  vector_mask <- vector_mask %>% sf::st_buffer(0) %>% 
    sf::st_collection_extract(type = "POLYGON", warn = FALSE) %>% 
    sf::st_cast("POLYGON", warn = FALSE)
  # %>% 
  #   rmapshaper::ms_clip(bbox = sf::st_bbox(sf::st_buffer(crw_sf,100000)), 
  #                       remove_slivers = TRUE)
  
  fix <- fix_segments(crw_sf = crw_sf, 
                      vector_mask = vector_mask, 
                      crwFit = crwFit, 
                      alpha = alpha, 
                      crwIS = crwIS)
  
  if (inherits(crw_object, "crwIS")) {
    alpha.sim <- crw_object$alpha.sim
    loc.types <- crw_object$locType
    
    for (i in 1:nrow(fix)) {
      start_idx <- fix$start_idx[i]
      end_idx <- fix$end_idx[i] - 1
      alpha.sim[start_idx:end_idx, ] <- as.matrix(fix$fixed_seg[[i]])
      loc.types[start_idx:end_idx] <- "f"
    }
    crw_object$alpha.sim <- alpha.sim
    crw_object$locType <- loc.types
    # crw_object$locType <- crw_object$locType[alpha$predicted_idx]
    # crw_object$Time <- crw_object$Time[alpha$predicted_idx]
    return(crw_object)
  }
  
  if (inherits(crw_object,"crwPredict")) {
    for (i in 1:nrow(fix)) {
      start_idx <- fix$start_idx[i]
      end_idx <- fix$end_idx[i] - 1
      crw_object[start_idx:end_idx,"mu.x"] <- fix$fixed_seg[[i]][,1]
      crw_object[start_idx:end_idx,"nu.x"] <- fix$fixed_seg[[i]][,2]
      crw_object[start_idx:end_idx,"mu.y"] <- fix$fixed_seg[[i]][,3]
      crw_object[start_idx:end_idx,"nu.x"] <- fix$fixed_seg[[i]][,4]
      crw_object[start_idx:end_idx,"locType"] <- "f"
    }
    return(crw_object)
  }
}

