#' @title Identify segments of a path that cross through a restricted area
#' @description This function is used to identify sections of a path that pass through 
#' a restricted area (e.g. for marine mammals or fish, a land mask). the CTCRW model in 
#' crawl cannot actively steer paths away 
#' from land. So, this function
#' will identify path segments from the unrestrained path that pass through these areas.
#' If the path/points end within the land area, those records will be removed.
#' The user can then use this information to adjust the path as desired. 
#' @param crw_object A \code{crwIS} object from the \code{crawl} package 
#' @param vector_mask A \code{sf} object from sf package that indicates 
#' restricted areas as a polygon feature.  
#' @return A data.frame with each row associated with each section of the path
#' that crosses a restricted area. The columns provide the start and end row indices of \code{xy} where
#' the section occurs and the previous and post locations that are in unrestricted space.
#' @author Josh M. London (josh.london@noaa.gov)
#' @export
#' 

get_mask_segments = function(crw_object, vector_mask) {
  # pull out the alpha.sim and Time columns from the crwIS object
  # filter to include only the "p" location types
  alpha.sim <- crw_object$alpha.sim[crw_object$locType == "p",]
  times.sim <- crw_object$Time[crw_object$locType == "p"]
  
  # convert crwIS to a POINT sf object
  xy <- crw_as_sf(crw_object,"POINT","p")
  
  # intersect xy with vector mask
  on_mask <- sf::st_intersects(xy, vector_mask) %>% 
    purrr::map_lgl(~ length(.x) > 0)
  if (any(is.na(on_mask))) {
    stop("points in xy fall outside the extent of vector_mask")
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
  xy <- xy[head_start:tail_end,]
  on_mask <- on_mask[head_start:tail_end]
  
  in.segment <- (on_mask == TRUE)
  
  start_idx <- which(c(FALSE, in.segment) == TRUE &
                       dplyr::lag(c(FALSE, in.segment) == FALSE)) - 2
  end_idx <- which(c(in.segment, FALSE) == TRUE & 
                     dplyr::lead(c(in.segment, FALSE) == FALSE)) + 1
  on_mask_segments <- data.frame(start_idx, end_idx) %>% 
    rowwise() %>% 
    dplyr::mutate(start_alpha = list(alpha.sim[start_idx,]),
                  end_alpha = list(alpha.sim[end_idx,]),
                  times = list(times.sim[start_idx:end_idx]))
  on_mask_segments <- list(
    on_mask_segments = on_mask_segments,
    fixed_range = c(head_start,tail_end)
  )
  
  return(on_mask_segments)
}

#' @title Simulate possible path points conditioned on fit
#' @description Simulates a set of possible points for a given time conditioned on fit
#' @param n integer specifying the number of points to return
#' @param t0 time value for the first location in the segment
#' @param t1 time value for the current location to be simulated
#' @param t2 time value for the last location in the segment
#' @param alpha0 coordinate and velocity values for t0
#' @param alpha2 coordinate and velocity values for t2
#' @param par par values from the crwFit object
#' @param active numeric 1 or 0 whether the animal is moving or not (should almost always = 1)
#' @param inf_fac Variance inflation factor to increase simulation area
#' @param bm boolean whether to draw from a Brownian process
#' @return matrix of coordinate and velocity values drawn
#' @export
#'

cond_sim = function(n=500, t0, alpha0, t2, alpha2, t1, par, active=1, inf_fac=1, bm=0){
  set.seed(42)
  if (inherits(t0,"POSIXct")) {
    t0 <- as.numeric(t0)
  }
  if (inherits(t1,"POSIXct")) {
    t1 <- as.numeric(t1)
  }
  if (inherits(t2, "POSIXct")) {
    t2 <- as.numeric(t2)
  }
  if (bm) {
    mu_cond = c(alpha0 + ((t1 - t0)/(t2 - t0))*(alpha2 - alpha0))[c(1,3)]
    sigma2 = exp(2*par[1])
    V_cond = sigma2*diag(rep((t1 - t0)*((t2 - t0) - (t1 - t0))/(t2 - t0), 2))
    sim = mvtnorm::rmvnorm(n, mu_cond, V_cond)
    smp = cbind(sim[,1], NA, sim[,2],NA)
    smp[,2] = (smp[,1] - alpha0[1])/(t1 - t0)
    smp[,4] = (smp[,3] - alpha0[3])/(t1 - t0)
  } else{
    beta = exp(par[2])
    sigma2 = exp(2*par[1])
    delta = diff(c(t0, t1, t2))
    T0 = makeT(b = beta, delta = delta[1], active = active)
    T1 = makeT(b = beta, delta = delta[2], active = active)
    Q0 = makeQ(b = beta, sig2 = sigma2, delta = delta[1], active = active)
    Q1 = makeQ(b = beta, sig2 = sigma2, delta = delta[2], active = active)
    V_inv = solve(Q0) + t(T1) %*% solve(Q1) %*% T1 
    v = solve(Q0) %*% T0 %*% alpha0 + t(T1) %*% solve(Q1) %*% alpha2
    mu_cond = solve(V_inv, v)
    V_cond = zapsmall(solve(V_inv))
    smp = mvtnorm::rmvnorm(n, mu_cond, inf_fac*V_cond)
  }
  return(smp)
}


#' @title Identify segments of a path that cross a restricted area
#' @description This function takes a crw_object (crwIS only for now) and an 'sf'
#' polygon object that defines the restricted area and identifies each segment
#' of the path that crosses the restricted area. Each segment begins and ends with
#' a coordinate that is outside the restricted area.
#' @param crw_object Coordinate locations for the path. Can be one of the following classes: 
#' (1) 'crwIS' object from the \code{crwPostIS} function
#' @param vector_mask an 'sf' polygon object that defines the restricted area
#' @param crwFit crwFit object that was used to generate the crw_object
#' @return a tibble with each record identifying the segments and pertinant values
# #' @importFrom rlang .data
#' @export

fix_segments <- function(crw_object, vector_mask, crwFit) {
  # get the epsg code of crw_object and par from crwFit
  crw_epsg <- attr(crw_object,"epsg")
  par = tail(crwFit$estPar, 2)
  
  # identify segments that are within the vector mask
  segments <- get_mask_segments(crw_object,vector_mask)
  segments <- segments$on_mask_segments
  # add an empty 'fixed_seg' column to segments
  segments[,"fixed_seg"] <- NA
  
  # begin loop through each segment that needs to be fixed
  for (i in 1:nrow(segments)) {
    message(paste("segment",i))
    # length of the segment is determined from the times column
    l = length(segments[i, ]$times[[1]])
    fixed_seg <- vector("list",l)
    
    # begin loop through each time point in segment
    for (j in 1:(l - 2)) {
      message(paste("pt", j))
      t0 = segments[i, ]$times[[1]][j]
      if (j == 1) {
        alpha0 = segments[i, ]$start_alpha[[1]] %>% as.numeric()
        fixed_seg[[1]] <- alpha0
        fixed_seg[[l]] <- segments[i, ]$end_alpha[[1]] %>% as.numeric()
      }
      t2 = segments[i, ]$times[[1]][l]
      alpha2 = segments[i, ]$end_alpha[[1]] %>% as.numeric()
      t1 = segments[i, ]$times[[1]][j + 1]
      
      # draw points and create lines to each from alpha0
      p_draw <- cond_sim(n = 500, t0, alpha0, t2, alpha2, t1, par)
      colnames(p_draw) <- c("mu.x", "nu.x", "mu.y", "nu.y")
      p_sf <- p_draw %>% as_tibble() %>%
        sf::st_as_sf(coords = c("mu.x", "mu.y")) %>% 
        sf::st_set_crs(crw_epsg)
      alpha0_pt <- matrix(alpha0,ncol = 4)
      colnames(alpha0_pt) <- c("mu.x", "nu.x", "mu.y", "nu.y")
      start_pt <- alpha0_pt %>% as_tibble() %>% 
        sf::st_as_sf(coords = c("mu.x", "mu.y")) %>% 
        sf::st_set_crs(crw_epsg)
      
      l_sf  <- vector("list", nrow(p_sf))
      for (k in 1:nrow(p_sf)) {
        l_sf[[k]] <- sf::st_linestring(rbind(sf::st_coordinates(start_pt),
                                             sf::st_coordinates(p_sf[k,])))
      }
      l_sfc <- sf::st_sfc(l_sf) %>% sf::st_set_crs(crw_epsg)
      
      # now that we have our lines, lets find which ones don't cross the vector mask
      cross_mask <- sf::st_intersects(l_sfc, vector_mask) %>% 
        purrr::map_lgl(~ length(.x) > 0)
      
      # if all lines cross the mask, then do a brownian draw
      if (length(cross_mask[!cross_mask]) == 0) {
        save(list = c("l_sfc", "p_sf", "start_pt", "cross_mask", "p_draw"), file = paste0("seg",i,"_pt",j,".rda"))
        # draw points and create lines to each from alpha0
        p_draw <- cond_sim(n = 500, t0, alpha0, t2, alpha2, t1, par, bm = 1)
        colnames(p_draw) <- c("mu.x", "nu.x", "mu.y", "nu.y")
        p_sf <- p_draw %>% as_tibble() %>%
          sf::st_as_sf(coords = c("mu.x", "mu.y")) %>% 
          sf::st_set_crs(crw_epsg)
        alpha0_pt <- matrix(alpha0,ncol = 4)
        colnames(alpha0_pt) <- c("mu.x", "nu.x", "mu.y", "nu.y")
        start_pt <- alpha0_pt %>% as_tibble() %>% 
          sf::st_as_sf(coords = c("mu.x", "mu.y")) %>% 
          sf::st_set_crs(crw_epsg)
        
        l_sf  <- vector("list", nrow(p_sf))
        for (k in 1:nrow(p_sf)) {
          l_sf[[k]] <- sf::st_linestring(rbind(sf::st_coordinates(start_pt),
                                               sf::st_coordinates(p_sf[k,])))
        }
        l_sfc <- sf::st_sfc(l_sf) %>% sf::st_set_crs(crw_epsg)
        
        # now that we have our brownian lines, lets find which 
        # ones don't cross the vector mask
        cross_mask <- sf::st_intersects(l_sfc, vector_mask) %>% 
          purrr::map_lgl(~ length(.x) > 0)
        
        save(list = c("l_sfc", "p_sf", "start_pt", "cross_mask", "p_draw"), file = paste0("seg",i,"_pt",j,"_brown.rda"))
        
        if (length(cross_mask[!cross_mask]) > 0) {
        idx <- seq_along(cross_mask)[!cross_mask]
        alpha0 <- p_draw[idx[1],] %>% as.numeric()
        
        } else {
          message(paste("brownian draw; segment = ", i,"; point = ",j))
          l_sfc <- sf::st_sf(id = seq_along(l_sfc), geometry = l_sfc)
          pts_intersect_mask <- sf::st_intersects(p_sf, vector_mask) %>% 
            purrr::map_lgl(~ length(.x) > 0)
          l_sfc <- l_sfc[!pts_intersect_mask,]
          
          cross_mask_length <- sf::st_intersection(l_sfc, vector_mask) %>%
            mutate(len = sf::st_length(.)) %>% 
            group_by(id) %>% summarise(total_len = sum(.data$len))
          
          idx <- which.min(cross_mask_length$total_len)
          alpha0 <- p_draw[idx[1],] %>% as.numeric()
        }
        
      } else {
        idx <- seq_along(cross_mask)[!cross_mask]
        alpha0 <- p_draw[idx[1],] %>% as.numeric()
      }
      fixed_seg[[j + 1]] <- alpha0
    }
    segments$fixed_seg[i] <- list(do.call(rbind,fixed_seg))
  }
  return(segments)
}


#' @title Project path away from restricted areas
#' @description Corrects a path so that it does not travel through a restricted area.
#' @param crw_object Coordinate locations for the path. Can be one of the following classes: 
#' (1) 'crwIS' object from the \code{crwPostIS} function
#' @param vector_mask an 'sf' polygon object that defines the restricted area
#' @param crwFit crwFit object that was used to generate the crw_object
#' @return a new crw_object (of type crwIS)
#' @export
#'

fix_path <- function(crw_object, vector_mask, crwFit) {
  vector_mask <- vector_mask %>% sf::st_buffer(0) %>% 
    rmapshaper::ms_clip(bbox = sf::st_bbox(crw_as_sf(crw_object,"POINT") %>% 
    sf::st_buffer(100000)), remove_slivers = TRUE)
  
  fix <- fix_segments(crw_object, vector_mask, crwFit)
  if (inherits(fix, "list")) {return(fix)}
  predicted_idx <- which(crw_object$locType == "p")
  alpha.sim <- crw_object$alpha.sim[predicted_idx,]
  
  for (i in 1:nrow(fix)) {
    start_idx <- fix$start_idx[i]
    end_idx <- fix$end_idx[i]
    alpha.sim[start_idx:end_idx,] <- fix$fixed_seg[[i]]
  }
  crw_object$alpha.sim <- alpha.sim
  crw_object$locType <- crw_object$locType[predicted_idx]
  crw_object$Time <- crw_object$Time[predicted_idx]
  return(crw_object)
  }

