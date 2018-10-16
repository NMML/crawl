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
  on_mask <- on_mask[head_start:tail_end]
  
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

fix_segments <- function(crw_sf, vector_mask, crwFit, alpha) {
  # get par from crwFit
  par <- tail(crwFit$estPar, 2)
  ts <- attr(crwFit,"time.scale")
  
  # identify segments that are within the vector mask
  segments <- get_mask_segments(crw_sf,vector_mask, alpha)
  segments <- segments$on_mask_segments
  # add an empty 'fixed_seg' column to segments
  segments[,"fixed_seg"] <- NA
  
  # begin loop through each segment that needs to be fixed
  for (i in 1:nrow(segments)) {
    start_idx <- segments[i,"start_idx"][[1]]
    end_idx <- segments[i,"end_idx"][[1]]
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
    
    fixed_seg <- vector("list",l)
    
    start_pt = sf::st_point(segments[i, ]$start_alpha[[1]][c("mu.x","mu.y")]) 
    end_pt = sf::st_point(segments[i, ]$end_alpha[[1]][c("mu.x","mu.y")])
    fix_line <- sf::st_linestring(rbind(start_pt,end_pt)) %>% 
      sf::st_sfc() %>% sf::st_set_crs(sf::st_crs(crw_sf))
    
    n_polys <- vector_mask %>% 
      sf::st_cast("POLYGON") %>% 
      st_intersects(fix_line,sparse=FALSE) %>% 
      sum()
    
    coast_line <- vector_mask %>% 
      sf::st_cast("POLYGON") %>% 
      dplyr::filter(lengths(st_intersects(., fix_line)) > 0) %>% 
      lwgeom::st_split(fix_line) %>% 
      sf::st_collection_extract("POLYGON") %>% 
      dplyr::slice(which.min(st_area(.))) %>% 
      sf::st_cast("LINESTRING") %>% 
      lwgeom::st_split(st_buffer(fix_line,dist = 1e-10)) %>% 
      sf::st_collection_extract("LINESTRING") %>% 
      dplyr::slice(which.max(st_length(.)))
    
    start_pt <- start_pt %>% 
      sf::st_sfc() %>%
      sf::st_set_crs(sf::st_crs(crw_sf))
    
    end_pt <- end_pt %>% 
      sf::st_sfc() %>% 
      sf::st_set_crs(sf::st_crs(crw_sf))
    
    coast_points <- coast_line %>% 
      sf::st_sample(l-2,type="regular",offset=runif(1)) %>% 
      sf::st_cast("POINT") %>% 
      sf::st_coordinates() %>%
      rbind(sf::st_coordinates(start_pt),.,st_coordinates(end_pt)) %>% 
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
      dplyr::full_join(data_times) %>% 
      dplyr::arrange(times) %>% 
      dplyr::select(type,mu.x,nu.x,mu.y,nu.y,times)
  }
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
  
  # convert crw_object to a POINT sf object
  crw_sf <- crawl::crw_as_sf(crw_object,"POINT")
  
  vector_mask <- vector_mask %>% sf::st_buffer(0) %>% 
    rmapshaper::ms_clip(bbox = sf::st_bbox(sf::st_buffer(crw_sf,100000)), 
                        remove_slivers = TRUE)
  
  fix <- fix_segments(crw_sf, vector_mask, crwFit, alpha)
  if (inherits(fix, "list")) {return(fix)}
  
  if (inherits(crw_object, "crwIS")) {
    alpha.sim <- crw_object$alpha.sim
    loc.types <- crw_object$locType
    for (i in 1:nrow(fix)) {
      start_idx <- fix$start_idx[i]
      end_idx <- fix$end_idx[i]
      alpha.sim[start_idx:end_idx, ] <- fix$fixed_seg[[i]]
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
      end_idx <- fix$end_idx[i]
      crw_object[start_idx:end_idx,"mu.x"] <- fix$fixed_seg[[i]][,1]
      crw_object[start_idx:end_idx,"nu.x"] <- fix$fixed_seg[[i]][,2]
      crw_object[start_idx:end_idx,"mu.y"] <- fix$fixed_seg[[i]][,3]
      crw_object[start_idx:end_idx,"nu.x"] <- fix$fixed_seg[[i]][,4]
      crw_object[start_idx:end_idx,"locType"] <- "f"
    }
    return(crw_object)
  }
  }

