
erf = function(x){
  return(2*pnorm(x*sqrt(2))-1)
}
erfinv = function(z){
  if(-1>=z | z>=1) stop("Correlation outside of (-1,1)\n")
  qnorm((z+1)/2)/sqrt(2)
}
                  
makeAvail <- function(i, Tmat, Qmat, predx, predy, vary, varx, driftMod, lonadj){
  .T <- matrix(0, 2+driftMod,2+driftMod)
  .T[1,1] <- 1
  .T[1,2] <- Tmat[i,1]
  .T[2,2] <- Tmat[i,2]
  if(driftMod){
    .T[1,3] <- Tmat[i,3]
    .T[3,3] <- Tmat[i,4]
  }
  Sy <- Qmat[i,1]
  Sx <- Qmat[i,1]/(lonadj[i]^2)
  ax <- as.double(predx[i,])
  ay <- as.double(predy[i,])
  c(c(.T%*%ax)[1], c(.T%*%ay)[1], (.T%*%varx[,,i]%*%t(.T))[1,1] + Sx, (.T%*%vary[,,i]%*%t(.T))[1,1] + Sy)
}

rmvtt <- function(mu, Sigma, df=Inf, lower, upper){
	p <- length(mu)
	out <- rep(NA,p)
	div <- ifelse(df==Inf, 1, sqrt(rchisq(1,df=df)/df))
	truncLo <- pnorm(div*(lower[1]-mu[1])/sqrt(Sigma[1,1]))
	truncUp <- pnorm((div*upper[1]-mu[1])/sqrt(Sigma[1,1]))
	out[1] <- mu[1] + sqrt(Sigma[1,1])*(qnorm(runif(1, truncLo, truncUp))/div)
	if(p>1){	
		for(i in 2:p){
			S12 <- Sigma[i,1:(i-1)]
			S22 <- Sigma[1:(i-1),1:(i-1)]
			res <- out[1:(i-1)]-mu[1:(i-1)]
			mu.c <- mu[i] + S12%*%solve(S22, res)
			S.c <- Sigma[i,i] - S12%*%solve(S22,S12)
			truncLo <- pnorm(div*(lower[i]-mu.c)/sqrt(S.c))
			truncUp <- pnorm((div*upper[i]-mu.c)/sqrt(S.c))
			out[i] <- mu.c + sqrt(S.c)*(qnorm(runif(1, truncLo, truncUp))/div)
		}
		return(out)
	}
	else return(out)
}

getSD <- function(x){
	d <- as.numeric(sapply(strsplit(as.character(x),"e-"), function(x) x[2]))
	if(any(!is.na(d))) return(max(d, na.rm=TRUE))
	else return(0)
}

check_fit <- function(mle) {
  checkMLE <- inherits(mle, 'try-error')
  checkConv <-
    ifelse(inherits(mle, 'try-error'), 1, mle$convergence > 0)
  
  C.tmp <- try(2 * solve(mle$hessian), silent = TRUE)
  if (inherits(C.tmp, "try-error")) {
    checkCovar <- 1
    checkDiag <- 1
  } else {
    checkCovar <- 0
    checkDiag <- ifelse(any(diag(C.tmp) <= 0), 1, 0)
  }
  return(sum(checkMLE, checkConv, checkCovar, checkDiag) > 0)
}

sfc_as_cols <- function(x, geometry, names = c("x","y")) {
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  } else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

# conf_elps <- function(x, y, V, prob=0.95){
#   require(sf, quietly = TRUE)
#   v <- eigen(V)$vectors
#   lambda <- eigen(V)$values
#   t <- seq(0,2*pi, length=50)
#   m <- qchisq(prob,2)
#   x_l <- x + sqrt(m*lambda[1])*v[1,1]*cos(t) + sqrt(m*lambda[2])*v[1,2]*sin(t)
#   y_l <- y + sqrt(m*lambda[1])*v[2,1]*cos(t) + sqrt(m*lambda[2])*v[2,2]*sin(t)
#   ell <- cbind(x=x_l, y=y_l)
#   ell <- st_linestring(ell) |> st_cast("POLYGON") 
#   return(st_sfc(list(ell)))
# }

