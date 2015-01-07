#' @export

crwMLE.argos = function(mov.model=~1, quality=NULL, diagnostic=NULL, activity=NULL, drift=FALSE, data, 
                        coord=c("x", "y"), Time.name, ...)
{
  if(drift) stop("At this time the drift model is not supported.")
  if(is.null(quality) & is.null(diagnostic)) stop("At least one of 'quality' or 'diagnostic' must be specified with a formula!")
  if(!is.null(diagnostic)){
    diaMat = model.matrix(diagnostic, data, na.action=na.pass)[,-1]
    if(inherits(data, "SpatialPoints")){
      data@data = cbind(data@data, argosDiag2Cov(diaMat[,1], diaMat[,2], diaMat[,3]))
    } else{
      data = cbind(data, argosDiag2Cov(diaMat[,1], diaMat[,2], diaMat[,3]))
    }
    x.err  = ~ln.sd.x-1
    y.err = ~ln.sd.y-1
    rho = ~error.corr
    npar = ncol(model.matrix(mov.model, data))
    if(!drift){
      fixPar = c(rep(NA,4),1,1,rep(NA,npar),rep(NA,npar))
      if(!is.null(activity)) fixPar = c(fixPar, NA)
    } else{
      fixPar = c(rep(NA,6),1,1,rep(NA,npar),rep(NA,npar),NA,NA)
      if(!is.null(activity)) fixPar = c(fixPar, NA)
    }
  } else{
    qual = model.frame(quality, data, na.action=na.pass)
  }
  out=crwMLE(
    mov.model=mov.model, err.model=list(x=x.err, y=y.err, rho=rho), drift=drift, activity=activity, data=data, 
    coord=coord, Time.name=Time.name, fixPar=fixPar,...)
  return(out)
}