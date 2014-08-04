#' @export 

argosDiag2Cov = function(Major, Minor, Orientation){
  a=Major
  b=Minor
  theta=Orientation
  if(any(theta<0 | theta>180)) stop("Argos diagnostic data orientation outside of [0,180]!")
  if(any(a < 0)) stop("Argos diagnostic data major axis < 0!")
  if(any(b<0)) stop("Argos diagnostic data minor axis < 0!")
  theta = pi*(theta/180)
  k=sqrt(2)
  v1 = (a/k)^2*sin(theta)^2 + (b/k)^2*cos(theta)^2
  v2 = (a/k)^2*cos(theta)^2 + (b/k)^2*sin(theta)^2
  c12 = ((a^2 - b^2)/k^2)*cos(theta)*sin(theta)
  rho = c12/(sqrt(v1)*sqrt(v2))
  if(any(rho > 1 | rho < -1)) stop("Faulty Argos error correlation calculated from 'argosDiag2Cov' function")
  return(data.frame(ln.sd.x=log(sqrt(v1)), ln.sd.y=log(sqrt(v2)), error.corr=rho))
}