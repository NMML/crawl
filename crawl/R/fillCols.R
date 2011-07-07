"fillCols" <- function(data) {
   nc <- ncol(data)
   getConst <- function(vec) {
      vals <- unique(vec)
      return(length(vals[!is.na(vals)])==1)
   }
   constCol <- apply(data, 2, getConst)
   data[,constCol] <- data[1,constCol]
   return(data)
}