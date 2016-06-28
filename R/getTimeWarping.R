#' A function to get the dynamic time warping between SHAPE trace magnitudes.
#' @title getTimeWarping
#' @aliases getTimeWarping
#' @keywords getTimeWarping change RNA
#' @usage getTimeWarping(sample, base=sample[1,], margin=1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @export
#' @import dtw
#' @details This function calculates the dynamic time warping between the base vector and each row (or column) in sample.
#' @return A numeric vector of time warping values.
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}}  
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get time warping
#' tw = getTimeWarping(samp_nreduce)
#'
getTimeWarping = function(sample, base=sample[1,], margin=1){
  
  #set optional paramater margin
  if(missing(margin)) {
    margin = 1
  } else {
    if(!(margin %in% c(1,2))){
      warning("Margin value not valid. Margin set to default.")
      margin = 1
    }
    if(margin==2){
      sample = t(sample)
    }
  }
  
  #set optional paramater base
  if(missing(base)) {
    base = sample[1,]
  } else {    
    base = base
  }
  
  #remove NA values
  s = sample
  b = base
  base[which(is.na(sample),arr.ind=T)[1]]=0
  sample[,is.na(base)]=0
  sample[is.na(sample)]=0
  base[is.na(base)]=0
  
  #calculate timewarp change
  timewarp = function(x, y){
    warp = dtw(x, y, distance.only=T)
    
    return(warp$normalizedDistance)
  }
  if(dim(sample)[1]==1){
    tw = dtw(as.vector(sample[1,]), as.vector(base))$normalizedDistance
  } else{ 
    tw = apply(sample, 1, timewarp, y=base)  
  }
  
  #return time warping
  return(tw)
}
