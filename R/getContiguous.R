#' A function to get the number of stretches of contiguous change.
#' @title getContiguous
#' @aliases getContiguous
#' @keywords trace getContiguous RNA
#' @usage getContiguous(sample, base=sample[1,], margin=1, tol=0.1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @export
#' @details This function calculates the number of stretches of contiguous change between the base vector and each row (or column) in sample.
#' @return A numeric vector of counts
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get trace difference
#' contig = getContiguous(samp_norm)
getContiguous = function(sample, base=sample[1,], margin=1, tol=0.1){
  
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
  
  #set optional paramater  tol
  if(missing(tol)){
    tol = 0.1
  } else {
    if(tol < 0){
      warning("Tol value not valid. Tol set to default.")
      tol = 0.1
    }
    tol = tol
  }
  
  #function to get the location of pattern changes
  contiguous = function(samp, base, tol){
    #initialize
    samp = as.numeric(samp)
    base = as.numeric(base)
    patterns = rep(0, length(samp)-1)
    patternb = rep(0, length(base)-1)
    
    #get samp pattern
    pat = (c(samp, NA) - c(NA, samp))[2:(length(samp))]
    patterns[pat>tol] = 1
    patterns[pat<(-tol)] = -1
    
    #get base pattern
    pat = (c(base, NA) - c(NA, base))[2:(length(base))]
    patternb[pat>tol] = 1
    patternb[pat<(-tol)] = -1
    
    #get change
    change = patterns-patternb
    
    #get range of change
    contChange = 0
    if(length(which(change!=0))>0){
      contChange = length(c(1, which(diff(change) != 1 & diff(change) != 0) + 1))
    }
    
    return(contChange)
  }
  
  #get contiguousness
  contig = t(apply(sample, 1, contiguous, base=base, tol=tol))
  contig = as.vector(contig)
  
  #return contiguousness
  return(contig)
}
