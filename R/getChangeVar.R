#' A function to get the spread or variance of change for a SHAPE trace.
#' @title getChangeVar
#' @aliases getChange Var
#' @keywords trace variance change RNA
#' @usage getChangeVar(sample, base=sample[1,], margin=1, tol=0.1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @export
#' @details This function calculates the variance of change positions for each row (or column) in sample.
#' @return A numeric vector of trace change variances.
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get change variance
#' cv = getChangeVar(samp_nreduce)
#'
getChangeVar = function(sample, base=sample[1,], margin=1, tol=0.1){
  
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
  
  #values cannot be negative (raise the minimum value to at least 0)
  if(sum(sample<0, na.rm=T)>0){
    sample = sample-min(sample, na.rm=T) 
    base = base-min(sample, na.rm=T) 
  }
  if(sum(base<0, na.rm=T)>0){
    base = base-min(base, na.rm=T)  
    sample = sample-min(base, na.rm=T) 
  }
  
  #function to get change variance
  changeSpread = function(samp, base, tol){
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
    
    #get change spread
    changeVar = 0
    if(length(which(change!=0))>0){
      changeVar = var(change)
    }
    
    return(changeVar)
  }
  
  #change variance
  chgvar = apply(sample, 1, changeSpread, base=base, tol=tol)
  
  #return change variance
  return(chgvar)
}
