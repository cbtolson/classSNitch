#' A between-sample normalization function for SHAPE traces.
#' @title normalize
#' @aliases normalize
#' @keywords normalize between-sample RNA
#' @usage normalize(sample, base=sample[1,], margin=1, outbase=FALSE, mean=1.5)
#' @param sample A numeric matrix containing values to be normalized (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the sample is to be normalized (e.g. a wild type SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if sample is organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param outbase An optional boolean indicating if the normalized base should be returned. Default is FALSE.
#' @param mean An optional number setting the mean SHAPE value. Default is 1.5.
#' @export
#' @details This function normalizes the average value of the base vector. Each row (or column) in sample is then normalized by minimizing the absolute difference between the base and the sample row (or column). 
#' @return 
#' \describe{
#'  \item{"samp_norm"}{A normalized numeric matrix with the same dimensions as sample.} 
#'  \item{"samp_norm"}{An optional list with two elements: normalized numeric matrix with the same dimensions as sample and a normalized vector the same length as base.}
#' }
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#'            
#' #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' base = sample(1:100, size=10)
#' #normalize
#' samp_norm = normalize(sample, base)
#'
normalize = function(sample, base=sample[1,], margin=1, outbase=FALSE, mean=1.5){
  
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
  
  #set optional paramater outbase
  if(missing(base)) {
    outbase = FALSE
  } else {
    outbase = outbase
  }
  
  #set optional paramater mean
  if(missing(mean)) {
    mean = 1.5
  } else {
    mean = mean
  }
  
  #remove large negative values
  base[base<(-0.5)] = NA
  sample[sample<(-0.5)] = NA
  
  #set average reactivity for wild-type 
  base = (mean*sum(!is.na(base))/sum(base, na.rm=T))*base
  base[base<(-0.5)] = NA
  
  #function to optimize difference between wild-type and sample
  optimizeNorm = function (x, samp, base){
    sum(abs(base-x*samp), na.rm=T)
  }
  
  #optimize each sample
  if(dim(sample)[1]==1){
    opt = optimize(f=optimizeNorm, interval=c(0,1), tol=0.0001, base=base, samp=sample)
    samp_norm = sample*opt$minimum
  } else{
    opt = apply(sample, 1, optimize, f=optimizeNorm, interval=c(0,1), tol=0.0001, base=base)
    samp_norm = do.call(rbind, lapply(1:length(opt), function(i){sample[i,]*opt[[i]]$minimum}))
  }
  
  #organize output data
  if(margin==2){
    samp_norm = t(samp_norm)
  }
  samp_norm = matrix(unlist(samp_norm), nrow=nrow(samp_norm))
  samp_norm[samp_norm<(-0.5)] = NA
  
  if(outbase==TRUE){
    samp_out = NULL
    samp_out[[1]] = samp_norm
    samp_out[[2]] = base
    samp_norm = samp_out
    names(samp_norm) = c("sample.normalized", "base.normalized")
  }
  
  #return normalized sample
  return(samp_norm)
}

