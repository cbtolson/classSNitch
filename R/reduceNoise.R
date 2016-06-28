#' A function to reduce noise in SHAPE data. This function removes peaks that are high in both the base and the comparison SHAPE traces.
#' @title reduceNoise
#' @aliases reduceNoise
#' @keywords reduce noise trim peak filter SHAPE RNA structure 
#' @usage reduceNoise(sample, base=sample[1,], margin=1, trim=0, high=boxplot(sample)$stats[4])
#' @param sample A numeric matrix containing reactivity scores to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the reactivity score to which the samples are to be compared (e.g. a wild type SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param trim An optional number indicating the number of nucleotides to be trimmed from the each end. Default is 0.
#' @param high An optional number indicating the value above which reactivities are considered high. Default is third quartile of sample.
#' @export
#' @details This function reduces the noise in SHAPE data. For positions where both the base vector and the sample row (or column) is above the high value, the position in the sample row (or column) is set equal to that position in the base vector. The function trims the data by setting the ends of sample equal to the ends of the base vector.
#' @return A noise reduced numeric matrix with the same dimensions as sample.
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#'
reduceNoise = function(sample, base=sample[1,], margin=1, trim=0, high=boxplot(sample)$stats[4]){
  
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
  
  #set optional paramater trim
  if(missing(trim)){
    trim = 0
  } else {
    if(trim < 0 || trim > dim(sample)[2]){
      warning("Trim value not valid. Trim set to default.")
      trim = 0
    } else if(trim > 0){
      trim = trim-1 
    }
  }
  
  #set optional paramater high
  if(missing(high)){
    high = boxplot(as.numeric(unlist(sample)), plot=F)$stats[4]
  } else {
    if(high < 0){
      warning("High value not valid. High set to default.")
      high = boxplot(as.numeric(unlist(sample)), plot=F)$stats[4]
    }
    high = high
  }
  
  #function to filter mutually high peaks from sample
  highPeakFilter = function(samp, base, high){
    inds = intersect(which(samp > high, arr.ind=T),  which(base > high, arr.ind=T))
    samp[inds] = base[inds]
    return(samp)
  }
  
  #high peak filter
  if(dim(sample)[1]==1 || dim(sample)[2]==1){
    samp_qual = highPeakFilter(sample, base=base, high=high)
  } else{
    samp_qual = do.call(rbind, apply(sample, 1, highPeakFilter, base=as.list(base), high=high))  
  }
  
  
  #low peak filter
  samp_qual[sample<(-0.5)] = 0
  base[base<(-0.5)] = 0
  
  #function to trim ends from sample
  trimEnds = function(samp, base, trim){
    samp[((length(samp)-trim):length(samp))] = base[((length(base)-trim):length(base))]
    samp[1:trim] = base[1:trim]
    return(samp)
  }
  
  #trim ends
  if(trim >= 0){
    if(dim(sample)[1]==1){
      samp_qual = trimEnds(samp_qual, base=base, trim=trim)
    } else{
      samp_qual = do.call(rbind, apply(samp_qual, 1, trimEnds, base=as.list(base), trim=trim))
    }
  }
  
  #organize output data
  if(margin==2){
    samp_qual = t(samp_qual)
  }
  samp_qual = matrix(unlist(samp_qual), nrow=nrow(samp_qual))
  
  #return noise reduced sample
  return(samp_qual)
}
