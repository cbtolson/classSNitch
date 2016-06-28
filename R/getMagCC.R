#' A function to get the correlation coefficient between SHAPE trace magnitudes.
#' @title getMagCC
#' @aliases getMagCC
#' @keywords magnitude correlation-coefficient RNA
#' @usage getMagCC(sample, base=sample[1,], margin=1)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wildtype SHAPE trace). Default is the first trace in sample.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @export
#' @details This function calculates the Pearson correlation coefficient between the base vector and each row (or column) in sample.
#' @return A numeric vector of correlation coefficients.
#' @author Chanin Tolson
#' @seealso  \code{\link{getFeatures}} 
#' @examples #sample data
#' sample = matrix(sample(1:100), ncol=10)
#' #normalize
#' samp_norm = normalize(sample)
#' #reduce noise
#' samp_nreduce = reduceNoise(samp_norm, trim=1, high=4)
#' #get magnitude correlation coefficient
#' mag = getMagCC(samp_nreduce)
#'
getMagCC = function(sample, base=sample[1,], margin=1){
  
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
  
  #calculate magnitude change
  if(dim(sample)[1]==1){
    mag = cor(as.vector(sample), base, method="pearson", use="pairwise.complete.obs")
  } else{
    mag = apply(sample, 1, cor, x=base, method="pearson", use="pairwise.complete.obs")
  }
  
  #return magnitude change
  return(mag)
}
