#' A function to get the features for describing RNA structure change. These features can be used in classification of RNA structure change.
#' @title getFeatures
#' @aliases getFeatures
#' @keywords change features RNA structure
#' @usage getFeatures(sample, base=NULL, margin=1,  norm=T, noise=T, mean=1.5, trim=0, high=NULL,
#'    tol=0.1, outfile=NULL, append=F)
#' @param sample A numeric matrix containing values to be compared (e.g. a set of mutant SHAPE traces).
#' @param base An optional numeric vector containing the values to which the samples are to be compared (e.g. a wild type SHAPE trace). Default is the first trace in each file.
#' @param margin An optional number indicating if the samples are organized by rows or columns, where 1 indicates rows and 2 indicates columns. Default is 1.
#' @param norm An optional boolean to normalize the sample. Default is TRUE. 
#' @param noise An optional boolean to reduce noise in the sample. Default is TRUE. 
#' @param mean An optional number setting the mean SHAPE value for normalization. Default is 1.5.
#' @param trim An optional number indicating the number of nucleotides to be trimed from the ends. Default is 0.
#' @param high An optional number indicating the reactivity above which reactivities are considered high. Default is third quartile of the sample in each file.
#' @param tol An optional number indicating the tolerance for the change. Default is 0.1.
#' @param outfile An optional string indicating the name of the output file. The output file will consist of feature columns. Default will not output a file.
#' @param append An optional boolean to append the file if an outfile is given. Default is FALSE. 
#' @export
#' @details This function calculates the  pattern correlation coefficient, dynamic time warping, contiguousness of change, magnitude correlation coefficient, change of variance, experimental structural disruption coefficient, change range and L2 norm. These are features used to describe structure change in SHAPE traces. 
#' @return 
#' \describe{
#'  \item{"outmat"}{A seven column numeric matrix for pattern correlation coefficient, dynamic time warping, contiguousness of change, magnitude correlation coefficient, change of variance, experimental structural disruption coefficient, change range, and L2 norm} 
#'  \item{"outfile"}{An optional output file for the matrix.}
#' }
#' @author Chanin Tolson
#' @seealso  \code{\link{normalize}} \code{\link{reduceNoise}} \code{\link{getPatternCC}} \code{\link{getTimeWarping}} \code{\link{getContiguous}} \code{\link{getMagCC}} \code{\link{getESDC}} \code{\link{getChangeRange}} \code{\link{getL2norm}}
#' @examples #input files
#' data("shape_ex")
#' #get features
#' params = getFeatures(shape_ex, trim=5, outfile="out.txt")
#'
getFeatures = function(sample, base=NULL, margin=1, norm=T, noise=T, mean=1.5, trim=0, high=NULL, tol=0.1, outfile=NULL, append=F){
  
  #set sample parameter
  sample = as.matrix(sample)
  if(dim(sample)[1]==1 || dim(sample)[2]==1){
    sample = t(sample)
  }
  
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
  
  #set optional paramater norm
  if(missing(norm)){
    norm = TRUE
  } else {
    norm = norm
  }
  
  #set optional paramater noise
  if(missing(noise)){
    noise = TRUE
  } else {
    noise = noise
  }
  
  #set optional paramater mean
  if(missing(mean)) {
    mean = 1.5
  } else {
    mean = mean
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
  
  #set optional paramater tol
  if(missing(tol)){
    tol = 0.1
  } else {
    if(tol < 0){
      warning("Tol value not valid. Tol set to default.")
      tol = 0.1
    }
    tol = tol
  }
  
  #set optional paramater append
  if(missing(append)){
    append = F
  } else {
    if(!(append %in% c(TRUE, FALSE))){
      warning("Append value not valid. Append set to default.")
      append = F
    }
    append = append
  }

  #set optional paramater normalize
  if(norm == TRUE){
    #normalize
    samp_norm = normalize(sample, base, mean=mean)
    base = (1.5*length(base)/sum(base, na.rm=T))*base
    base[base<(-0.5)] = 0  
  } else {
    samp_norm = sample
  }
  
  #set optional paramater high
  if(missing(high)){
    high = boxplot(as.numeric(unlist(samp_norm)), plot=F)$stats[4]
  } else {
    if(high < 0){
      warning("High value not valid. High set to default.")
      high = boxplot(as.numeric(unlist(samp_norm)), plot=F)$stats[4]
    }
    high = high
  }
  
  #reduce noise
  if(noise == TRUE){
    samp_qual = reduceNoise(samp_norm, base, trim=trim, high=high)
  } else{
    samp_qual = samp_norm
  }
      
  #pattern correlation coefficient
  pat = getPatternCC(samp_qual, base, tol=tol)
  
  #dynamic time warping
  tw = getTimeWarping(samp_qual, base)
  
  #change contiguousness
  contig = getContiguous(samp_qual, base=base, tol=tol)
  
  #magnitude correlation coefficient
  mag = getMagCC(samp_qual, base)
  
  #change variance
  var = getChangeVar(samp_qual, base=base, tol=tol)
  
  #experimental structural disruption coefficient
  eSDC = getESDC(samp_qual, base)
  
  #change range
  range = getChangeRange(samp_qual, base=base, tol=tol)
  
  #L2 norm
  l2norm = getL2norm(samp_qual, base)
  
  #combine features
  features = cbind(pat, tw, contig, mag, var, eSDC, range, l2norm)
  rownames(features) = rownames(sample)
  colnames(features) = c("pattern change", "time warping", "change contiguousness", "magnitude change", "change variance", "eSDC", "range", "l2 norm")
  
  #write parameter outfile
  if(missing(outfile)){
  } else{
    if(append==F){
      write.table(features, outfile, quote=F, sep="\t", row.names=T, col.names=T)
    } else{
      write.table(features, outfile, quote=F, sep="\t", row.names=T, col.names=F, append=T)
    }
  }
  
  #return features
  return(features)
}