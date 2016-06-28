#' A function to build a random forest classifier for RNA structure change
#' @title classifyRNA
#' @aliases classifyRNA
#' @keywords classifier RNA structure change random forest
#' @usage classifyRNA(data=NULL, cutoff=NULL)
#' @param data Optional data to build the classifier. Default is pre-loaded data.
#' @param cutoff An optional vector of length equal to number of classes. The winning class for an observation is the one with the maximum ratio of proportion of votes to cutoff. Default is 1/k where k is the number of classes (i.e., majority vote wins).
#' @export
#' @import randomForest
#' @details This function builds a random forest classifier for RNA structure change using the randomForest package.
#' @return A classifyRNA object, based on randomForest object (see randomForest package)
#' \describe{
#'  \item{call}{The original call to randomForest}
#'  \item{type}{One of regression, classification, or unsupervised.}
#'  \item{predicted}{The predicted values of the input data based on out-of-bag samples.}
#'  \item{importance}{A matrix with nclass + 2 columns. The first nclass columns are the class-specific measures computed as mean descrease in accuracy. The nclass + 1st column is the mean descrease in accuracy over all classes. The last column is the mean decrease in Gini index.}
#'  \item{importanceSD}{The standard errors of the permutation-based importance measure. A p by nclass + 1 matrix corresponding to the first nclass + 1 columns of the importance matrix.}
#'  \item{ntree}{Number of trees grown.}
#'  \item{mtry}{Number of predictors sampled for spliting at each node.}
#'  \item{forest}{A list that contains the entire forest}
#'  \item{err.rate}{Vector error rates of the prediction on the input data, the i-th element being the (OOB) error rate for all trees up to the i-th.}
#'  \item{confusion}{The confusion matrix of the prediction (based on OOB data).}
#'  \item{votes}{A matrix with one row for each input data point and one column for each class, giving the fraction or number of (OOB) votes from the random forest.}
#'  \item{oob.times}{Number of times cases are out-of-bag (and thus used in computing OOB error estimate)}
#'  \item{proximity}{A matrix of proximity measures among the input (based on the frequency that pairs of data points are in the same terminal nodes).}
#' }
#' @note Organization of the data file: header=TRUE, tab-delimited .txt file
#' \itemize{
#'  \item{"column 1"}{ class label} 
#'  \item{"column 2"}{ pattern correlation} 
#'  \item{"column 3"}{ dynamic time warping} 
#'  \item{"column 4"}{ contiguousness} 
#'  \item{"column 5"}{ magnitude correlation} 
#'  \item{"column 6"}{ change variance} 
#'  \item{"column 7"}{ eSDC} 
#'  \item{"column 8"}{ change range} 
#'  \item{"column 9"}{ L2 norm} 
#' }
#' The default data has been gathered from the RNA Mapping Database mutate and map experiments.
#' @author Chanin Tolson
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package) \cr\cr
#' \href{http://rmdb.stanford.edu/}{RNA Mapping Database}
#' @seealso  \code{\link{getFeatures}}
#' @examples
#' #build classifier
#' rf = classifyRNA()
#' #get confusion matrix
#' rf$confusion
#' 
classifyRNA = function(data=NULL, cutoff=NULL){
  
  #set optional paramater data
  if(missing(data)){
    data = classSNitch::classify_default
    responses = data[,1]
    input = data[,2:9]
  } else {
    data = data
    if(ncol(data) != 9){
      stop("Incorrect data file format.")
    }
    responses = data[,1]
    input = data[,2:9]
  }
  
  #set optional paramater cutoff
  if(!missing(cutoff)){
    if(sum(is.na(responses))>0){
      if(length(cutoff)!=length(unique(responses[-which(is.na(responses))]))){
        warning("cutoff set to default.")
        cutoff = rep(1/length(unique(responses[-which(is.na(responses))])), length(unique(responses[-which(is.na(responses))])))
      }
    } else {
      if(length(cutoff)!=length(unique(responses))){
        warning("cutoff set to default.")
        cutoff = rep(1/length(unique(responses)), length(unique(responses)))
      }
    }
    cutoff = cutoff
  } else {
    if(sum(is.na(responses))>0){
      cutoff = rep(1/length(unique(responses[-which(is.na(responses))])), length(unique(responses[-which(is.na(responses))])))
    } else{
      cutoff = rep(1/length(unique(responses)), length(unique(responses)))  
    }
  }

  #get features
  input = as.data.frame(cbind(responses, input))
  rownames(input) = rownames(data)
  colnames(input) = c("class", "pat", "tw", "contig", "mag", "var", "eSDC", "range", "l2norm") 
  if(sum(is.na(input[,1]), na.rm=T)>0){
    input = input[-which(is.na(input[,1]),arr.ind=T),]
  }
  input[,1] = factor(as.numeric(input[,1]))
  
  #random forest classification
  rf = randomForest(class~., data=input, importance=TRUE, proximity=TRUE, ntree=5001, cutoff=cutoff)
  
  #convert to a classifyRNA object
  cr = structure(rf, class = "classifyRNA")
  
  #return classifyRNA object
  return(cr)
}

