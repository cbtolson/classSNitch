#' A function to classify rna structure change from an existing classifier.
#' @title predict.classifyRNA
#' @aliases predict
#' @usage
#' \method{predict}{classifyRNA}(object, sample = NULL, resp="prob", ...)
#' @keywords predict prediction RNA structure change
#' @param object An object of classifyRNA (see classifyRNA function).
#' @param sample An optional matrix of predictors (e.g. output from getFeatures())
#' @param resp An optional string to determine type of return value. Default is "prob".
#' @param ... additional arguments for predict method
#' @export
#' @import randomForest
#' @details This function predicts RNA structure change in SHAPE data using a random forest classifier. 
#' @return A matrix of "response", "vote" or "prob" predictions.
#' \describe{
#'  \item{"response"}{ Predicted classes (classes with majority vote)}
#'  \item{"vote"}{ Vote count fraction (one column for each class and one row for each input)} 
#'  \item{"prob"}{ Class probabilities (one column for each class and one row for each input)}
#' }
#' @author Chanin Tolson
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package)
#' @seealso  \code{\link{getFeatures}} \code{\link{classifyRNA}}
#' @examples #input data
#' data("mutmap")
#' #build classifier
#' cr = classifyRNA()
#' #get prediction
#' cr_pred = predict(cr, mutmap[,2:9])
#'
predict.classifyRNA = function(object, sample=NULL, resp="prob", ...){
  
  #check argument object
  if(attr(object, "class") == "classifyRNA"){
    rf = structure(object, class = "randomForest")
  } else {
    stop("Object must be a classifyRNA object")
  }
  
  #get optional argument resp
  if(missing(resp)){
    resp = "prob"
  } else {
    if(!(resp %in% c("response", "vote", "prob"))){
      warning("resp set to default.")
    } 
    resp = resp
  }
  
  #get optional argument sample
  rf_pred = NULL
  if(!missing(sample)){
    sample = sample[,1:8]
    colnames(sample) = c("pat", "tw", "contig", "mag", "var", "eSDC", "range", "l2norm") 
    if(resp =="response"){
      rf_pred = predict(rf, sample, type="response")
    } else if(resp == "vote"){
      rf_pred = predict(rf, sample, type="vote", norm.votes=F)
    } else if(resp == "prob"){
      rf_pred = predict(rf, sample, type="prob")
    }
  } else{
    if(resp =="response"){
      rf_pred = predict(rf, type="response")
    } else if(resp == "vote"){
      rf_pred = predict(rf, type="vote", norm.votes=F)
    } else if(resp == "prob"){
      rf_pred = predict(rf, type="prob")
    }
  }

  #return predictions
  return(rf_pred)
}
