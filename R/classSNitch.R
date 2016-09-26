#' Package for the autonomous classification of RNA structure change. 
#' @docType package
#' @name classSNitch
#' @author Chanin Tolson
#' @import gplots
#' @import ROCR
#' @references A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News 2(3), 18--22 (randomForest package) \cr\cr
#' \href{http://rmdb.stanford.edu/}{RNA Mapping Database}
#' @seealso  \code{\link{getFeatures}} \code{\link{classifyRNA}}
#' @examples #get change features
#' library("ROCR")
#' library("gplots")
#' 
#' data("shape_ex")
#' sample = getFeatures(shape_ex[2:nrow(shape_ex),], base=shape_ex[1,], trim=5)
#' 
#' #predict change
#' data("mutmap")
#' cr = classifyRNA(mutmap)
#' cr_pred = predict(cr, sample, resp="response")
#' 
#' #plot ROC curve (no change v. local/global change)
#' data("mutmap")
#' predobj = prediction(cr$votes[,1], mutmap[,1]==1)
#' perfobj = performance(predobj, 'tpr', 'fpr')
#' aucobj = performance(predobj, 'auc')
#' plot(perfobj@@x.values[[1]], perfobj@@y.values[[1]], lwd=2, 
#'      type="l", xlab="Specificity", ylab="Sensitivity")
#' points(c(-1,2),c(-1,2), col="red", type="l")
#' text(0.8, 0.2, paste("AUC: ", format(aucobj@@y.values, digits=2), sep=""), cex=1)
#' 
NULL
