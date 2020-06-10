#' Quick computation of the AUC
#' 
#' @param label A factor vector with levels of c("Yes","No") in any order, or a character vector taking the same values.
#' @param scores A numeric vector of the same length above. Higher scores must correspond to the "Yes" level of the labels.
#' @return The AUC, defined as the area under the ROC curve computed using \code{label} and \code{scores}.
#' @examples 
#' \dontrun{
#' y=rbinom(100, size=1, prob=1/2)
#' label = factor(y, levels=c(0,1), labels=c("No", "Yes"))
#' scores = runif(100) # higher scores predict label == "Yes"
#' aucC::aucC_wrap(label=label, scores=scores)
#' }
#' @export
aucC_wrap <- function(label, scores){
  pos_predictor = scores[label == "Yes"]
  neg_predictor = scores[label != "Yes"]
  aucC(pos_predictor, neg_predictor)
}