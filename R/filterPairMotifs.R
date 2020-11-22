#' @title Filter significant motif interactions
#'
#' @description
#' Multiple hypothesis correction applied to filter for significant motif
#' interactions.
#'
#' @param interactionData TODO
#' @param method TODO
#' @param threshold TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom stats p.adjust
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats colMins
#' @export
filterPairMotifs <- function(interactionData, method=p.adjust.methods, threshold=0.05) {
  adjusted_p_interactions <- matrix(stats::p.adjust(as.vector(as.matrix(interactionData$pairMotifEnrich, method=method))), ncol=dim(interactionData$pairMotifEnrich)[2])
  rownames(adjusted_p_interactions) <- names(interactionData$anchorOneMotifIndices)
  colnames(adjusted_p_interactions) <- names(interactionData$anchorTwoMotifIndices)

  significant_anchor1 <- which(matrixStats::rowMins(adjusted_p_interactions) < threshold)
  adjusted_p_interactions_sig <- adjusted_p_interactions[significant_anchor1, ]
  significant_anchor2 <- which(matrixStats::colMins(adjusted_p_interactions) < threshold)
  adjusted_p_interactions_sig <- adjusted_p_interactions_sig[, significant_anchor2]

  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = interactionData$anchorOneMotifIndices,
                          anchorTwoMotifIndices = interactionData$anchorTwoMotifIndices,
                          pairMotifEnrich = adjusted_p_interactions,
                          pairMotifEnrich.sig = adjusted_p_interactions_sig,
                          is_multiple_hypothesis_corrected=FALSE)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
