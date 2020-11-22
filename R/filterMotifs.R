#' @title Filter motifs based on occurrence within interaction data
#'
#' @description
#' Select a subset of motifs that are in at least a threshold fraction of
#' regions. Motif subsets are selected separately for anchor one and anchor
#' two regions.
#'
#' @param interactionData an interactionData object of paired genomic regions
#' @param threshold fraction of interactions that should contain a motif for a
#' motif to be considered
#' @return an interactionData object where \code{obj$anchorOneMotifIndices}
#' and \code{obj$anchorTwoMotifIndices} have been filtered to motifs that are
#' present in a threshold fraction of interactions
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom SummarizedExperiment assays
#' @export
filter_motifs <- function(interactionData, threshold) {
  anchorOneIndices <- which(colMeans(as.matrix(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches)) > threshold)
  anchorTwoIndices <- which(colMeans(as.matrix(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches)) > threshold)
  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = anchorOneIndices,
                          anchorTwoMotifIndices = anchorTwoIndices)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
