#' @title Determine enriched motifs in anchors
#'
#' @description
#' Determine whether motifs between paired bed regions have a statistically
#' significant relationship. Options for significance are motif score
#' correlation, motif count correlation, or hypergeometric motif co-occurrence.
#'
#' @param interactionData an interactionData object of paired genomic regions
#' @param method choice of method for co-occurrence include
#' \code{countCorrelation}, \code{scoreCorrelation}, \code{countHypergeom}, or
#' \code{countFisher}
#' @return an interactionData object where \code{obj$pairMotifEnrich} contains
#' the p-values for significance of seeing a higher co-occurrence than
#' what we get by chance.
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom stats cor.test
#' @importFrom SummarizedExperiment assays
#' @importFrom stats phyper
#' @importFrom stats fisher.test
#' @export
anchorPairEnrich <- function(interactionData, method = c("countCorrelation", "scoreCorrelation", "countHypergeom", "countFisher")) {
  significance <- matrix(data = NA, nrow = length(interactionData$anchorOneMotifIndices),
                        ncol = length(interactionData$anchorTwoMotifIndices))
  indr <- 1
  for (i in interactionData$anchorOneMotifIndices) {
    indc <- 1
    for (j in interactionData$anchorTwoMotifIndices) {
      if (method == "countCorrelation") {
        significance[indr, indc] <- stats::cor.test(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifCounts[, i], SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifCounts[, j], alternative = 'greater', method = 'pearson')$p.value
      }
      if (method == "scoreCorrelation") {
        significance[indr, indc] <- stats::cor.test(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifScores[, i], SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifScores[, j], alternative = 'greater', method = 'pearson')$p.value
      }
      if (method == "countHypergeom") {
        significance[indr, indc] <- stats::phyper(sum((SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]) * (SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, j])),
               sum(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]),
               length(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]) - sum(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]),
               sum(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, j]), lower.tail = FALSE)
      }
      if (method == "countFisher") {
        dobpos <- sum((SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]) * (SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, j]))
        dobneg <- sum((!SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]) * (!SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, j]))
        fisher_mat <- matrix(c(dobpos,
                 sum(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, i]) - dobpos,
                 sum(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, j]) - dobpos,
                 dobneg), nrow = 2)
        significance[indr, indc] <- stats::fisher.test(fisher_mat, alternative = "greater")$p.value
      }
      indc <- indc + 1
    }
    indr <- indr + 1
  }
  rownames(significance) <- names(interactionData$anchorOneMotifIndices)
  colnames(significance) <- names(interactionData$anchorTwoMotifIndices)

  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = interactionData$anchorOneMotifIndices,
                          anchorTwoMotifIndices = interactionData$anchorTwoMotifIndices,
                          pairMotifEnrich = significance,
                          is_multiple_hypothesis_corrected = FALSE)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
