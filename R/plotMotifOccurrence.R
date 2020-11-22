#' @title Plot motif occurrence
#'
#' @description
#' Plots a histogram of motif values (either counts, instances, or scores)
#' for anchorOne and anchorTwo regions.
#'
#' @param interactionData TODO
#' @param method TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom BiocGenerics colMeans
#' @importFrom SummarizedExperiment assays
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_bar
#' @export
plotMotifOccurrence <- function(interactionData, method = c("counts", "matches", "scores")) {
  if (method == "counts") {
    anchor1_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifCounts)
    anchor2_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifCounts)
  }else if (method == "matches") {
    anchor1_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches)
    anchor2_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches)
  }else{
    anchor1_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifScores)
    anchor2_values <- BiocGenerics::colMeans(SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifScores)
  }
  plottingdata <- data.frame(id = c(names(anchor1_values), names(anchor2_values)),
                            value = c(anchor1_values, anchor2_values),
                            variable = c(rep("anchor1", length(anchor1_values)),
                                       rep("anchor2", length(anchor2_values))))
  ggplot2::ggplot(plottingdata, ggplot2::aes(x = factor(id), y = value)) +
    ggplot2::facet_wrap(~variable) +
    ggplot2::geom_bar(ggplot2::aes(fill = factor(id)), ylab = method)
}
