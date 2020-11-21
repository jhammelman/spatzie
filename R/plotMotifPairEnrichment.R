#' @title Plot motif pair occurrence
#'
#' @description
#' TODO
#'
#' @param pairEnrichmentMatrix TODO
#' @return TODO
#' @author Jennifer Hammelman
#' @export
plotMotifPairsHeatmap <- function(pairEnrichmentMatrix){
  pheatmap(pairEnrichmentMatrix,fontsize=6)
}
