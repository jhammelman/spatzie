#' @title Plot motif pair occurrence
#'
#' @description
#' TODO
#'
#' @param pairEnrichmentMatrix TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom pheatmap pheatmap
#' @export
plot_motif_pairs_heatmap <- function(pairEnrichmentMatrix) {
  pheatmap::pheatmap(pairEnrichmentMatrix, fontsize = 6)
}
