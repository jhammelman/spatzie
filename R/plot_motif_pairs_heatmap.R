#' @title Plot motif pair occurrence
#'
#' @description
#' TODO
#'
#' @param pair_enrichment_matrix TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom pheatmap pheatmap
#' @export
plot_motif_pairs_heatmap <- function(pair_enrichment_matrix) {
  return(pheatmap::pheatmap(pair_enrichment_matrix, fontsize = 6))
}
