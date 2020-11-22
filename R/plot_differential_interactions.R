#' @title Plot differential interactions
#'
#' @description
#' Plots clustered heatmap of log likelihood ration for motif interactions
#' between two datasets
#'
#' @param compared_motif_pairs TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom pheatmap pheatmap
#' @export
plot_differential_interactions <- function(compared_motif_pairs) {
  pheatmap::pheatmap(compared_motif_pairs, fontsize = 6)
}
