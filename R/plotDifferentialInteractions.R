#' @title Plot differential interactions
#'
#' @description
#' Plots clustered heatmap of log likelihood ration for motif interactions
#' between two datasets
#'
#' @param comparedMotifPairs TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom pheatmap pheatmap
#' @export
plotDifferentialInteractions <- function(comparedMotifPairs){
  pheatmap::pheatmap(comparedMotifPairs, fontsize = 6)
}
