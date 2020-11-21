#' @title Plot differential interactions
#'
#' @description
#' Plots clustered heatmap of log likelihood ration for motif interactions
#' between two datasets
#'
#' @param comparedMotifPairs TODO
#' @return TODO
#' @author Jennifer Hammelman
#' @export
plotDifferentialInteractions <- function(comparedMotifPairs){
  pheatmap(comparedMotifPairs,fontsize=6)
}
