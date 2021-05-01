#' @title Compare pairs of motifs between two interaction datasets
#'
#' @description
#' Compute the log-likelihood ratio that a motif pair is differential between
#' two interaction datasets. Note that motif pair significance should have
#' been computed using the same method for both datasets.
#'
#' @param interaction_data1 an interactionData object of paired genomic regions
#'                          that has been scanned for significant motif:motif
#'                          interactions
#' @param interaction_data2 an interactionData object of paired genomic regions
#'                          that has been scanned for significant motif:motif
#'                          interactions
#' @param differential_p threshold for significance of differential p-value
#' @return a matrix of the log likelihood ratio of motif pairs that are
#'         significantly differential between between two interactionData sets
#'
#' @examples
#' \dontrun{
#' pheatmap::pheatmap(compare_motif_pairs(spatzie:::interactionDataK562,
#'                                        spatzie:::interactionDataMSLCL),
#'                    fontsize = 6)
#' }
#' @author Jennifer Hammelman
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats colMaxs
#' @export
compare_motif_pairs <- function(interaction_data1, interaction_data2,
                                differential_p = 0.05) {
  data1_anchor1 <- (rownames(interaction_data1$pair_motif_enrich) %in%
                      rownames(interaction_data2$pair_motif_enrich))
  data1_anchor2 <- (colnames(interaction_data1$pair_motif_enrich) %in%
                      colnames(interaction_data2$pair_motif_enrich))
  data1_mat <- interaction_data1$pair_motif_enrich[data1_anchor1, ]
  data1_mat <- data1_mat[, data1_anchor2]

  data2_anchor1 <- (rownames(interaction_data2$pair_motif_enrich) %in%
                      rownames(interaction_data1$pair_motif_enrich))
  data2_anchor2 <- (colnames(interaction_data2$pair_motif_enrich) %in%
                      colnames(interaction_data1$pair_motif_enrich))
  data2_mat <- interaction_data2$pair_motif_enrich[data2_anchor1, ]
  data2_mat <- data2_mat[, data2_anchor2]

  pseudocount <- 1e-100
  differential <- (-(log2(data1_mat + pseudocount) -
                       log2(data2_mat + pseudocount)))
  differential <- differential[matrixStats::rowMaxs(abs(differential)) >
                                 -log2(differential_p), ]
  differential <- differential[, matrixStats::colMaxs(abs(differential)) >
                                 -log2(differential_p)]
  return(differential)
}
