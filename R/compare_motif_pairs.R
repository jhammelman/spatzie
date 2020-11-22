#' @title Compare pairs of motifs between two interaction datasets
#'
#' @description
#' Compute the log-likelihood ratio that a motif pair is differential between
#' two interaction datasets. Note that motif pair significance should have
#' been computed using the same method for both datasets.
#'
#' @param interaction_data1 TODO
#' @param interaction_data2 TODO
#' @param differential_p TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
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

  differential <- (-(log2(data1_mat + 1e-100) - log2(data2_mat + 1e-100)))
  differential <- differential[matrixStats::rowMaxs(abs(differential)) > -log2(differential_p), ]
  differential <- differential[, matrixStats::colMaxs(abs(differential)) > -log2(differential_p)]
  return(differential)
}
