#' @title Filter significant motif interactions
#'
#' @description
#' Multiple hypothesis correction applied to filter for significant motif
#' interactions.
#'
#' @param interaction_data TODO
#' @param method TODO
#' @param threshold TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom stats p.adjust
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats colMins
#' @export
filter_pair_motifs <- function(interaction_data, method = p.adjust.methods,
                               threshold = 0.05) {
  adjusted_p_interactions <- matrix(stats::p.adjust(
    as.vector(as.matrix(interaction_data$pair_motif_enrich, method = method))),
    ncol = dim(interaction_data$pair_motif_enrich)[2])
  rownames(adjusted_p_interactions) <- names(
    interaction_data$anchor1_motif_indices)
  colnames(adjusted_p_interactions) <- names(
    interaction_data$anchor2_motif_indices)

  anchor1_mask <- which(
    matrixStats::rowMins(adjusted_p_interactions) < threshold)
  adjusted_p_interactions_sig <- adjusted_p_interactions[anchor1_mask, ]
  anchor2_mask <- which(
    matrixStats::colMins(adjusted_p_interactions) < threshold)
  adjusted_p_interactions_sig <- adjusted_p_interactions_sig[, anchor2_mask]

  interaction_data <- list(
    interactions = interaction_data$interactions,
    anchor1_motifs = interaction_data$anchor1_motifs,
    anchor2_motifs = interaction_data$anchor2_motifs,
    anchor1_motif_indices = interaction_data$anchor1_motif_indices,
    anchor2_motif_indices = interaction_data$anchor2_motif_indices,
    pair_motif_enrich = adjusted_p_interactions,
    pair_motif_enrich_sig = adjusted_p_interactions_sig,
    is_multiple_hypothesis_corrected = FALSE)
  class(interaction_data) <- "interactionData"
  return(interaction_data)
}
