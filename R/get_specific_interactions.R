#' @title Get interactions that contain a specific motif pair
#'
#' @description
#' Select interactions that contain anchor1_motif within anchor 1 and
#' anchor2_motif within anchor 2.
#'
#' @param interaction_data TODO
#' @param anchor1_motif TODO
#' @param anchor2_motif TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom SummarizedExperiment assays
#' @export
get_specific_interactions <- function(interaction_data, anchor1_motif = "", anchor2_motif = "") {
  #TODO return a subset of interactions that are only containing
  #anchor1_motif or anchor2_motif
  if (anchor1_motif == "" && anchor2_motif == "") {
    return(interaction_data)
  }
  if (anchor1_motif == "") {
    motif_mask <- which(anchor2_motif == colnames(interaction_data$anchor2_motifs))
    interaction_mask <- (interaction_data$anchor2_motifs$motifInstances[, motif_mask] == TRUE)
    return(interaction_data$interactions[interaction_mask])
  }
  else if (anchor2_motif == "") {
    motif_mask <- which(anchor1_motif == colnames(interaction_data$anchor1_motifs))
    interaction_mask <- (interaction_data$anchor1_motifs$motifInstances[, motif_mask] == TRUE)
    return(interaction_data$interactions[interaction_mask])
  }
  else{
    motif_mask_anchor2 <- which(anchor2_motif == colnames(interaction_data$anchor2_motifs))
    interaction_mask_anchor2 <- (SummarizedExperiment::assays(interaction_data$anchor2_motifs)$motifMatches[, motif_mask_anchor2] == TRUE)

    motif_mask_anchor1 <- which(anchor1_motif == colnames(interaction_data$anchor1_motifs))
    interaction_mask_anchor1 <- (SummarizedExperiment::assays(interaction_data$anchor1_motifs)$motifMatches[, motif_mask_anchor1] == TRUE)
    return(interaction_data$interactions[interaction_mask_anchor2 & interaction_mask_anchor1])
  }
}
