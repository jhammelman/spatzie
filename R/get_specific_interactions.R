#' @title Get interactions that contain a specific motif pair
#'
#' @description
#' Select interactions that contain anchor1_motif within anchor 1 and
#' anchor2_motif within anchor 2.
#'
#' @param interaction_data an interactionData object of paired genomic regions
#' @param anchor1_motif Motif name from \code{interactionData$anchor1_motifs}
#' @param anchor2_motif Motif name from \code{interactionData$anchor2_motifs}
#' @return a \code{\link[GenomicInteractions]{GenomicInteractions}} object
#' containing a subset subset of interactions that contain an instance of
#' \code{anchor1_motif} in anchor 1 and \code{anchor2_motif} in anchor 2
#'
#' @examples
#' genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
#' motifs_file <- system.file("extdata/motifs_subset.txt.gz",
#'                            package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
#'
#' yy1_pd_interaction <- scan_motifs(spatzie::interactions_yy1, motifs, genome)
#' yy1_pd_interaction <- filter_motifs(yy1_pd_interaction, 0.4)
#' yy1_pd_count_corr <- anchor_pair_enrich(yy1_pd_interaction,
#'                                         method = "score")
#' yy1_yy1_interactions <- get_specific_interactions(
#'   yy1_pd_interaction,
#'   anchor1_motif = "YY1",
#'   anchor2_motif = "YY1")
#'
#' @author Jennifer Hammelman
#' @importFrom SummarizedExperiment assays
#' @export
get_specific_interactions <- function(interaction_data, anchor1_motif = NULL,
                                      anchor2_motif = NULL) {
  anchor1_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor1_motifs)
  anchor2_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor2_motifs)

  if (is.null(anchor1_motif) && is.null(anchor2_motif)) {
    return(interaction_data$interactions)
  } else if (is.null(anchor1_motif)) {
    motif_mask_anchor2 <- which(
      anchor2_motif == colnames(interaction_data$anchor2_motifs))

    if (motif_mask_anchor2 == 0) {
      stop(paste("No match found for motif ", anchor2_motif))
    }

    interaction_mask <- (anchor2_motifs$motifMatches[, motif_mask_anchor2] == TRUE)

    return(interaction_data$interactions[interaction_mask])
  } else if (is.null(anchor2_motif)) {
    motif_mask_anchor1 <- which(
      anchor1_motif == colnames(interaction_data$anchor1_motifs))

    if (motif_mask_anchor1 == 0) {
      stop(paste("No match found for motif ", anchor1_motif))
    }

    interaction_mask <- (anchor1_motifs$motifMatches[, motif_mask_anchor1] == TRUE)

    return(interaction_data$interactions[interaction_mask])
  } else {
    motif_mask_anchor1 <- which(
      anchor1_motif == colnames(interaction_data$anchor1_motifs))

    if (motif_mask_anchor1 == 0) {
      stop(paste("No match found for motif ", anchor1_motif))
    }

    interaction_mask_anchor1 <- (anchor1_motifs$motifMatches[, motif_mask_anchor1] == TRUE)

    motif_mask_anchor2 <- which(
      anchor2_motif == colnames(interaction_data$anchor2_motifs))

    if (motif_mask_anchor2 == 0) {
      stop(paste("No match found for motif ", anchor2_motif))
    }

    interaction_mask_anchor2 <- (anchor2_motifs$motifMatches[, motif_mask_anchor2] == TRUE)

    return(interaction_data$interactions[interaction_mask_anchor1 & interaction_mask_anchor2])
  }
}
