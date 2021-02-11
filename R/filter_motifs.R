#' @title Filter motifs based on occurrence within interaction data
#'
#' @description
#' Select a subset of motifs that are in at least a threshold fraction of
#' regions. Motif subsets are selected separately for anchor one and anchor
#' two regions.
#'
#' @param interaction_data an interactionData object of paired genomic regions
#' @param threshold fraction of interactions that should contain a motif for a
#' motif to be considered
#' @return an interactionData object where \code{obj$anchor1_motif_indices}
#' and \code{obj$anchor2_motif_indices} have been filtered to motifs that are
#' present in a threshold fraction of interactions
#'
#' @examples
#' genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
#' motif_file <- system.file(
#'   "extdata/consensus_HOCOMOCOv11_core_MOUSE-plus_YY1.piq",
#'   package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motif_file, matrixClass = "PFM")
#'
#' yy1_pd_interaction <- scan_motifs(spatzie:::yy1_interactions, motifs, genome)
#' yy1_pd_interaction <- filter_motifs(yy1_pd_interaction, 0.4)
#'
#' @author Jennifer Hammelman
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocGenerics colMeans
#' @export
filter_motifs <- function(interaction_data, threshold) {
  anchor1_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor1_motifs)
  anchor2_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor2_motifs)
  anchor1_indices <- which(
    BiocGenerics::colMeans(anchor1_motifs$motifMatches) > threshold)
  anchor2_indices <- which(
    BiocGenerics::colMeans(anchor2_motifs$motifMatches) > threshold)
  interaction_data <- list(interactions = interaction_data$interactions,
                           anchor1_motifs = interaction_data$anchor1_motifs,
                           anchor2_motifs = interaction_data$anchor2_motifs,
                           anchor1_motif_indices = anchor1_indices,
                           anchor2_motif_indices = anchor2_indices)
  class(interaction_data) <- "interactionData"
  return(interaction_data)
}
