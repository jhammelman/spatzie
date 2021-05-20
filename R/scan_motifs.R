#' @title Scans interaction file for motif instances
#'
#' @description
#' Uses motifmatchR to scan interaction regions for given motifs.
#'
#' @param interactions an interactionData object of paired genomic regions
#' @param motifs a TFBS tools matrix of DNA binding motifs
#' @param genome a Biostrings genome must match chromosomes from interaction
#' data file
#' @return an interaction data object where \code{obj$anchor1_motifs} and
#' \code{obj$anchor2_motifs} contain information about the scores and matches
#' to motifs from anchor one and anchor two of interaction data genomic regions
#'
#' @examples
#' genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm9")
#' motifs_file <- system.file("extdata/motifs_subset.txt.gz",
#'                            package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
#'
#' yy1_pd_interaction <- scan_motifs(spatzie::interactions_yy1, motifs, genome)
#'
#' @author Jennifer Hammelman
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicInteractions anchorTwo
#' @export
scan_motifs <- function(interactions, motifs, genome) {
  anchor1_matches <- motifmatchr::matchMotifs(
    motifs, GenomicInteractions::anchorOne(interactions),
    genome = genome, out = "scores", bg = "subject")
  anchor2_matches <- motifmatchr::matchMotifs(
    motifs, GenomicInteractions::anchorTwo(interactions),
    genome = genome, out = "scores", bg = "subject")
  interaction_data <- list(
    interactions = interactions,
    anchor1_motifs = anchor1_matches,
    anchor2_motifs = anchor2_matches,
    anchor1_motif_indices = seq(length(anchor1_matches$name)),
    anchor2_motif_indices = seq(length(anchor2_matches$name)))
  class(interaction_data) <- "interactionData"
  return(interaction_data)
}
