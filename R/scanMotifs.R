#' @title Scans interaction file for motif instances
#'
#' @description
#' Uses motifmatchR to scan interaction regions for given motifs.
#'
#' @param interactions an interactionData object of paired genomic regions
#' @param motifs a TFBS tools matrix of DNA binding motifs
#' @param genome a Biostrings genome must match chromosomes from interaction
#' data file
#' @return an interaction data object where \code{obj$anchorOneMotifs} and
#' \code{obj$anchorTwoMotifs} contain information about the scores and matches
#' to motifs from anchor one and anchor two of interaction data genomic regions
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicInteractions anchorTwo
#' @export
scanMotifs <- function(interactions, motifs, genome){
  anchorOneMatches <- motifmatchr::matchMotifs(motifs, GenomicInteractions::anchorOne(interactions),
                                  genome = genome, out='scores', bg='subject')
  anchorTwoMatches <- motifmatchr::matchMotifs(motifs, GenomicInteractions::anchorTwo(interactions),
                                  genome = genome, out='scores', bg='subject')
  interactionData <- list(interactions = interactions,
                          anchorOneMotifs = anchorOneMatches,
                          anchorTwoMotifs = anchorTwoMatches,
                          anchorOneMotifIndices = seq(length(anchorOneMatches$name)),
                          anchorTwoMotifIndices = seq(length(anchorTwoMatches$name)))
  class(interactionData) <- "interactionData"
  return(interactionData)
}
