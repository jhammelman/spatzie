#' @export
scanMotifs <- function(interactions,motifs,genome){
  anchorOneMatches <- matchMotifs(motifs, anchorOne(interactions),
                                  genome = genome, out='scores',bg='subject')
  anchorTwoMatches <- matchMotifs(motifs, anchorTwo(interactions),
                                  genome = genome, out='scores',bg='subject')
  interactionData <- list(interactions = interactions,
                          anchorOneMotifs = anchorOneMatches,
                          anchorTwoMotifs = anchorTwoMatches,
                          anchorOneMotifIndices = seq(length(anchorOneMatches$name)),
                          anchorTwoMotifIndices = seq(length(anchorTwoMatches$name)))
  class(interactionData) <- "interactionData"
  return(interactionData)
}
