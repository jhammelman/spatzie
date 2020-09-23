filterMotifs <- function(interactionData,threshold){
  anchorOneIndices <- which(colMeans(as.matrix(assays(interactionData$anchorOneMotifs)$motifMatches)) > threshold)
  anchorTwoIndices <- which(colMeans(as.matrix(assays(interactionData$anchorTwoMotifs)$motifMatches)) > threshold)
  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = anchorOneIndices,
                          anchorTwoMotifIndices = anchorTwoIndices)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
