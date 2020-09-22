anchorPairEnrich <- function(interactionData,method="countCorrelation"){
  significance = matrix(data=NA,nrow=length(interactionData$anchorOneMotifIndices),
                        ncol=length(interactionData$anchorTwoMotifIndices))
  indr=1
  for (i in interactionData$anchorOneMotifIndices){
    indc=1
    for (j in interactionData$anchorTwoMotifIndices){
      if (method == "countCorrelation"){
        significance[indr,indc] <- cor.test(assays(interactionData$anchorOneMotifs)$motifCounts[,i],assays(interactionData$anchorTwoMotifs)$motifCounts[,j],alternative='greater',method='pearson')$p.value
      }
      if (method == "scoreCorrelation"){
        significance[indr,indc] <- cor.test(assays(interactionData$anchorOneMotifs)$motifScores[,i],assays(interactionData$anchorTwoMotifs)$motifScores[,j],alternative='greater',method='pearson')$p.value
      }
      if (method == "countHypergeom"){
        significance[indr,indc] <- dhyper(sum((assays(interactionData$anchorOneMotifs)$motifMatches[,i])*(assays(interactionData$anchorTwoMotifs)$motifMatches[,j])),
               sum(assays(interactionData$anchorOneMotifs)$motifMatches[,i]),
               length(assays(interactionData$anchorOneMotifs)$motifMatches[,i])-sum(assays(interactionData$anchorOneMotifs)$motifMatches[,i]),
               sum(assays(interactionData$anchorTwoMotifs)$motifMatches[,j]))
      }
      indc= indc+1
    }
    indr= indr+1
  }
  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = interactionData$anchorOneMotifIndices,
                          anchorTwoMotifIndices = interactionData$anchorTwoMotifIndices,
                          pairMotifEnrich = significance)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
