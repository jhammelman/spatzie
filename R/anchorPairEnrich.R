anchorPairEnrich <- function(interactionData,method=c("countCorrelation","scoreCorrelation","countHypergeom","countFisher")){
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
        significance[indr,indc] <- phyper(sum((assays(interactionData$anchorOneMotifs)$motifMatches[,i])*(assays(interactionData$anchorTwoMotifs)$motifMatches[,j])),
               sum(assays(interactionData$anchorOneMotifs)$motifMatches[,i]),
               length(assays(interactionData$anchorOneMotifs)$motifMatches[,i])-sum(assays(interactionData$anchorOneMotifs)$motifMatches[,i]),
               sum(assays(interactionData$anchorTwoMotifs)$motifMatches[,j]),lower.tail=FALSE)
      }
      if (method == "countFisher"){
        dobpos <- sum((assays(interactionData$anchorOneMotifs)$motifMatches[,i])*(assays(interactionData$anchorTwoMotifs)$motifMatches[,j]))
        dobneg <- sum((!assays(interactionData$anchorOneMotifs)$motifMatches[,i])*(!assays(interactionData$anchorTwoMotifs)$motifMatches[,j]))
        fisher_mat <- matrix(c(dobpos,
                 sum(assays(interactionData$anchorOneMotifs)$motifMatches[,i])-dobpos,
                 sum(assays(interactionData$anchorTwoMotifs)$motifMatches[,j])-dobpos,
                 dobneg),nrow=2)
        significance[indr,indc] <- fisher.test(fisher_mat,alternative="greater")$p.value
      }
      indc= indc+1
    }
    indr= indr+1
  }
  rownames(significance) <- names(interactionData$anchorOneMotifIndices)
  colnames(significance) <- names(interactionData$anchorTwoMotifIndices)

  interactionData <- list(interactions = interactionData$interactions,
                          anchorOneMotifs = interactionData$anchorOneMotifs,
                          anchorTwoMotifs = interactionData$anchorTwoMotifs,
                          anchorOneMotifIndices = interactionData$anchorOneMotifIndices,
                          anchorTwoMotifIndices = interactionData$anchorTwoMotifIndices,
                          pairMotifEnrich = significance,
                          is_multiple_hypothesis_corrected=FALSE)
  class(interactionData) <- "interactionData"
  return(interactionData)
}
