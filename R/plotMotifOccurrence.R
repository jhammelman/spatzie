plotMotifOccurrence <- function(interactionData,method=c("counts","matches","scores")){
  if (method == "counts"){
    anchor1_values <- colMeans(assays(interactionData$anchorOneMotifs)$motifCounts)
    anchor2_values <- colMeans(assays(interactionData$anchorTwoMotifs)$motifCounts)
  }else if (method == "matches"){
    anchor1_values <- colMeans(assays(interactionData$anchorOneMotifs)$motifMatches)
    anchor2_values <- colMeans(assays(interactionData$anchorTwoMotifs)$motifMatches)
  }else{
    anchor1_values <- colMeans(assays(interactionData$anchorOneMotifs)$motifScores)
    anchor2_values <- colMeans(assays(interactionData$anchorTwoMotifs)$motifScores)
  }
  plottingdata <- data.frame(id=c(names(anchor1_values),names(anchor2_values)),
                            value=c(anchor1_values,anchor2_values),
                            variable=c(rep("anchor1",length(anchor1_values)),
                                       rep("anchor2",length(anchor2_values))))
  ggplot(plottingdata,aes(x=factor(id), y = value)) +
    facet_wrap(~variable) +
    geom_bar(aes(fill = factor(id)),ylab=method)
}
