compareMotifPairs <- function(interactionData1,interactionData2){
  data1_anchor1 <- (rownames(interactionData1$pairMotifEnrich.sig) %in% rownames(interactionData2$pairMotifEnrich.sig))
  data1_anchor2 <- (colnames(interactionData1$pairMotifEnrich.sig) %in% colnames(interactionData2$pairMotifEnrich.sig))
  data1_mat <- interactionData1$pairMotifEnrich.sig[data1_anchor1,]
  data1_mat <- data1_mat[,data1_anchor2]

  data2_anchor2 <- (names(interactionData2$anchorTwoMotifIndices) %in% names(interactionData1$anchorTwoMotifIndices))
  data2_anchor1 <- (names(interactionData2$anchorOneMotifIndices) %in% names(interactionData1$anchorOneMotifIndices))
  data2_mat <- interactionData2$pairMotifEnrich.sig[data2_anchor1,]
  data2_mat <- data2_mat[,data2_anchor2]

  differential <- (-(log2(data1_mat+1e-100)-log2(data2_mat+1e-100)))
  return(differential)

}
