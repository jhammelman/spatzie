#' @title Compare pairs of motifs between two interaction datasets
#'
#' @description
#' Compute the log-likelihood ratio that a motif pair is differential between
#' two interaction datasets. Note that motif pair significance should have
#' been computed using the same method for both datasets.
#'
#' @param interactionData1 TODO
#' @param interactionData2 TODO
#' @param differential_p TODO
#' @return TODO
#' @author Jennifer Hammelman
#' @export
compareMotifPairs <- function(interactionData1,interactionData2,differential_p=0.05){
  data1_anchor1 <- (rownames(interactionData1$pairMotifEnrich) %in% rownames(interactionData2$pairMotifEnrich))
  data1_anchor2 <- (colnames(interactionData1$pairMotifEnrich) %in% colnames(interactionData2$pairMotifEnrich))
  data1_mat <- interactionData1$pairMotifEnrich[data1_anchor1,]
  data1_mat <- data1_mat[,data1_anchor2]

  data2_anchor1 <- (rownames(interactionData2$pairMotifEnrich) %in% rownames(interactionData1$pairMotifEnrich))
  data2_anchor2 <- (colnames(interactionData2$pairMotifEnrich) %in% colnames(interactionData1$pairMotifEnrich))
  data2_mat <- interactionData2$pairMotifEnrich[data2_anchor1,]
  data2_mat <- data2_mat[,data2_anchor2]
  #print(dim(data1_mat))
  #print(dim(data2_mat))

  differential <- (-(log2(data1_mat+1e-100)-log2(data2_mat+1e-100)))
  differential <- differential[rowMaxs(abs(differential)) > -log2(0.05),]
  differential <- differential[,colMaxs(abs(differential)) > -log2(0.05)]
  return(differential)

}
