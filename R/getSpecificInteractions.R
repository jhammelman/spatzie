#' @title Get interactions that contain a specific motif pair
#'
#' @description
#' Select interactions that contain anchorOneMotif within anchorOne and
#' anchorTwoMotif within anchorTwo.
#'
#' @param interactionData TODO
#' @param anchorOneMotif TODO
#' @param anchorTwoMotif TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom SummarizedExperiment assays
#' @export
getSpecificInteractions <- function(interactionData, anchorOneMotif="", anchorTwoMotif=""){
  #TODO return a subset of interactions that are only containing
  #anchorOneMotif or anchorTwoMotif
  if (anchorOneMotif == "" && anchorTwoMotif == ""){
    return(interactionData)
  }
  if (anchorOneMotif == ""){
    findMotif <- which(anchorTwoMotif == colnames(interactionData$anchorTwoMotifs))
    interactionskeepMask <- (interactionData$anchorTwoMotifs$motifInstances[, findMotif] == TRUE)
    return(interactionData$interactions[interactionskeepMask])
  }
  else if (anchorTwoMotif == ""){
    findMotif <- which(anchorOneMotif == colnames(interactionData$anchorOneMotifs))
    interactionskeepMask <- (interactionData$anchorOneMotifs$motifInstances[, findMotif] == TRUE)
    return(interactionData$interactions[interactionskeepMask])
  }
  else{
    findTwoMotif <- which(anchorTwoMotif == colnames(interactionData$anchorTwoMotifs))
    interactionsTwokeepMask <- (SummarizedExperiment::assays(interactionData$anchorTwoMotifs)$motifMatches[, findTwoMotif] == TRUE)

    findOneMotif <- which(anchorOneMotif == colnames(interactionData$anchorOneMotifs))
    interactionsOnekeepMask <- (SummarizedExperiment::assays(interactionData$anchorOneMotifs)$motifMatches[, findOneMotif] == TRUE)
    return(interactionData$interactions[interactionsTwokeepMask & interactionsOnekeepMask])
  }
}
