#' @title TODO
#'
#' @description
#' TODO
#'
#' @param interactions TODO
#' @param allInteractions TODO
#' @param bedranges TODO
#' @param anchor TODO
#' @return TODO
#'
#' @examples
#' # TODO
#'
#' @author Jennifer Hammelman
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicRanges intersect
#' @importFrom BiocGenerics unique
#' @importFrom stats phyper
#' @importFrom GenomicInteractions anchorTwo
#' @export
overlapEnrich <- function(interactions, allInteractions, bedranges, anchor=c("anchorOne", "anchorTwo")){
  if (anchor=="anchorOne"){
    all <- unique(GenomicInteractions::anchorOne(allInteractions))
    allinterbed <- GenomicRanges::intersect(all, BiocGenerics::unique(bedranges))
    sel <- unique(GenomicInteractions::anchorOne(interactions))
    selinterbed <- GenomicRanges::intersect(sel, BiocGenerics::unique(bedranges))
    print(c(length(selinterbed),
            length(sel),
            length(all)-length(sel),
            length(allinterbed)))
    return(stats::phyper(length(selinterbed),
                  length(sel),
                  length(all)-length(sel),
                  length(allinterbed), lower.tail=FALSE))
  }
  else{
    all <- unique(GenomicInteractions::anchorTwo(allInteractions))
    allinterbed <- GenomicRanges::intersect(all, BiocGenerics::unique(bedranges))
    sel <- unique(GenomicInteractions::anchorTwo(interactions))
    selinterbed <- GenomicRanges::intersect(sel, BiocGenerics::unique(bedranges))
    print(c(length(selinterbed),
            length(sel),
            length(all)-length(sel),
            length(allinterbed)))
    return(stats::phyper(length(selinterbed),
                  length(sel),
                  length(all)-length(sel),
                  length(allinterbed), lower.tail=FALSE))
  }
}
