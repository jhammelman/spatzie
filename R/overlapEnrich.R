overlapEnrich <- function(interactions,allInteractions,bedranges,anchor=c("anchorOne","anchorTwo")){
  if (anchor=="anchorOne"){
    all <- unique(anchorOne(allInteractions))
    allinterbed <- intersect(all,unique(bedranges))
    sel <- unique(anchorOne(interactions))
    selinterbed <- intersect(sel,unique(bedranges))
    print(c(length(selinterbed),
            length(sel),
            length(all)-length(sel),
            length(allinterbed)))
    return(phyper(length(selinterbed),
                  length(sel),
                  length(all)-length(sel),
                  length(allinterbed),lower.tail=FALSE))
  }
  else{
    all <- unique(anchorTwo(allInteractions))
    allinterbed <- intersect(all,unique(bedranges))
    sel <- unique(anchorTwo(interactions))
    selinterbed <- intersect(sel,unique(bedranges))
    print(c(length(selinterbed),
            length(sel),
            length(all)-length(sel),
            length(allinterbed)))
    return(phyper(length(selinterbed),
                  length(sel),
                  length(all)-length(sel),
                  length(allinterbed),lower.tail=FALSE))
  }
}
