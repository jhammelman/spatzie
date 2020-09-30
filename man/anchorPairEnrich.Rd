\name{anchorPairEnrich}
\alias{anchorPairEnrich}
\title{
Determine enriched motifs in anchors
}
\description{
Determine whether motifs between paired bed regions have a statistically significant relationship. Options for significance are motif score correlation, motif count correlation, or hypergeometric motif co-occurrence.
}
\usage{
anchorPairEnrich(interactionData, method = "countCorrelation")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interactionData}{
  an interactionData object of paired genomic regions
}
  \item{method}{
  choice of method for co-occurrence include "countCorrelation", "scoreCorrelation", or "countHypergeom"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an interactionData object where obj$pairMotifEnrich contains the p-values for significance of seeing a higher co-occurrence than what we get by chance.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jennifer Hammelman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (interactionData, method = "countCorrelation")
{
    indr = 1
    for (i in interactionData$anchorOneMotifIndices) {
        indc = 1
        for (j in interactionData$anchorTwoMotifIndices) {
            if (method == "countCorrelation") {
                significance[indr, indc] <- cor.test(assays(interactionData$anchorOneMotifIndices)$motifCounts[,
                  i], assays(interactionData$anchorTwoMotifIndices)$motifCounts[,
                  j], alternative = "greater", method = "pearson")$p.value
            }
            if (method == "scoreCorrelation") {
                significance[indr, indc] <- cor.test(assays(interactionData$anchorOneMotifIndices)$motifScores[,
                  i], assays(interactionData$anchorTwoMotifIndices)$motifScores[,
                  j], alternative = "greater", method = "pearson")$p.value
            }
            if (method == "countHypergeom") {
                significance[indr, indc] <- dhyper(sum((assays(interactionData$anchorOneMotifIndices)$motifMatches[,
                  i]) * (assays(interactionData$anchorTwoMotifIndices)$motifMatches[,
                  j])), sum(assays(interactionData$anchorOneMotifIndices)$motifMatches[,
                  i]), length(assays(interactionData$anchorOneMotifIndices)$motifMatches[,
                  i]) - sum(assays(interactionData$anchorOneMotifIndices)$motifMatches[,
                  i]), sum(assays(interactionData$anchorTwoMotifIndices)$motifMatches[,
                  j]))
            }
            indc = indc + 1
        }
        indr = indr + 1
    }
    interactionData <- list(interactions = interactionData$interactions,
        anchorOneMotifs = interactionData$anchorOneMotifs, anchorTwoMotifs = interactionData$anchorTwoMotifs,
        anchorOneMotifIndices = interactionData$anchorOneMotifIndices,
        anchorTwoMotifIndices = interactionData$anchorTwoMotifIndices,
        pairMotifEnrich = significance)
    class(interactionData) <- "interactionData"
    return(interactionData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
