\name{filterMotifs}
\alias{filterMotifs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter motifs based on occurrence within interaction data
}
\description{
Select a subset of motifs that are in at least a threshold fraction of regions. Motif subsets are selected separately for anchor one and anchor two regions.
}
\usage{
filterMotifs(interactionData, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interactionData}{
  an interactionData object of paired genomic regions
}
  \item{threshold}{
  fraction of interactions that should contain a motif for a motif to be considered
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an interactionData object where obj$anchorOneMotifIndices and obj$anchorTwoMotifIndices have been filtered to motifs that are present in a threshold fraction of interactions
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
function (interactionData, threshold)
{
    anchorOneIndices <- which(colMeans(assays(interactionData$anchorOneMotifs)$motifMatches) >
        threshold)
    anchorTwoIndices <- which(colMeans(assays(interactionData$anchorTwoMotifs)$motifMatches) >
        threshold)
    interactionData <- list(interactions = interactionData$interactions,
        anchorOneMotifs = interactionData$anchorOneMotifs, anchorTwoMotifs = interactionData$anchorTwoMotifs,
        anchorOneMotifIndices = anchorOneIndices, anchorTwoMotifIndices = anchorTwoIndices)
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
