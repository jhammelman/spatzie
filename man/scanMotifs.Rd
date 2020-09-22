\name{scanMotifs}
\alias{scanMotifs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scans interactions for motifs
}
\description{
Uses motifmatchR to scan interaction regions for given motifs.
}
\usage{
scanMotifs(interactions, motifs, genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interactions}{
  an interactionData object of paired genomic regions
}
  \item{motifs}{
  a TFBS tools matrix of DNA binding motifs
}
  \item{genome}{
  a Biostrings genome must match chromosomes from interaction data file
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an interaction data object where obj$anchorOneMotifs and obj$anchorTwoMotifs contain information about the scores and matches to motifs from anchor one and anchor two of interaction data genomic regions
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
function (interactions, motifs, genome)
{
    anchorOneMatches <- matchMotifs(motifs, anchorOne(interactions),
        genome = genome, out = "scores")
    anchorTwoMatches <- matchMotifs(motifs, anchorTwo(interactions),
        genome = genome, out = "scores")
    interactionData <- list(interactions = interactions, anchorOneMotifs = anchorOneMatches,
        anchorTwoMotifs = anchorTwoMatches, anchorOneMotifIndices = seq(length(anchorOneMatches$name)),
        anchorTwoMotifIndices = seq(length(anchorTwoMatches$name)))
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
