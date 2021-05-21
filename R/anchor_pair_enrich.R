#' @title Determine enriched motifs in anchors
#'
#' @description
#' Determine whether motifs between paired bed regions have a statistically
#' significant relationship. Options for significance are motif score
#' correlation, motif count correlation, or hypergeometric motif co-occurrence.
#'
#' @param interaction_data an interactionData object of paired genomic regions
#' @param method method for co-occurrence, valid options include:
#' \tabular{rl}{
#'   \code{count}: \tab correlation between counts (for each anchor, tally
#'   positions where motif score > \eqn{5 * 10^{-5}})\cr
#'   \code{score}: \tab correlation between motif scores (for each anchor, use
#'   the maximum score over all positions)\cr
#'   \code{match}: \tab association between motif matches (for each anchor,
#'   a match is defined if the is at least one position with a motif score
#'   > \eqn{5 * 10^{-5}})
#' }
#' @return an interactionData object where \code{obj$pair_motif_enrich} contains
#' the p-values for significance of seeing a higher co-occurrence than
#' what we get by chance.
#'
#' @section Score-based correlation:
#'
#' We assume motif scores follow a normal distribution and are independent
#' between enhancers and promoters. We can therefore compute how correlated
#' scores of any two transcription factor motifs are between enhancer and
#' promoter regions using Pearson's product-moment correlation coefficient:
#' \deqn{r = \frac{\sum (x^{\prime}_i - \bar{x}^{\prime})(y^{\prime}_i -
#' \bar{y}^{\prime})}{\sqrt{\sum(x^{\prime}_i -
#' \bar{x}^{\prime})^2\sum(y^{\prime}_i - \bar{y}^{\prime})^2}}},
#' where the input vectors \eqn{\boldsymbol{x}} and \eqn{\boldsymbol{y}} from
#' above are transformed to vectors \eqn{\boldsymbol{x^{\prime}}} and
#' \eqn{\boldsymbol{y^{\prime}}} by replacing the set of scores with the
#' maximum score for each region:
#' \deqn{x^{\prime}_i = \max x_i}
#' \eqn{x^{\prime}_i} is then the maximum motif score of motif \eqn{a} in the
#' promoter region of interaction \eqn{i}, \eqn{y^{\prime}_i} is the maximum
#' motif score of motif \eqn{b} in the enhancer region of interaction \eqn{i},
#' and \eqn{\bar{x}^{\prime}} and \eqn{\bar{y}^{\prime}} are the sample means.
#'
#' Significance is then computed by transforming the correlation coefficient
#' \eqn{r} to test statistic \eqn{t}, which is Student \eqn{t}-distributed
#' with \eqn{n - 2} degrees of freedom.
#' \deqn{t = \frac{r\sqrt{n-2}}{\sqrt{1-r^2}}}
#'
#' All p-values are calculated as one-tailed p-values of the probability
#' that scores are greater than or equal to \eqn{r}.
#'
#' @section Count-based correlation:
#'
#' Instead of calculating the correlation of motif scores directly, the
#' count-based correlation metric first tallies the number of instances of a
#' given motif within an enhancer or a promoter region, which are defined as
#' all positions in those regions with motif score p-values of less than
#' \eqn{5 * 10^{-5}}. Formally, the input vectors \eqn{\boldsymbol{x}} and
#' \eqn{\boldsymbol{y}} are transformed to vectors
#' \eqn{\boldsymbol{x^{\prime\prime}}} and \eqn{\boldsymbol{y^{\prime\prime}}}
#' by replacing the set of scores with the cardinality of the set:
#' \deqn{x^{\prime\prime}_i = |x_i|}
#' And analogous for \eqn{y^{\prime\prime}_i}. Finally, the correlation
#' coefficient \eqn{r} between \eqn{\boldsymbol{x^{\prime\prime}}} and
#' \eqn{\boldsymbol{y^{\prime\prime}}} and its associated significance are
#' calculated as described above.
#'
#' @section Match-based association:
#'
#' Instance co-occurrence uses the presence or absence of a motif within an
#' enhancer or promoter to determine a statistically significant association,
#' thus \eqn{\boldsymbol{x^{\prime\prime\prime}}} and
#' \eqn{\boldsymbol{y^{\prime\prime\prime}}} are defined by:
#' \deqn{x^{\prime\prime\prime}_i = \boldsymbol{1}_{x^{\prime\prime}_i > 0}}
#'
#' Instance co-occurrence is computed using the hypergeometric test:
#' \deqn{p = \sum_{k=I_{ab}}^{P_a} \frac{\text{binom}(P_a, k)
#' \text{binom}(n - P_a, E_b - k)}{\text{binom}(n, E_b)},}
#' where \eqn{I_{ab}} is the number of interactions that contain a match for
#' motif \eqn{a} in the promoter and motif \eqn{b} in the enhancer, \eqn{P_a}
#' is the number of promoters that contain motif \eqn{a}
#' (\eqn{P_a = \sum^n_i x^{\prime\prime\prime}_i}), \eqn{E_b} is the number of
#' enhancers that contain motif \eqn{b}
#' (\eqn{E_b = \sum^n_i y^{\prime\prime\prime}_i}), and \eqn{n} is the total
#' number of interactions, which is equal to the number of promoters and to
#' the number of enhancers.
#'
#' @examples
#' \dontrun{
#' genome_id <- "BSgenome.Mmusculus.UCSC.mm9"
#' if (!(genome_id %in% rownames(utils::installed.packages()))) {
#'   BiocManager::install(genome_id, update = FALSE, ask = FALSE)
#' }
#' genome <- BSgenome::getBSgenome(genome_id)
#'
#' motifs_file <- system.file("extdata/motifs_subset.txt.gz",
#'                            package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
#'
#' yy1_pd_interaction <- scan_motifs(spatzie::interactions_yy1, motifs, genome)
#' yy1_pd_interaction <- filter_motifs(yy1_pd_interaction, 0.4)
#' yy1_pd_count_corr <- anchor_pair_enrich(yy1_pd_interaction, method = "count")
#' }
#'
#' res <- anchor_pair_enrich(spatzie::scan_interactions_example_filtered,
#'                           method = "score")
#'
#' @author Jennifer Hammelman
#' @author Konstantin Krismer
#' @importFrom stats cor.test
#' @importFrom SummarizedExperiment assays
#' @importFrom stats phyper
#' @export
anchor_pair_enrich <- function(interaction_data,
                               method = c("count", "score", "match")) {
  method <- match.arg(method, c("count", "score", "match"))
  significance <- matrix(data = NA,
                         nrow = length(interaction_data$anchor1_motif_indices),
                         ncol = length(interaction_data$anchor2_motif_indices))
  values <- matrix(data = NA,
                   nrow = length(interaction_data$anchor1_motif_indices),
                   ncol = length(interaction_data$anchor2_motif_indices))
  indr <- 1
  anchor1_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor1_motifs)
  anchor2_motifs <- SummarizedExperiment::assays(
    interaction_data$anchor2_motifs)
  for (i in interaction_data$anchor1_motif_indices) {
    indc <- 1
    for (j in interaction_data$anchor2_motif_indices) {
      if (method == "count") {
        significance[indr, indc] <- stats::cor.test(
          anchor1_motifs$motifCounts[, i],
          anchor2_motifs$motifCounts[, j],
          alternative = "greater", method = "pearson")$p.value
        values[indr, indc] <- stats::cor(
          anchor1_motifs$motifCounts[, i],
          anchor2_motifs$motifCounts[, j])
      } else if (method == "score") {
        significance[indr, indc] <- stats::cor.test(
          anchor1_motifs$motifScores[, i],
          anchor2_motifs$motifScores[, j],
          alternative = "greater", method = "pearson")$p.value
        values[indr, indc] <- stats::cor(
          anchor1_motifs$motifScores[, i],
          anchor2_motifs$motifScores[, j])
      } else if (method == "match") {
        significance[indr, indc] <- stats::phyper(
          sum((anchor1_motifs$motifMatches[, i]) *
                (anchor2_motifs$motifMatches[, j])),
          sum(anchor1_motifs$motifMatches[, i]),
          length(anchor1_motifs$motifMatches[, i]) -
            sum(anchor1_motifs$motifMatches[, i]),
          sum(anchor2_motifs$motifMatches[, j]), lower.tail = FALSE)
        values[indr, indc] <- sum((anchor1_motifs$motifMatches[, i]) *
                                    (anchor2_motifs$motifMatches[, j]))
        max_ep <- min(sum(anchor1_motifs$motifMatches[, i]),
                      sum(anchor2_motifs$motifMatches[, j]))
        values[indr, indc] <- values[indr, indc] / max_ep
      }
      indc <- indc + 1
    }
    indr <- indr + 1
  }
  rownames(significance) <- names(interaction_data$anchor1_motif_indices)
  colnames(significance) <- names(interaction_data$anchor2_motif_indices)
  rownames(values) <- names(interaction_data$anchor1_motif_indices)
  colnames(values) <- names(interaction_data$anchor2_motif_indices)

  interaction_data <- list(
    interactions = interaction_data$interactions,
    anchor1_motifs = interaction_data$anchor1_motifs,
    anchor2_motifs = interaction_data$anchor2_motifs,
    anchor1_motif_indices = interaction_data$anchor1_motif_indices,
    anchor2_motif_indices = interaction_data$anchor2_motif_indices,
    pair_motif_scores = values,
    pair_motif_enrich = significance,
    is_multiple_hypothesis_corrected = FALSE)
  class(interaction_data) <- "interactionData"
  return(interaction_data)
}
