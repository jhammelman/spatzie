#' @title Plot motif pair occurrence
#'
#' @description
#' Creates heatmap to plot values associated with interaction pairs
#'
#' @param pair_enrichment_matrix matrix of format \code{obj$pair_motif_enrich}
#' or \code{obj$pair_motif_enrich_sig}
#' @return heatmap of input data
#'
#' @examples
#' genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
#' motif_file <- system.file(
#'   "extdata/consensus_HOCOMOCOv11_core_MOUSE-plus_YY1.piq",
#'   package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motif_file, matrixClass = "PFM")
#'
#' yy1_pd_interaction <- scan_motifs(spatzie:::yy1_interactions, motifs, genome)
#' yy1_pd_interaction <- filter_motifs(yy1_pd_interaction, 0.4)
#' yy1_pd_score_corr <- anchor_pair_enrich(yy1_pd_interaction,
#'                                         method = "scoreCorrelation")
#' plot_motif_pairs_heatmap(-log2(yy1_pd_score_corr$pair_motif_enrich))
#'
#' @author Jennifer Hammelman
#' @importFrom pheatmap pheatmap
#' @export
plot_motif_pairs_heatmap <- function(pair_enrichment_matrix) {
  return(pheatmap::pheatmap(pair_enrichment_matrix, fontsize = 6))
}
