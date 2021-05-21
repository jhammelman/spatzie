#' @title Scans interaction file for motif instances
#'
#' @description
#' Uses motifmatchR to scan interaction regions for given motifs.
#'
#' @param interactions an interactionData object of paired genomic regions
#' @param motifs a TFBS tools matrix of DNA binding motifs
#' @param genome BSgenome object or DNAStringSet object, must match chromosomes
#' from interaction data file
#' @return an interaction data object where \code{obj$anchor1_motifs} and
#' \code{obj$anchor2_motifs} contain information about the scores and matches
#' to motifs from anchor one and anchor two of interaction data genomic regions
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
#' }
#'
#' motifs_file <- system.file("extdata/motifs_subset.txt.gz",
#'                            package = "spatzie")
#' motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
#' left <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr1", "chr1"),
#'   ranges = IRanges::IRanges(start = c(1, 15, 20),
#'                             end = c(10, 35, 31)))
#' right <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr2", "chr2"),
#'   ranges = IRanges::IRanges(start = c(17, 47, 41),
#'                             end = c(28, 54, 53)))
#' test_interactions <- GenomicInteractions::GenomicInteractions(left, right)
#'
#' # toy DNAStringSet to replace BSgenome object
#' seqs <- c("chr1" = "CCACTAGCCACGCGTCACTGGTTAGCGTGATTGAAACTAAATCGTATGAAAATCC",
#'           "chr2" = "CTACAAACTAGGAATTTAGGCAAACCTGTGTTAAAATCTTAGCTCATTCATTAAT")
#' toy_genome <- Biostrings::DNAStringSet(seqs, use.names = TRUE)
#'
#' res <- scan_motifs(test_interactions, motifs, toy_genome)
#'
#' @author Jennifer Hammelman
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicInteractions anchorTwo
#' @export
scan_motifs <- function(interactions, motifs, genome) {
  anchor1_matches <- motifmatchr::matchMotifs(
    motifs, GenomicInteractions::anchorOne(interactions),
    genome = genome, out = "scores", bg = "subject")
  anchor2_matches <- motifmatchr::matchMotifs(
    motifs, GenomicInteractions::anchorTwo(interactions),
    genome = genome, out = "scores", bg = "subject")
  interaction_data <- list(
    interactions = interactions,
    anchor1_motifs = anchor1_matches,
    anchor2_motifs = anchor2_matches,
    anchor1_motif_indices = seq(length(anchor1_matches$name)),
    anchor2_motif_indices = seq(length(anchor2_matches$name)))
  class(interaction_data) <- "interactionData"
  return(interaction_data)
}
