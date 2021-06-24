#' @title Find co-enriched motif pairs in enhancer-promoter interactions
#'
#' @description
#' Identifies co-enriched pairs of motifs in enhancer-promoter interactions
#' selected from a data frame of general genomic interactions.
#'
#' If \code{identify_ep}: Promoters and enhancers are identified
#' using genomic annotations, where anchors close to promoter annotations
#' (within 2500 base pairs) are considered promoters and all other anchors are
#' considered gene-distal enhancers. Only interactions in
#' \code{int_raw_data} between promoters and enhancers are used for motif
#' co-enrichment analysis.
#'
#' If \code{!identify_ep}: Instead of automatically identifying
#' promoters and enhancers based on genomic annotations, all interactions
#' in \code{int_raw_data} must be preprocessed in a way that anchor 1 contains
#' promoters and anchor 2 contains enhancers. Motif
#' co-enrichment analysis is performed under this assumption.
#'
#' Calls functions \code{\link{scan_motifs}}, \code{\link{filter_motifs}},
#' and \code{\link{anchor_pair_enrich}} internally.
#'
#' @param int_raw_data a \code{\link[GenomicInteractions]{GenomicInteractions}}
#' object or a data frame with at least six columns:
#' \tabular{rl}{
#'   column 1: \tab character; genomic location of interaction anchor 1 -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab integer; genomic location of interaction anchor 1 -
#'   start coordinate\cr
#'   column 3: \tab integer; genomic location of interaction anchor 1 -
#'   end coordinate\cr
#'   column 4: \tab character; genomic location of interaction anchor 2 -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5: \tab integer; genomic location of interaction anchor 2 -
#'   start coordinate\cr
#'   column 6: \tab integer; genomic location of interaction anchor 2 -
#'   end coordinate
#' }
#' @param motifs_file \href{http://jaspar.genereg.net/faq/}{JASPAR format}
#' matrix file containing multiple motifs to scan for, gz-zipped files allowed
#' @param motifs_file_matrix_format type of position-specific scoring matrices
#' in \code{motifs_file}, valid options include:
#' \tabular{rl}{
#'   \code{pfm}: \tab position frequency matrix, elements are absolute
#'   frequencies, i.e., counts (default)\cr
#'   \code{ppm}: \tab position probability matrix, elements are probabilities,
#'   i.e., Laplace smoothing corrected relative frequencies\cr
#'   \code{pwm}: \tab position weight matrix, elements are log likelihoods
#' }
#' @param genome_id ID of genome assembly interactions in \code{int_raw_data}
#' were aligned to, valid options include \code{hg19}, \code{hg38}, \code{mm9},
#' and \code{mm10}, defaults to \code{hg38}
#' @param identify_ep logical, set \code{FALSE} if enhancers and promoters
#' should not be identified based on genomic annotations, but instead
#' assumes anchor 1 contains promoters and anchor 2 contains enhancers,
#' for all interactions in \code{int_raw_data}, defaults to \code{TRUE}, i.e.,
#' do identify enhancers and promoters of interactions in \code{int_raw_data}
#' based on genomic interactions and filter all interactions which are not
#' between promoters and enhancers
#' @param cooccurrence_method method for co-occurrence, valid options include:
#' \tabular{rl}{
#'   \code{count}: \tab correlation between counts (for each anchor, tally
#'   positions where motif score > \eqn{5 * 10^{-5}})\cr
#'   \code{score}: \tab correlation between motif scores (for each anchor, use
#'   the maximum score over all positions)\cr
#'   \code{match}: \tab association between motif matches (for each anchor,
#'   a match is defined if the is at least one position with a motif score
#'   > \eqn{5 * 10^{-5}})
#' }
#' See \code{\link{anchor_pair_enrich}} for details.
#' @param filter_threshold fraction of interactions that should contain a
#' motif for a motif to be considered, see \code{\link{filter_motifs}},
#' defaults to \code{0.4}
#'
#' @return a list with the following items:
#' \tabular{rl}{
#'   \code{int_data} \tab
#'   \code{\link[GenomicInteractions]{GenomicInteractions}} object;
#'   promoter-enhancer interactions\cr
#'  \code{int_data_motifs}: \tab \code{interactionData} object; return value of
#'  \code{\link{scan_motifs}}\cr
#'   \code{filtered_int_data_motifs}: \tab \code{interactionData} object;
#'   return value of \code{\link{filter_motifs}}\cr
#'   \code{annotation_pie_chart}: \tab ggplot2 plot; return value of
#'   \code{\link[GenomicInteractions]{plotInteractionAnnotations}}\cr
#'   \code{motif_cooccurrence}: \tab \code{interactionData} object; return
#'   value of \code{\link{anchor_pair_enrich}}
#' }
#'
#' @examples
#' \dontrun{
#' interactions_file <- system.file("extdata/yy1_interactions.bedpe.gz",
#'                                  package = "spatzie")
#' motifs_file <- system.file("extdata/motifs_subset.txt.gz",
#'                            package = "spatzie")
#'
#' df <- read.table(gzfile(interactions_file), header = TRUE, sep = "\t")
#' res <- find_ep_coenrichment(df, motifs_file,
#'                             motifs_file_matrix_format = "pfm",
#'                             genome_id = "mm10")
#' }
#'
#' @author Jennifer Hammelman
#' @author Konstantin Krismer
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics unique
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicFeatures promoters
#' @importFrom GenomicInteractions annotateInteractions
#' @importFrom GenomicInteractions plotInteractionAnnotations
#' @importFrom GenomicInteractions isInteractionType
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicInteractions anchorTwo
#' @importFrom S4Vectors elementMetadata
#' @importFrom TFBSTools readJASPARMatrix
#' @importFrom S4Vectors elementMetadata
#' @importFrom S4Vectors elementMetadata
#' @importFrom BSgenome getBSgenome
#' @export
find_ep_coenrichment <- function(int_raw_data,
                                 motifs_file,
                                 motifs_file_matrix_format = c("pfm", "ppm",
                                                               "pwm"),
                                 genome_id = c("hg38", "hg19", "mm9", "mm10"),
                                 identify_ep = TRUE,
                                 cooccurrence_method = c("count",
                                                         "score",
                                                         "match"),
                                 filter_threshold = 0.4) {
  motifs_file_matrix_format <- match.arg(motifs_file_matrix_format,
                                         c("pfm", "ppm", "pwm"))
  genome_id <- match.arg(genome_id, c("hg19", "hg38", "mm9", "mm10"))

  if (inherits(int_raw_data, "GenomicInteractions")) {
    int_data <- int_raw_data
  } else if (inherits(int_raw_data, "data.frame")) {
    left_anchor <- GenomicRanges::GRanges(
      seqnames = int_raw_data[, 1],
      ranges = IRanges::IRanges(start = int_raw_data[, 2],
                                end = int_raw_data[, 3]),
      seqinfo = GenomeInfoDb::Seqinfo(genome = genome_id))
    right_anchor <- GenomicRanges::GRanges(
      seqnames = int_raw_data[, 4],
      ranges = IRanges::IRanges(start = int_raw_data[, 5],
                                end = int_raw_data[, 6]),
      seqinfo = GenomeInfoDb::Seqinfo(genome = genome_id))

    # trims out-of-bound ranges located on non-circular sequences
    left_anchor <- GenomicRanges::trim(left_anchor)
    right_anchor <- GenomicRanges::trim(right_anchor)

    int_data <- GenomicInteractions::GenomicInteractions(left_anchor,
                                                         right_anchor)
  } else {
    stop("'int_raw_data' data type unsupported: ", class(int_raw_data))
  }

  if (identify_ep) {
    if (genome_id == "hg19") {
      check_required_package("BSgenome.Hsapiens.UCSC.hg19")
      check_required_package("TxDb.Hsapiens.UCSC.hg19.knownGene")

      genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      ensembl_data_set <- "hsapiens_gene_ensembl"
      gene_symbol <- "hgnc_symbol"
    } else if (genome_id == "hg38") {
      check_required_package("BSgenome.Hsapiens.UCSC.hg38")
      check_required_package("TxDb.Hsapiens.UCSC.hg38.knownGene")

      genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      ensembl_data_set <- "hsapiens_gene_ensembl"
      gene_symbol <- "hgnc_symbol"
    } else if (genome_id == "mm9") {
      check_required_package("BSgenome.Mmusculus.UCSC.mm9")
      check_required_package("TxDb.Mmusculus.UCSC.mm9.knownGene")

      genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm9")
      txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
      ensembl_data_set <- "mmusculus_gene_ensembl"
      gene_symbol <- "mgi_symbol"
    } else if (genome_id == "mm10") {
      check_required_package("BSgenome.Mmusculus.UCSC.mm10")
      check_required_package("TxDb.Mmusculus.UCSC.mm10.knownGene")

      genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
      ensembl_data_set <- "mmusculus_gene_ensembl"
      gene_symbol <- "mgi_symbol"
    }

    promoter_ranges <- GenomicFeatures::promoters(
      txdb, upstream = 2500, downstream = 2500,
      columns = c("tx_name", "gene_id"))

    # trims out-of-bound ranges located on non-circular sequences
    promoter_ranges <- GenomicRanges::trim(promoter_ranges)

    # remove duplicate promoters from transcript isoforms
    promoter_ranges <- BiocGenerics::unique(promoter_ranges)

    promoters_df <- as.data.frame(promoter_ranges)
    promoters_df$gene_id <- as.character(promoters_df$gene_id)

    annotation_features <- list(promoter = promoter_ranges)

    GenomicInteractions::annotateInteractions(int_data, annotation_features)
    annotation_pie_chart <- GenomicInteractions::plotInteractionAnnotations(
      int_data)

    distal_promoter_idx <- GenomicInteractions::isInteractionType(
      int_data, "distal", "promoter")
    int_data <- int_data[distal_promoter_idx]
    anchor1 <- GenomicInteractions::anchorOne(int_data)
    anchor2 <- GenomicInteractions::anchorTwo(int_data)

    promoter_left <- S4Vectors::elementMetadata(anchor1)[, "node.class"] == "promoter"
    promoter_right <- S4Vectors::elementMetadata(anchor2)[, "node.class"] == "promoter"
    promoter_ranges <- c(anchor1[promoter_left],
                         anchor2[promoter_right])
    enhancer_ranges <- c(anchor2[promoter_left],
                         anchor1[promoter_right])
    int_data <- GenomicInteractions::GenomicInteractions(promoter_ranges,
                                                         enhancer_ranges)
  } else {
    if (genome_id == "hg19") {
      check_required_package("BSgenome.Hsapiens.UCSC.hg19")
      genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
    } else if (genome_id == "hg38") {
      check_required_package("BSgenome.Hsapiens.UCSC.hg38")
      genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
    } else if (genome_id == "mm9") {
      check_required_package("BSgenome.Mmusculus.UCSC.mm9")
      genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm9")
    } else if (genome_id == "mm10") {
      check_required_package("BSgenome.Mmusculus.UCSC.mm10")
      genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
    }
    annotation_pie_chart <- NULL
  }

  if (motifs_file_matrix_format == "pfm") {
    jaspar_matrix_class <- "PFM"
  } else if (motifs_file_matrix_format == "ppm") {
    jaspar_matrix_class <- "PWMProb"
  } else if (motifs_file_matrix_format == "pwm") {
    jaspar_matrix_class <- "PWM"
  }

  motifs <- TFBSTools::readJASPARMatrix(motifs_file,
                                        matrixClass = jaspar_matrix_class)
  int_data_motifs <- scan_motifs(int_data, motifs, genome)
  filtered_int_data_motifs <- filter_motifs(int_data_motifs, filter_threshold)

  motif_cooccurrence <- anchor_pair_enrich(
    filtered_int_data_motifs, method = cooccurrence_method)

  return(list(int_data = int_data,
              int_data_motifs = int_data_motifs,
              filtered_int_data_motifs = filtered_int_data_motifs,
              annotation_pie_chart = annotation_pie_chart,
              motif_cooccurrence = motif_cooccurrence))
}

#' @importFrom utils installed.packages
check_required_package <- function(package_name) {
  if (!(package_name %in% rownames(utils::installed.packages()))) {
    stop("Bioconductor package missing (", package_name,
         "), install with command: BiocManager::install(\"",
         package_name, "\")")
  }
}
