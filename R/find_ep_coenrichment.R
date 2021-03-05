#' @title Find co-enriched motif pairs in enhancer-promoter interactions
#'
#' @description
#' Identifies co-enriched pairs of motifs in enhancer-promoter interactions
#' selected from a data frame of general genomic interactions.
#'
#' Calls functions \code{\link{scan_motifs}}, \code{\link{filter_motifs}},
#' and \code{\link{anchor_pair_enrich}} internally.
#'
#' Uses \code{biomaRt} to retrieve annotations.
#'
#' @param int_raw_data a data frame with at least six columns:
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
#' @param piq_motif_file PIQ file containing motifs to scan for
#' @param genome_id ID of genome assembly interactions in \code{int_raw_data}
#' were aligned to, valid options include \code{hg19}, \code{hg38}, \code{mm9},
#' and \code{mm10}, defaults to \code{hg38}
#' @param cooccurrence_method choice of method for co-occurrence include
#' \code{countCorrelation}, \code{scoreCorrelation}, \code{countHypergeom}, or
#' \code{countFisher}, see \code{\link{anchor_pair_enrich}}, defaults to
#' \code{countCorrelation}
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
#' interactions_file <- system.file("extdata/yy1_interactions.bedpe",
#'                                  package = "spatzie")
#' motifs_file <- system.file(
#'   "extdata/consensus_HOCOMOCOv11_core_MOUSE-plus_YY1.piq",
#'   package = "spatzie")
#'
#' df <- read.table(interactions_file, header = FALSE, sep = "\t")
#' res <- find_ep_coenrichment(df, motifs_file, genome_id = "mm10")
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
#' @importFrom biomaRt useMart
#' @importFrom GenomicInteractions annotateInteractions
#' @importFrom GenomicInteractions plotInteractionAnnotations
#' @importFrom GenomicInteractions isInteractionType
#' @importFrom GenomicInteractions anchorOne
#' @importFrom GenomicInteractions anchorTwo
#' @importFrom S4Vectors elementMetadata
#' @importFrom TFBSTools readJASPARMatrix
#' @importFrom S4Vectors elementMetadata
#' @importFrom S4Vectors elementMetadata
#' @export
find_ep_coenrichment <- function(int_raw_data,
                                 piq_motif_file,
                                 genome_id = c("hg38", "hg19", "mm9", "mm10"),
                                 cooccurrence_method = c("countCorrelation",
                                                         "scoreCorrelation",
                                                         "countHypergeom",
                                                         "countFisher"),
                                 filter_threshold = 0.4) {
  genome_id <- match.arg(genome_id, c("hg19", "hg38", "mm9", "mm10"))

  if (genome_id == "hg19") {
    check_required_package("BSgenome.Hsapiens.UCSC.hg19")
    check_required_package("TxDb.Hsapiens.UCSC.hg19.knownGene")
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    ensembl_data_set <- "hsapiens_gene_ensembl"
    gene_symbol <- "hgnc_symbol"
  } else if (genome_id == "hg38") {
    check_required_package("BSgenome.Hsapiens.UCSC.hg38")
    check_required_package("TxDb.Hsapiens.UCSC.hg38.knownGene")
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    ensembl_data_set <- "hsapiens_gene_ensembl"
    gene_symbol <- "hgnc_symbol"
  } else if (genome_id == "mm9") {
    check_required_package("BSgenome.Mmusculus.UCSC.mm9")
    check_required_package("TxDb.Mmusculus.UCSC.mm9.knownGene")
    genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
    ensembl_data_set <- "mmusculus_gene_ensembl"
    gene_symbol <- "mgi_symbol"
  } else if (genome_id == "mm10") {
    check_required_package("BSgenome.Mmusculus.UCSC.mm10")
    check_required_package("TxDb.Mmusculus.UCSC.mm10.knownGene")
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    ensembl_data_set <- "mmusculus_gene_ensembl"
    gene_symbol <- "mgi_symbol"
  }

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

  promoter_ranges <- GenomicFeatures::promoters(
    txdb, upstream = 2500, downstream = 2500, columns = c("tx_name", "gene_id"))

  # trims out-of-bound ranges located on non-circular sequences
  promoter_ranges <- GenomicRanges::trim(promoter_ranges)

  # remove duplicate promoters from transcript isoforms
  promoter_ranges <- BiocGenerics::unique(promoter_ranges)

  promoters_df <- as.data.frame(promoter_ranges)
  promoters_df$gene_id <- as.character(promoters_df$gene_id)

  ensembl <- biomaRt::useMart("ensembl", dataset = ensembl_data_set)

  id_df <- biomaRt::getBM(attributes = c("entrezgene_id", gene_symbol),
                          filters = "entrezgene_id",
                          values = unique(promoters_df$gene_id),
                          mart = ensembl)

  names(promoter_ranges) <- id_df[match(promoters_df$gene_id,
                                        id_df$entrezgene_id), gene_symbol]
  missing_idx <- is.na(names(promoter_ranges)) | names(promoter_ranges) == ""
  names(promoter_ranges)[missing_idx] <- promoters_df$tx_name[missing_idx]

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

  motifs <- TFBSTools::readJASPARMatrix(piq_motif_file, matrixClass = "PFM")
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
    stop(paste0("Bioconductor package missing (", package_name,
                "), install with command: BiocManager::install(\"",
                package_name, "\")"))
  }
}
