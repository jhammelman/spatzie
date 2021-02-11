genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
ensembl_data_set <- "mmusculus_gene_ensembl"
gene_symbol <- "mgi_symbol"

yy1_interactions_file <- system.file("extdata/yy1_interactions.bedpe",
                                     package = "spatzie")
yy1_interactions <- GenomicInteractions::makeGenomicInteractionsFromFile(
  yy1_interactions_file,
  type = "bedpe",
  experiment_name = "yy1",
  description = "mESC yy1 chr1")

promoter_ranges <- GenomicFeatures::promoters(txdb,
                                              upstream = 2500,
                                              downstream = 2500,
                                              columns = c("tx_name", "gene_id"))
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
GenomicInteractions::annotateInteractions(yy1_interactions, annotation_features)

distal_promoter_idx <- GenomicInteractions::isInteractionType(
  yy1_interactions, "distal", "promoter")
yy1_pd <- yy1_interactions[distal_promoter_idx]
anchor1 <- GenomicInteractions::anchorOne(yy1_pd)
anchor2 <- GenomicInteractions::anchorTwo(yy1_pd)

promoter_left <- S4Vectors::elementMetadata(anchor1)[, "node.class"] == "promoter"
promoter_right <- S4Vectors::elementMetadata(anchor2)[, "node.class"] == "promoter"
promoter_ranges <- c(anchor1[promoter_left],
                     anchor2[promoter_right])
enhancer_ranges <- c(anchor2[promoter_left],
                     anchor1[promoter_right])
yy1_interactions <- GenomicInteractions::GenomicInteractions(promoter_ranges,
                                                             enhancer_ranges)
save(yy1_interactions, file = "data/yy1_interactions.rda", compress = "xz")
yy1p_interactions <- get_specific_interactions(
  yy1_interactions, anchor1_motif = "TYY1_MOUSE.H11MO.0.A")
save(yy1p_interactions, file = "data/yy1p_interactions.rda", compress = "xz")

yy1e_interactions <- get_specific_interactions(
  yy1_interactions, anchor2_motif = "TYY1_MOUSE.H11MO.0.A")
save(yy1e_interactions, file = "data/yy1e_interactions.rda", compress = "xz")

yy1p_yy1e_interactions <- get_specific_interactions(
  yy1_interactions,
  anchor1_motif = "TYY1_MOUSE.H11MO.0.A",
  anchor2_motif = "TYY1_MOUSE.H11MO.0.A")
save(yy1p_yy1e_interactions,
     file = "data/yy1p_yy1e_interactions.rda", compress = "xz")
