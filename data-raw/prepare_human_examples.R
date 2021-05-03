process_cell_type <- function(df, cell_type, motifs_file) {
  df <- dplyr::filter(df, !!rlang::sym(cell_type) == "Yes")
  df <- df[, 1:6]
  colnames(df) <- c("chr_a", "start_a", "end_a",
                    "chr_b", "start_b", "end_b")
  res <- spatzie::find_ep_coenrichment(df, motifs_file, genome_id = "hg19",
                                       cooccurrence_method = "score")
  int_data <- spatzie::filter_pair_motifs(res$motif_cooccurrence)
  int_data$interactions <- NA
  int_data$anchor1_motifs <- NA
  int_data$anchor2_motifs <- NA
  int_data$anchor1_motif_indices <- NA
  int_data$anchor2_motif_indices <- NA
  int_data$pair_motif_enrich_sig <- NA
  int_data$pair_motif_scores <- NA

  return(spatzie::filter_pair_motifs(res$motif_cooccurrence))
}

# data origin of 'interactions_file'
# publication DOI: https://doi.org/10.1038/s41586-020-2151-x
# item: Supplementary Table 4 | Pan-cell type cohesin-mediated chromatin loops
interactions_file <- "data-raw/41586_2020_2151_MOESM5_ESM.bedpe.txt.gz"

# data origin of 'motifs_file'
# website URL: https://hocomoco11.autosome.ru/downloads_v11
# item: H11_HUMAN_mono_jaspar_format.txt
motifs_file <- system.file("extdata/HOCOMOCOv11_core_HUMAN_mono.txt.gz",
                           package = "spatzie")

df <- read.table(gzfile(interactions_file), sep = "\t", header = TRUE)

int_data_mslcl <- process_cell_type(df, "MSLCL", motifs_file)
save(int_data_mslcl, file = "data/int_data_mslcl.rda", compress = "xz")

int_data_k562 <- process_cell_type(df, "K562", motifs_file)
save(int_data_k562, file = "data/int_data_k562.rda", compress = "xz")
