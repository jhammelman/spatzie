set.seed(17)

# data origin of 'interactions_file'
# publication DOI: https://doi.org/10.1016/j.cell.2017.11.008
# GEO ID: GSM2645440
# item: GSM2645440_YY1_rep1_rep2_merged.origami_v1.1.ints_all.csv.gz
# 40,000 interaction subset of original file
interactions_file <- system.file("extdata/yy1_interactions.bedpe.gz",
                                 package = "spatzie")

# data origin of 'motifs_file'
# website URL: https://hocomoco11.autosome.ru/downloads_v11
# item: H11_MOUSE_mono_jaspar_format.txt
# subset
motifs_file <- system.file("extdata/motifs_subset.txt.gz",
                           package = "spatzie")

df <- read.table(gzfile(interactions_file), sep = "\t", header = FALSE)

df <- df[, 1:6]
colnames(df) <- c("chr_a", "start_a", "end_a",
                  "chr_b", "start_b", "end_b")
res <- spatzie::find_ep_coenrichment(df, motifs_file, genome_id = "mm9",
                                     cooccurrence_method = "score")

interactions_yy1 <- res$int_data
save(interactions_yy1, file = "data/interactions_yy1.rda", compress = "xz")

df <- df[sample(nrow(df), 5000), ]

res <- spatzie::find_ep_coenrichment(df, motifs_file, genome_id = "mm9",
                                     cooccurrence_method = "score")

interactions_yy1_promoter <- spatzie::get_specific_interactions(
  res$int_data_motifs, anchor1_motif = "YY1")
save(interactions_yy1_promoter, file = "data/interactions_yy1_promoter.rda",
     compress = "xz")

interactions_yy1_enhancer <- spatzie::get_specific_interactions(
  res$int_data_motifs, anchor2_motif = "YY1")
save(interactions_yy1_enhancer, file = "data/interactions_yy1_enhancer.rda",
     compress = "xz")

interactions_yy1_ep <- spatzie::get_specific_interactions(
  res$int_data_motifs, anchor1_motif = "YY1", anchor2_motif = "YY1")
save(interactions_yy1_ep, file = "data/interactions_yy1_ep.rda", compress = "xz")

int_data_yy1 <- res$int_data_motifs
save(int_data_yy1, file = "data/int_data_yy1.rda", compress = "xz")
