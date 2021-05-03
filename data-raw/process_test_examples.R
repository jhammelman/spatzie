motif_file <- system.file("extdata/motifs_subset.txt.gz", package = "spatzie")

motifs <- TFBSTools::readJASPARMatrix(motif_file, matrixClass = "PFM")
left <- GenomicRanges::GRanges(seqnames = c("chr1", "chr1", "chr1"),
                ranges = IRanges::IRanges(start = c(100000, 100550, 101050),
                               end = c(100550, 101050, 101600)))
right <- GenomicRanges::GRanges(seqnames = c("chr1", "chr2", "chr2"),
                 ranges = IRanges::IRanges(start = c(200000, 200550, 201050),
                                end = c(200550, 201050, 201600)))
interactions <- GenomicInteractions::GenomicInteractions(left, right)
scan_interactions_example <- spatzie::scan_motifs(
  interactions, motifs,
  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
scan_interactions_example_filtered <- spatzie::filter_motifs(
  scan_interactions_example, threshold = 0.1)
anchor_pair_example_score <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "score")
anchor_pair_example_count <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "count")
anchor_pair_example_match <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "match")

save(scan_interactions_example,
     file = "data/scan_interactions_example.rda", compress = "xz")
save(scan_interactions_example_filtered,
     file = "data/scan_interactions_example_filtered.rda", compress = "xz")
save(anchor_pair_example_score,
     file = "data/anchor_pair_example_score.rda", compress = "xz")
save(anchor_pair_example_count,
     file = "data/anchor_pair_example_count.rda", compress = "xz")
save(anchor_pair_example_match,
     file = "data/anchor_pair_example_match.rda", compress = "xz")

compare_pairs_example <- spatzie::compare_motif_pairs(
  anchor_pair_example_score, anchor_pair_example_count)
save(compare_pairs_example,
     file = "data/compare_pairs_example.rda", compress = "xz")

filter_pairs_example <- spatzie::filter_pair_motifs(
  anchor_pair_example_count, threshold = 0.5)
save(filter_pairs_example,
     file = "data/filter_pairs_example.rda", compress = "xz")
