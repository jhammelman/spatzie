motif_file <- system.file(
  "extdata/motifs_subset.txt.gz", package = "spatzie")
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
anchor_pair_example_scorecorr <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "score")
anchor_pair_example_countcorr <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "count")
anchor_pair_example_counthyper <- spatzie::anchor_pair_enrich(
  scan_interactions_example_filtered, method = "match")

save(scan_interactions_example,
     file = "data/scan_interactions_example.rda", compress = "xz")
save(scan_interactions_example_filtered,
     file = "data/scan_interactions_example_filtered.rda", compress = "xz")
save(anchor_pair_example_scorecorr,
     file = "data/anchor_pair_example_scorecorr.rda", compress = "xz")
save(anchor_pair_example_countcorr,
     file = "data/anchor_pair_example_countcorr.rda", compress = "xz")
save(anchor_pair_example_counthyper,
     file = "data/anchor_pair_example_counthyper.rda", compress = "xz")

filter_pairs_example <- spatzie::filter_pair_motifs(
  anchor_pair_example_countcorr, threshold = 0.5)
save(filter_pairs_example,
     file = "data/filter_pairs_example.rda", compress = "xz")
