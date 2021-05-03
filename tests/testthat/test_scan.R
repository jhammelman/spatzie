context("scan_motifs")
data("scan_interactions_example", package = "spatzie")
test_that("can scan motifs with motifmatchr", {
  motifs_file <- system.file("extdata/motifs_subset.txt.gz",
                             package = "spatzie")
  motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
  left <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100000, 100550, 101050),
                              end = c(100550, 101050, 101600)))
  right <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2", "chr2"),
    ranges = IRanges::IRanges(start = c(200000, 200550, 201050),
                              end = c(200550, 201050, 201600)))
  test_interactions <- GenomicInteractions::GenomicInteractions(left, right)
  scanned_interactions <- scan_motifs(
    test_interactions, motifs,
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  expect_is(scanned_interactions, "interactionData")
  expect_is(scanned_interactions$anchor1_motifs, "RangedSummarizedExperiment")
  expect_is(scanned_interactions$anchor2_motifs, "RangedSummarizedExperiment")
  data("scan_interactions_example", package = "spatzie")
  expect_equal(scanned_interactions, scan_interactions_example)
})
