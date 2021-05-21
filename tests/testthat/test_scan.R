context("scan_motifs")
data("scan_interactions_example", package = "spatzie")
test_that("can scan motifs with motifmatchr", {
  motifs_file <- system.file("extdata/motifs_subset.txt.gz",
                             package = "spatzie")
  motifs <- TFBSTools::readJASPARMatrix(motifs_file, matrixClass = "PFM")
  left <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 15, 20),
                              end = c(10, 35, 31)))
  right <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2", "chr2"),
    ranges = IRanges::IRanges(start = c(17, 47, 41),
                              end = c(28, 54, 53)))
  test_interactions <- GenomicInteractions::GenomicInteractions(left, right)

  # toy DNAStringSet to replace BSgenome object
  seqs <- c("chr1" = "CCACTAGCCACGCGTCACTGGTTAGCGTGATTGAAACTAAATCGTATGAAAATCC",
            "chr2" = "CTACAAACTAGGAATTTAGGCAAACCTGTGTTAAAATCTTAGCTCATTCATTAAT")
  toy_genome <- Biostrings::DNAStringSet(seqs, use.names = TRUE)

  scanned_interactions <- scan_motifs(test_interactions, motifs, toy_genome)

  expect_is(scanned_interactions, "interactionData")
  expect_is(scanned_interactions$anchor1_motifs, "RangedSummarizedExperiment")
  expect_is(scanned_interactions$anchor2_motifs, "RangedSummarizedExperiment")
})
