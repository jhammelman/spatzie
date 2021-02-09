context("anchor_pair_enrich")
data("scan_interactions_example_filtered", package = "spatzie")
data("anchor_pair_example_scorecorr", package = "spatzie")

test_that("can compute significance with score correlation",{
  anchor_pair_scorecorr <- anchor_pair_enrich(scan_interactions_example_filtered,
                                              method="scoreCorrelation")
  expect_equal(anchor_pair_example_scorecorr$pair_motif_scores,
               anchor_pair_scorecorr$pair_motif_scores)
  expect_equal(anchor_pair_example_scorecorr$pair_motif_enrich,
               anchor_pair_scorecorr$pair_motif_enrich)
})

test_that("can compute significance with count correlation",{
  anchor_pair_countcorr <- anchor_pair_enrich(scan_interactions_example_filtered,
                                              method="countCorrelation")
  data("anchor_pair_example_countcorr", package = "spatzie")
  expect_equal(anchor_pair_example_countcorr$pair_motif_scores,
               anchor_pair_countcorr$pair_motif_scores)
  expect_equal(anchor_pair_example_countcorr$pair_motif_enrich,
               anchor_pair_countcorr$pair_motif_enrich)
})

test_that("can compute significance with hypergeometric test",{
  anchor_pair_counthyper <- anchor_pair_enrich(scan_interactions_example_filtered,
                                               method="countHypergeom")
  data("anchor_pair_example_counthyper", package = "spatzie")
  expect_equal(anchor_pair_example_counthyper$pair_motif_scores,
               anchor_pair_counthyper$pair_motif_scores)
  expect_equal(anchor_pair_example_counthyper$pair_motif_enrich,
               anchor_pair_counthyper$pair_motif_enrich)
})

test_that("can compute significance with fisher test",{
  anchor_pair_countfisher <- anchor_pair_enrich(scan_interactions_example_filtered,
                                                        method="countFisher")
  data("anchor_pair_example_countfisher", package = "spatzie")
  expect_equal(anchor_pair_example_countfisher$pair_motif_scores,
               anchor_pair_countfisher$pair_motif_scores)
  expect_equal(anchor_pair_example_countfisher$pair_motif_enrich,
               anchor_pair_countfisher$pair_motif_enrich)
})
