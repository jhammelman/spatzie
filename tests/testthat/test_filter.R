context("filter_motifs")
data("scan_interactions_example", package = "spatzie")
data("scan_interactions_example_filtered", package = "spatzie")

test_that("can filter motifs by presence in peaks", {
  scan_interactions_filtered <- filter_motifs(scan_interactions_example,
                                              threshold = 0.1)
  expect_is(scan_interactions_filtered, "interactionData")
  expect_equal(scan_interactions_filtered$anchor1_motif_indices,
               scan_interactions_example_filtered$anchor1_motif_indices)
  expect_equal(scan_interactions_filtered$anchor2_motif_indices,
               scan_interactions_example_filtered$anchor2_motif_indices)
  expect_equal(scan_interactions_filtered$anchor2_motif_indices,
               scan_interactions_example_filtered$anchor2_motif_indices)
})
