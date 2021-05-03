context("anchor_pair_enrich")
data("scan_interactions_example_filtered", package = "spatzie")
data("anchor_pair_example_score", package = "spatzie")

test_that("can compute significance with score correlation", {
  anchor_pair_score <- anchor_pair_enrich(
    scan_interactions_example_filtered, method = "score")
  expect_equal(anchor_pair_example_score$pair_motif_scores,
               anchor_pair_score$pair_motif_scores)
  expect_equal(anchor_pair_example_score$pair_motif_enrich,
               anchor_pair_score$pair_motif_enrich)
})

test_that("can compute significance with count correlation", {
  anchor_pair_count <- anchor_pair_enrich(
    scan_interactions_example_filtered, method = "count")
  data("anchor_pair_example_count", package = "spatzie")
  expect_equal(anchor_pair_example_count$pair_motif_scores,
               anchor_pair_count$pair_motif_scores)
  expect_equal(anchor_pair_example_count$pair_motif_enrich,
               anchor_pair_count$pair_motif_enrich)
})

test_that("can compute significance with hypergeometric test", {
  anchor_pair_match <- anchor_pair_enrich(
    scan_interactions_example_filtered, method = "match")
  data("anchor_pair_example_match", package = "spatzie")
  expect_equal(anchor_pair_example_match$pair_motif_scores,
               anchor_pair_match$pair_motif_scores)
  expect_equal(anchor_pair_example_match$pair_motif_enrich,
               anchor_pair_match$pair_motif_enrich)
})
