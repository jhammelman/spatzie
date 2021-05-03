context("helpers")
data("compare_pairs_example", package = "spatzie")
data("anchor_pair_example_score", package = "spatzie")
data("anchor_pair_example_count", package = "spatzie")
data("int_data_yy1", package = "spatzie")
data("interactions_yy1_promoter", package = "spatzie")
data("interactions_yy1_enhancer", package = "spatzie")
data("interactions_yy1_ep", package = "spatzie")

test_that("compare_motifs_heatmap produces output", {
  pairs_compared <- compare_motif_pairs(anchor_pair_example_score,
                                        anchor_pair_example_count)
  expect_equal(pairs_compared, compare_pairs_example)
})

test_that("get_specific_interactions.R", {
  yy1p_yy1e_test_interactions <- get_specific_interactions(
    int_data_yy1,
    anchor1_motif = "YY1",
    anchor2_motif = "YY1")
  expect_equal(interactions_yy1_ep, yy1p_yy1e_test_interactions)
  expect_is(yy1p_yy1e_test_interactions, "GenomicInteractions")

  yy1p_test_interactions <- get_specific_interactions(
    int_data_yy1, anchor1_motif = "YY1")
  expect_equal(interactions_yy1_promoter, yy1p_test_interactions)
  expect_is(yy1p_test_interactions, "GenomicInteractions")

  yy1e_test_interactions <- get_specific_interactions(
    int_data_yy1, anchor2_motif = "YY1")
  expect_equal(interactions_yy1_enhancer, yy1e_test_interactions)
  expect_is(yy1e_test_interactions, "GenomicInteractions")

  yy1_pd_test_interactions <- get_specific_interactions(int_data_yy1)
  expect_equal(yy1_pd_test_interactions, int_data_yy1$interactions)
  expect_is(yy1_pd_test_interactions, "GenomicInteractions")
})
