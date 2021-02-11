context("helpers")
data("compare_pairs_example", package = "spatzie")
data("yy1_interactions", package = "spatzie")
data("yy1P_yy1E_interactions", package = "spatzie")
data("yy1P_interactions", package = "spatzie")
data("yy1E_interactions", package = "spatzie")
data("yy1_pd_interaction", package = "spatzie")

test_that("compare_motifs_heatmap produces output",{
  pairs_compared <- compare_motif_pairs(anchor_pair_example_scorecorr,anchor_pair_example_countcorr)
  expect_equal(pairs_compared,compare_pairs_example)
})

test_that("get_specific_interactions.R",{
  yy1P_yy1E_test_interactions <- get_specific_interactions(yy1_pd_interaction,
                                                      anchor1_motif="TYY1_MOUSE.H11MO.0.A",
                                                      anchor2_motif="TYY1_MOUSE.H11MO.0.A")
  expect_equal(yy1P_yy1E_interactions,
               yy1P_yy1E_test_interactions)
  expect_is(yy1P_yy1E_test_interactions,"GenomicInteractions")

  yy1P_test_interactions <- get_specific_interactions(yy1_pd_interaction,
                                                      anchor1_motif="TYY1_MOUSE.H11MO.0.A")
  expect_equal(yy1P_interactions,
               yy1P_test_interactions)

  expect_is(yy1P_test_interactions,"GenomicInteractions")
  yy1E_test_interactions <- get_specific_interactions(yy1_pd_interaction,
                                                      anchor2_motif="TYY1_MOUSE.H11MO.0.A")
  expect_equal(yy1E_interactions,
               yy1E_test_interactions)
  expect_is(yy1E_test_interactions,"GenomicInteractions")

  yy1_pd_test_interactions <- get_specific_interactions(yy1_pd_interaction)
  expect_equal(yy1_pd_test_interactions,
               yy1_pd_interaction$interactions)
  expect_is(yy1_pd_test_interactions,"GenomicInteractions")

})
