context("visualizations")
data("anchor_pair_example_score", package = "spatzie")

test_that("plot_motif_occurrence produces output", {
  p <- plot_motif_occurrence(anchor_pair_example_score)
  expect_equal(class(p), c("gg", "ggplot"))
})
