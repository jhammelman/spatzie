context("visualizations")
data("anchor_pair_example_scorecorr", package = "spatzie")
data("anchor_pair_example_countcorr", package = "spatzie")

test_that("plot_motif_pairs_heatmap produces output", {
    p <- plot_motif_pairs_heatmap(
      -log2(anchor_pair_example_scorecorr$pair_motif_enrich + 1e-10))
    expect_equal(class(p), "pheatmap")
})

test_that("plot_motif_occurrence produces output", {
  p <- plot_motif_occurrence(anchor_pair_example_scorecorr)
  expect_equal(class(p), c("gg", "ggplot"))
})
