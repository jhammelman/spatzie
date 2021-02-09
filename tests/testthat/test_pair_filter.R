context('filter_pair_motifs')
data("anchor_pair_example_countcorr", package = "spatzie")
test_that('motif pairs filtered by multiple hypothesis corrected significance',{
  filter_pairs <- filter_pair_motifs(anchor_pair_example_countcorr,threshold=0.5)
  data("filter_pairs_example", package = "spatzie")
  expect_equal(filter_pairs_example$pair_motif_enrich_sig,
               filter_pairs$pair_motif_enrich_sig)
})
