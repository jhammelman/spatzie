library(GenomicFeatures)
library(GenomicInteractions)
library(spatzie)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(motifmatchr)
library(testthat)
motif_file <- system.file(  "extdata/consensus_HOCOMOCOv11_core_MOUSE-plus_YY1.piq", package = "spatzie")
motifs <- TFBSTools::readJASPARMatrix(motif_file, matrixClass = "PFM")
left <- GRanges(seqnames=c("chr1","chr1","chr1"),
                ranges=IRanges(start=c(100000,100550,101050),
                               end=c(100550,101050,101600)))
right <- GRanges(seqnames=c("chr1","chr2","chr2"),
                 ranges=IRanges(start=c(200000,200550,201050),
                                end=c(200550,201050,201600)))
interactions <- GenomicInteractions(left,right)
scan_interactions_example <- scan_motifs(interactions,motifs,BSgenome.Hsapiens.UCSC.hg19)
scan_interactions_example_filtered <- filter_motifs(scan_interactions_example,threshold=0.1)
anchor_pair_example_scorecorr <- anchor_pair_enrich(scan_interactions_example_filtered,method="scoreCorrelation")
anchor_pair_example_countcorr <- anchor_pair_enrich(scan_interactions_example_filtered,method="countCorrelation")
anchor_pair_example_counthyper <- anchor_pair_enrich(scan_interactions_example_filtered,method="countHypergeom")
anchor_pair_example_countfisher <- anchor_pair_enrich(scan_interactions_example_filtered,method="countFisher")

save(scan_interactions_example,file="data/scan_interactions_example.rda",compress='xz')
save(scan_interactions_example_filtered,file="data/scan_interactions_example_filtered.rda",compress='xz')
save(anchor_pair_example_scorecorr,file="data/anchor_pair_example_scorecorr.rda",compress='xz')
save(anchor_pair_example_countcorr,file="data/anchor_pair_example_countcorr.rda",compress='xz')
save(anchor_pair_example_counthyper,file="data/anchor_pair_example_counthyper.rda",compress='xz')
save(anchor_pair_example_countfisher,file="data/anchor_pair_example_countfisher.rda",compress='xz')
filter_pairs_example <- filter_pair_motifs(anchor_pair_example_countcorr,threshold=0.5)
save(filter_pairs_example,file="data/filter_pairs_example.rda",compress='xz')

