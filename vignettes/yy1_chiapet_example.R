library(motifmatchr)
library("biomaRt")
library('GenomicInteractions')
library('TFBSTools')
library('ggdendro')
library('ggplot2')
library('reshape2')
library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicFeatures)
library(gplots)
library(pheatmap)
library(spatzie)
library(topGO)

#------------------------------------------------------------------------#
#                             load data                                  #
#------------------------------------------------------------------------#
yy1_interactions_file <- system.file("extdata/yy1_interactions.bedpe",
                                     package = "spatzie")
yy1_interactions <- makeGenomicInteractionsFromFile(yy1_interactions_file,
                                type="bedpe",
                                experiment_name="yy1",
                                description="mESC yy1 chr1")
length(yy1_interactions)
#------------------------------------------------------------------------#
#                   load mouse promoter information                      #
#------------------------------------------------------------------------#
#todo genes -> gene names not working anymore
mm9.refseq.db <- makeTxDbFromUCSC(genome="mm9")
refseq.genes = GenomicFeatures::genes(mm9.refseq.db)
refseq.transcripts = transcriptsBy(mm9.refseq.db, by="gene")
refseq.transcripts = refseq.transcripts[names(refseq.transcripts)
                                         %in% unlist(refseq.genes$gene_id)]
mm9_refseq_promoters <- promoters(refseq.transcripts, 2500,2500)
# some duplicate promoters from different transcript isoforms
mm9_refseq_promoters <- unique(mm9_refseq_promoters)
mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")
tx_names_simple <- sapply(strsplit(mm9_refseq_promoters$tx_name,
                                   ".", fixed = TRUE),"[",1)
genes <- getBM(attributes = c("hgnc_symbol", "ucsc"), filter = "ucsc",
               values = tx_names_simple, mart = mart)
mm9_refseq_promoters$geneSymbol <- genes$hgnc_symbol[match(tx_names_simple,
                                                           genes$ucsc)]
names(mm9_refseq_promoters) <- mm9_refseq_promoters$geneSymbol
na.symbol <- is.na(names(mm9_refseq_promoters))
names(mm9_refseq_promoters)[na.symbol] <- mm9_refseq_promoters$tx_name[na.symbol]

#------------------------------------------------------------------------#
#              annotate interactions with promoters                      #
#------------------------------------------------------------------------#
annotation.features <- list(promoter = mm9_refseq_promoters)
annotateInteractions(yy1_interactions, annotation.features)
plotInteractionAnnotations(yy1_interactions)

#------------------------------------------------------------------------#
#         limit and sort interactions to enhancer:promoter               #
#------------------------------------------------------------------------#
yy1PD = yy1_interactions[isInteractionType(yy1_interactions,
                                           "distal", "promoter")]
promoter_left <- (elementMetadata(anchorOne(yy1PD))[,"node.class"] == "promoter")
promoter_right <- (elementMetadata(anchorTwo(yy1PD))[,"node.class"] == "promoter")
promoter_ranges <- c(anchorOne(yy1PD)[promoter_left],anchorTwo(yy1PD)[promoter_right])
enhancer_ranges <- c(anchorTwo(yy1PD)[promoter_left],anchorOne(yy1PD)[promoter_right])
yy1PD <- GenomicInteractions(promoter_ranges,enhancer_ranges)

#------------------------------------------------------------------------#
#                              scan motifs                               #
#------------------------------------------------------------------------#
motiffile <- system.file("extdata/consensus_HOCOMOCOv11_core_MOUSE-plus_YY1.piq",
                                     package = "spatzie")
motifs <- readJASPARMatrix(motiffile,
                           matrixClass="PFM")
yy1PDInteraction <- scanMotifs(yy1PD,motifs,BSgenome.Mmusculus.UCSC.mm9)
yy1PDInteraction <- filterMotifs(yy1PDInteraction,0.4)

#------------------------------------------------------------------------#
#                      compute motif co-significance                    #
#------------------------------------------------------------------------#
yy1PDcountCorr <- anchorPairEnrich(yy1PDInteraction,method='countCorrelation')
yy1PDscoreCorr <- anchorPairEnrich(yy1PDInteraction,method='scoreCorrelation')
yy1PDhyperCorr <- anchorPairEnrich(yy1PDInteraction,method='countHypergeom')

plotMotifPairsHeatmap(-log2(yy1PDcountCorr$pairMotifEnrich))
plotMotifPairsHeatmap(-log2(yy1PDscoreCorr$pairMotifEnrich))
plotMotifPairsHeatmap(-log2(yy1PDhyperCorr$pairMotifEnrich))

