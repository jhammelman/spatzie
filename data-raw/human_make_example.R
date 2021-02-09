library(motifmatchr)
library("biomaRt")
library('GenomicInteractions')
library('TFBSTools')
library('ggdendro')
library('ggplot2')
library('reshape2')
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicFeatures)
library(gplots)
library(pheatmap)
library(spatzie)
library(topGO)

setwd('../')
process_ct_data <- function(ct_data,annotation.features,motifs){
  left <- GRanges(seqnames=ct_data$Chromosome,
                  ranges=IRanges(start=as.vector(ct_data$Start),
                                 end=as.vector(ct_data$End)))
  right <- GRanges(seqnames=ct_data$Chromosome.1,
                   ranges=IRanges(start=as.vector(ct_data$Start.1),
                                  end=as.vector(ct_data$End.1)))
  interactions <- GenomicInteractions(left,right)
  annotateInteractions(interactions, annotation.features)
  length(interactions[isInteractionType(interactions, "promoter", "distal")])
  PD = interactions[isInteractionType(interactions, "distal", "promoter")]
  plotInteractionAnnotations(interactions)
  promoter_left <- (elementMetadata(anchorOne(PD))[,"node.class"] == "promoter")
  promoter_right <- (elementMetadata(anchorTwo(PD))[,"node.class"] == "promoter")
  promoter_ranges <- c(anchorOne(PD)[promoter_left],anchorTwo(PD)[promoter_right])
  enhancer_ranges <- c(anchorTwo(PD)[promoter_left],anchorOne(PD)[promoter_right])
  PD <- GenomicInteractions(promoter_ranges,enhancer_ranges)
  PDInteraction <- scan_motifs(PD,motifs,BSgenome.Hsapiens.UCSC.hg19)
  PDInteraction <- filter_motifs(PDInteraction,0.4)
  print(length(PDInteraction$anchorOneMotifIndices))
  print(length(PDInteraction$anchorTwoMotifIndices))

  PDscoreCorr <- anchor_pair_enrich(PDInteraction,method='scoreCorrelation')
  PDscoreCorr <- filter_pair_motifs(PDscoreCorr)
  return(PDscoreCorr)
}

hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19",tablename="ensGene")
refseq.genes = GenomicFeatures::genes(hg19.refseq.db)
refseq.transcripts = transcriptsBy(hg19.refseq.db, by="gene")
refseq.transcripts = refseq.transcripts[ names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) ]
hg19_refseq_promoters <- promoters(refseq.transcripts, 2500,2500)
hg19_refseq_promoters <- unlist(hg19_refseq_promoters)
hg19_refseq_promoters <- unique(hg19_refseq_promoters) # some duplicate promoters from different transcript isoforms

mart = useMart("ensembl", dataset="hgfemale_gene_ensembl")
tx_names_simple <- sapply(strsplit(hg19_refseq_promoters$tx_name, ".", fixed = TRUE),"[",1)
genes <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id"), filter = "ensembl_transcript_id",
               values = tx_names_simple, mart = mart)
hg19_refseq_promoters$geneSymbol <- genes$hgnc_symbol[match(tx_names_simple, genes$ensembl_transcript_id)]
names(hg19_refseq_promoters) <- hg19_refseq_promoters$geneSymbol
na.symbol <- is.na(names(hg19_refseq_promoters))
names(hg19_refseq_promoters)[na.symbol] <- hg19_refseq_promoters$tx_name[na.symbol]
annotation.features <- list(promoter = hg19_refseq_promoters)
motifs <- readJASPARMatrix("~/giffordlab/EPEnrich/HOCOMOCOv11_core_MOUSE_mono.piq", matrixClass="PFM")
all_data <- read.table('data/41586_2020_2151_MOESM5_ESM.bedpe.txt',sep="\t",header=TRUE)
i <- 1
ct_interaction_list <- list()
for (ct in colnames(all_data)[22:32]){
  ct_data <- all_data[all_data[ct]=="Yes",]
  ct_interaction_list[[i]] <- process_ct_data(ct_data,annotation.features,motifs)
  i <- i+1
}
colnames(all_data)[22:32]
lapply(ct_interaction_list,function(x){length(x$interactions)})
interactionDataMSLCL = ct_interaction_list[[9]]
interactionDataK562 = ct_interaction_list[[2]]

plot_motif_pairs_heatmap(compare_motif_pairs(mslcl_inter,k562_inter))
save(interactionDataMSLCL,file='spatzie/data/interactionDataMSLCL.rda',compress = "xz")
save(interactionDataK562,file='spatzie/data/interactionDataK562.rda',compress="xz")

