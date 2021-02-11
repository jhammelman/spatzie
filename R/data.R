#' Mouse YY1 Enhancer - Promoter Interactions Data Set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' and serves as example and unit test data. The same data set is used in
#' the vignette.
#'
#' @format A \code{\link[GenomicInteractions]{GenomicInteractions}} object
#' @usage data(yy1_interactions)
"yy1_interactions"

#' K562 Enhancer - Promoter Interactions Data Set
#'
#' This object contains genomic interactions obtained by human RAD21 ChIA-PET
#' from K562 cells and serves as unit test data.
#'
#' @format A \code{\link[GenomicInteractions]{GenomicInteractions}} object
#' @usage data(interactionDataK562)
"interactionDataK562"

#' MSLCL Enhancer - Promoter Interactions Data Set
#'
#' This object contains genomic interactions obtained by human RAD21 ChIA-PET
#' from MSLCL cells and serves as unit test data.
#'
#' @format A \code{\link[GenomicInteractions]{GenomicInteractions}} object
#' @usage data(interactionDataMSLCL)
"interactionDataMSLCL"

#' Interactions scanned for motifs - interactionData object
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs and serves as unit test data.
#' @format An interactionData object
#'
#' @usage data(scan_interactions_example)
"scan_interactions_example"

#' Interactions with motifs filtered for significance - interactionData object
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs and filtered for motifs present
#' in at least 10% of interactions. It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(scan_interactions_example_filtered)
"scan_interactions_example_filtered"


#' spatzie count correlation data set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs, filtered for motifs present
#' in at least 10% of interactions, and processed for significant motif:motif
#' interactions with count correlation. It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(anchor_pair_example_countcorr)
"anchor_pair_example_countcorr"

#' spatzie count fisher data set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs, filtered for motifs present
#' in at least 10% of interactions, and processed for significant motif:motif
#' interactions with using Fisher's test.  It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(anchor_pair_example_countfisher)
"anchor_pair_example_countfisher"

#' spatzie count hyper data set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs, filtered for motifs present
#' in at least 10% of interactions, and processed for significant motif:motif
#' interactions with using the Hypergeometric test.  It serves as unit test data.
#'
#' @format A interactionData object
#' @usage data(anchor_pair_example_counthyper)
"anchor_pair_example_counthyper"

#' spatzie score correlation data set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs, filtered for motifs present
#' in at least 10% of interactions, and processed for significant motif:motif
#' interactions with score correlation. It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(anchor_pair_example_scorecorr)
"anchor_pair_example_scorecorr"

#' spatzie score correlation filtered data set
#'
#' This object contains genomic interactions obtained by mouse YY1 ChIA-PET
#' scanned for mouse transcription factor motifs, filtered for motifs present
#' in at least 10% of interactions, processed for significant motif:motif
#' interactions with score correlation, and filtered for pairs with p < 0.5.
#' It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(filter_pairs_example)
"filter_pairs_example"

#' compare_motif_pairs example
#'
#' This is a matrix containing example result from compare_motif_pairs. It
#' serves as unit test data.
#'
#' @format A matrix
#' @usage data(compare_pairs_example)
#'
"compare_pairs_example"

#' yy1_pd_interaction
#'
#' This is a interactionData object containing processed results from YY1 ChIA-PET
#' interaction data. It serves as unit test data.
#'
#' @format An interactionData object
#' @usage data(yy1_pd_interaction)
"yy1_pd_interaction"

#' yy1E_interactions
#'
#' This is a GenomicInteractions object containing proccessed results from YY1 ChIA-PET
#' of interactions that contain a YY1 motif in the enhancer (anchor 2) region. It serves
#' as unit test data.
#'
#' @format A GenomicInteractions object
#' @usage data(yy1E_interactions)
"yy1E_interactions"

#' yy1P_interactions
#'
#' This is a GenomicInteractions object containing proccessed results from YY1 ChIA-PET
#' of interactions that contain a YY1 motif in the promoter (anchor 1) region. It serves
#' as unit test data.
#'
#' @format A GenomicInteractions object
#' @usage data(yy1P_interactions)
"yy1P_interactions"

#' yy1P_yy1E_interactions
#'
#' This is a GenomicInteractions object containing proccessed results from YY1 ChIA-PET
#' of interactions that contain a YY1 motif in the promoter (anchor 1) region and a YY1
#' motif in the enhancer (anchor 2) region. It serves as unit test data.
#'
#' @format A GenomicInteractions object
#' @usage data(yy1P_yy1E_interactions)
"yy1P_yy1E_interactions"
