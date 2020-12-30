#' spatzie
#'
#' Looks for motifs which are significantly co-enriched from enhancer-promoter
#' interaction data, derived from assays such as as HiC, ChIA-PET, etc. It can
#' also look for differentially enriched motif pairs between to interaction
#' experiments.
#'
#' @docType package
#' @author Jennifer Hammelman
#' @name spatzie
NULL

.onLoad <- function(libname = find.package("spatzie"), pkgname = "spatzie") {
  envir <- parent.env(environment())
  utils::data("interactions", package = pkgname, envir = envir)
}
