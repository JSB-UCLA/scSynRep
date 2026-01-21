#' scSynRep: Generate Synthetic Single-Cell RNA-seq Replicates
#'
#' @description
#' The scSynRep package provides a pipeline for generating synthetic
#' single-cell RNA-seq replicates by integrating bulk RNA-seq distribution
#' fitting with scDesign3-based single-cell modeling.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{scSynRep_prep}}: Preprocess and filter bulk and single-cell data
#'   \item \code{\link{fit_bulk}}: Fit multivariate normal distribution to bulk RNA-seq
#'   \item \code{\link{scDesign3_fit}}: Fit scDesign3 model to single-cell data
#'   \item \code{\link{scSynRep_gen_bulk}}: Generate synthetic bulk samples
#'   \item \code{\link{scSynRep_gen_sc}}: Generate synthetic single-cell replicates
#'   \item \code{\link{scSynRep_from_tissue}}: Complete pipeline wrapper
#' }
#'
#' @section Workflow:
#' The typical workflow is:
#' \enumerate{
#'   \item Use \code{scSynRep_prep()} to filter genes and compute PCA
#'   \item Use \code{fit_bulk()} to fit bulk distribution parameters
#'   \item Use \code{scDesign3_fit()} to model single-cell expression
#'   \item Use \code{scSynRep_gen_bulk()} to sample synthetic bulk profiles
#'   \item Use \code{scSynRep_gen_sc()} to generate synthetic single-cell replicates
#' }
#'
#' Alternatively, use \code{scSynRep_from_tissue()} to run the complete pipeline.
#'
#' @docType package
#' @name scSynRep-package
#' @aliases scSynRep
#' @keywords internal
"_PACKAGE"
