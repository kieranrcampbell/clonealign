
#' Assign scRNA-seq to clone-of-origin
#'
#' 
#' @param gene_expression_data A matrix of gene counts or a
#' \code{SingleCellExperiment}. See \link{details}.
#' @param copy_number_data A matrix of copy number calls for each clone or a \code{SummarizedExperiment}.
#' See \link{details}.
#' @param s Vector of cell size factors. If NULL computed using
#' \code{scran::computeSumFactors}
#' @param max_iter Maximum number of EM iterations before algorithm is terminated
#' @param rel_tol Relative tolerance (change in log-likelihood per EM iteration in percent) below which the 
#' EM algorithm is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param bp_param Parameters for parallel optimization of the Q function during EM. Default parameters
#' taken by call to \code{bpparam()}. See \code{?bpparam} for more details.
#' 
#' @details
#' \strong{Input format}
#' ABC
#' 
clonealign <- function(gene_expression_data,
                       copy_number_data,
                       s = NULL,
                       max_iter = 100, 
                       rel_tol = 1e-5,
                       gene_filter_threshold = 0,
                       bp_param = bpparam()) {
  
}