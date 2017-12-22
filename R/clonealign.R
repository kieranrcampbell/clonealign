
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
                       size_factors = NULL,
                       max_iter = 100, 
                       rel_tol = 1e-5,
                       gene_filter_threshold = 0,
                       bp_param = bpparam()) {
  
  N <- NA # Number of cells
  G <- NA # Number of genes
  C <- NA # Number of clones
  
  # Parse gene expression data first
  if(is(gene_expression_data, "SingleCellExperiment") || is(gene_expression_data, "SummarizedExperiment")) {
    assay_names <- names(assays(gene_expression_data))
    if(!("counts" %in% assay_names)) {
      stop(paste("counts not in assays(gene_expression_data). Available assays:", paste(assay_names, collapse = ",")))
    }
    Y <- t(assay(gene_expression_data, "counts"))
  } else if(is(gene_expression_data, "matrix")) {
    Y <- gene_expression_data
  } else {
    stop("Input gene_expression_data must be SingleCellExperiment, SummarizedExperiment, or matrix")
  }
  N <- nrow(Y)
  G <- ncol(Y)
  
  
  # Parse size factors
  if(is.null(size_factors) && if(is(gene_expression_data, "SingleCellExperiment"))) {
    if(!is.null(sizeFactors(gene_expression_data))) {
      size_factors <- sizeFactors(gene_expression_data)
    }
  }
  
  if(is.null(size_factors)) {
    size_factors <- scran::computeSumFactors(t(Y))
  }
  
  # Parse cnv data
  
  
  
  
}