
#' Assign scRNA-seq to clone-of-origin
#'
#' 
#' \code{clonealign} assigns single cells (measured with RNA-seq) to their clones of origin, where
#' the clones have been inferred from ultra-shallow scDNA-seq and collated into copy number profiles.
#' 
#' @param gene_expression_data A matrix of gene counts or a
#' \code{SingleCellExperiment}. This should contain raw counts. See \link{details}.
#' @param copy_number_data A matrix of copy number calls for each clone or a \code{SummarizedExperiment}.
#' See \link{details}.
#' @param s Vector of cell size factors. If NULL computed using
#' \code{scran::computeSumFactors}
#' @param max_iter Maximum number of EM iterations before algorithm is terminated
#' @param rel_tol Relative tolerance (change in log-likelihood per EM iteration in percent) below which the 
#' EM algorithm is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param verbose Should warnings and EM convergence information be printed? Default TRUE
#' @param bp_param Parameters for parallel optimization of the Q function during EM. Default parameters
#' taken by call to \code{bpparam()}. See \code{?bpparam} for more details.
#' 
#' @importFrom BiocParallel bpparam
#' 
#' @details
#' \strong{Input format}
#' ABC
#' 
#' \strong{EM convergence monitoring}
#' 
#' \strong{Multithreaded optimization}
#' 
#' 
#' @examples 
#' data(example_sce)
#' cnv_data <- rowData(example_sce)[,c("A", "B", "C")]
#' cal <- clonealign(example_sce, cnv_data)
#' print(cal)
#' clones <- cal$clone
clonealign <- function(gene_expression_data,
                       copy_number_data,
                       size_factors = NULL,
                       max_iter = 100, 
                       rel_tol = 1e-5,
                       gene_filter_threshold = 0,
                       verbose = TRUE,
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
  if(is.null(size_factors) && is(gene_expression_data, "SingleCellExperiment")) {
    if(!is.null(sizeFactors(gene_expression_data))) {
      size_factors <- sizeFactors(gene_expression_data)
    }
  }
  
  if(is.null(size_factors)) {
    size_factors <- scran::computeSumFactors(t(Y))
  }
  
  # Parse cnv data
  if(is(cnv_data, "data.frame") || is(cnv_data, "DataFrame")) {
    L <- as.matrix(cnv_data)
  } else if(is(cnv_data, "matrix")) {
    L <- cnv_data
  } else {
    stop(paste("cnv_data must be a matrix, data.frame or DataFrame. Current class:", class(cnv_data)))
  }
  
  if(nrow(L) != G) {
    stop("cnv_data must have same number of genes (rows) as gene_expression_data")
  }
  
  # Sanity checking done - time to call em algorithm
  em <- inference_em(
    Y,
    L,
    size_factors,
    max_iter,
    rel_tol,
    gene_filter_threshold,
    verbose,
    bp_param
  )
  
  rlist <- list(
    clone = clone_assignment(em)
  )
  
  ml_params <- list(
    clone_probs = em$gamma,
    mu = em$mu,
    phi = em$phi
  )
  
  rlist$ml_params <- ml_params
  
  # Finally map clone names back to fitted values
  clone_names <- colnames(L)
  if(!is.null(L)) {
    rlist$clone <- plyr::mapvalues(rlist$clone, 
                                   from = sort(unique(rlist$clone)),
                                   to = clone_names)
    colnames(rlist$ml_params$clone_probs) <- clone_names
  }
  
  rlist
  
  structure(rlist, class = "clonealign_fit")
  
}

#' Print a clonealign_fit
#' 
#' @param x An object of class \code{clonealign_fit}
#' @param ... Additional arguments
#' 
#' @method print clonealign_fit
#' @export
#' @return A character string representation of \code{x}
#' @importFrom glue glue
#' @examples 
#' data(example_clonealign_fit)
#' print(clonealign_fit)
print.clonealign_fit <- function(x, ...) {
  N <- length(x$clone)
  G <- length(x$ml_params$mu)
  C <- ncol(x$ml_params$clone_probs)
  
  msg <- glue("A clonealign_fit for {N} cells, {G} genes, and {C} clones")
  msg <- glue("{msg}\nTo access clone assignments, call x$clone")
  msg <- glue("{msg}\nTo access ML parameter estimates, call x$ml_params")
  cat(msg)
}



