
#' Assign scRNA-seq to clone-of-origin
#'
#'
#' \code{clonealign} assigns single cells (measured with RNA-seq) to their clones of origin, where
#' the clones have been inferred from ultra-shallow scDNA-seq and collated into copy number profiles.
#'
#' @param gene_expression_data A matrix of gene counts or a
#' \code{SingleCellExperiment}. This should contain raw counts. See \code{details}.
#' @param copy_number_data A matrix or data frame of copy number calls for each clone.
#' See \code{details}.
#' @param size_factors Vector of cell size factors. If NULL computed using
#' \code{scran::computeSumFactors}
#' @param max_iter Maximum number of EM iterations before algorithm is terminated
#' @param rel_tol Relative tolerance (change in log-likelihood per EM iteration in percent) below which the
#' EM algorithm is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param verbose Should warnings and EM convergence information be printed? Default TRUE
#' @param multithread Should the M-step be performed in parallel using \code{BiocParallel}? Default \code{TRUE}
#' @param bp_param Parameters for parallel optimization of the Q function during EM. Default parameters
#' taken by call to \code{bpparam()}. See \code{?bpparam} for more details.
#'
#' @importFrom BiocParallel bpparam
#'
#' @details
#' \strong{Input format}
#'
#' \code{gene_expression_data} must either be a \code{SingleCellExperiment} or \code{SummarizedExperiment}
#' with a \code{counts} assay
#' representing raw gene expression counts, or a cell by gene matrix of raw counts.
#'
#' \code{copy_number_data} must either be a \code{matrix}, \code{data.frame} or \code{DataFrame} with a
#' row for each gene in \code{gene_expression_data} and a column for each of the clones.
#' If \code{colnames(copy_number_data)} is not \code{NULL} then these names will be used for each of
#' the clones in the final output.
#'
#' \code{size_factors} should be a vector of size factors (one for each cell in \code{gene_expression_data}).
#' If \code{NULL} these are computed using \code{scran::computeSumFactors}.
#'
#' \strong{EM convergence monitoring}
#'
#' Inference is performed using the EM algorithm (\url{https://en.wikipedia.org/wiki/Expectation-maximization_algorithm})
#' which uses the log marginal likelihood to monitor convergence. This is controlled using the
#' \code{rel_tol} parameter. When the difference between the new and old log marginal likelihoods normalized
#' by the absolute value of the old falls below \code{rel_tol}, the EM algorithm is considered converged.
#' The maximum number of iterations to acheive this is set using the \code{max_iter} parameter.
#'
#' \strong{Multithreaded optimization}
#'
#' The M-step of the EM algorithm consists of \eqn{G} different optimizations for \eqn{G} genes. These optimizations
#' are performed in parallel using the \code{BiocParallel} package. The parallel
#' computations are controlled using the \code{bp_param} argument. For fine control over this see
#' the \code{\link[BiocParallel]{bpparam}} function. By default \code{clonealign} uses the parameters
#' that are returned by a call to \code{bpparam()}. To use only e.g. 1 core call \code{clonealign} with
#' \code{clonealign(..., bp_param = MulticoreParam(workers = 1))}.
#'
#' @return An object of class \code{clonealign_fit}. The maximum likelihood estimates of the
#' clone assignment paramters are in the \code{clone} slot. Maximum likelihood estimates of
#' all model parameters are in the \code{ml_params} slot.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assays
#'
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(example_sce)
#' copy_number_data <- rowData(example_sce)[,c("A", "B", "C")]
#' cal <- clonealign(example_sce, copy_number_data)
#' print(cal)
#' clones <- cal$clone
clonealign <- function(gene_expression_data,
                       copy_number_data,
                       size_factors = NULL,
                       max_iter = 100,
                       rel_tol = 1e-5,
                       gene_filter_threshold = 0,
                       verbose = TRUE,
                       multithread = TRUE,
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
    Y <- t(as.matrix(assay(gene_expression_data, "counts")))
  } else if(is(gene_expression_data, "matrix")) {
    Y <- gene_expression_data
  } else {
    stop("Input gene_expression_data must be SingleCellExperiment, SummarizedExperiment, or matrix")
  }
  N <- nrow(Y)
  G <- ncol(Y)


  # Parse size factors
  if(is.null(size_factors) && is(gene_expression_data, "SingleCellExperiment")) {
    if(!is.null(BiocGenerics::sizeFactors(gene_expression_data))) {
      size_factors <- BiocGenerics::sizeFactors(gene_expression_data)
    }
  }

  if(is.null(size_factors)) {
    size_factors <- scran::computeSumFactors(t(Y))
  }

  # Parse cnv data
  if(is(copy_number_data, "data.frame") || is(copy_number_data, "DataFrame")) {
    L <- as.matrix(copy_number_data)
  } else if(is(copy_number_data, "matrix")) {
    L <- copy_number_data
  } else {
    stop(paste("copy_number_data must be a matrix, data.frame or DataFrame. Current class:", class(copy_number_data)))
  }

  if(nrow(L) != G) {
    stop("copy_number_data must have same number of genes (rows) as gene_expression_data")
  }

  C <- ncol(L)

  # Sanity checking done - time to call em algorithm
  em <- inference_em(
    Y,
    L,
    size_factors,
    max_iter,
    rel_tol,
    gene_filter_threshold,
    verbose,
    multithread,
    bp_param
  )

  rlist <- list(
    clone = clone_assignment(em)
  )

  ml_params <- list(
    clone_probs = em$gamma,
    mu = em$mu,
    beta = em$beta,
    phi = em$phi
  )

  rlist$ml_params <- ml_params

  # Finally map clone names back to fitted values
  clone_names <- colnames(L)
  if(!is.null(L)) {
    rlist$clone <- plyr::mapvalues(rlist$clone,
                                   from = seq_len(C),
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
#' print(example_clonealign_fit)
print.clonealign_fit <- function(x, ...) {
  N <- length(x$clone)
  G <- length(x$ml_params$mu)
  C <- ncol(x$ml_params$clone_probs)

  msg <- glue("A clonealign_fit for {N} cells, {G} genes, and {C} clones")
  msg <- glue("{msg}\nTo access clone assignments, call x$clone")
  msg <- glue("{msg}\nTo access ML parameter estimates, call x$ml_params")
  cat(msg)
}


#' Example SingleCellExperiment
#'
#' An example \code{SingleCellExperiment} for 100 genes and 200 cells.
#' Copy number data is in \code{rowData(example_sce)[, c("A", "B", "C",]}
#'
#' @seealso example_clonealign_fit
#' @examples
#' data(example_sce)
"example_sce"


#' Example clonealign_fit
#'
#' An example \code{clonealign_fit} using the data in \code{example_sce}.
#'
#' @seealso example_sce
#' @examples
#' data(example_clonealign_fit)
"example_clonealign_fit"


