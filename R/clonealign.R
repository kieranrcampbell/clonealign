
#' Run clonealign across a range of initializations
#' 
#' Run \code{clonealign} across a range of initializations and select the fit 
#' that acheives the best evidence lower bound (ELBO).
#' 
#' @param gene_expression_data A matrix of gene counts or a
#' \code{SingleCellExperiment}. See \code{?clonealign}
#' @param copy_number_data A matrix or data frame of copy number calls for each clone.
#' See \code{?clonealign}.
#' @param initial_shrinks Initial shrinkages for the clone assignment variational parameters
#' @param n_repeats Number of fits to perform at each initial shrink
#' @param print_elbos Logical - should the ELBOs inferred be printed?
#' @param ... Additional arguments to pass to \code{clonealign(...)}
#' 
#' @details
#' This function essentially wraps \code{clonealign} and can be interacted with as such. The
#' parameter \code{initial_shrinks} controls how hard the variational parameters are initially
#' assigned, analagous to the E-step in EM. At 0, they are initialized evenly across clones,
#' while at 10 they are semi hard assigned to the most likely initial values.
#' 
#' 
#' @export
#' 
#' @return The \code{clonealign_fit} object for the fit that maximizes the ELBO. See \code{?clonealign} for details.
run_clonealign <- function(gene_expression_data,
                          copy_number_data,
                          initial_shrinks = c(0, 5, 10),
                          n_repeats = 3,
                          print_elbos = TRUE,
                          ...) {
  
  args <- list(...)
  args[['gene_expression_data']] <- gene_expression_data
  args[['copy_number_data']] <- copy_number_data
  
  
  fits <- lapply(initial_shrinks, function(is) {
    args[['initial_shrink']] <- is
    replicate(n_repeats, {
      do.call(clonealign, args)
    }, simplify = FALSE)
  })
  
  fits <- unlist(fits, recursive = FALSE)
  
  final_elbos <- sapply(fits, function(ca) ca$convergence_info$final_elbo)
  median_correlations <- sapply(fits, function(ca) median(ca$correlations))
  
  if(print_elbos) {
    message(paste("ELBOs: ", paste(final_elbos, collapse= " ")))
  }
  
  fit_to_return <- fits[[ which.max(final_elbos) ]]
  
  fit_to_return$multirun_info <- list()
  
  fit_to_return$multirun_info$clone_prevalences_at_different_shrinks <- lapply(fits, function(ca) table(ca$clone))
  
  fit_to_return$multirun_info$elbos <- final_elbos
  fit_to_return$multirun_info$median_correlations <- median_correlations
  
  fit_to_return
}


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
#' @param rel_tol Relative tolerance (change in ELBO per iteration in percent) below which the
#' inference is considered converged
#' @param max_iter Maximum number of Variational Bayes iterations to perform
#' @param gene_filter_threshold Genes with total counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param learning_rate The learning rate to be passed to the Adam optimizer
#' @param fix_alpha Should the underlying priors for clone frequencies be fixed? Default TRUE
#' (values are inferred from the data)
#' @param dtype The dtype for tensorflow useage, either "float32" or "float64"
#' @param saturate Should the CNV-expression relationship saturate above copy number = \code{saturation_threshold}? Default TRUE
#' @param saturation_threshold If \code{saturate} is true, copy numbers above this will be reduced to the threshold
#' @param verbose Should warnings and EM convergence information be printed? Default TRUE
#' @param x An optional vector of covariates, e.g. corresponding to batch or patient. Can be a vector of a single
#' covariate or a sample by covariate matrix. Note this should not contain an intercept.
#' @param K The dimensionality of the expression latent space. If left \code{NULL}, K is set to 1 if fewer than 100 genes
#' are present and 6 otherwise.
#' @param initial_shrink The strength with which the variational parameters for clone assignments are
#' initially shrunk towards the most likely assignments. See \code{?run_clonealign}.
#' @param clone_allele A clone-by-variant matrix of copy numbers for each variant
#' @param cov A cell-by-variant matrix of coverage counts
#' @param ref A cell-by-variant matrix of reference allele counts
#' @param mc_samples The number of Monte Carlo samples to use to estimate the ELBO
#' @param clone_call_probability The probability above which a cell is assigned to a clone. If no clone has probability greater
#' than this value, then the clone is "unassigned".
#' @param data_init_mu Should the mu parameters be initialized using the data? (This typically speeds up convergence)
#' @param seed The random seed. See \code{details}.
#'
#'
#'
#' @details
#'
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
#'
#' \strong{Size factors}
#'
#' If \code{size_factors == "fixed"}, the size factors will be set to the overall library size per cell
#' (total number of reads mapping to the cell).
#' If \code{size_factors == "infer"}, the size factors will be treated as a model paramter and jointly
#' optimized during inference.
#' Otherwise, \code{size_factors} can be a numeric vector of precomputed, custom size factors.
#'
#' \strong{Recommended parameter settings}
#'
#' As with any probabilistic model there are many parameters to set. Through comprehensive simulations regarding
#' the robustness of the model to mis-specification (ie what's the minimum proportion of genes for which the
#' CNV-expression relationship can be true and our inferences still valid) we have come up with the following
#' guidelines for parameter settings, reflected in the default values:
#'
#' \itemize{
#' \item{Number of ADAM iterations - if set to 1 we essentially perform gradient descent on the marginal log-likelihood
#' which empircally appears to have the best performance}
#' \item{Dispersions should be clone-specific with weak shrinkage (\code{sigma} = 1 appears best)}
#' \item{The generating probabilities should be fixed to be a priori equal (this corresponds to setting \code{alpha = TRUE})}
#' \item{The cell size factors are best fixed in advanced by multiplying the total counts of whatever genes are passed
#' to clonealign by the edgeR (TMM) normalization factors}
#' }
#'
#'
#' \strong{Controlling Variational inference}
#'
#' Inference is performed using reparametrization-gradient variational inference. Convergence is monitored via changes
#' to the evidence lower bound (ELBO) - this is controlled using the
#' \code{rel_tol} parameter. When the difference between the new and old ELBOs normalized
#' by the absolute value of the old falls below \code{rel_tol}, the  algorithm is considered converged.
#' The maximum number of iterations to acheive this is set using the \code{max_iter} parameter.
#'
#' In each step, maximization is performed using Adam, with learning rate given by \code{learning_rate}.
#'
#' \strong{Random seed}
#'
#' The random seed can be set using the \code{seed} parameter. However, note that this disables GPU computation
#' and parallelism. See \url{https://tensorflow.rstudio.com/tensorflow/reference/use_session_with_seed.html}
#' for details.
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
                       max_iter = 200,
                       rel_tol = 1e-6,
                       gene_filter_threshold = 0,
                       learning_rate = 0.1,
                       x = NULL,
                       clone_allele = NULL,
                       cov = NULL,
                       ref = NULL,
                       fix_alpha = FALSE,
                       dtype = "float32",
                       saturate = TRUE,
                       saturation_threshold = 6,
                       K = NULL,
                       mc_samples = 1,
                       verbose = TRUE,
                       seed = NULL,
                       initial_shrink = 5,
                       clone_call_probability = 0.95,
                       data_init_mu = TRUE) {
  

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

  if(is.null(K)) {
    if(G <= 100) {
      K <- 1
    } else {
      K <- 1
    }
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
  if(is.null(colnames(L))) {
    colnames(L) <- paste0("clone_", letters[seq_len(C)])
  }


  # Sanity checking done - time to call em algorithm
  tflow_res <- inference_tflow(Y,
                               L,
                               max_iter = max_iter,
                               rel_tol = rel_tol,
                               learning_rate = learning_rate,
                               gene_filter_threshold = gene_filter_threshold,
                               x = x,
                               clone_allele = clone_allele,
                               cov = cov,
                               ref = cov,
                               fix_alpha = fix_alpha,
                               dtype = dtype,
                               saturate = saturate,
                               saturation_threshold = saturation_threshold,
                               K = K,
                               mc_samples = mc_samples,
                               verbose = verbose,
                               seed = seed,
                               initial_shrink = initial_shrink,
                               data_init_mu = data_init_mu)
  



  tflow_res$clone = clone_assignment(tflow_res$ml_params$clone_probs, colnames(L), clone_call_probability)

  # Finally map clone names back to fitted values
  colnames(tflow_res$ml_params$clone_probs) <- colnames(L)

  if(!is.null(clone_allele)) {
      colnames(tflow_res$clone_probs_from_snv) <- colnames(L)
  }
  
  tflow_res$correlations <- compute_correlations(Y, L, tflow_res$clone)

  structure(tflow_res, class = "clonealign_fit")

}

#' Post-hoc compute correlations between copy number and expression
#' 
#' @param Y Cell by gene count matrix
#' @param L Gene by clone copy number matrix
#' @param clones Assigned clones
#' 
#' @keywords internal
#' 
#' @importFrom stats cor
#' 
#' @return A vector of correlations between the copy number and assigned gene expression
compute_correlations <- function(Y, L, clones) {
  unassigned <- clones == "unassigned"
  Y <- Y[!unassigned,]
  
  ## Scale Y expression
  Y <- scale(Y)
  
  clones <- clones[!unassigned]
  
  sapply(seq_len(ncol(Y)), function(i) {
    y <- Y[,i]
    x <- L[i, clones]
    suppressWarnings({
      cor(x,y)
    })
  })
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

#' Example raw copy number data
#'
#' An example data frame with copy number calls by
#' region for three clones (A, B, C)
#'
#' @examples
#' data(df_cnv)
"df_cnv"


#' Example clonealign_fit
#'
#' An example \code{clonealign_fit} using the data in \code{example_sce}.
#'
#' @seealso example_sce
#' @examples
#' data(example_clonealign_fit)
"example_clonealign_fit"

#' Saturate a copy number matrix above a certain threshold
#'
#' @keywords internal
#'
#' @return The input clipped (or saturated) below threshold
saturate <- function(x, threshold = 4) {
  x[ x > threshold ] <- threshold
  x
}

.onLoad <- function(libpath, pkgname) {
  if(!reticulate::py_module_available("tensorflow")) {
    msg <- "Tensorflow does not appear to be installed\n"
    msg <- c(msg, "To install run install.pacakges(\"tensorflow\") then tensorflow::install_tensorflow()\n")
    msg <- c(msg, "For more details see the clonealign vignette or https://tensorflow.rstudio.com/tensorflow/articles/installation.html")
    stop(msg)
  }
}


#' Compute mean squared error of a clonealign fit
#'
#' @keywords internal
#'
#' @return The mean square error of the clonealign fit on the given expression data using
#' the provided clones
compute_ca_fit_mse <- function(fit, Y, L,
                               model_mu = FALSE, random_clones = FALSE) {

  clones <- fit$clone
  if(random_clones) {
    distinct_clones <- unique(clones)
    clones <- sample(distinct_clones, nrow(Y), replace = TRUE)
  }
  predicted_expression <- L[,clones] # G by N
  #
  if(model_mu) {
    mu <- as.vector(fit$ml_params$mu)
    predicted_expression <- mu * predicted_expression#[fit$retained_genes,]
  }
  normalizer <- rowSums(Y) / colSums(predicted_expression)
  predicted_expression <- t(predicted_expression) * normalizer

  mse <- mean((predicted_expression - Y)^2)
  mse
}


