
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
#' @param B Number of basis functions for spline fitting
#' @param size_factors Either "fixed", "infer", or a numeric vector of size factors. See \code{details}.
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
                       max_iter = 100,
                       rel_tol = 1e-6,
                       gene_filter_threshold = 0,
                       learning_rate = 0.1,
                       x = NULL,
                       clone_allele = NULL,
                       cov = NULL,
                       ref = NULL,
                       fix_alpha = FALSE,
                       size_factors = "fixed",
                       dtype = "float64",
                       saturate = TRUE,
                       saturation_threshold = 6,
                       K = NULL,
                       B = 20,
                       verbose = TRUE,
                       seed = NULL,
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
      K <- 6
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
                               size_factors = size_factors,
                               dtype = dtype,
                               saturate = saturate,
                               saturation_threshold = saturation_threshold,
                               K = K,
                               B = B,
                               verbose = verbose,
                               seed = seed,
                               data_init_mu = data_init_mu)


  rlist <- list(
    clone = clone_assignment(tflow_res)
  )

  ml_params <- list(
    clone_probs = tflow_res$gamma,
    mu = tflow_res$mu,
    s = tflow_res$s,
    alpha = tflow_res$alpha,
    a = tflow_res$a,
    b = tflow_res$b
  )

  if("psi" %in% names(tflow_res)) {
    ml_params$psi <- tflow_res$psi
    ml_params$W <- tflow_res$W
    ml_params$chi <- tflow_res$chi
  }

  if("beta" %in% names(tflow_res)) {
    ml_params$beta <- tflow_res$beta
  }

  rlist$ml_params <- ml_params
  rlist$elbo <- tflow_res$elbo
  rlist$retained_genes <- tflow_res$retained_genes
  rlist$basis_means <- tflow_res$basis_means

  if("clone_probs_from_snv" %in% names(tflow_res)) {
    rlist$clone_probs_from_snv <- tflow_res$clone_probs_from_snv
  }

  # Finally map clone names back to fitted values
  clone_names <- colnames(L)
  if(!is.null(L)) {
    rlist$clone <- plyr::mapvalues(rlist$clone,
                                   from = seq_len(C),
                                   to = clone_names)
    colnames(rlist$ml_params$clone_probs) <- clone_names

    if("clone_probs_from_snv" %in% names(rlist)) {
      colnames(rlist$clone_probs_from_snv) <- clone_names
    }
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


#' Evaluate the quality of a clonealign fit
#'
#' Use mean squared error of predicting the expression of held-out genes to evaluate the quality of a clonealign fit.
#'
#' @param gene_expression_data Either a \code{SingleCellExperiment} or matrix of counts, same as input to \code{clonealign}
#' @param copy_number_data A gene-by-clone matrix of copy numbers, the same as the input to \code{clonealign}
#' @param clonealign_fit The fit returned by a call to \code{clonealign()} on the full geneset
#' @param prop_holdout The proportion of genes to hold out as a \emph{test} set for predicting gene expression
#' @param n_samples The number of random permutations to establish a null distribution of predictive performance
#' @param s Vector of cell size factors - defaults to the total counts per cell
#' @param ... Additional arguments to pass to the \code{clonealign} call
#'
#' @return A list object describing the evaluation results. See \code{details}
#'
#' @details
#'
#' This evaluation function is built around the idea of how good predicted expression is under the model given
#' the inferred (assigned) clones compared to how well you could predict expression given a random clonal assignment.
#'
#' \strong{Evaluations performed}
#'
#' \enumerate{
#' \item On the \emph{full} dataset, the mean square error is compared to the randomized error, and a message
#' printed about the ratio of the two errors and the proportion of time the observed error was less than the
#' error under a null distribution. If either the error under the null is smaller than the observed, or is
#' smaller some percentage of time, then the fit has real problems and should be abandoned as it suggests the
#' model is stuck in a local maxima (which could happen if the proposed clones aren't actually present in the
#' RNA-seq).
#' \item A certain proportion of genes (as decided by the \code{prop_holdout} parameter) are held out as a \emph{test}
#' set, and the clonealign fit is re-performed on the \code{1 - prop_holdout} proportion of genes. The function will then
#' print an agreement table of clone assignments between the full table and the reduced table. If these vastly disagree
#' for only a small change in gene set (ie \code{prop_holdout} is around 0.1 or 0.2), then the fit may be unreliable
#' as it is sensitive to the genes inputted.
#' \item The same metrics from (1) are then printed where the evaluations are performed on the heldout set. Again,
#' if the observed mean squared error given the clonealign fit isn't less than the null mean squared error
#' a large proportion of the time, the fit may be unreliable.
#'
#' }
#'
#' \strong{Object returned}
#'
#' Everything computed above is returned in a list object with the following entries:
#'
#' \itemize{
#' \item full_clonealign_fit - the original clonealign fit on the full gene set that was passed in as the
#' \code{clonealign_fit} argument
#' \item full_observed_mse - the observed mean square error using the full gene set
#' \item full_null_mse - A vector of mean square error under the randomized (permuted) distribution
#' \item reduced_clonealign_fit - a clonealign fit on only the reduced (train) set of genes
#' \item heldout_genes - the names of the genes held out (test set) for evaluating predictive performance
#' out-of-sample
#' \item kept_names - the names of retained genes as part of the reduced (train) set
#' \item heldout_observed_mse - the observed mean square error on the heldout (test) set of genes
#' \item heldout_null_mse - a vector of mean square errors under a null distribution of randomly permuting the clones
#'
#' }
#'
#'
#'
#' @importFrom glue glue
#' @export
#'
#' @examples
#' library(SingleCellExperiment)
#' data(example_clonealign_fit)
#' data(example_sce)
#' copy_number_data <- rowData(example_sce)[,c("A", "B", "C")]
#' evaluate_clonealign(example_sce, copy_number_data, example_clonealign_fit)
#'
evaluate_clonealign <- function(gene_expression_data,
                    copy_number_data,
                    clonealign_fit,
                    prop_holdout = 0.2,
                    n_samples = 2,
                    s = NULL,
                    ...) {

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

  rownames(L) <- colnames(Y)

  C <- ncol(L)

  # Compute mse on full data set
  observed_mse <- compute_ca_fit_mse(clonealign_fit, Y, L)
  null_mse <- replicate(n_samples, compute_ca_fit_mse(clonealign_fit,
                                                Y,
                                                L,
                                                random_clones = TRUE))

  cat(glue("On the full dataset, the observed MSE was on average {mean(null_mse) / observed_mse} times smaller than under a null model"))
  cat(glue(" and smaller {100 * mean(observed_mse < mean(null_mse))}% of the time"))
  cat("\n")

  # Fix which indices we're going to hold out
  genes <- colnames(Y)
  heldout_genes <- sample(genes, round(prop_holdout * G))
  kept_genes <- setdiff(genes, heldout_genes)

  # Fit reduced clonealign model
  message("Fitting reduced clonealign model...")
  ca <- clonealign(Y[, kept_genes], L[kept_genes,], verbose = FALSE, ...)

  tbl <- table(ca$clone, clonealign_fit$clone)


  cat("Agreement between original (rows) and reduced (columns):")
  print(tbl)

  # Compute mse on held out:
  observed_mse_ho <- compute_ca_fit_mse(ca, Y[, heldout_genes], L[heldout_genes,], model_mu = FALSE)
  null_mse_ho <- replicate(n_samples, compute_ca_fit_mse(ca, Y[, heldout_genes], L[heldout_genes,], model_mu = FALSE, random_clones = TRUE))

  cat(glue("On the held-out dataset, the observed MSE was on average {mean(null_mse_ho) / observed_mse_ho} times smaller than under a null model"))
  cat(glue(" and smaller {100 * mean(observed_mse_ho < mean(null_mse_ho))}% of the time"))
  cat("\n")

  list(
    full_clonealign_fit = clonealign_fit,
    full_observed_mse = observed_mse,
    full_null_mse = null_mse,
    reduced_clonealign_fit = ca,
    heldout_genes = heldout_genes,
    kept_genes = kept_genes,
    heldout_observed_mse = observed_mse_ho,
    heldout_null_mse = null_mse_ho
  )

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

#' Plot mean dispersion relationship
#'
#' @param em An object returned by a call to \code{clonealign}
#'
#' @return A \code{ggplot2} plot showing the dispersion as a function of mean
#' @examples
#' data(example_clonealign_fit)
#' plot_mean_dispersion(example_clonealign_fit)
#' @export
plot_mean_dispersion <- function(em) {
  basis_means <- em$basis_means
  a <- em$ml_params$a
  b <- em$ml_params$b
  x <- seq(from = min(basis_means), to = max(basis_means), length.out = 1000)

  y <- sapply(x, function(xx) {
    sum( a * exp(-b*(xx - basis_means)^2) )
  })

  ggplot2::qplot(x, y, geom = 'line') +
    labs(x = expression(mu), y = expression(phi))
}
