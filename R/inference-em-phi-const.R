#'
#' #' @keywords internal
#' likelihood_yn_phi_const <- function(y, L, s_n, pi, params, phi) {
#'   m <- L[, pi] * s_n * params[, 'mu']
#'   ll <- sum(dnbinom2(y, mu = m, size = phi))
#'   ll
#' }
#'
#' #' Computes gamma_{nc} = p(pi_n = c), returning
#' #' N by C matrix
#' #'
#' #' @importFrom matrixStats logSumExp
#' #' @param data Input data
#' #' @param params Model parameters
#' #'
#' #' @keywords internal
#' #'
#' #' @return The probability that each cell belongs to each clone, as a matrix
#' p_pi_phi_const <- function(data, params) {
#'   gamma <- matrix(NA, nrow = data$N, ncol = data$C)
#'   for(n in seq_len(data$N)) {
#'     for(c in seq_len(data$C)) {
#'       gamma[n,c] <- likelihood_yn_phi_const(y = data$Y[n,],
#'                                   L = data$L,
#'                                   s_n = data$s[n],
#'                                   pi = c,
#'                                   params = params,
#'                                   data$phi)
#'     }
#'     gamma[n,] <- exp(gamma[n,] - logSumExp(gamma[n,]))
#'   }
#'   gamma
#' }
#'
#' #' @keywords internal
#' #' @importFrom stats dnbinom
#' dnbinom2 <- function(x, mu, size) {
#'   dnbinom(x, size = size, mu = mu, log = TRUE)
#' }
#'
#'
#' #' Computes Q(theta|theta^(t))
#' #' (function to be optimised under EM) with constant phi
#' #' @keywords internal
#' #'
#' #' @param pars Parameters to optimise
#' #' @param y Gene expression for gene
#' #' @param l Copy number profiles for gene
#' #' @param gamma Expectation of clone assignments at current EM step
#' #' @param data Data used
#' #'
#' #' @return The g'th term in the expected complete data log likelihood
#' Q_g_phi_const <- function(pars, y, l, gamma, data) {
#'   mu <- pars[1]
#'   qq <- 0
#'   for(c in seq_len(data$C)) {
#'     m <- l[c] * data$s * mu # N length vector for given gene of means
#'     l_c <- dnbinom2(y, mu = m, size = data$phi) # p(y_g | pi)
#'     qq <- qq + sum(gamma[,c] * l_c )
#'   }
#'   -qq
#' }
#'
#' #' Computes Q(theta|theta^(t))
#' #' (function to be optimised under EM) with constant phi
#' #' @keywords internal
#' #'
#' #' @param pars Parameters to optimise
#' #' @param y Gene expression for gene
#' #' @param l Copy number profiles for gene
#' #' @param gamma Expectation of clone assignments at current EM step
#' #' @param data Data used
#' #'
#' #' @keywords internal
#' #'
#' #' @return The g'th term in the expected complete data log likelihood
#' Qgr_g_phi_const <- function(pars, y, l, gamma, data) {
#'   mu <- pars[1]
#'   phi <- data$phi
#'   gr <- c('mu' = 0)
#'   for(c in seq_len(data$C)) {
#'     mu_ng <- mu * data$s * l[c] # N-length vector
#'
#'     gr_1 <- (y / mu_ng - (y + phi) / (mu_ng + phi) ) * data$s * l[c]
#'     gr[1] <- gr[1] + sum(gamma[,c] * gr_1)
#'   }
#'   -gr
#' }
#'
#'
#' #' @keywords internal
#' log_likelihood_phi_const <- function(params, data) {
#'   ll <- 0
#'   mu <- params[,'mu']
#'
#'   for(n in seq_len(data$N)) {
#'     pnc <- sapply(seq_len(data$C), function(c) {
#'       sum(dnbinom2(data$Y[n,], mu * data$L[,c] * data$s[n], size = data$phi))
#'     })
#'     ll <- ll + logSumExp(pnc)
#'    }
#'   ll
#' }
#'
#' #' Expectation-maximization for assigning scRNA-seq data
#' #' to clone-of-origin
#' #'
#' #' @param Y cell-by-gene expression matrix of raw counts
#' #' @param L gene-by-clone matrix of copy number variation
#' #' @param s Vector of cell size factors. If NULL computed using
#' #' \code{scran::computeSumFactors}
#' #' @param max_iter Maximum number of EM iterations before algorithm is terminated
#' #' @param rel_tol Relative tolerance in percent below which the log-likelihood is considered converged
#' #' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' #' be filtered out (removes genes with no counts by default)
#' #' @param verbose Logical - should convergence information be printed?
#' #' @param multithread Should the M-step be performed in parallel using \code{BiocParallel}? Default \code{TRUE}
#' #' @param bp_param Parameters for multithreaded optimization of Q function. See \code{?bpparam()}
#' #'
#' #' @importFrom glue glue
#' #' @importFrom BiocParallel bplapply
#' #' @importFrom stats optim
#' #'
#' #' @keywords internal
#' #'
#' #' @return A list with ML estimates for each of the model parameters
#' inference_em_phi_const <- function(Y, L, s = NULL, max_iter = 100, rel_tol = 1e-5,
#'                           gene_filter_threshold = 0, verbose = TRUE,
#'                           multithread = TRUE,
#'                           bp_param = bpparam()) {
#'
#'   zero_gene_means <- colMeans(Y) <= gene_filter_threshold
#'
#'   if(verbose) {
#'     message(glue("Removing {sum(zero_gene_means)} genes with low counts"))
#'   }
#'
#'   Y <- Y[, !zero_gene_means]
#'   L <- L[!zero_gene_means,]
#'
#'
#'   N <- nrow(Y)
#'   G <- ncol(Y)
#'   C <- ncol(L)
#'
#'   # Sanity checks
#'   stopifnot(nrow(L) == G)
#'
#'
#'   if(is.null(s)) {
#'     s <- scran::computeSumFactors(t(Y))
#'   }
#'   stopifnot(length(s) == N)
#'   stopifnot(all(s > 0))
#'
#'
#'   # Always normalise L
#'   L <- t( t(L) / colMeans(L) )
#'
#'   # Initialise
#'   params <- cbind(
#'     colMeans(Y / s) + 0.01
#'   )
#'   colnames(params) <- c("mu")
#'
#'   data <- list(
#'     Y = Y,
#'     L = L,
#'     s = s,
#'     N = N,
#'     G = G,
#'     C = C,
#'     phi = compute_phi(Y)
#'   )
#'
#'   data$L[data$L == 0] <- 1
#'
#'   ll_old <- log_likelihood_phi_const(params, data)
#'
#'   lls <- ll_old
#'
#'   any_optim_errors <- FALSE
#'
#'   for(i in seq_len(max_iter)) {
#'
#'     # E step
#'     gamma <- p_pi_phi_const(data, params)
#'
#'     # M step
#'     if(multithread) {
#'       stop("Why is this multithreading")
#'       pnew <- bplapply(seq_len(data$G), function(g) {
#'         opt <- optim(par = params[g,],
#'                      fn = Q_g_phi_const,
#'                      gr = Qgr_g_phi_const,
#'                      y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
#'                      method = "L-BFGS-B",
#'                      lower = c(1e-10),
#'                      upper = c(max(data$Y)),
#'                      control = list())
#'         if(opt$convergence != 0) {
#'           warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
#'           any_optim_errors <- TRUE
#'         }
#'         c(opt$par, -opt$value)
#'       }, BPPARAM = bp_param)
#'     } else {
#'       pnew <- bplapply(seq_len(data$G), function(g) {
#'         opt <- optim(par = params[g,], # (mu,phi)
#'                      fn = Q_g_phi_const,
#'                      gr = Qgr_g_phi_const,
#'                      y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
#'                      method = "L-BFGS-B",
#'                      lower = c(1e-10),
#'                      upper = c(max(data$Y)),
#'                      control = list())
#'         if(opt$convergence != 0) {
#'           warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
#'           any_optim_errors <- TRUE
#'         }
#'         c('mu' = opt$par, 'value' = -opt$value)
#'       })
#'     }
#'
#'     pnew <- do.call(rbind, pnew)
#'     colnames(pnew) <- c("mu", "value")
#'     params <- pnew[,"mu",drop=FALSE] # just mu
#'     ll <- log_likelihood_phi_const(params, data)
#'
#'     ll_diff <- (ll - ll_old)  / abs(ll_old) * 100
#'
#'     lls <- c(lls, ll)
#'
#'     if(verbose) {
#'       message(glue("{i} Current: {ll_old}\tNew: {ll}\tChange: {ll_diff}"))
#'     }
#'
#'     if(!is.na(ll_diff)) {
#'       if(ll_diff < rel_tol) {
#'         if(verbose) {
#'           print(glue("EM converged after {i} iterations"))
#'         }
#'         break
#'       }
#'     }
#'     ll_old <- ll
#'   }
#'
#'   if(any_optim_errors) {
#'     message("There were errors in optimization of Q function. However, results may still be valid. See errors above.")
#'   }
#'
#'   gamma <- p_pi_phi_const(data, params)
#'   rlist <- list(
#'     gamma = gamma,
#'     mu = params[, 'mu'],
#'     phi = data$phi,
#'     lls = lls
#'   )
#'
#'   if(i == max_iter) {
#'     message("Maximum number of iterations reached; consider increasing max_iter")
#'   }
#'   return(rlist)
#' }
#'
#' #' Compute an emprical estimate of NB dispersion
#' #'
#' #' This function computes an empirical estimate of the negative binomial dispersion.
#' #'
#' #' @param Y A cell-by-gene matrix of counts
#' compute_phi <- function(Y, s = NULL) {
#'   if(is.null(s)) s <- scran::computeSumFactors(t(Y))
#'   means <- colMeans(Y / s)
#'   vars <- colVars(Y / s)
#'   mdf <- data_frame(means, vars)
#'   mdf <- dplyr::mutate(mdf,vars_over_means_minus_1 = vars / means - 1)
#'   fit <- lm(vars_over_means_minus_1 ~ 0 + means, data = mdf)
#'   dispersion <- 1 / coef(fit)
#'   names(dispersion) <- NULL
#'   dispersion
#' }
#'
