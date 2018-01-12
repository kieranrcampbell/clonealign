
#' @keywords internal
likelihood_yn <- function(y, L, s_n, pi, params) {
  m <- params[, 'mu'] * (1 + params[,'beta'] * L[, pi]) * s_n
  phi <- params[, 'phi']
  # m[m == 0] <- 0.1
  ll <- sum(dnbinom2(y, mu = m, size = phi))
  ll
}

#' Computes gamma_{nc} = p(pi_n = c), returning
#' N by C matrix
#'
#' @importFrom matrixStats logSumExp
#' @param data Input data
#' @param params Model parameters
#'
#' @keywords internal
#'
#' @return The probability that each cell belongs to each clone, as a matrix
p_pi <- function(data, params) {
  gamma <- matrix(NA, nrow = data$N, ncol = data$C)
  for(n in seq_len(data$N)) {
    for(c in seq_len(data$C)) {
      gamma[n,c] <- likelihood_yn(y = data$Y[n,],
                                  L = data$L,
                                  s_n = data$s[n],
                                  pi = c,
                                  params = params)
    }
    gamma[n,] <- exp(gamma[n,] - logSumExp(gamma[n,]))
  }
  gamma
}

#' @keywords internal
#' @importFrom stats dnbinom
dnbinom2 <- function(x, mu, size) {
  dnbinom(x, size = size, mu = mu, log = TRUE)
}


#' Computes Q(theta|theta^(t))
#' (function to be optimised under EM)
#' @keywords internal
#'
#' @param pars Parameters to optimise
#' @param y Gene expression for gene
#' @param l Copy number profiles for gene
#' @param gamma Expectation of clone assignments at current EM step
#' @param data Data used
#'
#' @return The g'th term in the expected complete data log likelihood
Q_g <- function(pars, y, l, gamma, data) {
  mu <- pars[1]
  beta <- pars[2]
  phi <- pars[3]
  qq <- 0
  for(c in seq_len(data$C)) {
    m <- mu * (1 + beta * l[c]) * data$s # N length vector for given gene of means
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(gamma[,c] * l_c )
  }
  -qq
}

#' Computes map clone assignment given EM object
#'
#' @param em List returned by \code{inference_em}
#' @return A vector of maximum likelihood clone assignments
#' @keywords internal
clone_assignment <- function(em) {
  apply(em$gamma, 1, which.max)
}



#' @keywords internal
log_likelihood <- function(params, data) {
  ll <- 0
  mu <- params[,'mu']
  beta <- params[,'beta']
  phi <- params[,'phi']

  for(n in seq_len(data$N)) {
    pnc <- sapply(seq_len(data$C), function(c) {
      m <- mu * (1 + beta * data$L[,c]) * data$s[n]
      sum(dnbinom2(data$Y[n,], m, size = phi))
    })
    ll <- ll + logSumExp(pnc)
   }
  ll
}

#' Expectation-maximization for assigning scRNA-seq data
#' to clone-of-origin
#'
#' @param Y cell-by-gene expression matrix of raw counts
#' @param L gene-by-clone matrix of copy number variation
#' @param s Vector of cell size factors. If NULL computed using
#' \code{scran::computeSumFactors}
#' @param max_iter Maximum number of EM iterations before algorithm is terminated
#' @param rel_tol Relative tolerance in percent below which the log-likelihood is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param verbose Logical - should convergence information be printed?
#' @param multithread Should the M-step be performed in parallel using \code{BiocParallel}? Default \code{TRUE}
#' @param bp_param Parameters for multithreaded optimization of Q function. See \code{?bpparam()}
#'
#' @importFrom glue glue
#' @importFrom BiocParallel bplapply
#' @importFrom stats optim
#'
#' @keywords internal
#'
#' @return A list with ML estimates for each of the model parameters
inference_em <- function(Y, L, s = NULL, max_iter = 100, rel_tol = 1e-5,
                          gene_filter_threshold = 0, verbose = TRUE,
                          multithread = TRUE,
                          bp_param = bpparam()) {

  zero_gene_means <- colMeans(Y) <= gene_filter_threshold

  if(verbose) {
    message(glue("Removing {sum(zero_gene_means)} genes with low counts"))
  }

  Y <- Y[, !zero_gene_means]
  L <- L[!zero_gene_means,]


  N <- nrow(Y)
  G <- ncol(Y)
  C <- ncol(L)

  # Sanity checks
  stopifnot(nrow(L) == G)


  if(is.null(s)) {
    s <- scran::computeSumFactors(t(Y))
  }
  stopifnot(length(s) == N)
  stopifnot(all(s > 0))


  # Always normalise L
  # L <- t( t(L) / colMeans(L) )

  # Initialise
  params <- cbind(
    colMeans(Y / s) + 0.01,
    rep(0.5, G),
    rep(1, G)
  )
  colnames(params) <- c("mu", "beta", "phi")

  data <- list(
    Y = Y,
    L = L,
    s = s,
    N = N,
    G = G,
    C = C
  )

  data$L[data$L == 0] <- 1

  ll_old <- log_likelihood(params, data)

  any_optim_errors <- FALSE

  for(i in seq_len(max_iter)) {

    # E step
    gamma <- p_pi(data, params)

    # M step
    if(multithread) {
      pnew <- bplapply(seq_len(data$G), function(g) {
        maxL <- max(data$L[g,])
        beta_lower_bound <- -1/maxL
        opt <- optim(par = params[g,], # (mu,phi)
                     fn = Q_g,
                     y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                     method = "L-BFGS-B",
                     lower = c(1e-10, 0, 1e-10),
                     upper = c(max(data$Y), 1e6, 1e6),
                     control = list())
        if(opt$convergence != 0) {
          warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
          any_optim_errors <- TRUE
        }
        c(opt$par, -opt$value)
      }, BPPARAM = bp_param)
    } else {
      pnew <- lapply(seq_len(data$G), function(g) {
        opt <- optim(par = params[g,], # (mu,phi)
                     fn = Q_g,
                     # gr = grad_g,
                     y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                     method = "L-BFGS-B",
                     lower = c(1e-10, 1e-10, 1e-10),
                     upper = c(max(data$Y), 1e6, 1e6),
                     control = list())
        if(opt$convergence != 0) {
          warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
          any_optim_errors <- TRUE
        }
        c(opt$par, -opt$value)
      })
    }

    pnew <- do.call(rbind, pnew)
    params <- pnew[,c('mu', 'beta', 'phi')]
    ll <- log_likelihood(params, data)

    ll_diff <- (ll - ll_old)  / abs(ll_old) * 100

    if(verbose) {
      message(glue("{i} Current: {ll_old}\tNew: {ll}\tChange: {ll_diff}"))
    }

    if(!is.na(ll_diff)) {
      if(ll_diff < rel_tol) {
        if(verbose) {
          print(glue("EM converged after {i} iterations"))
        }
        break
      }
    }
    ll_old <- ll
  }

  if(any_optim_errors) {
    message("There were errors in optimization of Q function. However, results may still be valid. See errors above.")
  }

  gamma <- p_pi(data, params)
  rlist <- list(
    gamma = gamma,
    mu = params[, 'mu'],
    beta = params[, 'beta'],
    phi = params[, 'phi']
  )

  if(i == max_iter) {
    message("Maximum number of iterations reached; consider increasing max_iter")
  }
  return(rlist)
}


