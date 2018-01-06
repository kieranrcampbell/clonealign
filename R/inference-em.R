
#' @keywords internal
likelihood_yn <- function(y, L, S_n, pi, params) {
  m <- L[, pi] * params[, 'mu'] # TODO - calculation of this can be moved outside the loop
  m <- m / sum(m) * S_n
  phi <- params[, 'phi']
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
                                  S_n = data$S[n],
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
  phi <- pars[2]
  qq <- 0
  for(c in seq_len(data$C)) {
    m <- l[c] * data$s * mu # N length vector for given gene of means
    # m[m == 0] <- 1
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(gamma[,c] * l_c )
  }
  -qq
}

Q <- function(pars, data, gamma) {
  G <- data$G
  mu <- c(1, pars[seq_len(G-1)])
  phi <- pars[seq(G, 2*G-1)]
  qq <- 0
  for(n in 1:data$N) {
    for(c in 1:data$C) {
      m <- data$L[,c] * mu
      m <- m / sum(m) * data$S[n]
      l_c <- dnbinom2(data$Y[n,], mu = m, size = phi)
      qq <- qq + gamma[n,c] * sum(l_c)
    }
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
  phi <- params[,'phi']

  for(n in 1:data$N) {
    for(c in 1:data$C) {
      m <- data$L[,c] * mu
      m <- m / sum(m) * data$S[n]
      ll <- ll + sum(dnbinom2(data$Y[n,], mu = m, size = phi))
    }
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
inference_em <- function(Y, L, S = NULL, max_iter = 100, rel_tol = 1e-5,
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


  if(is.null(S)) {
    S = rowSums(Y)
  }
  stopifnot(length(S) == N)
  stopifnot(all(s > 0))


  # Initialise
  params <- cbind(
    colMeans(Y / S) + 0.01, # TODO - how to initialise mu ?
    rep(1, G)
  )

  params[-1,1] <- params[-1,1] / params[1,1]
  params[1,1] <- 1

  colnames(params) <- c("mu", "phi")

  data <- list(
    Y = Y,
    L = L,
    S = S,
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

    # TODO delete
    # sourceCpp("../clonealign/Q.cpp")
    # system.time(Qcpp(par, gamma, data$Y, data$L, data$S))
    # system.time(Q(par, data = data, gamma = gamma))
    # Q(par, data = data, gamma = gamma)
    # Qcpp(par, gamma, data$Y, data$L, data$S)

    # M step
    par <- c(params[-1,1], params[,2]) # oh isn't this fun
    opt <- optim(par = par, # (mu[-1],phi)
                 fn = Qcpp,
                 gamma = gamma, Y = data$Y, L = data$L, S = data$S,
                 method = "L-BFGS-B",
                 lower = rep(1e-6, length(par)),
                 # upper = rep(1e6, length(par)),
                 control = list())

    s <- spg(par = par, # (mu[-1],phi)
             fn = Qcpp,
             gamma = gamma, Y = data$Y, L = data$L, S = data$S,
             lower = rep(1e-6, length(par)))

    params[,1] <- c(1, opt$par[1:(G-1)])
    params[,2] <- opt$par[G:(2*G-1)]


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
    phi = params[, 'phi']
  )

  if(i == max_iter) {
    message("Maximum number of iterations reached; consider increasing max_iter")
  }
  return(rlist)
}


