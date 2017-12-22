
#' @keywords internal
likelihood_yn <- function(y, L, s_n, pi, params) {
  m <- L[, pi] * s_n * params[, 'mu']
  phi <- params[, 'phi']
  ll <- sum(dnbinom2(y, mu = m, size = phi))
  ll
}

#' Computes gamma_{nc} = p(pi_n = c), returning
#' N by C matrix
#'
#' @importFrom matrixStats logSumExp
p_pi <- function(data, params) {
  gamma <- matrix(NA, nrow = data$N, ncol = data$C)
  for(n in seq_len(data$N)) {
    for(c in seq_len(data$C)) {
      gamma[n,c] <- likelihood_yn(y = data$Y[n,],
                                  L = data$L,
                                  s = data$s[n],
                                  pi = c,
                                  params = params)
    }
    gamma[n,] <- exp(gamma[n,] - logSumExp(gamma[n,]))
  }
  gamma
}

#' Computes Q(theta|theta^(t))
#' (function to be optimised under EM)
Q_g <- function(pars, y, l, gamma, data) {
  mu <- pars[1]
  phi <- pars[2]
  qq <- 0
  for(c in seq_len(data$C)) {
    m <- l[c] * data$s * mu # N length vector for given gene of means
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(gamma[,c] * l_c )
    # print(qq)
  }
  -qq
}

#' Computes map clone assignment given EM object
#' @export
clone_assignment <- function(em) {
  apply(em$gamma, 1, which.max)
}



#' @keywords internal
log_likelihood <- function(params, data) {
  ll <- 0
  mu <- params[,'mu']
  phi <- params[,'phi']

  for(n in seq_len(data$N)) {
    pnc <- sapply(seq_len(data$C), function(c) {
     sum(dnbinom2(data$Y[n,], mu * data$L[,c] * data$s[n], size = phi))
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
#' @param tol Relative tolerance in percent below which the log-likelihood is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#'
#' @importFrom glue glue
#' @importFrom BiocParallel bplapply
#'
#' @export
inference_em <- function(Y, L, s = NULL, max_iter = 100, tol = 1e-5,
                          gene_filter_threshold = 0,
                          bp_param = bpparam()) {

  zero_gene_means <- colMeans(Y) <= gene_filter_threshold
  message(glue("Removing {sum(zero_gene_means)} genes with low counts"))
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
  L <- t( t(L) / colMeans(L) )

  # Initialise
  params <- cbind(
    colMeans(Y / s) + 0.01,
    rep(1, G)
  )
  colnames(params) <- c("mu", "phi")

  data <- list(
    Y = Y,
    L = L,
    s = s,
    N = N,
    G = G,
    C = C
  )

  ll_old <- NA

  for(i in seq_len(max_iter)) {

    # E step
    gamma <- p_pi(data, params)

    # M step

    pnew <- bplapply(seq_len(data$G), function(g) {
      opt <- optim(par = params[g,], # (mu,phi)
                   fn = Q_g,
                   # gr = grad_g,
                   y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                   method = "L-BFGS-B",
                   lower = c(1e-10, 1e-10),
                   upper = c(max(data$Y), 1e6),
                   control = list())
      if(opt$convergence != 0) {
        # TODO - deal with nonconvergence somehow
        # print(opt$message)
      }
      c(opt$par, -opt$value)
    }, BPPARAM = bp_param)
    pnew <- do.call(rbind, pnew)
    params <- pnew[,c('mu', 'phi')]
    ll <- log_likelihood(params, data)

    ll_diff <- (ll - ll_old)  / abs(ll_old) * 100

    message(glue("Old: {ll_old}\tNew: {ll}\tChange: {ll_diff}"))
    # message(glue("New log-likelihood: {q_new}"))

    if(!is.na(ll_diff)) {
      if(ll_diff < tol) {
        print(glue("Converged after {i} iterations"))
        break
      }
    }
    ll_old <- ll
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


