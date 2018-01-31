
#' @keywords internal
likelihood_yn <- function(y, L, s_n, pi, params) {
  m <- L[, pi] * s_n * params[, 'mu']
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
#' @keywords internal
#'
#' @return The g'th term in the expected complete data log likelihood
Q_g <- function(pars, y, l, data, pi_traces, rho_traces, lambda, l_g_hat) {
  mu <- pars[1]
  beta <- pars[2]
  phi <- pars[3]
  qq <- 0
  for(it in seq_along(rho_traces)) {
    rho <- rho_traces[it]
    pi <- pi_traces[it,]
    m <- data$s * ((1 - rho) * mu + rho * beta * l[pi])
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(l_c )
  }
  -qq + lambda * (mu - l_g_hat * beta)^2
}

#' Gradient of Q(theta|theta^(t)) (function to be optimised under EM)
#'
#' @param pars Parameters to optimise
#' @param y Gene expression for gene
#' @param l Copy number profiles for gene
#' @param gamma Expectation of clone assignments at current EM step
#' @param data Data used
#'
#' @keywords internal
#'
#' @return The gradient g'th term in the expected complete data log likelihood
Qgr_g <- function(pars, y, l, data, pi_traces, rho_traces, lambda, l_g_hat) {
  mu <- pars[1]
  beta <- pars[2]
  phi <- pars[3]
  gr <- c('mu' = 0, 'beta' = 0, 'phi' = 0)

  for(it in seq_along(rho_traces)) {
    rho <- rho_traces[it]
    pi <- pi_traces[it,]

    mu_ng <- data$s * ((1 - rho) * mu +  rho * beta * l[pi])

    d_py_d_mung <- (y / mu_ng - (y + phi) / (mu_ng + phi))

    gr_1 <- d_py_d_mung * data$s * (1 - rho)
    gr[1] <- gr[1] + sum(gr_1)

    gr_2 <- d_py_d_mung * data$s * rho * l[pi]
    gr[2] <- gr[2] + sum(gr_2)

    gr_3 <- digamma(phi + y) - digamma(phi) - y / (phi + mu_ng) +
      log(phi) + 1 - log(phi + mu_ng) - phi / (phi + mu_ng)
    gr[3] <- gr[3] + sum(gr_3)
  }

  # Add on penalty - remember the minus sign!
  gr[1] <- gr[1] - 2 * lambda * (mu - l_g_hat * beta)
  gr[2] <- gr[2] - 2 * lambda * (l_g_hat * beta - mu) * l_g_hat

  -gr
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
     c(sum(dnbinom2(data$Y[n,], beta * data$L[,c] * data$s[n], size = phi)),
       sum(dnbinom2(data$Y[n,], mu * data$s[n], size = phi)))
    })
    ll <- ll + logSumExp(as.vector(pnc))
   }
  ll - sum((2 * mu - beta)^2)
}



# Main EM algorithm -------------------------------------------------------


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
inference_em <- function(Y,
                         L,
                         s = NULL,
                         max_iter = 100,
                         rel_tol = 1e-5,
                         gene_filter_threshold = 0,
                         verbose = TRUE,
                         multithread = FALSE,
                         mcmc_iteration_multiple = 5,
                         bp_param = bpparam(),
                         rho_init = NULL,
                         lambda = 1) {

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
  L <- t( t(L) / colMeans(L) )

  # Initialise
  params <- cbind(
    colMeans(Y / s) + 0.01,
    colMeans(Y / s) + 0.01,
    rep(1, G)
  )
  colnames(params) <- c("mu", "beta", "phi")

  data <- list(
    Y = Y,
    L = L,
    s = s,
    N = N,
    G = G,
    C = C,
    l_g_hat = rowMeans(L)
  )

  rho <- sample(c(0,1), data$G, replace = TRUE)
  if(!is.null(rho_init)) {
    rho <- rho_init
  }

  data$L[data$L == 0] <- 1

  ll_old <- log_likelihood(params, data)

  lls <- ll_old

  any_optim_errors <- FALSE

  for(i in seq_len(max_iter)) {

    # E step - now with added Gibbs sampling
    gibbs_samples <- gibbs_pi_rho(rho, data, params, n_iter = 10 + mcmc_iteration_multiple * i)
    rho <- gibbs_samples$rho_traces[1,]

    # M step
    n_optim_errors <- 0
    pnew <- lapply(seq_len(data$G), function(g) {
      opt <- optim(par = params[g,], # (mu,beta,phi)
                   fn = Q_g,
                   gr = Qgr_g,
                   y = data$Y[,g], l = data$L[g,], data = data,
                   pi_traces = gibbs_samples$pi_trace,
                   rho_traces = gibbs_samples$rho_trace[,g],
                   l_g_hat = data$l_g_hat[g],
                   lambda = 1,
                   method = "L-BFGS-B",
                   lower = c(1e-10, 1e-10, 1e-10),
                   upper = c(max(data$Y), 1e6, 1e6),
                   control = list())
      if(opt$convergence != 0) {
        n_optim_errors <<- n_optim_errors + 1
      }
      c(opt$par, -opt$value)
    })

    print(glue("{n_optim_errors} optimization errors"))


    pnew <- do.call(rbind, pnew)
    params <- pnew[,c('mu', 'beta', 'phi')]
    ll <- log_likelihood(params, data)

    ll_diff <- (ll - ll_old)  / abs(ll_old) * 100

    lls <- c(lls, ll)

    if(verbose) {
      message(glue("{i} Current: {ll_old}\tNew: {ll}\tChange: {ll_diff}"))
    }

    if(!is.na(ll_diff)) {
      if(abs(ll_diff) < rel_tol) { # TODO remove abs
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

  gibbs_samples <- gibbs_pi_rho(rho, data, params, n_mcmc_iter)
  rlist <- list(
    gamma = gamma,
    mu = params[, 'mu'],
    beta = params[,'beta'],
    phi = params[, 'phi'],
    lls = lls,
    gibbs_samples = gibbs_samples
  )

  if(i == max_iter) {
    message("Maximum number of iterations reached; consider increasing max_iter")
  }
  return(rlist)
}


