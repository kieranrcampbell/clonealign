
#' @keywords internal
likelihood_yn <- function(y, l, s_n, params) {
  m <- l * s_n * params[, 'mu']
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
                                  l = data$L[,c],
                                  s_n = data$s[n],
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
Qgr_g <- function(pars, y, l, gamma, data) {
  mu <- pars[1]
  phi <- pars[2]
  gr <- c('mu' = 0, 'phi' = 0)
  for(c in seq_len(data$C)) {
    mu_ng <- mu * data$s * l[c] # N-length vector

    gr_1 <- (y / mu_ng - (y + phi) / (mu_ng + phi) ) * data$s * l[c]
    gr[1] <- gr[1] + sum(gamma[,c] * gr_1)

    gr_2 <- digamma(phi + y) - digamma(phi) - y / (phi + mu_ng) +
      log(phi) + 1 - log(phi + mu_ng) - phi / (phi + mu_ng)
    gr[2] <- gr[2] + sum(gamma[,c] * gr_2)
  }
  -gr
}

#' Computes Q(theta|theta^(t)) for constant phi
#' (function to be optimised under EM)
#' @keywords internal
#'
#' @param pars Parameters to optimise
#' @param y Gene expression for gene
#' @param l Copy number profiles for gene
#' @param gamma Expectation of clone assignments at current EM step
#' @param data Data used
#' @param phi Dispersion value
#'
#' @keywords internal
#'
#' @return The g'th term in the expected complete data log likelihood
Q_g_phi_const <- function(pars, y, l, gamma, data, phi) {
  mu <- pars[1]
  qq <- 0
  for(c in seq_len(data$C)) {
    m <- l[c] * data$s * mu # N length vector for given gene of means
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(gamma[,c] * l_c )
  }
  -qq
}

#' Gradient of Q(theta|theta^(t)) (function to be optimised under EM) for constant phi
#'
#' @param pars Parameters to optimise
#' @param y Gene expression for gene
#' @param l Copy number profiles for gene
#' @param gamma Expectation of clone assignments at current EM step
#' @param data Data used
#' @param phi Dispersion value
#'
#' @keywords internal
#'
#' @return The gradient g'th term in the expected complete data log likelihood
Qgr_g_phi_const <- function(pars, y, l, gamma, data, phi) {
  mu <- pars[1]
  gr <- c('mu' = 0)
  for(c in seq_len(data$C)) {
    mu_ng <- mu * data$s * l[c] # N-length vector

    gr_1 <- (y / mu_ng - (y + phi) / (mu_ng + phi) ) * data$s * l[c]
    gr[1] <- gr[1] + sum(gamma[,c] * gr_1)

  }
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
#' @param rel_tol Relative tolerance in percent below which the log-likelihood is considered converged
#' @param gene_filter_threshold Genes with mean counts below or equal to this threshold will
#' be filtered out (removes genes with no counts by default)
#' @param verbose Logical - should convergence information be printed?
#' @param multithread Should the M-step be performed in parallel using \code{BiocParallel}? Default \code{TRUE}
#' @param bp_param Parameters for multithreaded optimization of Q function. See \code{?bpparam()}
#' @param phi_const Logical - should the dispersion be empirically estimated?
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
                         multithread = TRUE,
                         bp_param = bpparam(),
                         phi_const = FALSE) {

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
    colMeans(Y / s) + 0.01, # Mean
    rep(1, G) # Phi
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

  if(phi_const) {
    # Estimate empirical dispersion
    params[,'phi'] <- empirical_dispersion_estimate(data)
    params[params[,'phi'] < 0,'phi'] <- 0.1
  }


  data$L[data$L == 0] <- 1 # TODO remove

  ll_old <- log_likelihood(params, data)

  lls <- ll_old

  any_optim_errors <- FALSE

  for(i in seq_len(max_iter)) {

    # E step
    gamma <- p_pi(data, params)

    # M step
    if(!phi_const) {
      if(multithread) {
        pnew <- bplapply(seq_len(data$G), function(g) {
          opt <- optim(par = params[g,], # (mu,phi)
                       fn = Q_g,
                       gr = Qgr_g,
                       y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                       method = "L-BFGS-B",
                       lower = c(1e-10, 1e-10),
                       upper = c(max(data$Y), 1e6),
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
                       gr = Qgr_g,
                       y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                       method = "L-BFGS-B",
                       lower = c(1e-10, 1e-10),
                       upper = c(max(data$Y), 1e6),
                       control = list())
          if(opt$convergence != 0) {
            warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
            any_optim_errors <- TRUE
          }
          c(opt$par, -opt$value)
        })
      }

      pnew <- do.call(rbind, pnew)
      params <- pnew[,c('mu', 'phi')]

    } else {
      if(multithread) {
        pnew <- bplapply(seq_len(data$G), function(g) {
          opt <- optim(par = params[g,1], # (mu,)
                       fn = Q_g_phi_const,
                       gr = Qgr_g_phi_const,
                       y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                       phi = params[g,2],
                       method = "L-BFGS-B",
                       lower = c(1e-10),
                       upper = c(max(data$Y)),
                       control = list())
          if(opt$convergence != 0) {
            warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
            any_optim_errors <- TRUE
          }
          c(opt$par, -opt$value)
        }, BPPARAM = bp_param)
      } else {
        pnew <- lapply(seq_len(data$G), function(g) {
          opt <- optim(par = params[g,1], # (mu,)
                       fn = Q_g_phi_const,
                       gr = Qgr_g_phi_const,
                       y = data$Y[,g], l = data$L[g,], gamma = gamma, data = data,
                       phi = params[g,2],
                       method = "L-BFGS-B",
                       lower = c(1e-10),
                       upper = c(max(data$Y)),
                       control = list())
          if(opt$convergence != 0) {
            warning(glue("L-BFGS-B optimization of Q function warning: {opt$message}"))
            any_optim_errors <- TRUE
          }
          c(opt$par, -opt$value)
        })
      }
      pnew <- do.call(rbind, pnew)
      params <- cbind(pnew[,1], params[,'phi'])
      colnames(params) <- c("mu", "phi")
    }


    ll <- log_likelihood(params, data)

    ll_diff <- (ll - ll_old)  / abs(ll_old) * 100

    lls <- c(lls, ll)

    if(verbose) {
      message(glue("{i} Current: {ll_old}\tNew: {ll}\tChange: {ll_diff}"))
    }

    if(!is.na(ll_diff)) {
      if(ll_diff < rel_tol) {
        if(verbose) {
          message(glue("EM converged after {i} iterations"))
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
    phi = params[, 'phi'],
    lls = lls
  )

  if(i == max_iter) {
    message("Maximum number of iterations reached; consider increasing max_iter")
  }
  return(rlist)
}

#' @importFrom stats loess
#' @importFrom tibble data_frame
#'
empirical_dispersion_estimate <- function(data) {

  Y <- data$Y
  dsf <- data$s

  w <- colVars(Y / dsf)
  q <- colMeans(Y / dsf)

  z <- q * mean(1 / dsf)

  l <- loess(log(w) ~ log(q))

  df <- data_frame(q, w, z, w_q = exp(predict(l)))

  df <- dplyr::mutate(df, v = w_q - z, mu = sum(dsf) * q, var = mu + sum(dsf^2) * v) %>%
    dplyr::mutate(dispersion = 1 / ((var - mu) / (mu^2)))

  df$dispersion

}


