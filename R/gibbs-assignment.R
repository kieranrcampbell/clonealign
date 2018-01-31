

#' Computes p(pi_n = c | data, params, rho)
#' @importFrom matrixStats logSumExp
p_pi_n <- function(n, data, params, rho) {
  p_pi_n_equals_c <- sapply(seq_len(data$C), function(c) {
    m <- data$s[n] * ((1 - rho) * params[,'mu'] + rho * params[,'beta'] * data$L[,c])
    sum(dnbinom2(data$Y[n,], mu = m, size = params[,'phi']))
  })
  exp(p_pi_n_equals_c - logSumExp(p_pi_n_equals_c))
}

#' @import matrixStats logSumExp
p_rho_g <- function(g, data, params, pi) {
  p_rho_g_equals_i <- sapply(0:1, function(i) {
    m <- data$s * ((1 - i) * params[g,'mu'] + i * params[g,'beta'] * data$L[g, pi])
    sum(dnbinom2(data$Y[,g], m, params[g,'phi']))
  })
  exp(p_rho_g_equals_i - logSumExp(p_rho_g_equals_i))
}


#' @param n_iter Number of MCMC iterations
#' @param burn_prop Proportion of beginning traces to burn
gibbs_pi_rho <- function(rho, data, params, n_iter = 20, burn_prop = 0.5) {

  n_burn <- floor(0.5 * n_iter)
  n_iter_keep <- n_iter - n_burn

  pi_trace <- matrix(NA, nrow = n_iter_keep, ncol = data$N)
  rho_trace <- matrix(NA, nrow = n_iter_keep, ncol = data$G)

  for(it in seq_len(n_iter)) {
    p_pi <- t(sapply(seq_len(data$N), p_pi_n, data, params, rho))
    pi <- apply(p_pi, 1, function(p) sample(seq_len(data$C), size = 1, prob = p))

    p_rho <- t(sapply(seq_len(data$G), p_rho_g, data, params, pi))
    rho <- apply(p_rho, 1, function(p) sample(0:1, size = 1, prob = p))

    if(it > n_burn) {
      it_keep <- it - n_burn
      pi_trace[it_keep,] <- pi
      rho_trace[it_keep,] <- rho
    }
  }

  list(
    pi_trace = pi_trace,
    rho_trace = rho_trace
  )
}

clone_probs_from_gibbs <- function(pi_trace, C) {
  sapply(seq_len(C), function(c) {
    colMeans(pi_trace == c)
  })
}

rho_probs_from_gibbs <- function(rho_trace) {
  sapply(0:1, function(c) {
    colMeans(rho_trace == c)
  })
}

