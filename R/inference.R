

# -----
# List `params`` with components
# - mu
# - phi
# - pi
# - gamma
#
# List `hyperparams` with components
# - a
# - m_mu, s_mu
# - m_phi, s_phi
#





dnbinom2 <- function(x, mu, size) {
  dnbinom(x, size = size, mu = mu, log = TRUE)
}

calculated_means <- function(Y, L, params, s) {
  M <- matrix(NA, nrow = N, ncol = G)
  for(g in 1:G) {
    for(n in 1:N) {
      M[n, g] <- L[g, sim$pi[n] ] * params$mu[g] * s[n]
    }
  }
  M
}

likelihood_yn <- function(y, n, L, params, s) {
  ll <- 0

  G <- length(y)

  # Likelihood
  for(g in 1:G) {
    if(L[g, params$pi[n] ] > 0) {
      m <- L[g, params$pi[n] ] * s[n] * params$mu[g]
      ll <- ll + dnbinom2(y[g], mu = m, size = params$phi[g])
    } else {
      ll <- ll + dpois(y[g], 0.1, log = TRUE)
    }
  }


  ll
}


log_ratio_mug <- function(y, L, g, mu_g_prime, mu_g, params, hyperparams, s) {
  ll <- 0
  r <- params$phi[g]
  # Likelihood
  for(n in seq_along(y)) {
    if(L[g, params$pi[n] ] > 0) {
      m_prime <- L[g, params$pi[n] ] * s[n] * mu_g_prime
      m <- L[g, params$pi[n] ] * s[n] * mu_g
      ll <- ll + y[n] * log(m_prime) - y[n] * log(r + m_prime) - r * log(r + m_prime)
      ll <- ll + -y[n] * log(m) + y[n] * log(r + m) + r * log(r + m)
    } # TODO reimplement when L == 0
  }

  prior <- dnorm(log(mu_g_prime), params$nu_mu, 1 / sqrt(params$tau_mu), log = TRUE ) -
    dnorm(log(mu_g), params$nu_mu, 1 / sqrt(params$tau_mu), log = TRUE )
  jacobian <- log(mu_g) - log(mu_g_prime)
  ll + prior + jacobian
}

log_ratio_phig <- function(y, L, g, phi_g_prime, phi_g, params, hyperparams, s) {
  mu <- params$mu[g] * L[g, params$pi] * s # Mean of each cell for gene g
  ll <- sum(dnbinom2(y, mu, phi_g_prime)) - sum(dnbinom2(y, mu, phi_g))
  prior <- dnorm(log(phi_g_prime), params$nu_phi, 1 / sqrt(params$tau_phi), log = TRUE) -
    dnorm(log(phi_g), params$nu_phi, 1 / sqrt(params$tau_phi), log = TRUE)
  jacobian <- log(phi_g) - log(phi_g_prime)
  ll + prior + jacobian
}

joint_density <- function(Y, L, params, hyperparams, s) {
  C <- ncol(L)
  N <- nrow(Y)

  ll <- 0

  for(n in 1:N) {
    ll <- ll + likelihood_yn(Y[n,], n, L, params, s)
  }

  n <- tabulate(params$pi, C)
  ll <- ll + sum(params$gamma^n)

  # Prior
  prior <- 0
  prior <- prior + sum( dnorm(log(params$mu), hyperparams$m_mu,
                             hyperparams$s_mu, log = TRUE ) )
  prior <- prior + sum(dnorm(log(params$phi), hyperparams$m_phi,
                          hyperparams$s_phi, log = TRUE ))

  prior <- prior + ddirichlet(params$gamma, hyperparams$a)

  ll + prior
}

propose_normal <- function(x, delta) {
  rnorm(length(x), x, sqrt(delta * (1 + x^2)))
}

density_normal <- function(x_new, x, delta) {
  sum(dnorm(x_new, x, sqrt(delta * (1 + x^2)), log = TRUE))
}


#' @importFrom matrixStats colVars
#' @export
basic_dispersion_estimates <- function(Y, s) {
  variances <- colVars(Y / s)
  means <- colMeans(Y / s)
  phi <- means^2 / (variances - means)

  fit <- loess(phi ~ means, data.frame(phi, means))
  phi_hat <- predict(fit)
  stopifnot(all(phi_hat > 0))
  phi_hat
}



# profvis({

#' @importFrom glue glue
#' @importFrom MCMCpack rdirichlet
#' @export
infer_simple <- function(Y, L, s,
                  delta_phi = 1e-3,
                  delta_mu = 1e-2,
                  n_iter = 200,
                  burn = n_iter / 2,
                  phi_init = 1,
                  print_progress_every = n_iter / 10,
                  phi_hat = NULL) {

  print(glue("Testing glue"))
  N <- nrow(Y)
  G <- ncol(Y)
  C <- ncol(L)

  # Sanity checks
  stopifnot(nrow(L) == G)
  stopifnot(length(s) == N)
  stopifnot(all(s > 0))

  if(!is.null(phi_hat)) stopifnot(length(phi_hat) == G)

  # Initialise
  params <- list(
    mu = colMeans(Y / s),
    phi = phi_hat,
    pi = sample(1:C, N, replace = TRUE),
    gamma = rep(1, C) / C
  )


  hyperparams <- list(
    a = rep(1 / C, C),
    m_mu = mean(Y),
    s_mu = 0.1,
    m_phi = 0,
    s_phi = 1
  )


  pi_trace <- matrix(NA, nrow = n_iter, ncol = N)
  gamma_trace <- matrix(NA, nrow = n_iter, ncol = C)

  phi_trace <- rep(NA, n_iter)

  mu_trace <- matrix(NA, nrow = n_iter, ncol = G)

  phi_accept <- rep(0, n_iter)

  mu_accept <- rep(0, n_iter)


  for(i in 1:n_iter) {

    if(i %% print_progress_every == 0) {
      print(glue("Iteration {i}"))
    }

    # Sample pi
    pc <- matrix(NA, N, C)

    for(n in 1:N) {
      for(c in 1:C) {
        params_tmp <- params
        params_tmp$pi[n] <- c
        pc[n,c] <- likelihood_yn(Y[n,], n, L, params_tmp, s) + log(params$gamma[c])
      }
    }

    probs <- exp(pc - apply(pc, 1, logSumExp))
    pi_new <- sapply(1:N, function(n) sample(1:C, 1, prob = probs[n, ]))

    params$pi <- pi_new

    pi_trace[i,] <- pi_new

    # Sample Gamma
    tab <- tabulate(params$pi, C)
    gamma_new <- rdirichlet(1, tab + hyperparams$a)
    gamma_trace[i,] <- gamma_new

    params$gamma <- gamma_new

    ## MH proposal for phi
    if(FALSE) {
    phi_t_new <- propose_normal(log(params$phi), delta_phi)
    params_new <- params
    params_new$phi <- exp(phi_t_new)
    log_rat <- joint_density(Y, L, params_new, hyperparams, s) -
      joint_density(Y, L, params, hyperparams, s)
    log_rat <- log_rat + density_normal(log(params$phi), phi_t_new, delta_phi) -
      density_normal(phi_t_new, log(params$phi), delta_phi)
    log_rat <- log_rat + phi_t_new - log(params$phi) # Adjust Jacobian


    if(log_rat > log(runif(1))) {
      params$phi <- exp(phi_t_new)
      phi_accept[i] <- 1
    }
    phi_trace[i] <- params$phi
    }


    ## MH proposal for mu
    if(TRUE) {
    for(g in 1:G) {
      mu <- params$mu[g]
      mu_t_new <- propose_normal(log(mu), delta_mu)

      log_rat <- log_ratio_mug(Y[,g], L, g, exp(mu_t_new),
                               mu, params, hyperparams, s)
      log_rat <- log_rat + density_normal(log(mu), mu_t_new, delta_mu) -
        density_normal(mu_t_new, log(mu), delta_mu)
      log_rat <- log_rat + mu_t_new - log(mu) # Adjust Jacobian

      if(log_rat > log(runif(1))) {
        params$mu[g] <- exp(mu_t_new)
        mu_accept[i] <- mu_accept[i] + 1
      }
      mu_trace[i,] <- params$mu
    }
    mu_accept[i] <- mu_accept[i] / G
    }

  }


  print(glue("Mean phi accept: {mean(phi_accept[-seq_len(burn)])}"))
  print(glue("Mean mu accept: {mean(mu_accept[-seq_len(burn)])}"))

  list(
    traces = list(
      mu = mu_trace,
      phi = phi_trace,
      pi = pi_trace,
      gamma = gamma_trace
    ),
    accepts = list(
      mu = mu_accept,
      phi = phi_accept
    )

  )

}

#' Sample either nu_mu or nu_phi
gibbs_sample_mean <- function(mu_tilde, tau, tau_0, nu_0) {
  G <- length(mu_tilde)
  mu_tilde_sum <- sum(mu_tilde)
  precision <- tau * G + tau_0
  mean_unorm <- tau * mu_tilde_sum + tau_0 * nu_0
  rnorm(1, mean_unorm / precision, 1 / sqrt(precision))
}

#' Sample either tau_mu or tau_phi
gibbs_sample_precision <- function(a, b, mu_tilde, nu) {
  G <- length(mu_tilde)
  square_sum <- 0.5 * sum((mu_tilde - nu)^2)

  rgamma(1, shape = a + G / 2, rate = b + square_sum)
}


#' @importFrom glue glue
#' @importFrom MCMCpack rdirichlet
#' @export
infer_hierarchical <- function(Y, L, s,
                         delta_phi = 5e-1,
                         delta_mu = 1e-1,
                         n_iter = 200,
                         burn = n_iter / 2,
                         phi_init = 1,
                         print_progress_every = n_iter / 10,
                         phi_hat = NULL) {

  N <- nrow(Y)
  G <- ncol(Y)
  C <- ncol(L)

  # Sanity checks
  stopifnot(nrow(L) == G)
  stopifnot(length(s) == N)
  stopifnot(all(s > 0))

  if(!is.null(phi_hat)) {
    stopifnot(length(phi_hat) == G)
  } else {
    phi_hat <- rep(phi_init, G)
  }

  # Initialise
  params <- list(
    mu = colMeans(Y / s) + 0.01,
    phi = phi_hat,
    pi = sample(1:C, N, replace = TRUE),
    gamma = rep(1, C) / C,
    nu_mu = 0, # Hierarchical mean on mu
    tau_mu = 1, # Hierarchical precision on mu
    nu_phi = 0, # Hierarchical mean on phi
    tau_phi = 1 # Hierarchical precision on phi
  )


  hyperparams <- list(
    a = rep(1 / C, C),
    m_mu = mean(Y),
    s_mu = 0.1,
    m_phi = 0,
    s_phi = 1,
    a_mu = 2,
    b_mu = 1,
    a_phi = 2,
    b_phi = 1,
    nu_phi_0 = 0,
    tau_phi_0 = 0.01,
    nu_mu_0 = 0,
    tau_mu_0 = 0.01
  )

  # Always normalise L
  L <- t( t(L) / colMeans(L) )

  pi_trace <- matrix(NA, nrow = n_iter, ncol = N)
  gamma_trace <- matrix(NA, nrow = n_iter, ncol = C)
  phi_trace <- matrix(NA, nrow = n_iter, ncol = G)
  mu_trace <- matrix(NA, nrow = n_iter, ncol = G)
  nu_mu_trace <- rep(NA, n_iter)
  tau_mu_trace <- rep(NA, n_iter)
  nu_phi_trace <- rep(NA, n_iter)
  tau_phi_trace <- rep(NA, n_iter)

  phi_accept <- rep(0, n_iter)

  mu_accept <- rep(0, n_iter)


  for(i in 1:n_iter) {

    if(i %% print_progress_every == 0) {
      print(glue("Iteration {i}"))
    }

    # Sample pi
    pc <- matrix(NA, N, C)

    for(n in 1:N) {
      for(c in 1:C) {
        params_tmp <- params
        params_tmp$pi[n] <- c
        pc[n,c] <- likelihood_yn(Y[n,], n, L, params_tmp, s) + log(params$gamma[c])
      }
    }

    probs <- exp(pc - apply(pc, 1, logSumExp))
    pi_new <- sapply(1:N, function(n) sample(1:C, 1, prob = probs[n, ]))

    params$pi <- pi_new

    pi_trace[i,] <- pi_new

    # Sample Gamma
    tab <- tabulate(params$pi, C)
    gamma_new <- rdirichlet(1, tab + hyperparams$a)
    gamma_trace[i,] <- gamma_new

    params$gamma <- gamma_new

    ## MH proposal for mu
    for(g in 1:G) {
      mu <- params$mu[g]
      log_mu <- log(mu)
      mu_t_new <- propose_normal(log_mu, delta_mu)

      log_rat <- log_ratio_mug(Y[,g], L, g, exp(mu_t_new),
                               mu, params, hyperparams, s)
      log_rat <- log_rat + density_normal(log_mu, mu_t_new, delta_mu) -
        density_normal(mu_t_new, log_mu, delta_mu)
      log_rat <- log_rat + mu_t_new - log_mu # Adjust Jacobian

      if(is.na(log_rat)) {
        print(mu_t_new)
        print(log_mu)
        print(params$phi[g])
        stop("Waah")
      }
      if(log_rat > log(runif(1))) {
        params$mu[g] <- exp(mu_t_new)
        mu_accept[i] <- mu_accept[i] + 1
      }
    }
    mu_trace[i,] <- params$mu
    mu_accept[i] <- mu_accept[i] / G

    ## MH proposal for phi
    for(g in 1:G) {
      log_phi <- log(params$phi[g])
      phi_t_new <- propose_normal(log_phi, delta_phi)

      log_rat <- log_ratio_phig(Y[,g], L, g, exp(phi_t_new), params$phi[g],
                                params, hyperparams, s)
      log_rat <- log_rat + density_normal(log_phi, phi_t_new, delta_phi) -
        density_normal(phi_t_new, log_phi, delta_phi)
      log_rat <- log_rat + phi_t_new - log_phi # Adjust Jacobian

      if(log_rat > log(runif(1))) {
        params$phi[g] <- exp(phi_t_new)
        phi_accept[i] <- phi_accept[i] + 1
      }
    }
    phi_trace[i,] <- params$phi
    phi_accept[i] <- phi_accept[i] / G

    # plot(params$phi)
    # print(glue("Phi accept: {phi_accept[i]}"))

    # Gibbs updates for (nu, phi) "hyperparameters"
    params$nu_mu <- gibbs_sample_mean(log(params$mu), params$tau_mu,
                                      hyperparams$tau_mu_0, hyperparams$nu_mu_0)
   params$nu_phi <- gibbs_sample_mean(log(params$phi), params$tau_phi,
                                    hyperparams$tau_phi_0, hyperparams$nu_phi_0)

    params$tau_mu <- gibbs_sample_precision(hyperparams$a_mu, hyperparams$b_mu,
                                             log(params$mu), params$nu_mu)
    params$tau_phi <- gibbs_sample_precision(hyperparams$a_phi, hyperparams$b_phi,
                                             log(params$phi), params$nu_phi)

    nu_mu_trace[i] <- params$nu_mu
    nu_phi_trace[i] <- params$nu_phi

    tau_mu_trace[i] <- params$tau_mu
    tau_phi_trace[i] <- params$tau_phi

  } # End MCMC iteration for loop




  print(glue("Mean phi accept: {mean(phi_accept[-seq_len(burn)])}"))
  print(glue("Mean mu accept: {mean(mu_accept[-seq_len(burn)])}"))

  list(
    traces = list(
      mu = mu_trace,
      phi = phi_trace,
      pi = pi_trace,
      gamma = gamma_trace,
      nu_mu = nu_mu_trace,
      nu_phi = nu_phi_trace,
      tau_mu = tau_mu_trace,
      tau_phi = tau_phi_trace
    ),
    accepts = list(
      mu = mu_accept,
      phi = phi_accept
    )

  )

}


