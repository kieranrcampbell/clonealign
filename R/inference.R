

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
      ll <- ll + y[n] * log(m_prime) - (y[n] + r) * log(r + m_prime)
      ll <- ll + -y[n] * log(m) + (y[n] + r)* log(r + m)
    } else {
      # TODO reimplement when L == 0
      stop("Needs implemented")
    }
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



