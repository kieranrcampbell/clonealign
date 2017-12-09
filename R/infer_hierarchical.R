

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

#' To MCMC object
#' @export
#' @importFrom coda mcmc
to_mcmc <- function(mcmc_output) {
  traces <- mcmc_output$traces
  var_names <- names(traces)
  for(var_name in var_names) {
    if(!is.matrix(traces[[var_name]])) {
      traces[[var_name]] <- matrix(traces[[var_name]])
    }
    nparam <- ncol(traces[[var_name]])
    colnames(traces[[var_name]]) <- paste0(var_name, seq_len(nparam))
  }
  traces_mat <- do.call("cbind", traces)
  mcmc(traces_mat)
}
