
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
#' @param tflow_res List returned by \code{inference_tflow}
#' @return A vector of maximum likelihood clone assignments
#' @keywords internal
clone_assignment <- function(tflow_res) {
  apply(tflow_res$gamma, 1, which.max)
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

#' EM inference with tensorflow
#' @export
#' @importFrom glue glue
#' @import tensorflow
inference_tflow <- function(Y_dat,
                            L_dat,
                            max_em_iter = 100,
                            rel_em_tol = 1e-5,
                            max_adam_iter = 100,
                            rel_adam_tol = 1e-5,
                            gene_filter_threshold = 0,
                            verbose = TRUE) {

  zero_gene_means <- colMeans(Y_dat) <= gene_filter_threshold

  if(verbose) {
    message(glue("Removing {sum(zero_gene_means)} genes with low counts"))
  }

  Y_dat <- Y_dat[, !zero_gene_means]
  L_dat <- L_dat[!zero_gene_means,]

  N <- nrow(Y_dat)
  G <- ncol(Y_dat)
  C <- ncol(L_dat)

  # Sanity checks
  stopifnot(nrow(L_dat) == G)

  data <- list(
    Y = Y_dat,
    L = L_dat,
    S = rowSums(Y_dat),
    N = N,
    G = G,
    C = C
  )

  mu_guess <- colMeans(data$Y)
  mu_guess_1 <- mu_guess[-1] / mu_guess[1]

  beta_init <- rep(1, data$G)

  LOWER_BOUND <- 1e-6

  # Tensorflow variables
  # Data
  Y <- tf$placeholder(tf$float32, shape = c(N,G))
  L <- tf$placeholder(tf$float32, shape = c(G,C))
  S <- tf$placeholder(tf$float32, shape = N)

  # Unconstrained variables
  mu_log <- tf$Variable(tf$constant(log(mu_guess_1)))
  beta_log <- tf$Variable(tf$constant(log(beta_init)))
  phi_log <- tf$Variable(tf$zeros(G))

  # Constrained variables
  mu <- tf$exp(mu_log) + LOWER_BOUND
  beta <- tf$exp(beta_log) + LOWER_BOUND
  phi <- tf$exp(phi_log) + LOWER_BOUND

  mu_measured <- tf$concat(list(tf$reshape(tf$constant(1.0), shape(1L)), mu), axis = 0L)

  mu_gc <- tf$transpose(mu_measured * (1 + beta * tf$transpose(L)))
  mu_c <- tf$reduce_sum(mu_gc, 0L)
  mu_gc = mu_gc / mu_c

  mu_ncg <- tf$transpose(tf$einsum('gc,n -> gcn', mu_gc, S), c(2L,0L,1L))

  # Tflow uses p NB formulation
  p <- (mu_ncg) / (mu_ncg + phi)

  y_pdf <- tf$contrib$distributions$NegativeBinomial(probs = p, total_count = phi)

  Y_tiled <- tf$transpose(tf$tile(tf$reshape(Y, c(N, G, 1L)), c(1L,1L,C)), c(0L, 2L, 1L))

  y_log_prob <- y_pdf$log_prob(Y_tiled)

  p_y_on_c <- tf$reduce_sum(y_log_prob, 2L)

  p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c, 1L), c(1L, -1L))

  gamma <- tf$exp(tf$transpose(tf$transpose(p_y_on_c) - p_y_on_c_norm))

  # Q function
  gamma_fixed <- tf$placeholder(dtype = tf$float32, shape = c(N, C))

  Q_y <- tf$einsum('nc,ncg->g', gamma_fixed, y_log_prob)
  Q <- -(tf$reduce_sum(Q_y))

  optimizer <- tf$train$AdamOptimizer()
  train <- optimizer$minimize(Q, var_list = list(mu_log, beta_log, phi_log))

  eta_y <- tf$reduce_sum(y_log_prob, 2L)
  L_y <- tf$reduce_sum(tf$reduce_logsumexp(eta_y, 1L))

  # Inference
  mu_final <- phi_final <- NA
  l <- NA

  sess <- tf$Session()
  init <- tf$global_variables_initializer()

  sess$run(init)
  fd <- dict(Y = Y_dat, L = L_dat, S = data$S)
  LL <- sess$run(L_y, feed_dict = fd)

  for( i in seq_len(max_em_iter) ) {

    # E-step
    g <- sess$run(gamma, feed_dict = fd)

    # M-step
    gfd <- dict(Y = Y_dat, L = L_dat, S = data$S, gamma_fixed =  g)
    Q_old <- sess$run(Q, feed_dict = gfd)
    Q_diff <- rel_adam_tol + 1
    mi <- 0
    while(mi < max_adam_iter && Q_diff > rel_adam_tol) {
      mi <- mi + 1
      sess$run(train, feed_dict = dict(Y = data$Y, L = data$L, S = data$S, gamma_fixed = g))
      Q_new = sess$run(Q, feed_dict = dict(Y = data$Y, L = data$L, S = data$S, gamma_fixed = g))
      Q_diff = -(Q_new - Q_old) / abs(Q_old)

      if(mi %% 20 == 0) {
        print(glue("Gradient descent iteration {mi} % diff {100 * Q_diff}"))
      }
      Q_old <- Q_new
    }

    if(mi == max_adam_iter) {
      warning("Maxmim number of ADAM iterations reached reached")
    }

    L_new <- sess$run(L_y, feed_dict = fd)
    LL_diff <- (L_new - LL)/abs(LL)
    print(glue("{mi}\tL old: {LL}; L new: {L_new}; Difference (%): {LL_diff}"))
    LL <- L_new


    if(LL_diff < rel_em_tol)
      break

  }

  s <- sess$run(list(mu_measured, gamma, beta, phi), feed_dict = fd)
  names(s) <- c("mu", "gamma", "beta", "phi")

  return(s)
}





