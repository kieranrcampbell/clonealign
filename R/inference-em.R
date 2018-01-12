

#' Computes map clone assignment given EM object
#'
#' @param tflow_res List returned by \code{inference_tflow}
#' @return A vector of maximum likelihood clone assignments
#' @keywords internal
clone_assignment <- function(tflow_res) {
  apply(tflow_res$gamma, 1, which.max)
}



#' EM inference with tensorflow
#' @export
#' @importFrom glue glue
#' @import tensorflow
inference_tflow <- function(Y_dat,
                            L_dat,
                            max_em_iter = 100,
                            rel_em_tol = 1e-5,
                            max_adam_iter = 500,
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
    S = scran::computeSumFactors(t(Y_dat)),
    N = N,
    G = G,
    C = C
  )

  mu_guess <- colMeans(data$Y / data$S)

  beta_init <- rep(0.5, data$G)

  LOWER_BOUND <- 1e-10

  # Tensorflow variables
  # Data
  Y <- tf$placeholder(tf$float32, shape = c(N,G))
  L <- tf$placeholder(tf$float32, shape = c(G,C))
  S <- tf$placeholder(tf$float32, shape = N)

  # Unconstrained variables
  mu_log <- tf$Variable(tf$constant(log(mu_guess)))
  beta_log <- tf$Variable(tf$constant(log(beta_init)))
  phi_log <- tf$Variable(tf$zeros(G))

  # Constrained variables
  mu <- tf$exp(mu_log) + LOWER_BOUND
  beta <- tf$exp(beta_log) + LOWER_BOUND
  phi <- tf$exp(phi_log) + LOWER_BOUND


  mu_gc <- tf$transpose(mu * (1 + beta * tf$transpose(L)))
  #mu_c <- tf$reduce_sum(mu_gc, 0L)
  #mu_gc = mu_gc / mu_c

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

  Q_g <- tf$einsum('nc,ncg->g', gamma_fixed, y_log_prob)
  Q <- -(tf$reduce_sum(Q_y))

  optimizer <- tf$train$AdamOptimizer()
  trainers <- lapply(seq_len(data$G), function(g) {
    g <- as.integer(g)
    mu_g <- tf$slice(tf$reshape(mu_log, c(G, -1L)), c(g, 1L), c(g, 1L))
    beta_g <- tf$slice(tf$reshape(beta_log, c(G, -1L)), c(g, 1L), c(g, 1L))
    phi_g <- tf$slice(tf$reshape(phi_log, c(G, -1L)), c(g, 1L), c(g, 1L))
    optimizer$minimize(Q_g[g], var_list = list(mu_g, beta_g, phi_g))
  })

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
      for(gg in seq_len(g)) {
        print(gg)
        sess$run(trainers[[gg]], feed_dict = dict(Y = data$Y, L = data$L, S = data$S, gamma_fixed = g))
      }
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

  s <- sess$run(list(mu, gamma, beta, phi), feed_dict = fd)
  names(s) <- c("mu", "gamma", "beta", "phi")

  return(s)
}





