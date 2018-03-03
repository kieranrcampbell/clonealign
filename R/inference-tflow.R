

#' Computes map clone assignment given EM object
#'
#' @param tflow_res List returned by \code{inference_tflow}
#' @return A vector of maximum likelihood clone assignments
#' @keywords internal
clone_assignment <- function(tflow_res) {
  apply(tflow_res$gamma, 1, which.max)
}

#' Round a number to two significant figures
#' @keywords internal
round2 <- function(x) formatC(signif(x,digits=2), digits=2,format="fg", flag="#")

#' EM inference with tensorflow
#' @export
#' @importFrom glue glue
#' @import tensorflow
#' @importFrom progress progress_bar
inference_tflow <- function(Y_dat,
                            L_dat,
                            max_em_iter = 50,
                            rel_em_tol = 1e-5,
                            max_adam_iter = 200,
                            rel_adam_tol = 1e-6,
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
    N = N,
    G = G,
    C = C
  )

  mu_guess <- colMeans(data$Y / rowMeans(data$Y))
  mu_guess <- mu_guess[-1] / mu_guess[1]

  beta_init <- rep(0.5, data$G)

  LOWER_BOUND <- 1e-10

  s_init = rowMeans(data$Y)

  if(verbose) {
    message("Creating Tensorflow graph...")
  }


  # Tensorflow variables
  # Data
  Y <- tf$placeholder(tf$float32, shape = c(N,G))
  L <- tf$placeholder(tf$float32, shape = c(G,C))

  # Unconstrained variables
  mu_log <- tf$Variable(tf$constant(log(mu_guess)))
  s_log <- tf$Variable(tf$constant(s_init))
  phi_log <- tf$Variable(tf$zeros(G))

  # Constrained variables
  mu <- tf$concat( list(tf$constant(matrix(1.0), dtype = tf$float32),
                        tf$reshape(tf$exp(mu_log) + LOWER_BOUND, c(1L, -1L))),
                   axis = 1L)
  mu <- tf$squeeze(mu)

  s <- tf$exp(s_log) + LOWER_BOUND
  phi <- tf$exp(phi_log) + LOWER_BOUND


  mu_gc <- tf$einsum('g,gc->gc', mu, L)

  mu_gc_norm_fctr <- 1 / tf$reduce_sum(mu_gc, 0L)

  mu_gc_norm = tf$einsum('gc,c->gc', mu_gc, mu_gc_norm_fctr)

  mu_ncg <- tf$einsum('n,gc->ncg', s, mu_gc_norm)

  p <- (mu_ncg) / (mu_ncg + phi)

  y_pdf <- tf$contrib$distributions$NegativeBinomial(probs = p, total_count = phi)


  Y_tiled <- tf$stack(rep(list(Y), data$C), axis = 1)

  y_log_prob <- y_pdf$log_prob(Y_tiled)

  p_y_on_c <- tf$reduce_sum(y_log_prob, 2L)

  p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c, 1L), c(1L, -1L))

  gamma <- tf$exp(tf$transpose(tf$transpose(p_y_on_c) - p_y_on_c_norm))

  # Q function
  gamma_fixed <- tf$placeholder(dtype = tf$float32, shape = c(N, C))

  Q_g <- tf$einsum('nc,ncg->g', gamma_fixed, y_log_prob)
  Q <- -(tf$reduce_sum(Q_g))

  optimizer <- tf$train$AdamOptimizer(learning_rate = 0.1)
  train <- optimizer$minimize(Q)

  eta_y <- tf$reduce_sum(y_log_prob, 2L)
  L_y <- tf$reduce_sum(tf$reduce_logsumexp(eta_y, 1L))

  # Inference
  mu_final <- s_final <- phi_final <- NA
  l <- NA

  sess <- tf$Session()
  init <- tf$global_variables_initializer()

  sess$run(init)
  fd <- dict(Y = Y_dat, L = L_dat)
  LL <- sess$run(L_y, feed_dict = fd)

  LL_diff <- Inf

  pb <- progress_bar$new(total = max_em_iter,
                         format = "  running EM [:bar] :percent | change in log-lik :change")
  pb$tick(0,tokens = list(change = glue("{LL_diff}%")))

  for( i in seq_len(max_em_iter) ) {

    # E-step
    g <- sess$run(gamma, feed_dict = fd)

    # M-step
    gfd <- dict(Y = Y_dat, L = L_dat, gamma_fixed =  g)
    l <- Q_old <- sess$run(Q, feed_dict = gfd)
    Q_diff <- Inf
    pb$tick(tokens = list(change = glue("{round2(LL_diff)}%")))

    mi <- 0

    while(mi < max_adam_iter && Q_diff > rel_adam_tol) {
      mi <- mi + 1
      sess$run(train, feed_dict = gfd)
      Q_new = sess$run(Q, feed_dict = dict(Y = data$Y, L = data$L, gamma_fixed = g))
      Q_diff = -(Q_new - Q_old) / abs(Q_old)

      Q_old <- Q_new
    }

    if(mi == max_adam_iter) {
      warning("Maximim number of ADAM iterations reached reached")
    }

    L_new <- sess$run(L_y, feed_dict = fd)
    LL_diff <- (L_new - LL)/abs(LL)
    # print(glue("EM iteration {i}\t New log-likelihood: {L_new}; Difference (%): {LL_diff}"))
    l <- c(l, LL)
    LL <- L_new


    if(LL_diff < rel_em_tol)
      break

  }

  rlist <- sess$run(list(mu, gamma, s, phi, l), feed_dict = fd)
  names(rlist) <- c("mu", "gamma", "s", "phi", "log_lik")

  return(rlist)
}





