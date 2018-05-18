

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
#' @return  The rounded number
round2 <- function(x) formatC(signif(x,digits=2), digits=2,format="fg", flag="#")

#' EM inference with tensorflow
#'
#' @param Y_dat Cell by gene matrix of counts
#' @param L_dat Gene by clone matrix of copy number
#'
#' @importFrom glue glue
#' @import tensorflow
#' @importFrom progress progress_bar
#' @importFrom stats prcomp
#'
#' @keywords internal
#'
#' @return A list of maximum likelihood parameter estimates and log likelihood values
inference_tflow <- function(Y_dat,
                            L_dat,
                            max_iter_em = 50,
                            max_iter_adam = 200,
                            rel_tol_em = 1e-5,
                            rel_tol_adam = 1e-6,
                            learning_rate = 1e-3,
                            gene_filter_threshold = 0,
                            fix_alpha = FALSE,
                            clone_specific_phi = TRUE,
                            fix_s = NULL,
                            sigma_hyper = 1,
                            saturate = TRUE,
                            verbose = TRUE) {

  # Do a first check that we actually have tensorflow support
  if(!reticulate::py_module_available("tensorflow")) {
    msg <- "Tensorflow does not appear to be installed\n"
    msg <- c(msg, "To install run install.pacakges(\"tensorflow\") then tensorflow::install_tensorflow()\n")
    msg <- c(msg, "For more details see the clonealign vignette or https://tensorflow.rstudio.com/tensorflow/articles/installation.html")
    stop(msg)
  }

  zero_gene_means <- colSums(Y_dat) <= gene_filter_threshold

  if(verbose) {
    message(glue("Removing {sum(zero_gene_means)} genes with low counts"))
  }

  Y_dat <- Y_dat[, !zero_gene_means]
  L_dat <- L_dat[!zero_gene_means,]

  retained_genes <- NULL
  if(!is.null(colnames(Y_dat))) {
    retained_genes <- colnames(Y_dat)
  } else {
    retained_genes <- which(!zero_gene_means)
  }

  N <- nrow(Y_dat)
  G <- ncol(Y_dat)
  C <- ncol(L_dat)

  # Sanity checks
  stopifnot(nrow(L_dat) == G)

  # Saturate
  if(saturate) {
    L_dat <- saturate(L_dat)
  }

  data <- list(
    Y = Y_dat,
    L = L_dat,
    N = N,
    G = G,
    C = C
  )

  #
  pca <- prcomp(Y_dat, center = TRUE, scale = TRUE)
  pc1 <- pca$x[,1]
  pc1 <- (pc1 - mean(pc1)) / sd(pc1)

  W <- tf$Variable(tf$zeros(shape = c(1, G), dtype = tf$float64))
  psi <- tf$Variable(tf$reshape(tf$constant(pc1, dtype = tf$float64), shape(-1,1)))
  psi_times_W <- tf$matmul(psi,W)


  mu_guess <- colMeans(data$Y / rowMeans(data$Y)) / rowMeans(data$L)
  mu_guess <- mu_guess[-1] / mu_guess[1]

  LOWER_BOUND <- 1e-10

  s_init = rowSums(data$Y)

  if(verbose) {
    message("Creating Tensorflow graph...")
  }


  # Tensorflow variables
  # Data
  Y <- tf$placeholder(shape = c(N,G), dtype = tf$float64)
  L <- tf$placeholder(shape = c(G,C), dtype = tf$float64)

  # Unconstrained variables
  mu_log <- tf$Variable(tf$constant(log(mu_guess), dtype = tf$float64))

  s <- NULL
  if(!is.null(fix_s)) {
    s <- tf$constant(s_init, dtype = tf$float64)
  } else {
    s_log <- tf$Variable(tf$constant(log(s_init), dtype = tf$float64))
    s <- tf$exp(s_log) + LOWER_BOUND
  }

  phi_log <- NULL
  if(clone_specific_phi) {
    phi_log <- tf$Variable(tf$constant(value = 0, shape = shape(C,G), dtype = tf$float64))
  } else {
    phi_log <- tf$Variable(tf$zeros(shape(G), dtype = tf$float64))
  }

  # Variance shrinkage variables
  phi_bar <- tf$Variable(tf$ones(G, dtype = tf$float64))
  # sigma_log <- tf$Variable(tf$ones(1, dtype = tf$float64))
  # TODO CHANGE ME
  sigma_log <- tf$constant(log(sigma_hyper), dtype = tf$float64)


  log_alpha <- NULL

  if(!fix_alpha) {
    # alpha_unconstr_0 <- tf$Variable(tf$zeros(C-1))
    # alpha_unconstr <- tf$concat(list(alpha_unconstr_0, tf$constant(0, dtype = tf$float32, shape = shape(1))), axis = 0L)
    alpha_unconstr <- tf$Variable(tf$zeros(C, dtype = tf$float64))
    log_alpha <- tf$nn$log_softmax(alpha_unconstr)
  } else {
    log_alpha <- tf$constant(rep(-log(C), C), dtype = tf$float64)
  }


  # Constrained variables
  mu <- tf$concat( list(tf$constant(matrix(1.0), dtype = tf$float64),
                        tf$reshape(tf$exp(mu_log) + LOWER_BOUND, c(1L, -1L))),
                   axis = 1L)
  mu <- tf$squeeze(mu)

  phi <- tf$exp(phi_log) + LOWER_BOUND


  mu_gc <- tf$einsum('g,gc->gc', mu, L)

  mu_gcn <- tf$einsum('gc,ng->gcn', mu_gc, exp(psi_times_W))



  mu_gc_norm_fctr <- tf$constant(1, dtype = tf$float64) / tf$reduce_sum(mu_gcn, 0L)


  mu_gc_norm = tf$einsum('gcn,cn->gcn', mu_gcn, mu_gc_norm_fctr)

  mu_ncg <- tf$einsum('n,gcn->ncg', s, mu_gc_norm)

  p <- (mu_ncg) / (mu_ncg + phi)

  y_pdf <- tf$contrib$distributions$NegativeBinomial(probs = p, total_count = phi)


  Y_tiled <- tf$stack(rep(list(Y), data$C), axis = 1)

  y_log_prob <- y_pdf$log_prob(Y_tiled)

  p_y_on_c <- tf$reduce_sum(y_log_prob, 2L)
  p_y_on_c <- p_y_on_c + log_alpha


  p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c, 1L), c(1L, -1L))

  gamma <- tf$exp(tf$transpose(tf$transpose(p_y_on_c) - p_y_on_c_norm)) #+ tf$constant(1e-10)

  # Prior on phi
  phi_pdf <- tf$contrib$distributions$Normal(loc = phi_bar, scale = tf$exp(sigma_log))
  p_phi <- phi_pdf$log_prob(tf$log(phi))

  # Prior on psi
  psi_pdf <- tf$contrib$distributions$Normal(loc = tf$zeros(1, dtype = tf$float64), scale = tf$ones(1, dtype = tf$float64))
  p_psi <- psi_pdf$log_prob(psi)

  # Q function
  gamma_fixed <- tf$placeholder(dtype = tf$float64, shape = c(N, C))

  y_log_prob_g_sum <- tf$reduce_sum(y_log_prob, 2L) + log_alpha

  Q <- -tf$einsum('nc,nc->', gamma_fixed, y_log_prob_g_sum) - tf$reduce_sum(p_phi) - tf$reduce_sum(p_psi)

  optimizer <- tf$train$AdamOptimizer(learning_rate = learning_rate)
  train <- optimizer$minimize(Q)

  eta_y <- y_log_prob_g_sum#tf$reduce_sum(y_log_prob, 2L)
  L_y <- tf$reduce_sum(tf$reduce_logsumexp(eta_y, 1L)) + tf$reduce_sum(p_phi)

  # Inference
  mu_final <- s_final <- phi_final <- NA

  sess <- tf$Session()
  init <- tf$global_variables_initializer()

  sess$run(init)

  initial_mu <- sess$run(mu)

  fd <- dict(Y = Y_dat, L = L_dat)
  l <- LL <- sess$run(L_y, feed_dict = fd)

  LL_diff <- Inf

  pb <- progress_bar$new(total = max_iter_em,
                         format = "  running EM [:bar] :percent | change in log-lik :change")
  pb$tick(0,tokens = list(change = glue("{LL_diff}%")))



  for( i in seq_len(max_iter_em) ) {

    # E-step
    g <- sess$run(gamma, feed_dict = fd)

    # M-step
    gfd <- dict(Y = Y_dat, L = L_dat, gamma_fixed =  g)
    Q_old <- sess$run(Q, feed_dict = gfd)
    Q_diff <- Inf
    pb$tick(tokens = list(change = glue("{round2(LL_diff)}%")))

    mi <- 0

    while(mi < max_iter_adam && Q_diff > rel_tol_adam) {
      mi <- mi + 1
      sess$run(train, feed_dict = gfd)
      Q_new = sess$run(Q, feed_dict = dict(Y = data$Y, L = data$L, gamma_fixed = g))
      Q_diff = -(Q_new - Q_old) / abs(Q_old)

      # print(Q_new)
      # print(Q_diff)

      Q_old <- Q_new
    }

    if(mi == max_iter_adam) {
      warning("Maximim number of ADAM iterations reached reached")
    }

    L_new <- sess$run(L_y, feed_dict = fd)
    LL_diff <- (L_new - LL)/abs(LL)
    # print(glue("EM iteration {i}\t New log-likelihood: {L_new}; Difference (%): {LL_diff}"))
    l <- c(l, L_new)
    LL <- L_new


    if(abs(LL_diff) < rel_tol_em)
      break

  }

  pb$tick(100, tokens = list(change = glue("{round2(LL_diff)}%")))

  if(verbose) {
    message("clonealign inference complete")
  }

  rlist <- sess$run(list(mu, gamma, s, phi, tf$exp(log_alpha), phi_bar, psi, W), feed_dict = fd)


  probs_eval_init_mu <- sess$run(gamma,
                                 feed_dict = dict(Y = Y_dat, L = L_dat, mu = initial_mu))

  # Close the tensorflow session
  sess$close()

  rlist$l <- l

  names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "phi_bar", "psi", "W", "log_lik")

  rlist$initial_mu <- initial_mu
  rlist$retained_genes <- retained_genes
  # rlist$probs_eval_init_mu <- probs_eval_init_mu

  return(rlist)
}





