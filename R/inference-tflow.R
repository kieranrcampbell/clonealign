

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
                            x = NULL,
                            fix_alpha = FALSE,
                            clone_specific_phi = TRUE,
                            fix_s = NULL,
                            sigma_hyper = 1,
                            dtype = c("float32", "float64"),
                            saturate = TRUE,
                            saturation_threshold = 4,
                            K = 1,
                            B = 10,
                            verbose = TRUE) {

  # Do a first check that we actually have tensorflow support
  if(!reticulate::py_module_available("tensorflow")) {
    msg <- "Tensorflow does not appear to be installed\n"
    msg <- c(msg, "To install run install.pacakges(\"tensorflow\") then tensorflow::install_tensorflow()\n")
    msg <- c(msg, "For more details see the clonealign vignette or https://tensorflow.rstudio.com/tensorflow/articles/installation.html")
    stop(msg)
  }

  # Reset graph
  tf$reset_default_graph()

  # Sort out the dtype
  dtype <- match.arg(dtype)
  dtype <- switch(dtype,
                  float32 = tf$float32,
                  float64 = tf$float64)

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
    L_dat <- saturate(L_dat, saturation_threshold)
  }

  # Get covariate information
  P <- 0
  if(!is.null(x)) {
    if(is.vector(x)) x <- matrix(x, ncol = 1)
    stopifnot(is.matrix(x))
    P <- ncol(x)
    stopifnot(nrow(x) == N)
  }


  data <- list(
    Y = Y_dat,
    L = L_dat,
    N = N,
    G = G,
    x = x,
    C = C
  )

  ## Added 8-8-2018
  B <- as.integer(B)

  # basis_means_fixed <- quantile(y_means, probs = seq(0, 1, length.out = B))
  basis_means_fixed <- seq(from = min(Y_dat), to = max(Y_dat), length.out = B)
  basis_means <- tf$constant(basis_means_fixed, dtype = dtype)

  b_init <- 2 * (basis_means_fixed[2] - basis_means_fixed[1])^2

  if(verbose) {
    message("Creating Tensorflow graph...")
  }


  # Data
  Y <- tf$placeholder(shape = shape(N,G), dtype = dtype)
  L <- tf$placeholder(shape = shape(G,C), dtype = dtype)
  if(P > 0) {
    X <- tf$placeholder(shape = shape(N, P), dtype = dtype)
  }
  LOWER_BOUND <- 1e-10

  # Initializations
  pca <- prcomp(Y_dat, center = TRUE, scale = TRUE)
  pcs <- pca$x[,seq_len(K),drop=FALSE]
  pcs <- scale(pcs)

  s_init = rowSums(data$Y)

  mu_guess <- colMeans(data$Y / rowMeans(data$Y)) / rowMeans(data$L)
  mu_guess <- mu_guess[-1] / mu_guess[1]

  # Spline variables
  a <- tf$exp(tf$Variable(tf$zeros(shape = B, dtype = dtype)))
  # b <- tf$exp(tf$Variable(tf$constant(rep(-log(b_init), B), dtype = dtype)))
  b <- tf$exp(tf$constant(rep(-log(b_init), B), dtype = dtype))

  # Variables ------
  W <- tf$Variable(tf$zeros(shape = c(G, K), dtype = dtype))
  psi <- tf$Variable(tf$constant(pcs, dtype = dtype))
  psi_times_W <- tf$matmul(psi, W, transpose_b = TRUE)
  if(P > 0) {
    beta <- tf$get_variable("beta", initializer = tf$zeros(shape = shape(G,P), dtype = dtype))
    X_times_beta = tf$matmul(X, beta, transpose_b = TRUE)
  }

  mu_log <- tf$Variable(tf$constant(log(mu_guess), dtype = dtype))

  s <- NULL
  if(!is.null(fix_s)) {
    s <- tf$constant(s_init, dtype = dtype)
  } else {
    s <- tf$exp(tf$Variable(tf$constant(log(s_init), dtype = dtype))) + LOWER_BOUND
  }

  # phi <- NULL
  # if(clone_specific_phi) {
  #   phi <- tf$exp(tf$Variable(tf$constant(value = 0, shape = shape(C,G), dtype = dtype))) + LOWER_BOUND
  # } else {
  #   phi <- tf$exp(tf$Variable(tf$zeros(shape(G), dtype = dtype))) + LOWER_BOUND
  # }
  #
  # phi_bar <- tf$Variable(tf$ones(G, dtype = dtype))
  # sigma_log <- tf$constant(log(sigma_hyper), dtype = dtype)


  log_alpha <- NULL

  if(!fix_alpha) {
    alpha_unconstr <- tf$Variable(tf$zeros(C, dtype = dtype))
    log_alpha <- tf$nn$log_softmax(alpha_unconstr)
  } else {
    log_alpha <- tf$constant(rep(-log(C), C), dtype = dtype)
  }


  # Build likelihood
  mu <- tf$concat( list(tf$constant(matrix(1.0), dtype = dtype),
                        tf$reshape(tf$exp(mu_log) + LOWER_BOUND, c(1L, -1L))),
                   axis = 1L)
  mu <- tf$squeeze(mu)

  random_fixed_effects <- NULL
  if(P == 0) {
    random_fixed_effects <- tf$exp(psi_times_W)
  } else {
    random_fixed_effects <- tf$exp(psi_times_W + X_times_beta)
  }

  mu_cg <- tf$transpose(L) * mu
  mu_gcn <- tf$einsum('cg,ng->gcn', mu_cg, random_fixed_effects)
  mu_gc_norm_fctr <- s / tf$reduce_sum(mu_gcn, 0L) # C by N
  mu_gcn_norm <- mu_gcn * mu_gc_norm_fctr
  mu_ncg <- tf$transpose(mu_gcn_norm, perm = c(2L, 1L, 0L))

  mu_ncgb <- tf$tile(tf$expand_dims(mu_ncg, axis = 3L), c(1L, 1L, 1L, B))
  phi <- tf$reduce_sum(a * tf$exp(-b * tf$square(mu_ncgb - basis_means)), 3L) + LOWER_BOUND # n by c by g


  p <- (mu_ncg) / (mu_ncg + phi)
  y_pdf <- tf$contrib$distributions$NegativeBinomial(probs = p, total_count = phi)
  Y_tiled <- tf$stack(rep(list(Y), data$C), axis = 1)
  y_log_prob <- y_pdf$log_prob(Y_tiled)

  p_y_on_c <- tf$reduce_sum(y_log_prob, 2L)
  p_y_on_c <- p_y_on_c + log_alpha
  p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c, 1L), c(1L, -1L))

  gamma <- tf$exp(tf$transpose(tf$transpose(p_y_on_c) - p_y_on_c_norm))

  # Prior on phi
  # phi_pdf <- tf$contrib$distributions$Normal(loc = phi_bar, scale = tf$exp(sigma_log))
  # p_phi <- phi_pdf$log_prob(tf$log(phi))

  # Prior on psi
  psi_pdf <- tf$contrib$distributions$Normal(loc = tf$zeros(1, dtype = dtype), scale = tf$ones(1, dtype = dtype))
  p_psi <- psi_pdf$log_prob(psi)

  # Q function
  gamma_fixed <- tf$placeholder(dtype = dtype, shape = c(N, C))

  y_log_prob_g_sum <- tf$reduce_sum(y_log_prob, 2L) + log_alpha

  Q <- -tf$einsum('nc,nc->', gamma_fixed, y_log_prob_g_sum) - tf$reduce_sum(p_psi)

  optimizer <- tf$train$AdamOptimizer(learning_rate = learning_rate)
  train <- optimizer$minimize(Q)

  eta_y <- y_log_prob_g_sum
  L_y <- tf$reduce_sum(tf$reduce_logsumexp(eta_y, 1L)) # + tf$reduce_sum(p_phi)

  # Inference
  mu_final <- s_final <- phi_final <- NA

  sess <- tf$Session()
  init <- tf$global_variables_initializer()

  sess$run(init)

  initial_mu <- sess$run(mu)

  fd <- NULL
  if(P == 0) {
    fd <- dict(Y = Y_dat, L = L_dat)
  } else {
    fd <- dict(Y = Y_dat, L = L_dat, X = x)
  }

  l <- LL <- sess$run(L_y, feed_dict = fd)

  if(is.na(l)) {
    stop("Initial log likelihood is NA")
  }

  LL_diff <- Inf

  pb <- progress_bar$new(total = max_iter_em+1,
                         format = "  running EM [:bar] :percent | change in log-lik :change")
  pb$tick(0,tokens = list(change = glue("{LL_diff}%")))



  for( i in seq_len(max_iter_em) ) {
    pb$tick(tokens = list(change = glue("{round2(LL_diff)}%")))

    # E-step
    g <- sess$run(gamma, feed_dict = fd)

    # M-step
    if(P == 0) {
      gfd <- dict(Y = Y_dat, L = L_dat, gamma_fixed =  g)
    } else {
      gfd <- dict(Y = Y_dat, L = L_dat, X = x, gamma_fixed = g)
    }

    for(mi in seq_len(max_iter_adam)) {
      sess$run(train, feed_dict = gfd)
    }

    L_new <- sess$run(L_y, feed_dict = fd)
    LL_diff <- (L_new - LL)/abs(LL)
    l <- c(l, L_new)
    LL <- L_new

    # print(l)
    if(abs(LL_diff) < rel_tol_em)
      break

  }

  if(verbose) {
    message("clonealign inference complete")
  }

  rlist <- sess$run(list(mu, gamma, s, phi, tf$exp(log_alpha), a, b, psi, W), feed_dict = fd)
  if(P > 0) {
    rlist$beta <- sess$run(beta, feed_dict = fd)
  }

  # Close the tensorflow session
  sess$close()

  rlist$l <- l
  if(P > 0) {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "psi", "W", "beta", "log_lik")
  } else {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "psi", "W", "log_lik")
  }

  rlist$initial_mu <- initial_mu
  rlist$retained_genes <- retained_genes
  rlist$basis_means <- basis_means_fixed

  return(rlist)
}





