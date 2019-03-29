
inverse_softplus <- function(x) {
  exp(log(x) - 1)
}

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
                            max_iter = 100,
                            rel_tol = 1e-5,
                            learning_rate = 0.1,
                            gene_filter_threshold = 0,
                            x = NULL,
                            clone_allele = NULL,
                            cov = NULL,
                            ref = NULL,
                            fix_alpha = FALSE,
                            size_factors = "fixed",
                            dtype = c("float32", "float64"),
                            saturate = TRUE,
                            saturation_threshold = 6,
                            K = 6,
                            B = 20,
                            mc_samples = 1,
                            verbose = TRUE,
                            seed = NULL,
                            data_init_mu = data_init_mu) {

  # Do a first check that we actually have tensorflow support
  if(!reticulate::py_module_available("tensorflow")) {
    msg <- "Tensorflow does not appear to be installed\n"
    msg <- c(msg, "To install run install.pacakges(\"tensorflow\") then tensorflow::install_tensorflow()\n")
    msg <- c(msg, "For more details see the clonealign vignette or https://tensorflow.rstudio.com/tensorflow/articles/installation.html")
    stop(msg)
  }





  # Get distributions
  tfp <- reticulate::import("tensorflow_probability")
  tfd <- tfp$distributions
  tfb <- tfp$bijectors

  # Reset graph
  tf$reset_default_graph()

  if(!is.null(seed)) {
    use_session_with_seed(seed)
    ss <- tfd$SeedStream(seed, salt = 'qmu')
  }

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


  if(verbose) {
    message("Creating Tensorflow graph...")
  }

  # Allelic imbalance setup ----------
  use_allele <- !is.null(clone_allele) && !is.null(ref) && !is.null(cov)
  v_log_prob <- NULL
  if(verbose && use_allele) {
    message("Using allelic imbalance info")
  }

  if(use_allele) { # Use allele-specific info
    V <- nrow(clone_allele) # Number of variants
    sanitize_allele_info(V, clone_allele, cov, ref, N, C)

    cov <- t(cov) # Legacy
    ref <- t(ref)

    alt <- cov - ref

    clone_allele_ph <- tf$placeholder(shape = shape(V,C), dtype = dtype)
    alt_ph <- tf$placeholder(shape = shape(V,N), dtype = dtype)
    cov_ph <- tf$placeholder(shape = shape(V,N), dtype = dtype)

    v_log_prob <- construct_ai_likelihood(clone_allele_ph, alt_ph, cov_ph, dtype, C, N)
  }

  # Data specification ----------
  Y <- tf$placeholder(shape = shape(N,G), dtype = dtype)
  L <- tf$placeholder(shape = shape(G,C), dtype = dtype)

  if(P > 0) {
    X <- tf$placeholder(shape = shape(N, P), dtype = dtype)
  }

  # Global lower bound for positive variables ----------
  LOWER_BOUND <- 1e-10

  # Initialization for ----------
  # (i) latent space (PCA)
  # (ii) size factors
  # (iii) mean parameters
  pca <- prcomp(log2(Y_dat+1), center = TRUE, scale = TRUE)
  pcs <- pca$x[,seq_len(K),drop=FALSE]
  pcs <- scale(pcs)

  s_init <- rowSums(data$Y)

  mu_guess <- colMeans(data$Y / rowMeans(data$Y)) / rowMeans(data$L)
  mu_guess <- mu_guess[-1] / mu_guess[1]

  if(!data_init_mu) {
    mu_guess <- rep(1, length(mu_guess))
  }


  # Tensors for dispersion basis functions ----------
  B <- as.integer(B)
  basis_means_fixed <- seq(from = min(Y_dat), to = max(Y_dat+1)/3, length.out = B)
  basis_means <- tf$constant(basis_means_fixed, dtype = dtype)
  b_init <- 2 * (basis_means_fixed[2] - basis_means_fixed[1])^2
  a <- tf$exp(tf$Variable(tf$constant(rep(log(1), B), dtype = dtype)))
  b <- tf$exp(tf$constant(rep(-log(b_init), B), dtype = dtype))

  # Variables to be optimized ----------
  W <- tf$Variable(tf$zeros(shape = c(G, K), dtype = dtype))
  chi <- tf$exp(tf$get_variable("chi", initializer = tf$zeros(K, dtype = dtype))) # Prior variance on W
  psi <- tf$Variable(tf$constant(pcs, dtype = dtype))
  psi_times_W <- tf$matmul(psi, W, transpose_b = TRUE)
  if(P > 0) {
    beta <- tf$get_variable("beta", initializer = tf$zeros(shape = shape(G,P), dtype = dtype))
    X_times_beta = tf$matmul(X, beta, transpose_b = TRUE)
  }

  s <- NULL

  # Sort out size factors
  if(length(size_factors) == 1 && N > 1) {
    if(!is.character(size_factors)) {
      stop("If size factors does not contain sizes per cell it must be either 'fixed' or 'infer'. See ?clonealign")
    }
    if(!(size_factors %in% c("fixed", "infer"))) {
      stop("Size factors must be either 'fixed' or 'infer'. See ?clonealign")
    }

    if(size_factors == "fixed") {
      s <- tf$constant(s_init, dtype = dtype)
    } else {
      s <- tf$exp(tf$Variable(tf$constant(log(s_init), dtype = dtype))) + LOWER_BOUND
    }

  } else if(length(size_factors) > 1) {
    if(length(size_factors) != N || !is.numeric(size_factors) || !all(size_factors > 0)) {
      stop("If size factors is not specified as 'fixed' or 'infer' a positive numeric vector of length N_cells must be supplied. See ?clonealign")
    }
    s <- tf$constant(size_factors, dtype = dtype)
  }


  log_alpha <- NULL

  if(!fix_alpha) {
    alpha_unconstr <- tf$Variable(tf$zeros(C, dtype = dtype))
    log_alpha <- tf$nn$log_softmax(alpha_unconstr)
  } else {
    log_alpha <- tf$constant(rep(-log(C), C), dtype = dtype)
  }

  # Variational variables for mu and z ----------
  sdinit <- 1e-1
  qmu <- tfd$TransformedDistribution(
    bijector = tfb$Softplus(),
    distribution = tfd$Normal(loc = tf$Variable(tf$constant(inverse_softplus(mu_guess),dtype=dtype)),
                              scale = tf$exp(tf$Variable(tf$constant(rep(log(sdinit), G-1), dtype = dtype)))),

    name = "qmu"
  )

  S <- as.integer(mc_samples)
  if(!is.null(seed)) {
    mmu_samples <- qmu$sample(S, seed = ss())
  } else {
    mmu_samples <- qmu$sample(S)
  }
  mu_one_sample <- tf$ones(shape = shape(S, 1L), dtype = dtype)
  mu_samples <- tf$concat(list(mu_one_sample, mmu_samples), axis = 1L)

  gamma_logits <- tf$Variable(tf$zeros(shape(N,C) , dtype = dtype))
  gamma <- tf$nn$softmax(gamma_logits)


  # Build likelihood model  ----------

  random_fixed_effects <- NULL
  if(P == 0 && K > 0) {
    random_fixed_effects <- tf$exp(psi_times_W)
  } else if(P > 0 && K > 0) {
    random_fixed_effects <- tf$exp(psi_times_W + X_times_beta)
  } else {
    random_fixed_effects <- tf$ones(shape = shape(N,G), dtype = dtype)
  }


  mu_scg <- tf$einsum('sg,gc->scg', mu_samples, L)
  mu_sgcn <- tf$einsum('scg,ng->sgcn', mu_scg, random_fixed_effects)
  mu_scn_norm_fctr <- s * tf$ones(1, dtype=dtype) / ( tf$reduce_sum(mu_sgcn, 1L) )  # REMOVE ME!!! # S by C by N
  mu_sgcn_norm <- tf$einsum('sgcn,scn->sgcn', mu_sgcn, mu_scn_norm_fctr)
  mu_sncg <- tf$transpose(mu_sgcn_norm, perm = c(0L, 3L, 2L, 1L))
  mu_sncgb <- tf$tile(tf$expand_dims(mu_sncg, axis = 4L), c(1L, 1L, 1L, 1L, B))
  phi_samples <- tf$reduce_sum(a * tf$exp(-b * tf$square(mu_sncgb - basis_means)), 4L) + LOWER_BOUND # s by n by c by g

  p <- (mu_sncg) / (mu_sncg + phi_samples)
  y_pdf <- tfd$NegativeBinomial(probs = tf$transpose(p, c(0L, 2L, 1L, 3L)), # S C N G
                                total_count = tf$transpose(phi_samples, c(0L, 2L, 1L, 3L)))

  y_log_prob <- y_pdf$log_prob(Y)

  # Build ELBO ----------

  # (i) E_q[log p(y | z, theta)]
  p_y_on_c <- tf$reduce_sum(y_log_prob, 3L) # Reduce over genes

  if(use_allele) {
    p_y_on_c <- p_y_on_c + tf$transpose(v_log_prob)
  }

  E_p_y_on_c <- tf$reduce_mean(p_y_on_c, 0L) # Reduce over MC samples

  EE_p_y <- tf$reduce_sum(gamma * tf$transpose(E_p_y_on_c)) # E_q[p(y | theta)] (reduce over cells and clones)
  EE_p_y <- tf$reduce_mean(tf$einsum('nc,scn->s', gamma, p_y_on_c))



  # Prior on w, psi and chi
  if(K > 0) {
    W_log_prob <- tf$reduce_sum(tfd$Normal(loc = tf$zeros(1, dtype = dtype),
                                           scale = tf$sqrt(tf$ones(1, dtype = dtype) / chi))$log_prob(W))

    chi_log_prob <- tf$reduce_sum(tfd$Gamma(concentration = tf$constant(2, dtype = dtype),
                                              rate = tf$ones(1, dtype = dtype))$log_prob(chi))

    psi_pdf <- tfd$Normal(loc = tf$zeros(1, dtype = dtype), scale = tf$ones(1, dtype = dtype))
    p_psi <- psi_pdf$log_prob(psi)
  }
  # (ii) E_q[log p(theta)]
  E_log_p_p <- tf$reduce_sum(log_alpha * gamma) +
    tf$reduce_sum(tfd$Normal(loc = tf$zeros(1, dtype = dtype), scale = tf$ones(1, dtype = dtype))$log_prob(tf$log(mmu_samples))) / tf$constant(S, dtype = dtype) +
    tf$reduce_sum(tfd$Dirichlet(tf$constant(rep(1, C), dtype = dtype))$log_prob(tf$exp(log_alpha)))

  if(K>0) {
    E_log_p_p <- E_log_p_p + W_log_prob + chi_log_prob + tf$reduce_sum(p_psi)
  }


  # (iii) E_q[log q]
  E_log_q <- tf$reduce_sum(tf$reduce_mean(qmu$log_prob(mmu_samples), 0L)) +
    # tf$reduce_sum(tf$reduce_mean(qs$log_prob(s_samples), 0L)) +
    tf$reduce_sum(gamma * tf$log(gamma))


  elbo <- EE_p_y + E_log_p_p - E_log_q

  # Extra tensors required to initialize gamma
  p_y_on_c_ <- p_y_on_c
  p_y_on_c_norm <- tf$reshape(tf$reduce_logsumexp(p_y_on_c, 1L), c(S, 1L, -1L))
  gamma_init <- tf$transpose(tf$reduce_mean(tf$exp(p_y_on_c - p_y_on_c_norm), 0L))
  gamma_init_ph <- tf$placeholder(shape = shape(N,C), dtype=dtype)
  init_gamma <- tf$assign(gamma_logits, tf$log(gamma_init_ph) )

  # Optimizer for ELBO ----------
  optimizer <- tf$train$AdamOptimizer(learning_rate = learning_rate)
  train <- optimizer$minimize(-elbo)

  # Inference ----------
  mu_final <- s_final <- phi_final <- NA

  sess <- tf$Session()
  init <- tf$global_variables_initializer()
  sess$run(init)

  fd <- dict(Y = Y_dat, L = L_dat)

  if(P > 0) {
    fd$update(dict(X = x))
  }

  if(use_allele) {
    fd$update(dict(clone_allele_ph = clone_allele))
    fd$update(dict(alt_ph = alt))
    fd$update(dict(cov_ph = cov))
  }

  # Initialize gamma correctly
  gi <- sess$run(gamma_init, feed_dict = fd)
  sess$run(init_gamma, feed_dict = dict(gamma_init_ph = gi))


  elbo_val <- sess$run(elbo, feed_dict = fd)

  if(is.na(elbo_val)) {
    stop("Initial elbo is NA")
  }

  elbo_diff <- Inf
  elbos <- NULL

  pb <- progress_bar$new(total = max_iter+1,
                         format = "  running VB [:bar] :percent | change in elbo :change")

  pb$tick(0,tokens = list(change = glue("{elbo_diff}%")))



  for( i in seq_len(max_iter) ) {
    pb$tick(tokens = list(change = glue("{round2(elbo_diff)}%")))

    sess$run(train, feed_dict = fd)

    elbo_new <- sess$run(elbo, feed_dict = fd)
    elbo_diff <- (elbo_new - elbo_val)/abs(elbo_val)
    elbos <- c(elbos, elbo_new)
    elbo_val <- elbo_new

    if(abs(elbo_diff) < rel_tol)
      break

  }

  if(verbose) {
    message("clonealign inference complete")
  }

  rlist <- sess$run(list( tf$reduce_mean(qmu$sample(50L), 0L), gamma, s, tf$reduce_mean(phi_samples, 0L), tf$exp(log_alpha), a, b), feed_dict = fd)
  if(P > 0) {
    rlist$beta <- sess$run(beta, feed_dict = fd)
  }

  if(K > 0) {
    k_params <- sess$run(list(psi, W, chi), feed_dict = fd)
    rlist$psi <- k_params[[1]]
    rlist$W <- k_params[[2]]
    rlist$chi <- k_params[[3]]
  }

  clone_probs_from_snv <- NULL
  if(use_allele) {
    clone_probs_from_snv <- tf$exp(tf$transpose(v_log_prob) - tf$reduce_logsumexp(v_log_prob, 1L))
    clone_probs_from_snv <- t(sess$run(clone_probs_from_snv, feed_dict = fd))
  }

  # Close the session ----------
  sess$close()

  rlist$elbos <- elbos

  # Correctly name the return vector
  if(P > 0 && K > 0) {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "beta", "psi", "W", "chi", "elbo")
  } else if (K > 0 && P == 0) {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "psi", "W", "chi", "elbo")
  } else if (K == 0 && P == 0) {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "elbo")
  } else {
    names(rlist) <- c("mu", "gamma", "s", "phi", "alpha", "a", "b", "elbo")
  }

  rlist$mu <- c(1, rlist$mu) # Add on the constant at the end

  rlist$retained_genes <- retained_genes
  rlist$basis_means <- basis_means_fixed
  rlist$clone_probs_from_snv <- clone_probs_from_snv

  return(rlist)
}
