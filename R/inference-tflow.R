
inverse_softplus <- function(x) {
  log(exp(x) - 1)
}

safe_inverse_softplus <- function(x) {
  if(any(x < 0)) {
    stop("Inverse softplus only takes positive values")
  }
  log(1 - exp(-abs(x))) + pmax(x, 0)
}

softplus <- function(x) {
  log(1 + exp(x))
}

#' Computes map clone assignment given EM object
#'
#' @param tflow_res List returned by \code{inference_tflow}
#' @return A vector of maximum likelihood clone assignments
#' @keywords internal
clone_assignment <- function(gamma, clone_names, clone_assignment_probability = 0.95) {
  apply(gamma, 1, function(r) {
    if(max(r) < clone_assignment_probability) {
      return("unassigned")
    }
    return(clone_names[which.max(r)])
  })
}

#' Re-computes map clone assignment at a given threshold
#'
#' @param ca A clonealign fit
#' @return A clonealign fit with clones re-called at different probabilities
#' @keywords internal
recompute_clone_assignment <- function(ca, clone_assignment_probability = 0.95) {
  clone_names <- colnames(ca$ml_params$clone_probs)
  clones <- apply(ca$ml_params$clone_probs, 1, function(r) {
    if(max(r) < clone_assignment_probability) {
      return("unassigned")
    }
    return(clone_names[which.max(r)])
  })
  ca$clone <- clones
  ca
}

#' @keywords internal
get_next_seed <- function() {
  sample(.Machine$integer.max - 1, 1)
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
                            dtype = c("float32", "float64"),
                            saturate = TRUE,
                            saturation_threshold = 6,
                            K = 1,
                            mc_samples = 1,
                            verbose = TRUE,
                            initial_shrink = 5,
                            data_init_mu = data_init_mu) {


  tf <- tf$compat$v1
  tf$disable_v2_behavior()
  
  tfp <- reticulate::import('tensorflow_probability')
  tfd <- tfp$distributions
  tfb <- tfp$bijectors
  
  tf$reset_default_graph()

  
  if(verbose) {
    message("Constructing tensorflow graph")
  }

  # if(!is.null(seed)) {
  #   use_session_with_seed(seed)
  #   ss <- tfd$SeedStream(seed, salt = 'qmu')
  # }

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
  K <- as.integer(K)

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
  
  pcs <- pcs + matrix(rnorm(nrow(pcs) * ncol(pcs), mean=0, sd = .05), nrow = nrow(pcs))

  s_init <- rowSums(data$Y)
  
  if(any(s_init == 0)) {
    stop("Some cells have no counts mapping")
  }


  
  beta_init <- matrix(0, nrow = G, ncol = P)

  if(is.logical(data_init_mu)) {
    if(data_init_mu) {
      mu_guess <- colMeans(data$Y / rowMeans(data$Y)) # old: colMeans(data$Y / rowMeans(data$Y)) / rowMeans(data$L)
    } else {
      mu_guess <- rep(1, data$G)
    }
  } else {
    if(is.numeric(data_init_mu)) {
      if(length(data_init_mu == data$G)) {
        if(verbose) {
          message("Using user-provided mu values to start")
        }
        mu_guess <- data_init_mu / mean(data_init_mu)
      }
    }
  }


  # Variables to be optimized ----------
  # browser()
  W <- tf$Variable(tf$zeros(shape = c(G, K), dtype = dtype))
  chi <- tf$exp(tf$get_variable("chi", initializer = tf$zeros(shape = K, dtype = dtype))) # Prior variance on W
  psi <- tf$Variable(tf$constant(pcs, dtype = dtype))
  psi_times_W <- tf$matmul(psi, W, transpose_b = TRUE)
  if(P > 0) {
    beta <- tf$get_variable("beta", initializer = tf$constant(beta_init, dtype = dtype))
    X_times_beta = tf$matmul(X, beta, transpose_b = TRUE)
  }


  s <- tf$constant(s_init, dtype = dtype)

  log_alpha <- NULL

  alpha_unconstr <- tf$Variable(tf$zeros(C, dtype = dtype))
  log_alpha <- tf$nn$log_softmax(alpha_unconstr)

  # Variational variables for mu and z ----------
  sdinit <- 1
  
  qmu <- tfd$TransformedDistribution(
    bijector = tfb$Softplus(),
    distribution = tfd$Normal(loc = tf$Variable(tf$constant(safe_inverse_softplus(mu_guess),dtype=dtype)),
                              scale = tf$exp(tf$Variable(tf$constant(rep(log(sdinit), G), dtype = dtype)))),

    name = "qmu"
  )

  S <- as.integer(mc_samples) 
  mu_samples <- qmu$sample(S, seed = get_next_seed())


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
  mu_scn_norm_fctr <- tf$ones(shape(1), dtype=dtype) / ( tf$reduce_sum(mu_sgcn, 1L) ) # 
  mu_sgcn_norm <- tf$einsum('sgcn,scn->sgcn', mu_sgcn, mu_scn_norm_fctr)
  mu_scng <- tf$transpose(mu_sgcn_norm, perm = c(0L, 2L, 3L, 1L))
  
  y_pdf <- tfd$Multinomial(total_count = s_init, probs = mu_scng)
  
  p_y_on_c <- y_pdf$log_prob(Y)

  # Build ELBO ----------

  # (i) E_q[log p(y | z, theta)]

  if(use_allele) {
    p_y_on_c <- p_y_on_c + tf$transpose(v_log_prob)
  }

  E_p_y_on_c <- tf$reduce_mean(p_y_on_c, 0L) # Reduce over MC samples

  EE_p_y <- tf$reduce_sum(gamma * tf$transpose(E_p_y_on_c)) # E_q[p(y | theta)] (reduce over cells and clones)

  # Prior on w, psi and chi
  if(K > 0) {
    W_log_prob <- tf$reduce_sum(tfd$Normal(loc = tf$zeros(shape(1), dtype = dtype),
                                           scale = tf$sqrt(tf$ones(shape(1), dtype = dtype) / chi))$log_prob(W))

    chi_log_prob <- tf$reduce_sum(tfd$Gamma(concentration = tf$constant(2, dtype = dtype),
                                              rate = tf$ones(shape(1), dtype = dtype))$log_prob(chi))

    psi_pdf <- tfd$Normal(loc = tf$zeros(shape(1), dtype = dtype), scale = tf$ones(shape(1), dtype = dtype))
    p_psi <- psi_pdf$log_prob(psi)
  }
  # (ii) E_q[log p(theta)]
  E_log_p_p <- tf$reduce_sum(log_alpha * gamma) +
    tf$reduce_sum(tfd$Normal(loc = tf$zeros(shape(1), dtype = dtype), scale = tf$ones(shape(1), dtype = dtype))$log_prob(tf$log(mu_samples))) / tf$to_float(S) +
    tf$reduce_sum(tfd$Dirichlet(tf$constant(rep(1/C, C), dtype = dtype))$log_prob(tf$exp(log_alpha) + tf$constant(1e-3, dtype = dtype)))

  if(K > 0) {
    E_log_p_p <- E_log_p_p + W_log_prob + chi_log_prob + tf$reduce_sum(p_psi)
  }


  # (iii) E_q[log q]
  E_log_q <- tf$reduce_sum(tf$reduce_mean(qmu$log_prob(mu_samples), 0L)) +
    tf$reduce_sum(tf$where(gamma == 0, tf$zeros(shape = gamma$shape, dtype = dtype), gamma * tf$nn$log_softmax(gamma_logits)))


  elbo <- EE_p_y + E_log_p_p - E_log_q

  gamma_init <- tf$reduce_sum(p_y_on_c, axis=0L)
  gamma_init <- gamma_init - tf$reduce_logsumexp(gamma_init, 0L)
  gamma_init <- tf$transpose(gamma_init)
  gamma_init_ph <- tf$placeholder(shape = shape(N,C), dtype=dtype)
  init_gamma <- tf$assign(gamma_logits, gamma_init_ph)

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
  elbo_diffs <- rep(1e3, 10)
  elbos <- elbo_val
  
  if(verbose) {
    message("Optimizing ELBO")
    pb <- progress_bar$new(total = max_iter+1,
                           format = "  running VB [:bar] :percent | elbo :elbo | chg elbo :meanchange")
    
    pb$tick(0,tokens = list(change = glue("{elbo_diff}%"),
                            elbo = elbo_val,
                            meanchange = NA))
  }



  for( i in seq_len(max_iter) ) {
    if(verbose) {
      pb$tick(tokens = list(change = glue("{round2(elbo_diff)}%"),
                            elbo = elbo_val,
                            meanchange = glue("{100*signif(mean(abs(elbo_diffs)), 2)}%")))
    }

    sess$run(train, feed_dict = fd)

    elbo_new <- sess$run(elbo, feed_dict = fd)
    elbo_diff <- (elbo_new - elbo_val)/abs(elbo_val)
    elbo_diffs <- c(elbo_diffs[-1], elbo_diff)
    elbos <- c(elbos, elbo_new)
    elbo_val <- elbo_new
    

    
    if(is.na(mean(abs(elbo_diffs)))) {
    }
    
    if(mean(abs(elbo_diffs)) < rel_tol)
      break

  }

  if(verbose) {
    message("\nELBO converged or reached max iterations")
  }
  

  rlist <- sess$run(list( tf$nn$softplus(qmu$distribution$loc),gamma, s, tf$exp(log_alpha)), feed_dict = fd)
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
  
  if(verbose) {
    message("Computing final ELBO")
  }
  
  # Get the final ELBO
  final_elbo <- replicate(20, {
    elbo_val <- sess$run(elbo, feed_dict = fd)
  })
  
  convergence_info <- list()
  
  convergence_info$final_elbo <- mean(final_elbo)
  convergence_info$sd_final_elbo <- sd(final_elbo)

  # Close the session ----------
  sess$close()

  convergence_info$elbos <- elbos
  
  names(convergence_info) <- c("final_elbo", "sd_final_elbo", "elbo")
  

  # Correctly name the return vector
  if(P > 0 && K > 0) {
    names(rlist) <- c("mu", "clone_probs", "s", "alpha", "beta", "psi", "W", "chi")
  } else if (K > 0 && P == 0) {
    names(rlist) <- c("mu", "clone_probs", "s", "alpha", "psi", "W", "chi")
  } else if (K == 0 && P == 0) {
    names(rlist) <- c("mu", "clone_probs", "s", "alpha")
  } else {
    names(rlist) <- c("mu", "clone_probs", "s", "alpha")
  }

  list(
    ml_params = rlist,
    convergence_info = convergence_info,
    retained_genes = retained_genes,
    clone_probs_from_snv = clone_probs_from_snv
  )
}
