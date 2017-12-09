

#' @export
likelihood_yn <- function(y, L, s, pi, params) {
  ll <- 0

  # Likelihood
  # for(g in seq_along(y)) {
    # if(L[g, pi ] > 0) {
      m <- L[, pi] * s * params[, 'mu']
      ll <- ll + sum(dnbinom2(y, mu = m, size = params[, 'phi']))
    #} else {
    #  stop("Zero L's not supported")
    #  ll <- ll + dpois(y[g], 0.1, log = TRUE)
    #}
  #}
  ll
}

#' Computes gamma_{nc} = p(pi_n = c), returning
#' N by C matrix
#'
#' @importFrom matrixStats logSumExp
p_pi <- function(data, params) {
  gamma <- matrix(NA, nrow = data$N, ncol = data$C)
  for(n in seq_len(data$N)) {
    for(c in seq_len(data$C)) {
      gamma[n,c] <- likelihood_yn(y = data$Y[n,],
                                  L = data$L,
                                  s = data$s[n],
                                  pi = c,
                                  params = params)
    }
    gamma[n,] <- exp(gamma[n,] - logSumExp(gamma[n,]))
  }
  gamma
}

#' Computes Q(theta|theta^(t))
#' (function to be optimised under EM)
#' @export
Q_g <- function(pars, y, l, gamma, data) {
  mu <- pars[1]
  phi <- pars[2]
  qq <- 0
  for(c in seq_len(data$C)) {
    m <- l[c] * data$s * mu # N length vector for given gene of means
    l_c <- dnbinom2(y, mu = m, size = phi) # p(y_g | pi)
    qq <- qq + sum(gamma[,c] * l_c )
  }
  -qq
}

#' @export
Q <- function(pars, Y, L, gamma, data) {
  mu <- pars[1:data$G]
  phi <- pars[(data$G+1):(2*data$G)]
  qq <- 0
  for(g in seq_len(data$G)) {
    for(c in seq_len(data$C)) {
      m <- L[g,c] * data$s * mu[g] # N length vector for given gene of means
      l_c <- dnbinom2(Y[,g], mu = m, size = phi[g]) # p(y_g | pi)
      qq <- qq + sum(gamma[,c] * l_c )
    }
  }
  -qq
}

#' importFrom glue glue
#' @export
clonealign_em <- function(Y, L, s = NULL, n_iter, tol = 1e-3,
                          parallelise = TRUE) {

  zero_gene_means <- colMeans(Y) == 0
  message(glue("Removing {sum(zero_gene_means)} genes with zero counts"))
  Y <- Y[, !zero_gene_means]
  L <- L[!zero_gene_means,]


  N <- nrow(Y)
  G <- ncol(Y)
  C <- ncol(L)

  # Sanity checks
  stopifnot(nrow(L) == G)
  stopifnot(length(s) == N)
  stopifnot(all(s > 0))


  if(is.null(s)) {
    s <- scran::computeSumFactors(t(Y))
  }

  # Always normalise L
  L <- t( t(L) / colMeans(L) )

  # Initialise
  params <- cbind(
    colMeans(Y / s) + 0.01,
    rep(1, G)
  )
  colnames(params) <- c("mu", "phi")

  data <- list(
    Y = Y,
    L = L,
    s = s,
    N = N,
    G = G,
    C = C
  )

  q_old <- NA
  qvals <- matrix(NA, nrow = n_iter, ncol = G)

  for(i in seq_len(n_iter)) {

    # E step
    gamma <- p_pi(data, params)

    # print(c(sqrt(mean(gamma^2)), sqrt(mean(params[,'mu']^2)), sqrt(mean(params[,'phi']^2))))
    # Q_g(params[g,], Y[,g], L[g,], gamma, data)
    # Q(as.vector(params), Y, L, gamma, data)

    # M step
    if(parallelise) {
      pnew <- sapply(seq_len(data$G), function(g) {
        opt <- optim(par = params[g, ], # (mu,phi)
                     fn = Q_g,
                     y = Y[,g], l = L[g,], gamma = gamma, data = data,
                     method = "L-BFGS-B",
                     lower = c(1e-9,1e-9))
        c(opt$par, -opt$value)
      })
      params <- t(pnew[c('mu', 'phi'),])
      # print(c(dim(qvals), dim(pnew)))
      qvals[i,] <- pnew[3,]
      q_new <- sum(pnew[3,])
    } else {
      par <- as.vector(params)
      opt <- optim(par = par, # (mu,phi)
                   fn = Q,
                   Y = Y, L = L, gamma = gamma, data = data,
                   method = "L-BFGS-B",
                   lower = rep(0.00001, 2 * data$G),
                   control = list(maxit = 1e3))
      params <- matrix(opt$par, ncol = 2)
      colnames(params) <- c("mu", "phi")
      q_new <- -opt$value
    }
    q_diff <- ((q_new - q_old)) #  / q_old * 100)

    message(glue("Old: {q_old}\tNew: {q_new}\tChange: {q_diff}"))
    # message(glue("New log-likelihood: {q_new}"))

    # if(!is.na(q_diff)) {
    #   if(q_diff < tol) {
    #     print(glue("Converged after {i} iterations"))
    #     break
    #   }
    # }
    q_old <- q_new
  }

  gamma <- p_pi(data, params)
  rlist <- list(
    gamma = gamma,
    mu = params[, 'mu'],
    phi = params[, 'phi'],
    qvals = qvals
  )
}


## Testing time

test <- function() {
  data2 <- list(
    Y = Y[1:3, 1:2],
    L = matrix(c(0,.5,1,0.1), ncol = 2),
    s = s[1:3],
    N = 3,
    G = 2,
    C = 2
  )
  params2 <- matrix(1, nrow = data2$G, ncol = 2)
  colnames(params2) <- c('mu', 'phi')
  p_pi(data2, params2)
}

