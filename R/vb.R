

#' @export
clonealign_vb <- function(Y, L, s, maxiter = 200,
                          normalise_L = TRUE) {
  Lp <- L
  if(normalise_L) {
    Lp <- t(t(L) / colMeans(L))
  }

  N <- nrow(Y)
  C <- ncol(L)
  G <- ncol(Y)
  stopifnot(nrow(L) == G)

  mu <- colMeans(Y)

  # Variational parameters
  a <- mu
  b <- rep(1, G)
  phi <- matrix(1/3, nrow = N, ncol = C)
  epsilon <- rep(1/3, C)

  hyperparams <- list(
    alpha = mean(Y),
    beta = 1,
    kappa = rep(1 / C, C)
  )

  phi1_trace <- NULL
  mu1_trace <- NULL

  for(i in seq_len(maxiter)) {
    # Update epsilon
    epsilon <- hyperparams$kappa + colSums(phi)

    # Update phi
    digamma_a <- digamma(a)
    log_b <- log(b)
    digamma_sum <- digamma(sum(epsilon))
    for(c in 1:C) {
      dg <- digamma(epsilon[c]) - digamma_sum
      for(n in 1:N) {
        p <- dg + sum(Y[n,] * (digamma_a - log_b + log(s[n] * Lp[,c])) - s[n] * Lp[,c] * a / b)
        phi[n,c] <- p
      }
    }
    for(n in 1:N) { # Bring this outside to avoid computing digamma twice
      phi[n,] <- phi[n,] - logSumExp(phi[n,])
    }
    phi <- exp(phi)

    # phi1_trace <- c(phi1_trace, phi[1,1])
    # mu1_trace <- c(mu1_trace, a[1] / b[1])

    # Update mu
    a <- hyperparams$alpha + colSums(Y)
    b <- hyperparams$beta + colSums(s * (phi %*% t(Lp)))


  }

  l <- list(a = a, b = b,
            phi = phi, epsilon = epsilon)

}

#' @export
clone_assignments <- function(phi) {
  apply(phi, 1, which.max)
}
