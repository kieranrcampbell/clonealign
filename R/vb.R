


clonealign_vb <- function(Y, L, s, maxiter = 200) {
  Lp <- t(t(L) / colMeans(L))
  
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
    alpha = 2,
    beta = 1,
    kappa = rep(1, C)
  )
  
  for(i in seq_len(maxiter)) {
    print(i)
    # Update epsilon
    epsilon <- hyperparams$kappa + colSums(phi)
    
    # Update phi
    for(n in 1:N) {
      for(c in 1:C) {
        p <- digamma(epsilon[c]) - digamma(sum(epsilon))
        p <- p + sum(sapply(seq_len(G), function(g) {
          Y[n,g] * (digamma(a[g]) - log(b[g]) + log(s[n] * Lp[g,c])) - s[n] * Lp[g,c] * a[g] / b[g]
        }))
        phi[n,c] <- p
      }
      phi[n,] <- phi[n,] - logSumExp(phi[n,])
    }
    phi <- exp(phi)
    
    # Update mu
    a <- hyperparams$alpha + colSums(Y)
    b <- hyperparams$beta + colSums(s * (phi %*% t(Lp)))
    
  }
}
