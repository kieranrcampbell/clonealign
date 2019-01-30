library(ggplot2)

set.seed(2345234)
G <- 200
N <- 500
C <- 3

rho <- sample(c(0,1), G, replace = TRUE, prob = c(0.2, 0.9))

pi <- sample(seq_len(C), N, replace = TRUE)

mu <- runif(G, 1, 2)
beta <- mu
phi <- rgamma(G, 4, 1)

L <- matrix(sample(seq_len(C), C * G, replace = TRUE), ncol = C)
Lp <- t(t(L) / colMeans(L))

s <- runif(N, 500, 10000)


Y <- matrix(NA, nrow = N, ncol = G)

for(n in seq_len(N)) {
  for(g in seq_len(G)) {
    m <- s[n] * ((1 - rho[g]) * mu[g] + rho[g] * beta[g] * Lp[g,pi[n]])
    Y[n,g] <- rnbinom(1, mu = m, size = phi[g])
  }
}

data <- list(
  N = N,
  G = G,
  C = C,
  Y = Y,
  L = Lp,
  s = s,
  l_g_hat = rowMeans(Lp)
)

params <- cbind(mu, beta, phi)
names(params) <- c("mu", "beta", "phi")

devtools::load_all(".")

traces <- clonealign::gibbs_pi_rho(rho, data, params, n_iter = 20)

pi_traces <- traces$pi_trace
cp <- clone_probs_from_gibbs(pi_traces, C)
clone <- apply(cp, 1, which.max)
table(clone, pi)

rh <- rho_probs_from_gibbs(traces$rho_trace)
rh_mle <- apply(rh, 1, which.max)
table(rho, rh_mle)



# So that works - what about Q ? ------------------------------------------


devtools::load_all(".")
n_convergence_failures <- 0
pnew <- lapply(seq_len(data$G), function(g) {
  opt <- optim(par = params[g,], # (mu,beta,phi)
               fn = Q_g,
               gr = Qgr_g,
               y = data$Y[,g], l = data$L[g,], data = data,
               pi_traces = matrix(pi, nrow = 1),
               rho_traces = rho[g],
               lambda = 10,
               l_g_hat = data$l_g_hat[g],
               method = "L-BFGS-B",
               lower = c(1e-10, 1e-10, 1e-10),
               upper = c(max(data$Y), 1e6, 1e6),
               control = list())
  if(opt$convergence != 0) {
    n_convergence_failures <<- n_convergence_failures + 1
  }
  c(opt$par, -opt$value)
})
print(n_convergence_failures)
pnew <- do.call(rbind, pnew)
params <- pnew[,c('mu', 'beta', 'phi')]

qplot(mu, params[,'mu'], color = factor(rho))
qplot(beta, params[,'beta'], color = factor(rho))
qplot(params[,'mu'], params[,'beta'], color = factor(rho))



# Try it all together -----------------------------------------------------


devtools::load_all(".")
em <- inference_em(Y, Lp, rel_tol = 1e-9, multithread = FALSE, max_iter = 10,
                   lambda = 1, rho_init = rho)

cp <- clone_probs_from_gibbs(em$gibbs_samples$pi_trace, ncol(Lp))
rp <- rho_probs_from_gibbs(em$gibbs_samples$rho_trace)

em_clone <- apply(cp, 1, which.max)
table(em_clone, pi)

rho_mle <- apply(rp, 1, which.max) - 1
table(rho, rho_mle)

qplot(mu, em$mu, color = factor(rho))
qplot(beta, em$beta, color = factor(rho))



# Test gradients ----------------------------------------------------------

library(numDeriv)

g <- 1

Q_g2 <- function(p) {
  Q_g(p, y = data$Y[,g], l = data$L[g,], data = data,
     pi_traces = matrix(pi, nrow = 1),
     rho_traces = rho[g],
     lambda = 10,
     l_g_hat = data$l_g_hat[g])
}

g <- 1
print(rho[g])
Qgr_g(params[g,], y = data$Y[,g], l = data$L[g,], data = data,
    pi_traces = matrix(pi, nrow = 1),
    rho_traces = rho[g],
    lambda = 10,
    l_g_hat = data$l_g_hat[g])

grad(Q_g2, params[g,])

