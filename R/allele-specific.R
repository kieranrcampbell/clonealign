

#' Construct allele specific complete data likelihood
#'
#' Returns cell-by-clone complete likelihood
#'
#' @export
construct_ai_likelihood <- function(clone_allele_ph, alt_ph, cov_ph, dtype, C, N) {
  params_1 <- list(low = c('alpha'=0.1, 'beta'=1.9),
                   high = c('alpha'=1.9, 'beta'=0.1))
  params_2 <- c('alpha' = 2, 'beta' = 2)

  p1_low <- tf$log(tf$constant(0.5, dtype = dtype)) + beta_binomial_log_prob(alt_ph, cov_ph,
                                                                             tf$constant(params_1$low['alpha'], dtype = dtype),
                                                                             tf$constant(params_1$low['beta'], dtype = dtype))

  p1_high <- tf$log(tf$constant(0.5, dtype = dtype)) + beta_binomial_log_prob(alt_ph, cov_ph,
                                                                              tf$constant(params_1$high['alpha'], dtype = dtype),
                                                                              tf$constant(params_1$high['beta'], dtype = dtype))

  p1 <- tf$reduce_logsumexp(tf$stack(list(p1_low, p1_high)), axis = 0L)

  p2 <- beta_binomial_log_prob(alt_ph, cov_ph,
                               tf$constant(params_2['alpha'], dtype = dtype),
                               tf$constant(params_2['beta'], dtype = dtype))


  p1c <- tf$stack(rep(list(p1), C)) # clone by variant by cell
  p2c <- tf$stack(rep(list(p2), C))

  is_copy_num_2 <- clone_allele_ph == 2
  is_copy_num_2 <- tf$stack(rep(list(is_copy_num_2), N)) # cell by variant by clone
  is_copy_num_2 <- tf$transpose(is_copy_num_2, c(2L, 1L, 0L)) # clone by variant by cell

  L <- tf$where(is_copy_num_2, p2c, p1c) #

  L <- tf$transpose(tf$reduce_sum(L, axis = 1L))
  L
}



beta_binomial_log_prob <- function(k, n, alpha, beta) {

  ll <- tf$lgamma(n+1) - tf$lgamma(k+1) - tf$lgamma(n-k+1)
  ll <- ll + tf$lgamma(k+alpha) + tf$lgamma(n-k+beta) - tf$lgamma(alpha+beta+n)
  ll <- ll - tf$lgamma(alpha) - tf$lgamma(beta) + tf$lgamma(alpha+beta)
  ll
}

#' @keywords internal
sanitize_allele_info <- function(V, clone_allele, cov, ref, N, C) {


  stopifnot(ncol(clone_allele) == C) # Number of clones
  stopifnot(nrow(cov) == N)
  stopifnot(nrow(ref) == N)
  stopifnot(ncol(ref) == V)
  stopifnot(ncol(cov) == V)


}



