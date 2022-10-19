# Private function to calculate the density of differences between two
# independent beta distributions by simulation.
#
# The differences are computed as f1 - f0 where f0 is beta(a0, b0) and f1 is
# beta(a1, b1). f0 parameters are based on counts s0 (successes) and n0 (total),
# while f1 parameters are based on counts s1 and n1. Both distributions use the
# same beta_prior parameters.
#
#
ddiffbeta_sim <- function(x, s0, n0, s1, n1, beta_prior = c(1,1), nsim = 1e5, seed = NULL) {
  if (!is.null(seed[1])) set.seed(seed[1])

  a0 <- s0 + beta_prior[1]
  b0 <- n0 - s0 + beta_prior[2]

  a1 <- s1 + beta_prior[1]
  b1 <- n1 - s1 + beta_prior[2]

  rr <- rbeta(nsim, a1, b1) - rbeta(nsim, a0, b0)
  d <- density(rr, from=-1, to = 1)

  dfit <- approx(d, xout = x)
  dfit$y
}


# Private helper function to calculate probability distribution of differences
# between two independent beta distributions by simulation.
#
# The differences are computed as f1 - f0 where f0 is beta(a0, b0) and f1 is
# beta(a1, b1). f0 parameters are based on counts s0 (successes) and n0 (total),
# while f1 parameters are based on counts s1 and n1. Both distributions use the
# same beta_prior parameters.
#
pdiffbeta_sim <- function(x, s0, n0, s1, n1, beta_prior = c(1,1), nsim = 1e5, seed = NULL) {
  if (!is.null(seed[1])) set.seed(seed[1])

  a0 <- s0 + beta_prior[1]
  b0 <- n0 - s0 + beta_prior[2]

  a1 <- s1 + beta_prior[1]
  b1 <- n1 - s1 + beta_prior[2]

  rr <- rbeta(nsim, a1, b1) - rbeta(nsim, a0, b0)
  ecdf(rr)(x)
}


# Private helper function to calculate quantiles of differences between two
# independent beta distributions by simulation.
#
# The differences are computed as f1 - f0 where f0 is beta(a0, b0) and f1 is
# beta(a1, b1). f0 parameters are based on counts s0 (successes) and n0 (total),
# while f1 parameters are based on counts s1 and n1. Both distributions use the
# same beta_prior parameters.
#
qdiffbeta_sim <- function(probs, s0, n0, s1, n1, beta_prior = c(1,1), nsim = 1e5, seed = NULL) {
  if (!is.null(seed[1])) set.seed(seed[1])

  a0 <- s0 + beta_prior[1]
  b0 <- n0 - s0 + beta_prior[2]

  a1 <- s1 + beta_prior[1]
  b1 <- n1 - s1 + beta_prior[2]

  rr <- rbeta(nsim, a1, b1) - rbeta(nsim, a0, b0)
  quantile(rr, probs, type=4)
}

