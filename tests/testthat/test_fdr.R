library(fdr)
library(LaplacesDemon)
context("Output")
test_that("rmultn/dmultn", {
  set.seed(10237)
  mu <- matrix(rnorm(100), 100)
  Sigma <- crossprod(matrix(rnorm(100*100), 100, 100))

  set.seed(1)
  ld_out <- LaplacesDemon::rmvn(n = 1, mu = c(mu), Sigma = Sigma)
  set.seed(1)
  fdr_out <- rmultn(mu, Sigma)
  expect_equal(c(ld_out), c(fdr_out))

  ld_out <- LaplacesDemon::dmvn(ld_out, c(mu), Sigma, log = TRUE)
  fdr_out <- dmultn(t(fdr_out), c(mu), Sigma, logd = TRUE)
  expect_equal(c(ld_out), c(fdr_out))
})

test_that("rmatn/dmatn", {
  mu <- matrix(rnorm(100), 20, 5)
  P <- crossprod(matrix(rnorm(5*5), 5, 5))
  Q <- crossprod(matrix(rnorm(20*20), 20, 20))

  L <- chol(P)
  C <- t(chol(Q))
  set.seed(1)
  r <- mu + C %*% matrix(rnorm(100), 20, 5) %*% L
  set.seed(1)
  fdr_out <- rmatn(mu, P, Q)
  expect_equal(r, fdr_out)

  ld_out <- LaplacesDemon::dmatrixnorm(fdr_out, mu, Q, P, log = TRUE)
  fdr_out <- dmatn(fdr_out, mu, P, Q, logd = TRUE)
  expect_equal(c(ld_out), c(fdr_out))
})


test_that("rinvwish/dinvwish", {
  set.seed(1)
  S <- crossprod(matrix(rnorm(10*10), 10, 10))
  nu <- 20
  set.seed(1)
  ld_out <- LaplacesDemon::rinvwishart(nu, S)
  set.seed(1)
  fdr_out <- rinvwish(nu, S)
  expect_equal(ld_out, fdr_out)

  ld_out <- LaplacesDemon::dinvwishart(fdr_out, nu, S, log = TRUE)
  fdr_out <- dinvwish(fdr_out, nu, S, logd = TRUE)
  expect_equal(ld_out, fdr_out)
})
