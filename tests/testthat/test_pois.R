context("ashr with Poisson likelihoods")

test_that("lik_pois (identity link) fitted g is close to true g",{
  set.seed(1)
  trueg = unimix(c(0.5,0.5),c(1,1),c(1,5))
  lambda = c(rep(1,500), runif(500,1,5))
  x = rpois(1000,lambda)
  ash.pois.out = ash_pois(x, link="identity", g=trueg)
  expect_equal(ash.pois.out$fitted_g$pi, trueg$pi, tolerance = 0.05)
})

test_that("lik_pois (identity link) fitted mode is close to true mode",{
  set.seed(1)
  truemode = 50
  lambda = c(rnorm(1000,truemode,5))
  x = rpois(1000,lambda)
  ash.pois.out = ash_pois(x)
  expect_equal(ash.pois.out$fitted_g$a[1], truemode, tolerance = 0.05, scale=truemode)
})

test_that("lik_pois (identity link) with high intensities gives similar answers to normal likelihood",{
  set.seed(1)
  lambda = c(rnorm(1000,200,5))
  x = rpois(1000,lambda)
  ash.pois.out = ash_pois(x)
  # For large lambda, Poisson(lambda) is approximately N(lambda, lambda)
  ash.norm.out = ash(x, sqrt(lambda), mode="estimate", prior="uniform")
  expect_equal(ash.norm.out$result$PosteriorMean,
               ash.pois.out$result$PosteriorMean,
               tolerance = 0.05)
})

test_that("lik_pois (log link) fitted g is close to true g",{
  set.seed(1)
  trueg = unimix(c(0.8,0.2),c(0,-3),c(0,3)) 
  loglambda = c(rep(0,800), runif(200,-3,3)) 
  x = rpois(1000, exp(loglambda))
  ash.pois.out = ash_pois(x, link="log", g=trueg)
  expect_equal(ash.pois.out$fitted_g$pi, trueg$pi, tolerance = 0.05)
})

test_that("lik_pois (log link) fitted mode is close to true mode",{
  set.seed(1)
  truemode = 4
  N = 500
  # Set the SD in such a way that exp(loglambda) does not become too big.
  # Otherwise, expint::gammainc will fail
  loglambda = rnorm(N, 4, .1)
  x = rpois(N, exp(loglambda))
  ash.pois.out = ash_pois(x, link="log", mode="estimate")
  expect_equal(ash.pois.out$fitted_g$a[1], truemode, tolerance = 1e-2, scale=1)
})

test_that("Mode estimation for lik_pois finds an acceptable solution", {
  set.seed(1)
  ## Load example 10X Genomics data
  dat = readRDS("test_pois_data.Rds")
  m0 = ashr::ash(rep(0, nrow(dat)), 1, lik=ashr::lik_pois(dat$x, scale=dat$scale, link="identity"), mode="estimate")
  lam = dat$x / dat$scale
  m1 = ashr::ash(rep(0, nrow(dat)), 1, lik=ashr::lik_pois(dat$x, scale=dat$scale, link="identity"), mode=c(min(lam), max(lam)))
  expect_equal(m0$loglik, m1$loglik, tolerance=1, scale=1)
})

test_that("Mode estimation for lik_pois gives same answer under identity and log link", {
  set.seed(1)
  ## Typical values for scRNA-seq data from Sarkar et al. 2019
  s = 1e5
  log_mu = runif(n=1, min=-12, max=-8)
  log_phi = runif(n=1, min=-6, max=0)
  N = 1000
  lam = rgamma(n=N, shape=exp(-log_phi), scale=exp(log_mu + log_phi))
  x = rpois(n=N, lambda=s * lam)
  dat = data.frame(cbind(x, s))
  fit0 = ashr::ash_pois(dat$x, dat$s, link="identity", mixcompdist="halfuniform")
  fit1 = ashr::ash_pois(dat$x, dat$s, link="log", mixcompdist="halfuniform")
  expect_equal(fit0$loglik, fit1$loglik, tolerance=1e-2, scale=1)
  expect_equal(fit0$fitted_g$a[1], exp(fit1$fitted_g$a[1]), tolerance=1e-6, scale=1)
})
