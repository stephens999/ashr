context("ashr with other likelihoods")

test_that("general likelihood with multiple df works", {
  df = c(rep(100,50),rep(2,50))
  s = rgamma(100,1,1)
  betahat = s*rt(n=100,df=df)
  data = set_data(betahat,s,lik = lik_normal(),alpha=0)
  expect_equal(is_normal(data$lik),TRUE)
  expect_equal(is_const(data$lik),TRUE)
  
  data =set_data(betahat,s,lik = lik_t(df),alpha=0)
  expect_equal(is_normal(data$lik),FALSE)
  expect_equal(is_const(data$lik),FALSE)
  
  data =set_data(betahat,s,lik = lik_t(1),alpha=0)
  expect_equal(is_normal(data$lik),FALSE)
  expect_equal(is_const(data$lik),TRUE)
})

test_that("general likelihood with multiple df works", {
  set.seed(10) 
  df = c(rep(100,50),rep(2,50))
  s = rgamma(100,1,1)
  betahat = s*rt(n=100,df=df)
  
  #calc_null_loglik(data)
  
  data =set_data(betahat,s,lik = lik_normal(),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dnorm(betahat,sd=s,log=TRUE)))
  
  data =set_data(betahat,s,lik = lik_t(df),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dt(betahat/s,df=df,log=TRUE)-log(s)))
})

test_that("Poisson (identity link) marginal PMF is correct", {
  y = seq(0, 100)
  g = unimix(1, .1, .2)
  lik = lik_pois(y=y, link="identity")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function (z) {integrate(function(lam) {dpois(z, lam)}, g$a[1], g$b[1])$value / (g$b[1] - g$a[1])}))
  expect_equal(py, true_py)
})

test_that("Poisson (identity link) marginal PMF for point mass is correct", {
  y = seq(0, 100)
  g = unimix(1, .1, .1)
  lik = lik_pois(y=y, link="identity")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function(x) {dpois(x, g$a[1])}))
  expect_equal(py, true_py)
})

test_that("Poisson (identity link) marginal PMF is correct with scale factors", {
  y = seq(0, 100)
  s = 100
  g = unimix(1, 1e-3, 2e-3)
  lik = lik_pois(y=y, scale=s, link="identity")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function (z) {integrate(function(lam) {dpois(z, s * lam)}, g$a[1], g$b[1])$value / (g$b[1] - g$a[1])}))
  expect_equal(py, true_py, tolerance=1e-5, scale=1)
})

test_that("Poisson (log link) marginal PMF is correct", {
  y = seq(0, 100)
  g = unimix(1, log(.1), log(.2))
  lik = lik_pois(y=y, link="log")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function (z) {integrate(function(log_lam) {dpois(z, exp(log_lam))}, g$a[1], g$b[1])$value / (g$b[1] - g$a[1])}))
  expect_equal(py, true_py, tolerance=1e-5, scale=1)
})

test_that("Poisson (log link) marginal PMF for point mass is correct", {
  y = seq(0, 100)
  g = unimix(1, log(.1), log(.1))
  lik = lik_pois(y=y, link="log")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function(x) {dpois(x, exp(g$a[1]))}))
  expect_equal(py, true_py)
})

test_that("Poisson (log link) marginal PMF is correct with scale factors", {
  y = seq(0, 100)
  s = 100
  g = unimix(1, log(1e-3), log(2e-3))
  lik = lik_pois(y=y, scale=s, link="log")
  data = set_data(rep(0, length(y)), 1, lik)
  py = drop(comp_dens_conv(g, data))
  true_py = unlist(lapply(y, function (z) {integrate(function(log_lam) {dpois(z, s * exp(log_lam))}, g$a[1], g$b[1])$value / (g$b[1] - g$a[1])}))
  expect_equal(py, true_py, tolerance=1e-5, scale=1)
})

test_that("ln p(x = 0 | s) is correct for Poisson likelihood (log link)", {
  s = 10 ^ seq(log10(1e3), log10(1e6), .1)
  y = rep(0, length(s))
  g = unimix(1, log(1e-3), log(2e-3))
  lik = lik_pois(y=y, scale=s, link="log")
  data = set_data(rep(0, length(y)), 1, lik)
  log_py = drop(log(comp_dens_conv(g, data)))
  true_log_py = log(unlist(lapply(s, function (z) {integrate(function(log_lam) {dpois(0, z * exp(log_lam))}, g$a[1], g$b[1])$value / (g$b[1] - g$a[1])})))
  ## Work around a recent bug comparing Inf (called from expect_equal)
  ##
  ## https://github.com/r-lib/testthat/issues/789
  ## https://github.com/r-lib/testthat/commit/992ddd82fd7b6f1fdc5bb66c31db94277f3df126
  expect_true(all((true_log_py == log_py | abs(true_log_py - log_py) < 1e-5)))
})

test_that("Poisson (identity link) returns sensible marginal PMF on real data", {
  skip("save time")
  dat = readRDS("test_pois_data.Rds")
  fit <- ash_pois(dat$x, dat$scale, link="identity", mixcompdist="halfuniform")
  F = comp_dens_conv(fit$fitted_g, fit$data)
  expect_true(all(F < 1))
})

test_that("Poisson (log link) returns sensible marginal PMF on real data", {
  dat = readRDS("test_pois_data.Rds")
  lam <- dat$x / dat$scale
  eps = 1 / mean(dat$scale)
  log_lam <- log(lam + eps)
  se_log_lam <- sqrt(var(lam) / (lam + eps)^2)
  mixsd <- seq(.1 * min(se_log_lam), max(2 * sqrt(log_lam^2 - se_log_lam^2)), by=.5 * log(2))
  fit <- ash_pois(dat$x, dat$scale, link="log", mixcompdist="halfuniform", mixsd=mixsd)
  F = comp_dens_conv(fit$fitted_g, fit$data)
  expect_true(all(F < 1))
})
