context("ashr with Poisson likelihoods")

test_that("lik_pois (identity link) fitted g is close to true g",{
  # Simulate a Poisson dataset
  set.seed(1)
  trueg = unimix(c(0.5,0.5),c(1,1),c(1,5)) # true prior g: 0.5*U(1,5)+0.5*delta(1)
  lambda = c(rep(1,500), runif(500,1,5)) # generate lambda from g
  x = rpois(1000,lambda) # Poisson observations
  out <- capture.output(
    ash.pois.out <- ash(rep(0,length(x)),1,lik=lik_pois(x),g=trueg,
                        control = list(verbose = TRUE)))
  
  # Check if the estimated mixture proportion for components delta(0.5) and U(0.1,0.9)
  # is close to the true mixture proportion (0.5,0.5)
  expect_equal(ash.pois.out$fitted_g$pi, c(0.5,0.5), tolerance = 0.05)
})

test_that("lik_pois (identity link) fitted mode is close to true mode",{
  # Simulate a Poisson dataset
  set.seed(1)
  truemode = 50 # set mode of prior g
  lambda = c(rnorm(1000,truemode,5)) # generate lambda from g
  x = rpois(1000,lambda) # Poisson observations
  ash.pois.out = ash(rep(0,length(x)),1,lik=lik_pois(x),mode="estimate")
  
  # Check if the estimated mode is close to the true mode 50
  expect_equal(ash.pois.out$fitted_g$a[1], truemode, tolerance = 0.05, scale=truemode)
})

test_that("lik_pois (identity link) with high intensities gives similar answers to normal likelihood",{
  # Simulate a Poisson data set with relatively high intensities
  set.seed(1)
  lambda = c(rnorm(1000,200,5)) # simulate intensities around 200
  x = rpois(1000,lambda)
  
  # Fit the ash model with two different likelihood densities: (1) the
  # normal distribution with (s.e.) to be match the standard deviations of 
  # Poisson distribution sqrt(lambda), and (2) the Poisson distribution
  ash.pois.out = ash(rep(0,length(x)),1,lik=lik_pois(x))
  ash.norm.out = ash(x, sqrt(lambda), mode="estimate", prior="uniform")
  
  # Compare the posterior mean estimates from ash using the two
  # different likelihood densities. We expect that the difference
  # between the two estimates should always be small (relative error
  # at most 5%).
  expect_equal(ash.norm.out$result$PosteriorMean,
               ash.pois.out$result$PosteriorMean,
               tolerance = 0.05)
})

test_that("lik_pois (log link) fitted g is close to true g",{
  # Simulate a Poisson dataset
  set.seed(1)
  trueg = unimix(c(0.8,0.2),c(0,-3),c(0,3)) 
  loglambda = c(rep(0,800), runif(200,-3,3)) 
  lambda = exp(loglambda)
  x = rpois(1000,lambda) # Poisson observations
  out <- capture.output(
    ash.pois.out <- ash(rep(0,length(x)),1,lik = lik_pois(x,link="log"),
                        g = trueg,control = list(verbose = TRUE)))
  
  # Check if the estimated mixture proportion for components delta(0)
  # and U(-3,3) is close to the true mixture proportion (0.8,0.2)
  expect_equal(ash.pois.out$fitted_g$pi, c(0.8,0.2), tolerance = 0.05)
})

test_that("lik_pois (log link) fitted mode is close to true mode",{
  # Simulate a Poisson dataset
  set.seed(1)
  truemode = 4
  loglambda = c(rep(4,500), rnorm(500,4,1)) # simulate log(lambda) from distn w/ mode at 4
  lambda = exp(loglambda)
  x = rpois(1000,lambda) # Poisson observations
  ash.pois.out = ash(rep(0,length(x)),1,lik=lik_pois(x,link="log"),mode="estimate")
  
  # Check if the estimated mode is close to the true mode 50
  expect_equal(ash.pois.out$fitted_g$a[1], truemode, tolerance = 0.05, scale=truemode)
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
