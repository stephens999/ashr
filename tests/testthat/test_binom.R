context("ashr with Binomial likelihoods")

test_that("lik_binom (identity link) fitted g is close to true g",{
  set.seed(1)

  # true prior g: 0.4*U(0.1,0.9)+0.6*delta(0.5)
  trueg = unimix(c(0.4,0.6),c(0.5,0.1),c(0.5,0.9)) 
  p = c(rep(0.5,400), runif(600,0.1,0.9)) # generate p from g
  n = rep(100,1000)
  x = rbinom(1000,n,p) # Binomial observations
  out <- capture.output(
    ash.binom.out <- ash(rep(0,length(x)),1,lik = lik_binom(x,n),g = trueg,
                         control = list(maxiter.sqp = 40,verbose = TRUE)))
  
  # Check if the estimated mixture proportion for components delta(0.5) and U(0.1,0.9)
  # is close to the true mixture proportion (0.4,0.6)
  expect_equal(ash.binom.out$fitted_g$pi, c(0.4,0.6), tolerance = 0.025)
})

test_that("lik_binom (identity link) fitted g is close to true g",{
  # Simulate a Binomial dataset
  set.seed(1)
  truemode = 0.3
  trueg = unimix(c(0.5,0.5),c(0.3,0.1),c(0.3,0.5)) 
  p = c(rep(0.3,500), runif(500,0.1,0.5)) # generate p from g
  n = rep(100,1000)
  x = rbinom(1000,n,p) # Binomial observations
  ash.binom.out = ash(rep(0,length(x)),1,lik=lik_binom(x,n),mode="estimate")
  
  # Check if the estimated mode is close to the true mode 0.3
  expect_equal(ash.binom.out$fitted_g$a[1], truemode, tolerance = 0.05, scale=truemode)
})

test_that("lik_binom (identity link) with big n gives similar answers to normal likelihood",{
  # Simulate a Binomial data set with n=200
  set.seed(1)
  p = c(rep(0.3,500), runif(500,0.1,0.5)) # generate p
  n = rep(200,1000) 
  x = rbinom(1000,n,p) # Binomial observations
  
  # Fit the ash model with two different likelihood densities: (1) the
  # normal distribution with (s.e.) to be match the standard deviations of 
  # Binomial distribution, and (2) the Binomial distribution
  ash.binom.out = ash(rep(0,length(x)),1,lik=lik_binom(x,n))
  ash.norm.out = ash(x/n, sqrt(p*(1-p)/n), mode="estimate", prior="uniform")
  
  # Compare the posterior mean estimates from ash using the two
  # different likelihood densities. We expect that the difference
  # between the two estimates should always be small (relative error
  # at most 5%).
  expect_equal(ash.norm.out$result$PosteriorMean,
               ash.binom.out$result$PosteriorMean,
               tolerance = 0.05)
})

test_that("lik_binom (logit link) fitted g is close to true g",{
    
  # Simulate a Binomial dataset
  set.seed(1)
  trueg = unimix(c(0.5,0.5),c(0,-3),c(0,3)) 
  logitp = c(rep(0,500), runif(500,-3,3))
  p = 1/(1+exp(-logitp))
  n = rep(1000,1000)
  x = rbinom(1000,n,p) # Binomial observations
  out <- capture.output(
    ash.binom.out <- ash(rep(0,length(x)),1,
                         lik = lik_binom(x,n,link = "logit"),
                         g = trueg,prior = "uniform",
                         control = list(verbose = TRUE)))
  
  # Check if the estimated mixture proportion for components
  # delta(0.5) and U(-3,3) is close to the true mixture proportion
  # (0.5,0.5).
  expect_equal(ash.binom.out$fitted_g$pi, c(0.5,0.5), tolerance = 0.05)
})

test_that("lik_binom (logit link) fitted g is close to true g",{
  # Simulate a Binomial dataset
  set.seed(1)
  truemode = 0
  logitp = c(rep(0,800), runif(200,-3,3))
  p = 1/(1+exp(-logitp))
  n = rep(100,1000)
  x = rbinom(1000,n,p) # Binomial observations
  ash.binom.out = ash(rep(0,length(x)),1,lik=lik_binom(x,n,link="logit"),
                      mode = "estimate")
  
  # Check if the estimated mode is close to the true mode
  expect_equal(ash.binom.out$fitted_g$a[1], truemode, tolerance = 0.05)
})
