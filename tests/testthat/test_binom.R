test_that("binom_lik fitted g is close to true g",{
  # Simulate a Binomial dataset
  set.seed(1)
  trueg = unimix(c(0.5,0.5),c(0.5,0.1),c(0.5,0.9)) # true prior g: 0.5*U(0.1,0.9)+0.5*delta(0.5)
  p = c(rep(0.5,500), runif(500,0.1,0.9)) # generate p from g
  n = rep(100,1000)
  x = rbinom(1000,n,p) # Binomial observations
  ash.binom.out = ash(rep(0,length(x)),1,lik=binom_lik(x,n),g=trueg)
  
  # Check if the estimated mixture proportion for components delta(0.5) and U(0.1,0.9)
  # is close to the true mixture proportion (0.5,0.5)
  expect_equal(ash.pois.out$fitted_g$pi, c(0.5,0.5), tolerance = 0.05)
})

test_that("binom_lik fitted g is close to true g",{
  # Simulate a Binomial dataset
  set.seed(1)
  truemode = 0.3
  trueg = unimix(c(0.5,0.5),c(0.3,0.1),c(0.3,0.5)) 
  p = c(rep(0.3,500), runif(500,0.1,0.5)) # generate p from g
  n = rep(100,1000)
  x = rbinom(1000,n,p) # Binomial observations
  ash.binom.out = ash(rep(0,length(x)),1,lik=binom_lik(x,n),mode="estimate")
  
  # Check if the estimated mixture proportion for components delta(0.5) and U(0.1,0.9)
  # is close to the true mixture proportion (0.5,0.5)
  expect_equal(ash.binom.out$fitted_g$a[1], truemode, tolerance = 0.05, scale=x)
})

test_that("binom_lik with big n gives similar answers to normal likelihood",{
  # Simulate a Binomial data set with n=200
  set.seed(1)
  p = c(rep(0.3,500), runif(500,0.1,0.5)) # generate p
  n = rep(200,1000) 
  x = rbinom(1000,n,p) # Binomial observations
  
  # Fit the ash model with two different likelihood densities: (1) the
  # normal distribution with (s.e.) to be match the standard deviations of 
  # Binomial distribution, and (2) the Poisson distribution
  ash.binom.out = ash(rep(0,length(x)),1,lik=binom_lik(x,n))
  ash.norm.out = ash(x/n, sqrt(p*(1-p)/n), mode="estimate", prior="uniform")
  
  # Compare the posterior mean estimates from ash using the two
  # different likelihood densities. We expect that the difference
  # between the two estimates should always be small (relative error
  # at most 5%).
  expect_equal(ash.norm.out$result$PosteriorMean,
               ash.binom.out$result$PosteriorMean,
               tolerance = 0.05, scale=x)
})
