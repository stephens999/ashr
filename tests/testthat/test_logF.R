context("ashr with logF likelihood")

test_that("logF error with df=(10,10) gives similar answers to normal error",{

  # Simulate a data set.
  set.seed(1)
  x <- rnorm(100) + log(rf(100,df1 = 10,df2 = 10))

  # Fit the ash model with two different likelihood densities: (1) the
  # normal distribution, and (2) the log-F distribution with degrees
  # of freedom (10,10). In the first case, we set the standard errors
  # (s.e.) to be match the empirical standard deviation of random
  # draws from the log-F distribution.
  s            <- sd(log(rf(10000,df1 = 10,df2 = 10)))
  ash.norm.out <- ash(x,s)
  ash.logF.out <- ash.workhorse(x,1,lik = lik_logF(df1 = 10,df2 = 10),
                                optmethod = "mixEM",control = list(tol = 1e-4))

  # Compare the posterior mean estimates from ash using the two
  # different likelihood densities. We expect that the difference
  # between the two estimates should always be small (relative error
  # at most 5%).
  expect_equal(ash.norm.out$result$PosteriorMean,
               ash.logF.out$result$PosteriorMean,
               tolerance = 0.05)
})
