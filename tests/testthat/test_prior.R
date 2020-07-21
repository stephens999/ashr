context("ashr \"prior\" argument")

test_that("numeric and (partial) string arguments both work", {
  betahat <- c(1.01636974224394, -2.05686254738995, -0.7135781676358,
               -1.16906745227838, -0.917039991627176)
  
  sebetahat <- c(1.02572223086898, 0.499285201440522, 0.476520330150983,
                 0.624576594477857, 0.198152636610839)
  
  aout1 <- ash.workhorse(betahat = betahat[1:5], sebetahat = sebetahat[1:5],
                         g = normalmix(rep(0, 5), rep(0, 5), 0:4), fixg = FALSE,
                         prior = "uniform")
  aout2 <- ash.workhorse(betahat = betahat[1:5], sebetahat = sebetahat[1:5],
                         g = normalmix(rep(0, 5), rep(0, 5), 0:4), fixg = FALSE,
                         prior = rep(1, 5))
  expect_identical(aout1$result$PosteriorMean, aout2$result$PosteriorMean)
  
  aout3 <- ash.workhorse(betahat = betahat[1:5], sebetahat = sebetahat[1:5],
                         g = normalmix(rep(0, 5), rep(0, 5), 0:4), fixg = FALSE,
                         prior = "null")
  expect_false(identical(aout1$result$PosteriorMean, aout3$result$PosteriorMean))
})

test_that("pi is nonzero for mixture components where prior > 1", {
  x <- 10:20
  s <- rep(1, 11)
  g <- unimix(rep(0, 3), c(0, -1, -20), c(0, 1, 20))
  aout1 <- ash(x, s, g = g, fixg = FALSE, prior = c(1, 1, 1), mixcompdist = "uniform")
  expect_false(all(aout1$fitted_g$pi > 0))
  
  aout2 <- ash(x, s, g = g, fixg = FALSE, prior = c(10, 10, 10), mixcompdist = "uniform")
  expect_true(all(aout2$fitted_g$pi > 0))
  
  aout3 <- ash(x, s, g = g, fixg = FALSE, prior = c(10, 1, 10), mixcompdist = "uniform")
  expect_true(aout3$fitted_g$pi[1] > 0)
  expect_false(aout3$fitted_g$pi[2] > 0)
})