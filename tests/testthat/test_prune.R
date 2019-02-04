context("pruning mixture components in ashr")

test_that("test pruning", {
  g= normalmix(c(0.5,0.4,0.1),c(0,0,0),c(1,2,3))
  g.pruned = prune(g,0.2)
  expect_equal(g.pruned$mean,c(0,0))
  g = unimix(c(0.5,0.4,0.1),c(0,0,0),c(1,2,3))
  g.pruned = prune(g,0.2)
  expect_equal(g.pruned$a,c(0,0))
})