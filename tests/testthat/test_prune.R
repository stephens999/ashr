test_that("test pruning ", {
  g= normalmix(c(0.5,0.4,0.1),c(0,0,0),c(1,2,3))
  g.pruned = prune(g,0.2)
  expect_equal(g.pruned$mean,c(0,0))
  g = unimix(c(0.5,0.4,0.1),c(0,0,0),c(1,2,3))
  g.pruned = prune(g,0.2)
  expect_equal(g.pruned$a,c(0,0))
}
)

test_that("test pruning with and without flash ", {
  set.seed(1)
    z = rnorm(100,0,2)
  ash.z <- ash(z,1,outputlevel=4)
  expect_equal(ncomp(get_fitted_g(ash.z)), nrow(ash.z$flash_data$comp_postmean))
}
)