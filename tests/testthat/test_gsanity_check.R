context("sanity check for initial g")

test_that("The sanity check on g passes and fails as expected", {
  set.seed(666)
  
  # Should pass when g is null and the data comes from the null model.
  g = normalmix(pi = 1, mean = 0, sd = 0)
  data = list(x = rnorm(1000), s = rep(1, 1000), lik = lik_normal())
  expect_true(gsanity_check(data, g))
  
  # Shouldn't pass when the data clearly does not come from the null model.
  data$x = rnorm(1000) + rnorm(1000)
  expect_false(gsanity_check(data, g))
  
  # Passing the results of ashr back in should always work.
  ash.u = ash(1:100, 1, mixcompdist = "uniform")
  expect_true(gsanity_check(ash.u$data, ash.u$fitted_g))
  
  # But it shouldn't work if g is constrained beyond a reasonable range.
  constrained.g = constrain_mix(ash.u$fitted_g, rep(1, ncomp(ash.u$fitted_g)), 
                                c(0, 90), "uniform")$g
  expect_false(gsanity_check(ash.u$data, constrained.g))
  
  # Check other mixcompdists.
  ash.hu = ash(1:10, 1, mixcompdist = "halfuniform")
  expect_true(gsanity_check(ash.hu$data, ash.hu$fitted_g))
  ash.hn = ash(-3:10, 1, mixcompdist = "halfuniform")
  expect_true(gsanity_check(ash.hn$data, ash.hn$fitted_g))
  ash.uplus = ash(-5:10, 1, mixcompdist = "+uniform")
  expect_true(gsanity_check(ash.uplus$data, ash.uplus$fitted_g))
  ash.uminus = ash(-10:1, 1, mixcompdist = "-uniform")
  expect_true(gsanity_check(ash.uminus$data, ash.uminus$fitted_g))
})