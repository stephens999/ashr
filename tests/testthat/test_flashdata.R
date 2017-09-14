test_that("PosteriorMean is equal to weighted component-specific posterior means", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1, outputlevel=4); expect_equal(z.ash$res$PosteriorMean, colSums(z.ash$flash_data$comp_postprob * z.ash$flash_data$comp_postmean)) 
})

test_that("PosteriorSD squared is equal to weighted component mean2 minus posterior mean all squared", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1, outputlevel=4,mixcompdist="normal")
  expect_equal(z.ash$res$PosteriorSD^2, colSums(z.ash$flash_data$comp_postprob * z.ash$flash_data$comp_postmean2)-z.ash$res$PosteriorMean^2) 
})

test_that("comp_postmean for null component is 0", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1, outputlevel=4,mixcompdist="normal")
  expect_equal(z.ash$flash_data$comp_postmean[1,],rep(0,100)) 
})

test_that("comp_postmean for missing data is 0", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(c(NA,z),1, outputlevel=4)
  expect_true(all(z.ash$flash_data$comp_postmean[,1]==0)) 
})

test_that("comp_postprob for missing data is fitted pi", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(c(NA,z),1, outputlevel=4)
  expect_true(all(z.ash$flash_data$comp_postprob[,1]==z.ash$fitted_g$pi)) 
})

test_that("prior is returned with correct components for flash data", {
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1,method="fdr",outputlevel=5,pi_thresh=1e-6)
  expect_true(length(z.ash$flash_data$prior)==ncomp(z.ash$flash_data$fitted_g)) 
})
