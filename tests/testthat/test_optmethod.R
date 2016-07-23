test_that("optmethod switches to EM when given negative probabilities", {
	set.seed(777)
        betahat <- c(0.7775, 0.6091, 0.6880, 0.6897, 0.6565, 0.7505,
                     0.7125, 0.7201, 0.7498, 0.7553)
        sebetahat <- rep(0.01, length = length(betahat))
        ash_out <- ash(betahat = betahat, sebetahat = sebetahat,
                             mixcompdist = "uniform", outputlevel=4)
        expect_true(min(ash_out$fitted.g$pi) > -10 ^ -12) 
        expect_true(ash_out$fit_details$optmethod == "mixEM")
})

test_that("control is passed to optmethod correctly when method is mixIP", {
  set.seed(1); z=rnorm(10,0,2) 
  z.ash= ash(z,1,control=list(rtol=1e-1),outputlevel=4)    
  expect_true(z.ash$fit_details$optreturn$control$rtol==1e-1)
  z.ash= ash(z,1,outputlevel=4)
  expect_true(z.ash$fit_details$optreturn$control$rtol==1e-6)
})
