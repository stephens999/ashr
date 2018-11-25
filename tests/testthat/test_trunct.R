context("ashr, context unclear")

test_that("my_e2truncnorm matches simulations", {
  set.seed(1); x = rnorm(1000000)
  expect_equal(mean(x[abs(x)<1]^2),my_e2truncnorm(-1,1),tolerance=0.01)
})

test_that("comp_postmean2 is 0 for null", {
  set.seed(1)
  z=rnorm(10)
  expect_equal(comp_postmean2(unimix(1,0,0),set_data(z,rep(1,10))),matrix(0,10,nrow=1))
  expect_equal(comp_postmean2(unimix(1,0,0),set_data(z,rep(1,10),lik_t(4))), matrix(0,10,nrow=1))
  expect_equal(comp_postmean2(normalmix(1,0,0),set_data(z,rep(1,10))),matrix(0,10,nrow=1))
  z.ash = ash(z,1,df=4,g=unimix(1,0,0),fixg=TRUE,outputlevel=3)
  expect_equal(z.ash$res$PosteriorSD,rep(0,10))
  expect_equal(z.ash$res$PosteriorMean,rep(0,10))
})

test_that("comp_postmean2.unimix matches simulations", {
  bhat = 3
  s = 2
  x = bhat+s*rt(100000,df=3)
  m = c(mean(x[x<2 & x>0]),mean(x[x<2 & x>1]),mean(x[x<0 & x>(-2)]))
  m2 = cbind(mean(x[x<2 & x>0]^2),mean(x[x<2 & x>1]^2),mean(x[x<0 & x>(-2)]^2))
  g= unimix(c(0.5,0.2,0.3),c(0,1,-2),c(2,2,0))
  temp2=as.vector(comp_postmean2(g,set_data(3,2,lik_t(3))))
  temp = as.vector(comp_postmean(g,set_data(3,2,lik_t(3))))
  expect_equal(mean((temp2-m2)^2), 0,tolerance=0.01)
  expect_equal(mean((temp-m)^2), 0,tolerance=0.01)
})

test_that("posterior means and sds computed for unimix from very flat prior are correct", {
  set.seed(1); z = rnorm(10,0,2); s=rgamma(10,10,10)
  #fit under t likelihood
  z.ash=ash(z,s,df=5,g=unimix(c(0.5,0.5),c(-100,-20),c(100,20)),fixg=TRUE,outputlevel=3)
  expect_equal(z.ash$res$PosteriorSD,s*sd(rt(1000000,df=5)),tolerance=0.01)
  #now do normal version
  z.ash=ash(z,s,df=NULL,g=unimix(c(0.5,0.5),c(-100,-20),c(100,20)),fixg=TRUE)
  expect_equal(z.ash$res$PosteriorSD,s,tolerance=0.01)
})

test_that("posterior means and sds computed for unimix from NAs match prior mean and sd", {
  set.seed(1); z = c(NA,rnorm(10,0,2)); s=c(rgamma(10,10,10),NA)
  z.ash=ash(z,s,df=5,g=unimix(c(0.5,0.5),c(-100,-20),c(100,20)),fixg=TRUE,outputlevel=3)
  priorsd = sd(c(runif(1000000,-100,100),runif(1000000,-20,20)))
  expect_equal(z.ash$res$PosteriorMean[1],0)
  expect_equal(z.ash$res$PosteriorMean[11],0)
  expect_equal(z.ash$res$PosteriorSD[1],priorsd,tolerance=0.01)
  expect_equal(z.ash$res$PosteriorSD[11],priorsd,tolerance=0.01)
})
