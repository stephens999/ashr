context("ashr with mixture-of-truncated-normal priors")

test_that("compdens_conv for tnormal gives same results as compdens_conv for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  x=c(-10,2)
  s = c(1,2)
  data = list(x=x,s=s,lik=lik_normal(),alpha=0)
  expect_equal(comp_dens_conv(gn, data),comp_dens_conv(gtn,data))
})

test_that("log_compdens_conv for tnormal gives reasonable results", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  x=c(-10,2)
  s = c(1,2)
  data = list(x=x,s=s,lik=lik_normal(),alpha=0)
  expect_equal(exp(log_comp_dens_conv(gtn,data)),comp_dens_conv(gtn,data))
})
          
test_that("comp_cdf for tnormal gives same results as comp_cdf for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  y=c(-1,-5,0.5,1)
  expect_equal(comp_cdf(gn, y),comp_cdf(gtn,y))
})

test_that("compcdf_post for tnormal gives same results as compcdf_post for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,-2,1,2)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  c = 0.5
  expect_equal(comp_cdf_post(gn,c,data),comp_cdf_post(gtn,c,data))
})

test_that("comp_postmean for tnormal gives same results as comp_postmean for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,1),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,1),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,0.4,1.3,5)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  expect_equal(comp_postmean(gn,data),comp_postmean(gtn,data))
})

test_that("comp_postsd for tnormal gives same results as comp_postsd for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,1),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,1),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,0.4,1.3,5)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  expect_equal(comp_postsd(gn,data),comp_postsd(gtn,data))
  expect_equal(comp_postmean2(gn,data),comp_postmean2(gtn,data))
})
