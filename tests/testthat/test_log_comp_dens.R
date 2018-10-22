context("ashr convolved density computations")

test_that("normalmix functions behave as expected", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  data = set_data(c(1,2,3),c(1,10,100))
  expect_equal(log(comp_dens_conv(gn,data)),log_comp_dens_conv(gn,data))
})


test_that("exp(log_comp_dens_conv) gives same results as comp_dens_conv", {
  g = unimix(c(0.1,0.45,0.45),c(0,0,0),c(0,1,2))
  gn = normalmix(c(0.1,0.45,0.45),c(0,0,0),c(0,0.1,1))
#  gig = igmix(c(0.5,0.5),c(1,2),c(3,4))
  x=c(-10,2)
  s = c(1,2)
  data = set_data(x,s)
  data2 = set_data(x,s,lik_t(df=2))
  expect_equal(comp_dens_conv(g, data), exp(log_comp_dens_conv(g,data)))
  expect_equal(comp_dens_conv(g, data2), exp(log_comp_dens_conv(g,data2)))
  expect_equal(comp_dens_conv(gn, data), exp(log_comp_dens_conv(gn,data)))
  
  data = set_data(x,s,alpha = 1)
  expect_equal(comp_dens_conv(g, data), exp(log_comp_dens_conv(g,data)))
  expect_equal(comp_dens_conv(gn, data), exp(log_comp_dens_conv(gn,data)))
})

test_that("comp_postprob is numerically stable", {
  g = unimix(c(0.5,0.5),c(1,2),c(0,0))
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gig = igmix(c(0.5,0.5),c(1,2),c(3,4))
  x=c(-10,2)
  s = c(1,2)
  data = set_data(x,s)
  data2 = set_data(x,s,lik_t(2))
  expect_equal(comp_postprob(g,set_data(-10,0.5)),cbind(c(2/3,1/3)))
  expect_equal(comp_postprob(g,set_data(-20,0.5)),cbind(c(2/3,1/3)))
})
