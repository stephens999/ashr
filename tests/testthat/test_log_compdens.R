test_that("normalmix functions behave as expected", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  data = list(x=c(1,2,3),s=c(1,10,100))
  expect_equal(log(compdens_conv(gn,data)),log_compdens_conv(gn,data))
  #comp_postmean(gn,data)
  #comppostprob(gn,data)
  #compcdf_post(gn,1,data)
  #cdf_post(gn,1,data)
}
)


test_that("exp(log_compdens_conv) gives same results as compdens_conv", {
  g = unimix(c(0.1,0.45,0.45),c(0,1,2),c(0,0,0))
  gn = normalmix(c(0.1,0.45,0.45),c(0,0,0),c(0,0.1,1))
#  gig = igmix(c(0.5,0.5),c(1,2),c(3,4))
  x=c(-10,2)
  s = c(1,2)
  data = set_data(x,s,df=NULL)
  data2 = set_data(x,s,df=rep(2,2))
  expect_equal(compdens_conv(g, data), exp(log_compdens_conv(g,data)))
  expect_equal(compdens_conv(g, data2), exp(log_compdens_conv(g,data2)))
  expect_equal(compdens_conv(gn, data), exp(log_compdens_conv(gn,data)))
  
  data = set_data(x,s,df=NULL,alpha = 1)
  expect_equal(compdens_conv(g, data), exp(log_compdens_conv(g,data)))
  expect_equal(compdens_conv(gn, data), exp(log_compdens_conv(gn,data)))
#  expect_equal(compdens_conv(gig, data2), exp(log_compdens_conv(gig,data2)))
})

test_that("comppostprob is numerically stable", {
  g = unimix(c(0.5,0.5),c(1,2),c(0,0))
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gig = igmix(c(0.5,0.5),c(1,2),c(3,4))
  x=c(-10,2)
  s = c(1,2)
  data = set_data(x,s,df=NULL)
  data2 = set_data(x,s,df=rep(2,2))
#  expect_equal(comppostprob(g,data),old.comppostprob.default(g,x,s,NULL))
#  expect_equal(comppostprob(g,data2),old.comppostprob.default(g,x,s,2))
#  expect_equal(comppostprob(gn,data),old.comppostprob.default(gn,x,s,NULL))
#  expect_equal(comppostprob(gig,data2),old.comppostprob.default(gig,x,s,2))
  
  expect_equal(comppostprob(g,set_data(-10,0.5,NULL)),cbind(c(2/3,1/3)))
  expect_equal(comppostprob(g,set_data(-20,0.5,NULL)),cbind(c(2/3,1/3)))
})