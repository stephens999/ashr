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

test_that("log1mexp matches log(1 - exp(z)) on safe range and stays finite near 0", {
  # On the safe range, log1mexp should match the naive log(1 - exp(z))
  # to machine precision.
  z = c(-10, -5, -2, -1, -log(2), -0.5, -0.1, -1e-3, -1e-6, -1e-9)
  expect_equal(ashr:::log1mexp(z), log(1 - exp(z)))

  # When |z| is below ~1e-16, the naive form returns -Inf because exp(z)
  # rounds to 1. log1mexp should stay finite and close to log(|z|).
  z_small = c(-1e-17, -1e-18, -1e-20)
  out = ashr:::log1mexp(z_small)
  expect_true(all(is.finite(out)))
  expect_equal(out, log(-z_small), tolerance = 1e-12)

  # NaN / NA inputs must propagate as NaN / NA respectively (not get
  # silently converted to NA by an internal comparison).
  expect_true(is.nan(ashr:::log1mexp(NaN)))
  expect_true(is.na(ashr:::log1mexp(NA_real_)))
})

test_that("logscale_sub matches naive log(exp - exp) on safe range and stays finite when logx ~ logy", {
  # On the safe range (logx > logy by a comfortable margin), logscale_sub
  # should match the naive log(exp(logx) - exp(logy)) to high precision.
  logx = c(  0,  -5,  10, -50, -100)
  logy = c(-1,  -7,   8, -55, -200)
  ref  = log(exp(logx) - exp(logy))
  expect_equal(ashr:::logscale_sub(logx, logy), ref, tolerance = 1e-12)

  # When logx and logy are very close, the naive form returns -Inf
  # because exp(logy - logx) rounds to 1. logscale_sub should stay
  # finite and equal logx + log(logx - logy) to first order.
  # We anchor at logx = 0 so that very small (logx - logy) values
  # remain representable in double precision.
  logx_close = c(0, 0, 0)
  logy_close = c(-1e-12, -1e-15, -1e-17)
  out = ashr:::logscale_sub(logx_close, logy_close)
  expect_true(all(is.finite(out)))
  expect_equal(out, logx_close + log(logx_close - logy_close), tolerance = 1e-6)

  # Corner case used by my_etruncnorm with (a = -Inf, b = +Inf):
  # dnorm(+-Inf, log = TRUE) = -Inf, so logscale_sub is called with
  # both arguments equal to -Inf. The mathematically correct answer
  # is log(0 - 0) = -Inf, not NaN.
  expect_equal(ashr:::logscale_sub(-Inf, -Inf), -Inf)
  expect_equal(ashr:::logscale_sub(c(-Inf, 0), c(-Inf, -1)),
               c(-Inf, log(1 - exp(-1))))
})

test_that("log_comp_dens_conv.unimix limit as b -> a matches normal density / s", {
  # As (b - a) -> 0 with a normal likelihood, the convolution density
  # approaches (1/s) * dnorm((x - a)/s). Check this in the tail, where
  # the old log(1 - exp(.)) form loses precision.
  # NOTE: log_comp_dens_conv.unimix is not registered as an S3 method
  # in NAMESPACE, so we call the method directly instead of going
  # through UseMethod dispatch.
  x = 0
  s = 1
  a = 5
  ref = dnorm(x - a, 0, s, log = TRUE) - log(s)
  for (eps in c(1e-3, 1e-6, 1e-9)) {
    g = unimix(1, a, a + eps)
    data = set_data(x, s)
    lcd = as.numeric(ashr:::log_comp_dens_conv.unimix(g, data))
    expect_true(is.finite(lcd))
    expect_equal(lcd, ref, tolerance = 1e-3)
  }
})
