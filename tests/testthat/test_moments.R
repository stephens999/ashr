context("postmean2 and postsd works")

test_that("postsd.default and postmean2.default don't return NaN's or negative values", {
    temp <- readRDS("error_dat.Rds")
    expect_false(any(is.nan(postsd.default(m = temp$m, betahat = temp$betahat,
                                           sebetahat = temp$sebetahat, v = temp$v))))
    expect_false(any(postmean2(m = temp$m, betahat = temp$betahat,
                               sebetahat = temp$sebetahat, v = temp$v) < 0))
}
)
