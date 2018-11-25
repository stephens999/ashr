context("ashr moments")

test_that("postsd.default and postmean2.default don't return NaN's or negative values", {
    temp <- readRDS("error_dat.Rds")
    data = set_data(temp$betahat,temp$sebetahat,lik_normal(),alpha=1)
    expect_false(any(is.nan(postsd.default(m = temp$m, data))))
    expect_false(any(postmean2(m = temp$m, data) < 0))
}
)
