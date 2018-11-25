context("ashr \"df\" argument")

test_that("df switches", {
    betahat <- c(1.01636974224394, -2.05686254738995, -0.7135781676358,
                 -1.16906745227838, -0.917039991627176)

    sebetahat <- c(1.02572223086898, 0.499285201440522, 0.476520330150983,
                   0.624576594477857, 0.198152636610839)

    aout <- ash.workhorse(betahat = betahat[1:5], sebetahat = sebetahat[1:5], df = Inf)
    expect_true(all(!is.nan(aout$result$PosteriorMean)))
}
)
