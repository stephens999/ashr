\name{ash}
\alias{ash}
\title{Main Adaptive SHrinkage function}
\usage{
  ash(betahat, sebetahat, method = c("shrink", "fdr"),
    mixcompdist = c("normal", "uniform", "halfuniform"),
    lambda1 = 1, lambda2 = 0, nullcheck = TRUE, df = NULL,
    randomstart = FALSE, pointmass = FALSE,
    onlylogLR = FALSE, prior = c("uniform", "nullbiased"),
    mixsd = NULL, VB = FALSE, gridmult = sqrt(2),
    minimaloutput = FALSE, g = NULL)
}
\arguments{
  \item{betahat,}{a p vector of estimates}

  \item{sebetahat,}{a p vector of corresponding standard
  errors}

  \item{method:}{specifies how ash is to be run. Can be
  "shrinkage" (if main aim is shrinkage) or "fdr" (if main
  aim is to assess fdr or fsr) This is simply a convenient
  way to specify certain combinations of parameters:
  "shrinkage" sets pointmass=FALSE and prior="uniform";
  "fdr" sets pointmass=TRUE and prior="nullbiased".}

  \item{mixcompdist:}{distribution of components in mixture
  ("normal", "uniform" or "halfuniform")}

  \item{lambda1:}{multiplicative "inflation factor" for
  standard errors (like Genomic Control)}

  \item{lambda2:}{additive "inflation factor" for standard
  errors (like Genomic Control)}

  \item{nullcheck:}{whether to check that any fitted model
  exceeds the "null" likelihood in which all weight is on
  the first component}

  \item{df:}{appropriate degrees of freedom for (t)
  distribution of betahat/sebetahat}

  \item{randomstart:}{bool, indicating whether to
  initialize EM randomly. If FALSE, then initializes to
  prior mean (for EM algorithm) or prior (for VBEM)}

  \item{pointmass:}{bool, indicating whether to use a point
  mass at zero as one of components for a mixture
  distribution}

  \item{onlylogLR:}{bool, indicating whether to use this
  function to get logLR. Skip posterior prob, posterior
  mean, lfdr...}

  \item{prior:}{string, or numeric vector indicating
  Dirichlet prior on mixture proportions (defaults to
  "uniform", or 1,1...,1; also can be "nullbiased"
  1,1/k-1,...,1/k-1 to put more weight on first component)}

  \item{mixsd:}{vector of sds for underlying mixture
  components}

  \item{VB:}{whether to use Variational Bayes to estimate
  mixture proportions (instead of EM to find MAP estimate)}

  \item{gridmult:}{the multiplier by which the default grid
  values for mixsd differ by one another. (Smaller values
  produce finer grids)}

  \item{minimal_output:}{if TRUE, just outputs the fitted g
  and the lfsr (useful for very big data sets where memory
  is an issue)}

  \item{g:}{the prior distribution for beta (usually
  estimated from the data; this is used primarily in
  simulated data to do computations with the "true" g)}
}
\value{
  a list with elements fitted.g is fitted mixture logLR :
  logP(D|mle(pi)) - logP(D|null)
}
\description{
  Takes vectors of estimates (betahat) and their standard
  errors (sebetahat), and applies shrinkage to them, using
  Empirical Bayes methods, to compute shrunk estimates for
  beta.
}
\details{
  See readme for more details
}
\examples{
beta = c(rep(0,100),rnorm(100))
sebetahat = abs(rnorm(200,0,1))
betahat = rnorm(200,beta,sebetahat)
beta.ash = ash(betahat, sebetahat)
summary(beta.ash)
plot(betahat,beta.ash$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
}

