\name{ash}
\alias{ash}
\title{Main Adaptive SHrinkage function}
\usage{
  ash(betahat, sebetahat, mixcompdist = "normal",
    nullcheck = TRUE, df = NULL, randomstart = FALSE,
    usePointMass = FALSE, onlylogLR = FALSE,
    localfdr = TRUE, localfsr = TRUE, prior = NULL,
    sigmaavec = NULL, auto = FALSE, sigma.est = FALSE,
    nc = NULL, VB = FALSE)
}
\arguments{
  \item{betahat}{(p vector);}

  \item{sebetahat}{(p vector of standard errors)}

  \item{mixcompdist:}{distribution of components in mixture
  ("normal", "uniform" or "halfuniform")}

  \item{df:}{degrees of freedome used to compute sebetahat}

  \item{randomstart:}{bool, indicating whether to
  initialize EM randomly}

  \item{usePointMass:}{bool, indicating whether to use a
  point mass at zero as one of components for a mixture
  distribution}

  \item{onlylogLR}{(= FALSE) : bool, indicating whether to
  use this function to get logLR. Skip posterior prob,
  posterior mean, localfdr...}

  \item{localfdr}{(=TRUE) : bool, indicating whether to
  compute localfdr and q-value}

  \item{auto}{(=FALSE): bool, whether to try to select the
  sigmaavec vector automatically (beta functionality)}

  \item{sigma.est:}{bool, whether to estimate sigma rather
  than fixing (Beta version)}

  \item{nc:}{number of components to use (only relevant
  when sigma.est=TRUE)}
}
\value{
  a list with elements fitted.g is fitted mixture logLR :
  logP(D|mle(pi)) - logP(D|null)
}
\description{
  takes a vector of betahats and ses; fits a mixture of
  normals to it and returns posteriors
}
\details{
  See readme for more details
}

