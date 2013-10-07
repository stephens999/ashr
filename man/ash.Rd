\name{ash}
\alias{ash}
\title{Main Adaptive SHrinkage function}
\usage{
  ash(betahat, sebetahat, mixcompdist = "normal",
    nullcheck = TRUE, df = NULL, randomstart = FALSE,
    usePointMass = FALSE, onlylogLR = FALSE,
    localfdr = TRUE, prior = NULL, sigmaavec = NULL,
    VB = FALSE)
}
\arguments{
  \item{betahat,}{a p vector of estimates}

  \item{sebetahat,}{a p vector of corresponding standard
  errors}

  \item{mixcompdist:}{distribution of components in mixture
  ("normal", "uniform" or "halfuniform")}

  \item{nullcheck:}{whether to check that any fitted model
  exceeds the "null" likelihood in which all weight is on
  the first component}

  \item{df:}{appropriate degrees of freedom for (t)
  distribution of betahat/sebetahat}

  \item{randomstart:}{bool, indicating whether to
  initialize EM randomly}

  \item{usePointMass:}{bool, indicating whether to use a
  point mass at zero as one of components for a mixture
  distribution}

  \item{onlylogLR}{(= FALSE) : bool, indicating whether to
  use this function to get logLR. Skip posterior prob,
  posterior mean, localfdr...}

  \item{localfdr}{(=TRUE) : bool, indicating whether to
  compute localfdr, localfsr, and q-value}

  \item{prior=NULL:}{vector of Dirichlet prior on mixture
  proportions (defaults to 1,1/k-1,...,1/k-1 to put more
  weight on first component)}

  \item{sigmaavec=NULL:}{vector of parameters for
  underlying mixture components}

  \item{VB=FALSE:}{whether to use Variational Bayes to
  estimate mixture proportions (instead of EM to find MAP
  estimate)}
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

