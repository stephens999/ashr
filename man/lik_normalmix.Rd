% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lik.R
\name{lik_normalmix}
\alias{lik_normalmix}
\title{Likelihood object for normal mixture error distribution}
\usage{
lik_normalmix(pilik, sdlik)
}
\arguments{
\item{pilik}{a k vector of mixture proportions (k is the number of mixture components), 
or an n*k matrix that the j'th row the is mixture proportions for betahat_j}

\item{sdlik}{a k vector of component-wise standard deviations, 
or an n*k matrix that the j'th row the is component-wise standard deviations for betahat_j}
}
\description{
Creates a likelihood object for ash for use with normal mixture error distribution
}
\examples{
   e = rnorm(100,0,0.8) 
   e[seq(1,100,by=2)] = rnorm(50,0,1.5) # generate e~0.5*N(0,0.8^2)+0.5*N(0,1.5^2)
   betahat = rnorm(100)+e
   ash(betahat, 1, lik=lik_normalmix(c(0.5,0.5),c(0.8,1.5)))
}
