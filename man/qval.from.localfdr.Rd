\name{qval.from.localfdr}
\alias{qval.from.localfdr}
\title{Function to compute q values from local false discovery rates}
\usage{
  qval.from.localfdr(localfdr)
}
\arguments{
  \item{localfdr,}{a vector of local fdr estimates}
}
\value{
  vector of q values
}
\description{
  Computes q values from a vector of local fdr estimates
}
\details{
  The q value for a given localfdr is an estimate of the
  (tail) False Discovery Rate for all findings with a
  smaller localfdr, and is found by the average of the
  localfdr for all more significant findings. See Storey
  (2003), Annals of Statistics, for definition of q value.
}

