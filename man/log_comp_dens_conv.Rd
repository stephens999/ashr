% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mix.R
\name{log_comp_dens_conv}
\alias{log_comp_dens_conv}
\title{log_comp_dens_conv}
\usage{
log_comp_dens_conv(m, data)
}
\arguments{
\item{m}{mixture distribution with k components}

\item{data}{details depend on the model}
}
\value{
a k by n matrix of log densities
}
\description{
compute the log density of the components of the
    mixture m when convoluted with a normal with standard deviation
    s or a scaled (se) student.t with df v, the density is
    evaluated at x
}
