[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ashr)](https://cran.r-project.org/package=ashr)
[![Build Status](https://travis-ci.org/stephens999/ashr.svg)](https://travis-ci.org/stephens999/ashr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/stephens999/ashr?branch=master&svg=true)](https://ci.appveyor.com/project/stephens999/ashr)
[![Coverage Status](https://coveralls.io/repos/github/stephens999/ashr/badge.svg?branch=master)](https://coveralls.io/github/stephens999/ashr?branch=master)
[![Coverage Status](https://img.shields.io/codecov/c/github/stephens999/ashr/master.svg)](https://codecov.io/github/stephens999/ashr?branch=master)

This repository contains an R package for performing "Adaptive Shrinkage."

To install the ashr package first you need to install devtools:

```R
install.packages("devtools")
library(devtools)
install_github("stephens999/ashr")
```


## Running Adaptive Shrinkage

The main function in the ashr package is `ash`. To get minimal help:

```R
library(ashr)
?ash
```

## More background

The ashr ("Adaptive SHrinkage") package aims to provide simple,
generic, and flexible methods to derive "shrinkage-based" estimates
and credible intervals for unknown quantities
$\beta=(\beta_1,\dots,\beta_J)$, given only estimates of those
quantities ($\hat\beta=(\hat\beta_1,\dots, \hat\beta_J)$) and their
corresponding estimated standard errors ($s=(s_1,\dots,s_J)$).

The "adaptive" nature of the shrinkage is two-fold. First, the
appropriate amount of shrinkage is determined from the data, rather
than being pre-specified. Second, the amount of shrinkage undergone by
each $\hat\beta_j$ will depend on the standard error $s_j$:
measurements with high standard error will undergo more shrinkage than
measurements with low standard error.

### Methods Outline

The methods are based on treating the vectors $\hat\beta$ and $s$ as
"observed data", and then performing inference for $\beta$ from these
observed data, using a standard hierarchical modelling framework to
combine information across $j=1,\dots,J$.

Specifically, we assume that the true $\beta_j$ values are independent
and identically distributed from some unimodal distribution $g$.  By
default we assume $g$ is unimodal about zero and symmetric.  You can
specify or estimate a different mode using the `mode` parameter.  You
can allow for asymmetric $g$ by specifying
`mixcompdist="halfuniform"`.

Then, we assume that the observations $\hat\beta_j \sim
N(\beta_j,s_j)$, or alternatively the normal assumption can be
replaced by a $t$ distribution by specifying `df`, the number of
degrees of freedom used to estimate $s_j$.  Actually this is
important: do be sure to specify `df` if you can.
