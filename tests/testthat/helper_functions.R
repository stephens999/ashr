# The Rmosek package on CRAN will not work with REBayes. This function
# is used for some of the tests to check whether the correct Rmosek
# package (the one downloaded from mosek.com) is installed.
#
#' @importFrom testthat skip_if_not_installed
#' @importFrom testthat skip_if
#' 
skip_if_mixkwdual_doesnt_work <- function() {
  skip_if_not_installed("REBayes")
  skip_if_not_installed("Rmosek")
  skip_if(!is.element("mosek_lptoprob",getNamespaceExports("Rmosek")))
}
