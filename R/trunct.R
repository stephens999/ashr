
#' @title my_e2trunct
#' @description Compute second moment of the truncated t. Uses results from O'Hagan, Biometrika, 1973
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param df degree of freedom of error distribution
#' @export
my_e2trunct = function(a,b,df){
  etrunct::e_trunct(a,b,df,r=2)
}

#' @title my_etrunct
#' @description Compute second moment of the truncated t. Uses results from O'Hagan, Biometrika, 1973
#'
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param df degree of freedom of error distribution
#' @export
my_etrunct = function(a,b,df){
  etrunct::e_trunct(a,b,df,r=1)
}

