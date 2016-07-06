test_that("matrix list functions work",{
  flist = list(function(x){pt(x,df=4)},pnorm)
  x = rbind(rnorm(10),rnorm(10))
  expect_equal(apply_functions_to_rows_of_matrix(flist,x),
               rbind(pt(x[1,],df=4),pnorm(x[2,])))
  expect_equal(apply_functions_to_rows_of_matrix(list(pnorm),x),
               rbind(pnorm(x[1,]),pnorm(x[2,])))
  expect_equal(apply_functions_to_elements_of_vector(flist, x[,1]),
               c(pt(x[1,1],df=4),pnorm(x[2,1])))
  etruncFUN = function(a,b){my_etrunct(a,b,df=4)}
  flist = list(etruncFUN,etruncFUN)
  x=c(1,2); a = 1:10; b=11:20
  alpha = outer(-x,a,"+")
  beta = outer(-x,b,"+")
  expect_equal(apply_functions_to_rows_of_two_matrices(flist,alpha,beta),
               etruncFUN(alpha,beta))
}
)