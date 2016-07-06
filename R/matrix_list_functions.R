# split a matrix into a list, by rows
matrix2list_byrows =function(x){
  lapply(1:nrow(x),function(i){x[i,]})
}

#given a list of functions, applies the ith function to the ith row of x
#if the matrix is n by k, and the function returns a k vector when given 
#a k vector then returned matrix is n by k
#(If just a single function, or a list with one function, applies it to the matrix)
apply_functions_to_rows_of_matrix = function(flist,x){
  if(!is.list(flist)){flist = list(flist)} #if a single funciton, put in list
  if(length(flist)==1){return(flist[[1]](x))}
  l = matrix2list_byrows(x)
  t(mapply(function(f,x){f(x)},flist,l))
}

#given a list of functions, applies the ith function to the ith row of x and ith row of y
#where x and y are concordant matrices
# if x is n by k, and the function maps R^k -> R^k then returned matrix is n by k
#(If just a single function, or a list with one function, applies it to the matrix)
apply_functions_to_rows_of_two_matrices = function(flist,x,y){
  if(!is.list(flist)){flist = list(flist)} #if a single funciton, put in list
  if(length(flist)==1){return(flist[[1]](x,y))}
  lx = matrix2list_byrows(x)
  ly = matrix2list_byrows(y)
  t(mapply(function(f,x,y){f(x,y)},flist,lx,ly))
}

#given a list of functions, applies the ith function to the ith element of x
#(if just a single function, or a list with one function, applies it to the vector)
apply_functions_to_elements_of_vector = function(flist,x){
  if(!is.list(flist)){flist=list(flist)}
  if(length(flist)==1){return(flist[[1]](x))} #if a single function, just apply to x
  mapply(function(f,x){f(x)},flist,as.list(x))
}
