# a wrapper function that estimates the mode, using optim
# called by ash if mode="estimate"
ash.estmode = function(betahat,...){
  test.op = function(c){return(-ash(betahat=betahat,mode=c,outputlevel="loglik",...)$loglik)}
  opt = stats::optimize(test.op,interval=c(min(betahat),max(betahat)))
  return(opt$minimum)
}


