# a wrapper function that estimates the mode, using optim
# called by ash if mode="estimate"
ash.estmode = function(betahat,modemin,modemax,...){
  test.op = function(c){return(-ash(betahat=betahat,mode=c,outputlevel="loglik",...)$loglik)}
  #opt = stats::optimize(test.op,interval=c(min(betahat),max(betahat)))
  opt = stats::optimize(test.op,interval=c(modemin,modemax))
  return(opt$minimum)
}


