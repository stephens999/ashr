#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

List fixptfn(NumericVector pi_est,NumericMatrix matrix_lik, NumericVector prior){
    int n=matrix_lik.nrow(), k=matrix_lik.ncol();
    NumericVector pi_new(k);
    for (int i=0;i<k;i++){
        pi_est[i]=std::max(0.0,pi_est[i]);
    }
    pi_est = pi_est/sum(pi_est); //normalize pi
    
    double loglik,lpriordens=0.0;
    NumericMatrix m(n,k);
    NumericMatrix classprob(m);
    NumericVector m_rowsum(n);
    
    for (int i=0;i<k;i++){
        m.column(i)=pi_est[i]*matrix_lik.column(i);
        m_rowsum=m_rowsum+m.column(i);
    }
    for (int i=0;i<k;i++){
        classprob.column(i)=classprob.column(i)/m_rowsum;
    }
    //calculating objective value--probability
    loglik=sum(log(m_rowsum));
    
    for (int i=0;i<k;i++){
        if(prior[i]!=1.0){
            lpriordens +=(prior[i]-1.0)*log(pi_est[i]);
        }
        
    }
    //generating new pi
    for (int i=0;i<k;i++){//set any estimates that are less than zero, which can happen with prior<1, to 0
        pi_new[i]=std::max(0.0,sum(classprob.column(i))+prior[i]-1.0);
    }
    
    pi_new=pi_new/sum(pi_new); //normalize pi
    
    return(List::create(Named("fixedpointvector")=pi_new,
                        Named("objfn")=-loglik-lpriordens));
}



List squarem1(NumericVector par,NumericMatrix matrix_lik,NumericVector prior,List control){
    //control variables
    int maxiter=control["maxiter"];
    int method_=control["method"];
    bool trace=control["trace"];
    double stepmin0=control["step.min0"];
    double stepmax0=control["step.max0"];
    //double kr=control["kr"];
    double objfninc=control["objfn.inc"];
    double tol=control["tol"];
    double mstep=control["mstep"];
    
    List plist,p1list,p2list,p3list;
    NumericVector loldcpp,lnewcpp,pcpp,p1cpp,p2cpp,pnew,p;
    NumericVector q1,q2,sr2,sq2,sv2,srv;
    NumericVector lold,lnew;
    
    double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;
    int iter,feval;
    bool conv,extrap;
    stepmin=stepmin0;
    stepmax=stepmax0;
    
    if(trace){Rcout<<"Squarem-1"<<std::endl;}
    iter=1;
    p=par;
    
    try{p1list=fixptfn(p,matrix_lik,prior);feval=1;}
    catch(...){
        Rcout<<"Error in fixptfn function evaluation";
        return 1;
    }
    p=p1list["fixedpointvector"];
    lold=p1list["objfn"];
    if(trace){Rcout<<"Objective fn: "<<lold[0]<<std::endl;}
    conv=true;
    
    p=par;
    
    while(feval<maxiter){
        //Step 1
        extrap = true;
        pcpp=p;
        try{p1list=fixptfn(p,matrix_lik,prior); feval++;}
        catch(...){
            Rcout<<"Error in fixptfn function evaluation";
            return 1;
        }
        p1cpp=p1list["fixedpointvector"];
        q1=p1cpp-pcpp;//
        sr2=q1*q1;
        sr2_scalar=0.0;
        for (int i=0;i<sr2.length();i++){sr2_scalar+=sr2[i];}
        if(sqrt(sr2_scalar)<tol){break;}
        
        
        //Step 2
        try{p2list=fixptfn(p1cpp,matrix_lik,prior);feval++;}
        catch(...){
            Rcout<<"Error in fixptfn function evaluation";
            return 1;
        }
        p2cpp=p2list["fixedpointvector"];
        lold=p2list["objfn"];
        q2=p2cpp-p1cpp;
        sq2=q2*q2;
        sq2_scalar=0;
        for (int i=0;i<q2.length();i++){sq2_scalar+=sq2[i];}
        sq2_scalar=sqrt(sq2_scalar);
        
        
        if (sq2_scalar<tol){break;}
        sv2=q2-q1;
        sv2_scalar=0;
        for (int i=0;i<sv2.length();i++){sv2_scalar+=sv2[i]*sv2[i];}
        srv_scalar=0;
        for (int i=0;i<q2.length();i++){srv_scalar+=sv2[i]*q1[i];}
        
        
        //Step 3 Proposing new value
        switch (method_){
            case 1: alpha= -srv_scalar/sv2_scalar;
            case 2: alpha= -sr2_scalar/srv_scalar;
            case 3: alpha= sqrt(sr2_scalar/sv2_scalar);
                //default: {
                //    Rcout<<"Misspecification in method, when K=1, method should be either 1, 2 or 3!";
                //    break;}
        }
        
        alpha=std::max(stepmin,std::min(stepmax,alpha));
        pnew = pcpp + 2.0*alpha*q1 + alpha*alpha*(q2-q1);
        
        //Step 4 stabilization
        if(std::abs(alpha-1)>0.01){
            try{p3list=fixptfn(pnew,matrix_lik,prior);feval++;}
            catch(...){
                pnew=p2cpp;
                lnew=p2list["objfn"];
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
                if(alpha==stepmax){stepmax=mstep*stepmax;}
                if(stepmin<0.0 && alpha==stepmin){stepmin=mstep*stepmin;}
                p=pnew;
                lnewcpp=lnew;
                if(!std::isnan(lnewcpp[0])){lold=lnew;}
                if(trace){Rcout<<"Objective fn: "<<lnewcpp[0]<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
                iter++;
                continue;//next round in while loop
            }
            pnew=p3list["fixedpointvector"];
            lnew=p3list["objfn"];
            lnewcpp=lnew;
            if (lnewcpp[0]>loldcpp[0]+objfninc) {
                pnew=p2list["fixedpointvector"];
                lnew=p2list["objfn"];
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
            }
        }else{//same as above, when stablization is not performed.
            lnew=lold;
            lnewcpp=lnew;
            if (lnewcpp[0]>loldcpp[0]+objfninc) {
                pnew=p2list["fixedpointvector"];
                lnew=p2list["objfn"];
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
            }
        }
        
        if(alpha==stepmax){stepmax=mstep*stepmax;}
        if(stepmin<0 && alpha==stepmin){stepmin=mstep*stepmin;}
        
        p=pnew;
        lnewcpp=lnew;
        if(!std::isnan(lnewcpp[0])){lold=lnew;}
        loldcpp=lold;
        if(trace){Rcout<<"Objective fn: "<<lnewcpp[0]<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
        
        iter++;
    }
    
    if (feval >= maxiter){conv=false;}
    
    return(List::create(Named("par")=p,
                        Named("value.objfn")=lold,
                        Named("iter")=iter,
                        Named("fpevals")=feval,
                        Named("objfevals")=feval,
                        Named("convergence")=conv));
}

//' @title Brief description of function.
//' @description Explain here what this function does.
//' @param matrix_lik Description of argument goes here.
//' @param prior Description of argument goes here.
//' @param pi_init Description of argument goes shere.
//' @param control Description of argument goes here.
//' @export
// [[Rcpp::export]]
List cxxMixSquarem (NumericMatrix matrix_lik,
		    NumericVector prior, NumericVector pi_init,
		    List control){ //note: no default pi_init=NULL
    int  k=matrix_lik.ncol(),niter;
    bool converged=NA_LOGICAL;
    double loglik;
    List res;
    NumericVector pi(k);
    
    if(Rf_isNull(pi_init))
        std::fill(pi.begin(), pi.end(), 1./(double)k);
    else{
        pi=clone(pi_init);
        for (int i=0;i<k;i++)//set any estimates that are very small to be very small
            pi[i]=std::max(0.0, pi[i]);
        pi=pi/sum(pi); //normalize pi
    }
    
    res=squarem1(pi,matrix_lik,prior,control);
    pi=res["par"];
    loglik=res["value.objfn"];
    niter=res["iter"];
    converged=res["convergence"];
    pi=pi/sum(pi); //normalize pi again
    return(List::create(Named("pihat")=pi,
                        Named("B")=loglik,
                        Named("niter")=niter,
                        Named("converged")=wrap(converged)));
}
