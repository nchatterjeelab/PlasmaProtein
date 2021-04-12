#include <cstdlib>
#include <time.h>
#include <string>
#include <iostream>// std::cout
#include <cmath>
#include <algorithm>// std::min
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;

//--------------------------------
  //--------------------------------
  // Construct path

//' Construct lambda path for single ethnic group
//'
//' @param summ summary statistics
//' @param s tuning parameter for L1/L2
//' @param nlambda number of lambda in the path
//' @param lambda_min_ratio the ratio of the maximum and minimum of lambdas
//' @return lambda path
//' @keywords internal
// [[Rcpp::export]]
vec l_path(vec summ, double alpha=0.5, int nlambda=20, double lambda_min_ratio = 0.1){
  nlambda += 1;
  double loglambda_max = log( max(abs(summ))/alpha);
  double loglambda_min = log( max(abs(summ))/alpha * lambda_min_ratio);
  vec lambdapath(nlambda-1);
  for (int i=1; i < nlambda; i++) {
    lambdapath(i-1) = exp(loglambda_max - (loglambda_max-loglambda_min)/(nlambda-1)  * i);
  }
  return(lambdapath);
}

//' Construct c path for multi ethnic group
//'
//' @param lambda_max1 maximum lambda1
//' @param lambda_max2 maximum lambda2
//' @param nlambda number of lambda3 in the path
//' @return lambda path for lambda3 (smalles is 0)
//' @keywords internal
//'
// [[Rcpp::export]]
vec c_path(double maxc, int nlambda=10){
  
  vec lambdapath(nlambda);
  for (int i=0; i < nlambda; i++) {
    lambdapath(i) = maxc/nlambda * (i+1);
  }
  return(lambdapath);
  
}


//--------------------------------
  //--------------------------------
  // Single-ethic analysis

//' Compute beta
//'
//' @param summ summary statistics
//' @param TaggingSNPinx tagging SNPs
//' @param Correlation LD of tagging SNPs
//' @param shared shared SNPs of ethnics
//' @return beta vector and computation details
//' @keywords internal
// [[Rcpp::export]]
List enet_singlethnic_bl(vec summ,
                         mat R,
                         double lambda,
                         double alpha=0.5,
                         double thresh=1e-04, int maxiter=1000){
  
  int p1 = summ.n_elem;
  double dlx, del, tmp, ui;
  vec b(p1); b.fill(0.0);
  
  double denom0 = 1.0 + lambda*(1.0-alpha);
  double Lambda = lambda*alpha;
  
  int conv=0;
  int niter;
  for(int k=0;k<maxiter ;k++) {
    dlx=0.0;
    
    for (int i=0; i<p1; i++) {
      
      tmp = b(i);
      
      ui = summ(i) - (dot(R.col(i), b)-b(i));
      b(i) = sign(ui) * std::max(0.0, abs(ui)-Lambda)/denom0;
      
      if(b(i)==tmp){
        continue;
      }else{
        del=b(i)-tmp;
        dlx=std::max(dlx,abs(del));
      }
    }
    
    if(dlx < thresh) {
      conv=1;
      niter=k;
      break;
    }
    
  }
  if(conv==0){
    niter=maxiter;
    b.fill(0.0);
  }

  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda") = lambda,
                      Named("alpha") = alpha,
                      Named("b") = b);

}



//' Compute beta joint tuning
//'
//' @param summ summary statistics
//' @param TaggingSNPinx tagging SNPs
//' @param Correlation LD of tagging SNPs
//' @param shared shared SNPs of ethnics
//' @return beta vector and computation details
//' @keywords internal
// [[Rcpp::export]]
List enet_singlethnic(vec summ,
                      mat R,
                      int nlambda,
                      double alpha=0.5,
                      double lambda_min_ratio = 0.1,
                      double thresh=1e-04, int maxiter=1000){

  vec Lambda = l_path(summ, alpha, nlambda, lambda_min_ratio);

  int M = summ.n_elem;

  int alltuning = nlambda;

  vec conv(alltuning);
  vec niter(alltuning);
  mat b(M, alltuning);
  vec lambda(alltuning);

  int ii=0;
  double li;

  for(int l=0; l<nlambda ;l++) {

    li = Lambda[l];

    List res_fixtuning = enet_singlethnic_bl( summ, R, li, alpha, thresh,  maxiter);

    lambda(ii) = li;

    conv(ii) = as<int>(res_fixtuning["conv"]);
    niter(ii) = as<int>(res_fixtuning["niter"]);
    b.col(ii) = as<vec>(res_fixtuning["b"]);

    Rcout << "*** Completed " << ii << "/" << alltuning << " tuning; Niter = " << niter(ii) << std::endl;

    ii += 1;


  }



  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda") = lambda,
                      Named("b") = b);
}



//--------------------------------
  //--------------------------------
  // Multi-ethic analysis


//' Compute beta joint
//'
//' @param summ summary statistics
//' @param TaggingSNPinx tagging SNPs
//' @param Correlation LD of tagging SNPs
//' @param shared shared SNPs of ethnics
//' @return beta vector and computation details
//' @keywords internal
// [[Rcpp::export]]
List enet_multiethnic_bl(vec summ1, vec summ2,
                         mat R1, mat R2,
                         vec shared1, vec shared2,
                         double lambda1, double lambda2, double c,
                         double alpha1=0.5, double alpha2=0.5,
                         double thresh=1e-04, int maxiter=1000){

  int p1 = summ1.n_elem;
  int p2 = summ2.n_elem;
  double dlx, del, tmp, ui;
  int tmp0;
  vec b1(p1); b1.fill(0.0);
  vec b2(p2); b2.fill(0.0);


  double denom10 = 1.0 + lambda1*(1.0-alpha1);
  double denom20 = 1.0 + lambda2*(1.0-alpha2);

  double denom1 = denom10 + c;
  double denom2 = denom20 + c;

  double Lambda1 = lambda1*alpha1;
  double Lambda2 = lambda2*alpha2;

  int conv=0;
  int niter;
  for(int k=0; k<maxiter ;k++) {
    dlx=0.0;

    for (int i=0; i<p1; i++) {

      tmp = b1(i);
      if(shared1(i)!=0){
        tmp0 = shared1(i)-1;
        ui = summ1(i) - (dot(R1.col(i), b1)-b1(i)) + c*b2(tmp0);
        b1(i) = sign(ui) * std::max(0.0, abs(ui)-Lambda1) / denom1;
      }else{
        ui = summ1(i) - (dot(R1.col(i), b1)-b1(i));
        b1(i) = sign(ui) * std::max(0.0, abs(ui)-Lambda1) / denom10;
      }

      if(b1(i)==tmp){
        continue;
      }else{
        del=b1(i)-tmp;
        dlx=std::max(dlx,abs(del));
      }
    }

    for (int i=0; i<p2; i++) {

      tmp = b2(i);
      if(shared2(i)!=0){
        tmp0 = shared2(i)-1;
        ui = summ2(i) - (dot(R2.col(i), b2)-b2(i)) + c*b1(tmp0);
        b2(i) = sign(ui) * std::max(0.0, abs(ui)-Lambda2) / denom2;
      }else{
        ui = summ2(i) - (dot(R2.col(i), b2)-b2(i));
        b2(i) = sign(ui) * std::max(0.0, abs(ui)-Lambda2) / denom20;
      }

      if(b2(i)==tmp){
        continue;
      }else{
        del=b2(i)-tmp;
        dlx=std::max(dlx,abs(del));
      }
    }

    if(dlx < thresh) {
      conv=1;
      niter=k;
      break;
    }

  }
  if(conv==0){
    niter=maxiter;
    b1.fill(0.0);
    b2.fill(0.0);
  }
  
  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("b1") = b1,
                      Named("b2") = b2);
}







//' Compute beta joint tuning
//'
//' @param summ summary statistics
//' @param TaggingSNPinx tagging SNPs
//' @param Correlation LD of tagging SNPs
//' @param shared shared SNPs of ethnics
//' @return beta vector and computation details
//' @keywords internal
// [[Rcpp::export]]
List enet_multiethnic(vec summ1, vec summ2,
                      mat R1, mat R2,
                      vec shared1, vec shared2,
                      int nlambda1, int nlambda2, int nc,
                      double alpha1=0.5, double alpha2=0.5,
                      double lambda_min_ratio1 = 0.1, double lambda_min_ratio2 = 0.1,
                      double thresh=1e-04, int maxiter=1000){
  
  vec Lambda1 = l_path(summ1, alpha1, nlambda1, lambda_min_ratio1);
  vec Lambda2 = l_path(summ2, alpha2, nlambda2, lambda_min_ratio2);
  vec C = c_path(nc, nc);
  
  int M1 = summ1.n_elem;
  int M2 = summ2.n_elem;
  
  int alltuning = nlambda1*nlambda2*nc;
  
  vec conv(alltuning);
  vec niter(alltuning);
  
  mat b1(M1, alltuning);
  mat b2(M2, alltuning);
  
  vec lambda1(alltuning);
  vec lambda2(alltuning);
  vec c(alltuning);
  
  int ii=0;
  double l1i,l2i,ci;
  
  for(int l1=0; l1<nlambda1 ;l1++) {
    for(int l2=0; l2<nlambda2 ;l2++) {
      for(int cc=0; cc<nc ;cc++) {
        
        l1i = Lambda1[l1];
        l2i = Lambda2[l2];
        ci = C[cc];
        
        List res_fixtuning = enet_multiethnic_bl( summ1,  summ2,
                                                  R1,  R2,
                                                  shared1,  shared2,
                                                  l1i,  l2i,  ci,
                                                  alpha1,  alpha2,
                                                  thresh,  maxiter);
        
        lambda1(ii) = l1i;
        lambda2(ii) = l2i;
        c(ii) = ci;
        conv(ii) = as<int>(res_fixtuning["conv"]);
        niter(ii) = as<int>(res_fixtuning["niter"]);
        b1.col(ii) = as<vec>(res_fixtuning["b1"]);
        b2.col(ii) = as<vec>(res_fixtuning["b2"]);
        
        Rcout << "*** Completed " << ii << "/" << alltuning << " tuning; Niter = " << niter(ii) << std::endl;
        
        ii += 1;
        
      }
      
    }
    
    
  }
  
  
  
  return List::create(Named("conv") = conv,
                      Named("niter") = niter,
                      Named("lambda1") = lambda1,
                      Named("lambda2") = lambda2,
                      Named("c") = c,
                      Named("b1") = b1,
                      Named("b2") = b2);
}
