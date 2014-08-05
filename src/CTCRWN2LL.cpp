// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//const double log2pi = std::log(2.0 * M_PI);
//double ln_dmvnrm(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) { 
//  int xdim = x.n_elem;
//  double out;
//  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
//  double rootisum = arma::sum(log(rooti.diag()));
//  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
//  arma::vec z = rooti * ( x - mean) ;    
//  out =  - 0.5 * arma::sum(z%z) + rootisum;      
//  return(out);
//} 

// [[Rcpp::export]]
arma::mat makeT(const double& b, const double& delta){
  arma::mat T(4,4, fill::zeros);
  T(0,0) = 1; T(2,2) = 1;
  T(0,1) = exp(R::pexp(delta,1/b,1,1) - log(b));
  T(1,1) = exp(-b*delta);
  T(2,3) = exp(R::pexp(delta,1/b,1,1) - log(b));
  T(3,3) = exp(-b*delta);
  return T;
}

// [[Rcpp::export]]
arma::mat makeQ(const double& b, const double& sig2, const double& delta){
  arma::mat Q(4,4, fill::zeros);
  Q(0,0) =  sig2*(delta - 2*exp(R::pexp(delta,1/b,1,1)-log(b))  + exp(R::pexp(delta,1/(2*b),1,1)-log(2*b)));
  Q(1,1) = sig2*exp(log(b) + R::pexp(delta,1/(2*b),1,1))/2;
  Q(0,1) = sig2*(1-2*exp(-b*delta)+exp(-2*b*delta))/2;
  Q(1,0) = Q(0,1);
  Q.submat(2,2,3,3) = Q.submat(0,0,1,1);
  return Q;
}


// [[Rcpp::export]]
Rcpp::List CTCRWNLL(const arma::mat& y, const  arma::mat& Hmat, 
const arma::mat& Qmat, const arma::mat& Tmat, 
const arma::vec& noObs,const arma::vec& active, const arma::colvec& a,
const arma::mat& P)
{
  // Define fixed matrices
  int N = y.n_rows;
  double detF;
  arma::mat Z(2,4, fill::zeros);
  Z(0,0) = 1; Z(1,2) = 1;
  arma::mat T(4,4, fill::zeros);
  T(0,0) = 1; T(2,2) = 1;
  arma::mat Q(4,4, fill::zeros);
  arma::mat F(2,2, fill::zeros);
  arma::mat H(2,2, fill::zeros);
  arma::mat K(4,2, fill::zeros);
  arma::mat L(4,4, fill::zeros);
  arma::colvec v(2, fill::zeros);
  arma::colvec aest(4);
  aest=a;
  arma::mat Pest(4,4);
  Pest=P;
  double ll=0;
  //Begin Kalman filter
  for(int i=0; i<N; i++){
    if(active(i)==0){
      T(0,1) = 0; T(1,1)=0; T(2,3)=0; T(3,3)=0;
      Q = 0.0*Q;
    } else{
      Q(0,0) = Qmat(i,0);
      Q(1,0) = Qmat(i,1);
      Q(0,1) = Qmat(i,1);
      Q(1,1) = Qmat(i,2);     
      Q(2,2) = Qmat(i,0);
      Q(3,2) = Qmat(i,1);
      Q(2,3) = Qmat(i,1);
      Q(3,3) = Qmat(i,2);
      T(0,1) = Tmat(i,0);
      T(1,1) = Tmat(i,1);
      T(2,3) = Tmat(i,0);
      T(3,3) = Tmat(i,1);
    }
    // prediction
    if(noObs(i)==1){
      aest = T*aest;
      Pest = T*Pest*T.t() + Q;
    } else {
      H(0,0) = Hmat(i,0); 
      H(1,1) = Hmat(i,1); 
      H(0,1) = Hmat(i,2); 
      H(1,0) = Hmat(i,2);
      v = y.row(i).t()-Z*aest;
      F = Z*Pest*Z.t() + H; 
      detF = F(0,0)*F(1,1) - F(1,0)*F(0,1);
      if(detF<=0){
        aest = T*aest;
        Pest = T*Pest*T.t() + Q;
      } else{
        ll += - (log(detF) + dot(v,solve(F,v)))/2; 
        K = T*Pest*Z.t()*F.i();     
        L = T - K*Z;
        aest = T*aest + K*v;
        Pest = T*Pest*L.t() + Q;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("ll") = ll);
}
