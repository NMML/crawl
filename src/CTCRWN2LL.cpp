// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function prototypes
arma::mat makeQ(const double& b, const double& sig2, const double& delta, const double& active);
arma::mat makeT(const double& b, const double& delta, const double& active);

// [[Rcpp::export]]
Rcpp::List CTCRWNLL(const arma::mat& y, const  arma::mat& Hmat, 
                    const arma::vec& beta, const arma::vec& sig2, const arma::vec& delta,
                    const arma::vec& noObs,const arma::vec& active, const arma::colvec& a,
                    const arma::mat& P)
{
  // Define fixed matrices
  int N = y.n_rows;
  double detF;
  arma::mat Z(2,4, fill::zeros);
  Z(0,0) = 1; Z(1,2) = 1;
  arma::mat T(4,4);
  arma::mat Q(4,4);
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
    T = makeT(beta(i), delta(i), active(i));
    Q = makeQ(beta(i), sig2(i), delta(i), active(i));
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
        // K = T*Pest*Z.t()*inv_sympd(F);  
        K = T*Pest*Z.t()*F.i();
        L = T - K*Z;
        aest = T*aest + K*v;
        Pest = T*Pest*L.t() + Q;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("ll") = ll);
}
