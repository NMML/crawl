// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function prototypes
arma::mat makeT_drift(const double& b, const double& b_drift, const double& delta, const double& active);
arma::mat makeQ_drift(const double& b, const double& b_drift, const double& sig2, const double& sig2_drift, 
                      const double& delta, const double& active);


// [[Rcpp::export]]
Rcpp::List CTCRWNLL_DRIFT(const arma::mat& y, const  arma::mat& Hmat, 
                    const arma::vec& beta, const arma::vec& beta_drift, const arma::vec& sig2, 
                    const arma::vec& sig2_drift, const arma::vec& delta,
                    const arma::vec& noObs,const arma::vec& active, const arma::colvec& a,
                    const arma::mat& P)
{
  // Define fixed matrices
  int N = y.n_rows;
  double detF;
  arma::mat Z(2,6, fill::zeros);
  Z(0,0) = 1; Z(1,3) = 1;
  arma::mat T(6,6);
  arma::mat Q(6,6);
  arma::mat F(2,2, fill::zeros);
  arma::mat H(2,2, fill::zeros);
  arma::mat K(6,2, fill::zeros);
  arma::mat L(6,6, fill::zeros);
  arma::colvec v(2, fill::zeros);
  arma::colvec aest(6);
  aest=a;
  arma::mat Pest(6,6);
  Pest=P;
  double ll=0;
  //Begin Kalman filter
  for(int i=0; i<N; i++){
    T = makeT_drift(beta(i), beta_drift(i), delta(i), active(i));
    Q = makeQ_drift(beta(i), beta_drift(i), sig2(i), sig2_drift(i), delta(i), active(i));
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
