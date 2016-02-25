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
Rcpp::List CTCRWPREDICT_DRIFT(const arma::mat& y, const  arma::mat& Hmat, 
                              const arma::vec& beta, const arma::vec& beta_drift, const arma::vec& sig2, 
                              const arma::vec& sig2_drift, const arma::vec& delta,
                              const arma::vec& noObs,const arma::vec& active, const arma::colvec& a,
                              const arma::mat& P)
{
  int I = y.n_rows;
  arma::mat u(2,I, fill::zeros);
  arma::mat jk(2,I, fill::zeros);
  arma::cube M(2,2,I, fill::zeros);
  arma::mat pred(6,I, fill::zeros);
  arma::cube predVar(6,6,I, fill::zeros);
  arma::mat Z(2,6, fill::zeros); Z(0,0) = 1; Z(1,3) = 1;
  arma::mat T(6,6, fill::zeros);
  arma::mat Q(6,6, fill::zeros);
  arma::cube F(2,2,I, fill::zeros);
  arma::mat H(2,2, fill::zeros);
  arma::cube K(6,2,I, fill::zeros);
  arma::cube L(6,6,I, fill::zeros);
  arma::mat v(2,I, fill::zeros);
  arma::mat aest(6,I+1, fill::zeros);
  aest.col(0)=a;
  arma::cube Pest(6,6,I+1, fill::zeros);
  Pest.slice(0)=P;
  arma::colvec r(6, fill::zeros);
  arma::mat N(6,6, fill::zeros);
  arma::vec chisq(I, fill::zeros);
  
  double ll=0;
  //Forward filter
  for(int i=0; i<I; i++){
    Q = makeQ_drift(beta(i), beta_drift(i), sig2(i), sig2_drift(i), delta(i), active(i));
    T = makeT_drift(beta(i), beta_drift(i), delta(i), active(i));
    if(noObs(i)==1){
      aest.col(i+1) = T*aest.col(i);
      Pest.slice(i+1) = T*Pest.slice(i)*T.t() + Q;
      L.slice(i) = T;
    } else {
      H(0,0) = Hmat(i,0); 
      H(1,1) = Hmat(i,1); 
      H(0,1) = Hmat(i,2); 
      H(1,0) = Hmat(i,2);
      v.col(i) = y.row(i).t()-Z*aest.col(i);
      F.slice(i) = Z*Pest.slice(i)*Z.t() + H; 
      ll += - (log(det(F.slice(i))) + dot(v.col(i),solve(F.slice(i),v.col(i))))/2; 
      K.slice(i) = T*Pest.slice(i)*Z.t()*F.slice(i).i();     
      L.slice(i) = T - K.slice(i)*Z;
      aest.col(i+1) = T*aest.col(i) + K.slice(i)*v.col(i);
      Pest.slice(i+1) = T*Pest.slice(i)*L.slice(i).t() + Q;
    }
  }
  
  // Backward smooth
  for(int j=I; j>0; j--){
    if(noObs(j-1)==1 || F.slice(j-1)(0,0)*F.slice(j-1)(1,1)==0){
      r = L.slice(j-1).t() * r;
      N = L.slice(j-1).t() * N * L.slice(j-1);
    } else{
      u.col(j-1) = solve(F.slice(j-1),v.col(j-1))-K.slice(j-1).t()*r;
      M.slice(j-1) = F.slice(j-1).i() + K.slice(j-1).t()*N*K.slice(j-1);
      chisq(j-1) = dot(u.col(j-1),solve(M.slice(j-1),u.col(j-1)));
      jk.col(j-1) = y.row(j-1).t() - solve(M.slice(j-1),u.col(j-1));
      r = Z.t()*solve(F.slice(j-1),v.col(j-1)) + L.slice(j-1).t() * r;
      N = Z.t() * solve(F.slice(j-1),Z) + L.slice(j-1).t()*N*L.slice(j-1); 
    }
    pred.col(j-1) = aest.col(j-1) + Pest.slice(j-1)*r;
    predVar.slice(j-1) = Pest.slice(j-1) - Pest.slice(j-1)*N*Pest.slice(j-1); 
  }
  return Rcpp::List::create(
    Rcpp::Named("ll") = ll, 
    Rcpp::Named("pred") = pred, Rcpp::Named("predVar")=predVar,
    Rcpp::Named("chisq")=chisq, Rcpp::Named("predObs")=jk
    );
}
