// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function prototypes
arma::mat makeQ(const double& beta, const double& sig2, const double& delta, const double& active);
arma::mat makeT(const double& beta, const double& delta, const double& active);
arma::vec mvn(const arma::vec& mu, const arma::mat& V);



// [[Rcpp::export]]
Rcpp::List CTCRWSAMPLE(const arma::mat& y, const  arma::mat& Hmat, 
const arma::vec& beta, const arma::vec& sig2, const arma::vec& delta, 
const arma::vec& noObs,const arma::vec& active, const arma::colvec& a,
const arma::mat& P)
{
  int I = y.n_rows;
  // SIMULATION RELATED MATRICES
  arma::mat y_plus(I, 2);
  arma::mat alpha_plus(4,I+1, fill::zeros);
  alpha_plus.col(0) = mvn(a,P);
  arma::mat alpha_plus_hat(4,I+1, fill::zeros);
  alpha_plus_hat.col(0) = a;
  arma::mat v_plus(2,I, fill::zeros);
  arma::colvec r_plus(4, fill::zeros);
  arma::mat sim(I,4);
  // KFS MATRICES FOR DATA INFORMED PART
  arma::mat Z(2,4, fill::zeros); Z(0,0) = 1; Z(1,2) = 1;
  arma::mat T(4,4, fill::zeros); T(0,0) = 1; T(2,2) = 1;
  arma::mat Q(4,4, fill::zeros);
  arma::cube F(2,2,I, fill::zeros);
  arma::mat H(2,2, fill::zeros);
  arma::cube K(4,2,I, fill::zeros);
  arma::cube L(4,4,I, fill::zeros);
  arma::mat v(2,I, fill::zeros);
  arma::mat alpha_hat(4,I+1, fill::zeros);
  alpha_hat.col(0) = a;
  arma::cube P_hat(4,4,I+1, fill::zeros);
  P_hat.slice(0)=P;
  arma::colvec r(4, fill::zeros);
  
  double ll=0;
  //Forward filter and simulation
  for(int i=0; i<I; i++){
    Q = makeQ(beta(i), sig2(i), delta(i), active(i));
    T = makeT(beta(i), delta(i), active(i));
    if(active(i)==0){
      alpha_plus.col(i+1) = alpha_plus.col(i);
    } else{
      alpha_plus.col(i+1) = mvn(T*alpha_plus.col(i), Q);
    }
    if(noObs(i)==1){
      alpha_hat.col(i+1) = T*alpha_hat.col(i);
      P_hat.slice(i+1) = T*P_hat.slice(i)*T.t() + Q;
      L.slice(i) = T;
      alpha_plus_hat.col(i+1) = T*alpha_plus_hat.col(i);
    } else {
      H(0,0) = Hmat(i,0); 
      H(1,1) = Hmat(i,1); 
      H(0,1) = Hmat(i,2); 
      H(1,0) = Hmat(i,2);
      // KF simulation for data section
      v.col(i) = y.row(i).t()-Z*alpha_hat.col(i);
      F.slice(i) = Z*P_hat.slice(i)*Z.t() + H; 
      ll += - (log(det(F.slice(i))) + dot(v.col(i),solve(F.slice(i),v.col(i))))/2; 
      K.slice(i) = T*P_hat.slice(i)*Z.t()*F.slice(i).i();     
      L.slice(i) = T - K.slice(i)*Z;
      alpha_hat.col(i+1) = T*alpha_hat.col(i) + K.slice(i)*v.col(i);
      P_hat.slice(i+1) = T*P_hat.slice(i)*L.slice(i).t() + Q;
      // Simulation component
      y_plus.row(i) = mvn(Z*alpha_plus.col(i), H).t();
      v_plus.col(i) = y_plus.row(i).t()-Z*alpha_plus_hat.col(i);
      alpha_plus_hat.col(i+1) = T*alpha_plus_hat.col(i) + K.slice(i)*v_plus.col(i);
      
    }
  }
  
  // Backward smooth
  for(int j=I; j>0; j--){
    if(noObs(j-1)==1 || F.slice(j-1)(0,0)*F.slice(j-1)(1,1)==0){
      r = L.slice(j-1).t() * r;
      r_plus = L.slice(j-1).t() * r_plus;
    } else{
      r = Z.t()*solve(F.slice(j-1),v.col(j-1)) + L.slice(j-1).t() * r;
      r_plus = Z.t()*solve(F.slice(j-1),v_plus.col(j-1)) + L.slice(j-1).t() * r_plus; 
    }
    alpha_hat.col(j-1) += P_hat.slice(j-1)*r;
    alpha_plus_hat.col(j-1) += P_hat.slice(j-1)*r_plus;
  }
  sim = (alpha_hat.cols(0,I-1) - (alpha_plus.cols(0,I-1)-alpha_plus_hat.cols(0,I-1))).t();
  return Rcpp::List::create(
    Rcpp::Named("ll") = ll, 
    Rcpp::Named("sim") = sim
    );
}
