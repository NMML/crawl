// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Random normal draws for posterior sampling 
arma::vec armaNorm(int n){
  NumericVector x = Rcpp::rnorm(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}

arma::vec mvn(const arma::vec& mu, const arma::mat& Sig){
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd(U, s, V, Sig);
  arma::mat out = mu + U*diagmat(sqrt(s))*armaNorm(mu.n_elem);
  return out; 
}

// [[Rcpp::export]]
arma::mat makeT(const double& b, const double& delta, const double& active){
  arma::mat T(4,4, fill::zeros);
  T(0,0) = 1; 
  T(2,2) = 1;
  if(active > 0){
    T(0,1) = exp(R::pexp(delta,1/b,1,1) - log(b));
    T(1,1) = exp(-b*delta);
    T(2,3) = exp(R::pexp(delta,1/b,1,1) - log(b));
    T(3,3) = exp(-b*delta);
  }
  return T;
}

// [[Rcpp::export]]
arma::mat makeQ(const double& b, const double& sig2, const double& delta, const double& active){
  arma::mat Q(4,4, fill::zeros);
  if(active > 0){
    Q(0,0) =  sig2*(delta - 2*exp(R::pexp(delta,1/b,1,1)-log(b))  + exp(R::pexp(delta,1/(2*b),1,1)-log(2*b)));
    Q(1,1) = sig2*exp(log(b) + R::pexp(delta,1/(2*b),1,1))/2;
    Q(0,1) = sig2*(1-2*exp(-b*delta)+exp(-2*b*delta))/2;
    Q(1,0) = Q(0,1);
    Q.submat(2,2,3,3) = Q.submat(0,0,1,1);
  }
  return Q;
}

// [[Rcpp::export]]
arma::mat makeT_drift(const double& b, const double& b_drift, const double& delta, const double& active){
  arma::mat T(6,6, fill::zeros);
  T(0,0) = 1; 
  if(active > 0){
    T(0,1) = exp(R::pexp(delta,1/b,1,1) - log(b));
    T(0,2) = exp(R::pexp(delta,1/b_drift,1,1) - log(b_drift));
    T(1,1) = exp(-b*delta);
    T(2,2) = exp(-b_drift*delta);
  }
  T.submat(3,3,5,5) = T.submat(0,0,2,2);
  return T;
}

// [[Rcpp::export]]
arma::mat makeQ_drift(const double& b, const double& b_drift, const double& sig2, const double& sig2_drift, 
                      const double& delta, const double& active){
  arma::mat Q(6,6, fill::zeros);
  if(active > 0){
    Q(0,0) =  sig2*(delta - 2*exp(R::pexp(delta,1/b,1,1)-log(b))  + exp(R::pexp(delta,1/(2*b),1,1)-log(2*b))) + 
      sig2_drift*(delta - 2*exp(R::pexp(delta,1/b_drift,1,1)-log(b_drift))  + exp(R::pexp(delta,1/(2*b_drift),1,1)-log(2*b_drift)));
    Q(1,1) = sig2*exp(log(b) + R::pexp(delta,1/(2*b),1,1))/2;
    Q(2,2) = sig2_drift*exp(log(b_drift) + R::pexp(delta,1/(2*b_drift),1,1))/2;
    Q(0,1) = sig2*(1-2*exp(-b*delta)+exp(-2*b*delta))/2;
    Q(1,0) = Q(0,1);
    Q(0,2) = sig2_drift*(1-2*exp(-b_drift*delta)+exp(-2*b_drift*delta))/2;
    Q(2,0) = Q(0,2);
    Q.submat(3,3,5,5) = Q.submat(0,0,2,2);
  }
  return Q;
}
