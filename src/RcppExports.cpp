// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// makeT
arma::mat makeT(const double& b, const double& delta);
RcppExport SEXP crawl_makeT(SEXP bSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const double& >::type b(bSEXP );
        Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP );
        arma::mat __result = makeT(b, delta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// makeQ
arma::mat makeQ(const double& b, const double& sig2, const double& delta);
RcppExport SEXP crawl_makeQ(SEXP bSEXP, SEXP sig2SEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const double& >::type b(bSEXP );
        Rcpp::traits::input_parameter< const double& >::type sig2(sig2SEXP );
        Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP );
        arma::mat __result = makeQ(b, sig2, delta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CTCRWNLL
Rcpp::List CTCRWNLL(const arma::mat& y, const arma::mat& Hmat, const arma::mat& Qmat, const arma::mat& Tmat, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP crawl_CTCRWNLL(SEXP ySEXP, SEXP HmatSEXP, SEXP QmatSEXP, SEXP TmatSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Qmat(QmatSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Tmat(TmatSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP );
        Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP );
        Rcpp::List __result = CTCRWNLL(y, Hmat, Qmat, Tmat, noObs, active, a, P);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// CTCRWPREDICT
Rcpp::List CTCRWPREDICT(const arma::mat& y, const arma::mat& Hmat, const arma::mat& Qmat, const arma::mat& Tmat, const arma::vec& noObs, const arma::vec& active, const arma::colvec& a, const arma::mat& P);
RcppExport SEXP crawl_CTCRWPREDICT(SEXP ySEXP, SEXP HmatSEXP, SEXP QmatSEXP, SEXP TmatSEXP, SEXP noObsSEXP, SEXP activeSEXP, SEXP aSEXP, SEXP PSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Hmat(HmatSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Qmat(QmatSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Tmat(TmatSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type noObs(noObsSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type active(activeSEXP );
        Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP );
        Rcpp::List __result = CTCRWPREDICT(y, Hmat, Qmat, Tmat, noObs, active, a, P);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
