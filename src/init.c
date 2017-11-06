#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _crawl_CTCRWNLL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_CTCRWNLL_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_CTCRWPREDICT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_CTCRWPREDICT_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_CTCRWSAMPLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_CTCRWSAMPLE_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_makeQ(SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_makeQ_drift(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _crawl_makeT(SEXP, SEXP, SEXP);
extern SEXP _crawl_makeT_drift(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_crawl_CTCRWNLL",           (DL_FUNC) &_crawl_CTCRWNLL,            9},
    {"_crawl_CTCRWNLL_DRIFT",     (DL_FUNC) &_crawl_CTCRWNLL_DRIFT,     11},
    {"_crawl_CTCRWPREDICT",       (DL_FUNC) &_crawl_CTCRWPREDICT,        9},
    {"_crawl_CTCRWPREDICT_DRIFT", (DL_FUNC) &_crawl_CTCRWPREDICT_DRIFT, 11},
    {"_crawl_CTCRWSAMPLE",        (DL_FUNC) &_crawl_CTCRWSAMPLE,         9},
    {"_crawl_CTCRWSAMPLE_DRIFT",  (DL_FUNC) &_crawl_CTCRWSAMPLE_DRIFT,  11},
    {"_crawl_makeQ",              (DL_FUNC) &_crawl_makeQ,               4},
    {"_crawl_makeQ_drift",        (DL_FUNC) &_crawl_makeQ_drift,         6},
    {"_crawl_makeT",              (DL_FUNC) &_crawl_makeT,               3},
    {"_crawl_makeT_drift",        (DL_FUNC) &_crawl_makeT_drift,         4},
    {NULL, NULL, 0}
};

void R_init_crawl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}