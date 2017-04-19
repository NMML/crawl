#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP crawl_CTCRWNLL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_CTCRWNLL_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_CTCRWPREDICT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_CTCRWPREDICT_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_CTCRWSAMPLE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_CTCRWSAMPLE_DRIFT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_makeQ(SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_makeQ_drift(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crawl_makeT(SEXP, SEXP, SEXP);
extern SEXP crawl_makeT_drift(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"crawl_CTCRWNLL",           (DL_FUNC) &crawl_CTCRWNLL,            9},
    {"crawl_CTCRWNLL_DRIFT",     (DL_FUNC) &crawl_CTCRWNLL_DRIFT,     11},
    {"crawl_CTCRWPREDICT",       (DL_FUNC) &crawl_CTCRWPREDICT,        9},
    {"crawl_CTCRWPREDICT_DRIFT", (DL_FUNC) &crawl_CTCRWPREDICT_DRIFT, 11},
    {"crawl_CTCRWSAMPLE",        (DL_FUNC) &crawl_CTCRWSAMPLE,         9},
    {"crawl_CTCRWSAMPLE_DRIFT",  (DL_FUNC) &crawl_CTCRWSAMPLE_DRIFT,  11},
    {"crawl_makeQ",              (DL_FUNC) &crawl_makeQ,               4},
    {"crawl_makeQ_drift",        (DL_FUNC) &crawl_makeQ_drift,         6},
    {"crawl_makeT",              (DL_FUNC) &crawl_makeT,               3},
    {"crawl_makeT_drift",        (DL_FUNC) &crawl_makeT_drift,         4},
    {NULL, NULL, 0}
};

void R_init_crawl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}