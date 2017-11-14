#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP systemicrisk_cloneMatrix(SEXP);
extern SEXP systemicrisk_ERE_step_cycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP systemicrisk_findFeasibleMatrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP systemicrisk_GibbsSteps_kcycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"systemicrisk_cloneMatrix",        (DL_FUNC) &systemicrisk_cloneMatrix,        1},
    {"systemicrisk_ERE_step_cycle",     (DL_FUNC) &systemicrisk_ERE_step_cycle,     6},
    {"systemicrisk_findFeasibleMatrix", (DL_FUNC) &systemicrisk_findFeasibleMatrix, 4},
    {"systemicrisk_GibbsSteps_kcycle",  (DL_FUNC) &systemicrisk_GibbsSteps_kcycle,  6},
    {NULL, NULL, 0}
};

void R_init_systemicrisk(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
