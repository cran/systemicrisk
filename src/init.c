#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _systemicrisk_cloneMatrix(SEXP);
extern SEXP _systemicrisk_ERE_step_cycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _systemicrisk_findFeasibleMatrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP _systemicrisk_GibbsSteps_kcycle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_systemicrisk_cloneMatrix",        (DL_FUNC) &_systemicrisk_cloneMatrix,        1},
    {"_systemicrisk_ERE_step_cycle",     (DL_FUNC) &_systemicrisk_ERE_step_cycle,     6},
    {"_systemicrisk_findFeasibleMatrix", (DL_FUNC) &_systemicrisk_findFeasibleMatrix, 4},
    {"_systemicrisk_GibbsSteps_kcycle",  (DL_FUNC) &_systemicrisk_GibbsSteps_kcycle,  6},
    {NULL, NULL, 0}
};

void R_init_systemicrisk(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
