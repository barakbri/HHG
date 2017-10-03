#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP HHG_R_C(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
                    SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"HHG_R_C", (DL_FUNC) &HHG_R_C, 16},
  
  {NULL, NULL, 0}
};

void R_init_repfdr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}