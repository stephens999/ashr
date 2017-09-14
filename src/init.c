#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// Declarations for .Call entry points.
SEXP _ashr_cxxMixSquarem (SEXP matrix_likSEXP, SEXP priorSEXP,
			 SEXP pi_initSEXP, SEXP controlSEXP);

// See "Registering native routines" in "Writing R Extensions" manual
// for an explanation of what these lines of code do.
#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

const static R_CallMethodDef R_CallDef[] = {
  CALLDEF(_ashr_cxxMixSquarem,4),
  {NULL, NULL, 0}
};

void attribute_visible R_init_varbvs(DllInfo *dll)
{
  R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
  R_useDynamicSymbols(dll,FALSE);
  R_forceSymbols(dll,TRUE);
}
