#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void solve_coxpath(int *n, double *x, double *lo, double *hi, int *solved, double *z, double *mz);
extern void solve_glmpath(int *n, double *x, double *lo, double *hi, int *solved, double *z, double *mz);

static const R_CMethodDef CEntries[] = {
    {"solve_coxpath", (DL_FUNC) &solve_coxpath, 7},
    {"solve_glmpath", (DL_FUNC) &solve_glmpath, 7},
    {NULL, NULL, 0}
};

void R_init_glmpath(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
