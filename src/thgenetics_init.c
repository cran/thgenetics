// tools::package_native_routine_registration_skeleton("thgenetics") # give the path to the package source
// creates all of the code listed below for this

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void zstat_pathway_perm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zstat_perm(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"zstat_pathway_perm", (DL_FUNC) &zstat_pathway_perm, 13},
    {"zstat_perm",         (DL_FUNC) &zstat_perm,         13},
    {NULL, NULL, 0}
};

void R_init_thgenetics(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
