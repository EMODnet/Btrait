#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif

void F77_NAME(extendtraitf)(int*, int*, int*, int*, int*,  
       double*, int*, int*, int*, double*, int*, int*, int*);

void F77_NAME(long2wide)(int*, int*, int*, int*, int*,  
      int*, double*, double*);
     
void F77_NAME(long2long)(int*, int*, int*, int*, int*,  
      int*, double*, double*, int*);
          
R_FortranMethodDef fortranMethods[] = {
 {"long2wide", (DL_FUNC) &F77_SUB(long2wide), 8},
 {"long2long", (DL_FUNC) &F77_SUB(long2long), 9},
 {"extendtraitf", (DL_FUNC) &F77_SUB(extendtraitf), 13},
 {NULL, NULL, 0}
};
void R_init_btrait(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
