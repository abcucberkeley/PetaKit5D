#ifndef __MEXUTILS_H__
#define __MEXUTILS_H__

#include "mex.h"

/*************************************************************
 wrappers to MATLAB function calls through MEX interface

 implemented in mexCalls.c
************************************************************/

void mexCallCholSol(double*, double*,double*,int);
void mexCallGaussSol(double*, double*,double*,int);
void mexCallCholInv(double*,int);
void mexCallGaussInv(double*,int);
void mexCallSymEig(double*,int,double*,double*,int);
void mexCallDisp(double* a, int nRows, int nCols);


/*********************************************************
  dataconverters between c and MATLAB arrays, vectors, etc.

  implemented in mexDataConverters.c
***********************************************************/

int mat2dblMtx(double**, int* , int* , const mxArray*);
mxArray* dblMtx2mat(double*, int, int);
int mat2uCharMtx(unsigned char**, int* , int* , const mxArray*);
mxArray* uCharMtx2mat(unsigned char*, int, int);
int mat2dblVec(double**, int*, const mxArray*);
int mat2dblScalar(double*, const mxArray*);
mxArray* dblScalar2mat(double);
int mat2intVec(int**, int*, mxArray*);
int mat2intScalar(int*,mxArray*);
mxArray* intScalar2mat(int);
void string2charBuf(char**,int*,mxArray*);

#endif
