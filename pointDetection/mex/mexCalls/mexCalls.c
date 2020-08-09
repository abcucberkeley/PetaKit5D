#include <memory.h>
#include "mex.h"
#include "mexUtils.h"

void mexCallCholSol(double* a, double* b,double* x,int n)
/*
  solves the equation system a*x = b with Cholesky decomposition, where 
  a is an nxn positive definite matrix (only upper triangle stored)
  b is an nx1 vector

  can only be called from a MATLAB mex function;

  otherwise use the MATLAB engine
  
  the vector x must be allocated from outside

*/

{
  double*   dPtr;
  double   *inA,*inFullA;

  mxArray*  fullA;
  mxArray*  matB;
  mxArray*  matX;
  mxArray  *aux1,*aux2,*aux3;
  mxArray  *inputArray[2];

  int i,j;
  
  /* pack the matrix a */
  fullA = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(fullA);
  inFullA = dPtr;
  inA = a;
  for(i=0; i<n; i++){
    memcpy(inFullA+i,inA,(n-i)*sizeof(double));
    for(j=0;j<i;j++){
      *(inFullA+j) = *(dPtr + j*n + i);
    }
    inFullA += n;
    inA += n-i;
  }

  /*   mexCallMATLAB(0,&aux1,1,&fullA,"disp"); */
  
  /* pack vector b */
  matB = mxCreateDoubleMatrix(n,1,mxREAL);
  dPtr = mxGetPr(matB);
  memcpy(dPtr,b,n*sizeof(double));

  /*  mexCallMATLAB(0,&aux1,1,&matB,"disp"); */

  /* do computation */
  mexCallMATLAB(1,&aux1,1,&fullA,"chol");
  mexCallMATLAB(1,&aux2,1,&aux1,"transpose");

  inputArray[0] = aux2;
  inputArray[1] = matB;
  mexCallMATLAB(1,&aux3,2,inputArray,"mldivide");
  
  inputArray[0] = aux1;
  inputArray[1] = aux3;
  mexCallMATLAB(1,&matX,2,inputArray,"mldivide");

  /* unpack vector X */
  dPtr =  mxGetPr(matX);
  memcpy(x,dPtr,n*sizeof(double));
  
  /*  mexCallMATLAB(0,&aux1,1,&matX,"disp"); */

  mxDestroyArray(aux1);
  mxDestroyArray(aux2);
  mxDestroyArray(aux3);
  mxDestroyArray(fullA);
  mxDestroyArray(matB);
  mxDestroyArray(matX);

  return;
}


void mexCallGaussSol(double* a, double* b,double* x,int n)
/*
  solves the equation system a*x = b with Gauss elimination, where 
  a is an nxn regular, symmetric matrix (only upper triangle stored)
  b is an nx1 vector

  can only be called from a MATLAB mex function;

  otherwise use the MATLAB engine
  
  the vector x must be allocated from outside

*/

{
  double*   dPtr;
  double   *inA,*inFullA;

  mxArray*  fullA;
  mxArray*  matB;
  mxArray*  matX;
  mxArray  *inputArray[2];

  int i,j;
  
  /* pack the matrix a */
  fullA = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(fullA);
  inFullA = dPtr;
  inA = a;
  for(i=0; i<n; i++){
    memcpy(inFullA+i,inA,(n-i)*sizeof(double));
    for(j=0;j<i;j++){
      *(inFullA+j) = *(dPtr + j*n + i);
    }
    inFullA += n;
    inA += n-i;
  }

  /*  mexCallMATLAB(0,NULL,1,&fullA,"disp");*/
  
  /* pack vector b */
  matB = mxCreateDoubleMatrix(n,1,mxREAL);
  dPtr = mxGetPr(matB);
  memcpy(dPtr,b,n*sizeof(double));

  /*  mexCallMATLAB(0,NULL,1,&matB,"disp");*/

  /* do computation */
  inputArray[0] = fullA;
  inputArray[1] = matB;
  mexCallMATLAB(1,&matX,2,inputArray,"mldivide");

  /* unpack vector X */
  dPtr =  mxGetPr(matX);
  memcpy(x,dPtr,n*sizeof(double));
  
  /*  mexCallMATLAB(0,&aux1,1,&matX,"disp"); */
  mxDestroyArray(fullA);
  mxDestroyArray(matB);
  mxDestroyArray(matX);

  return;
}


void mexCallCholInv(double* a,int n)
/*
  overwrites the positive definite nxn matrix a with a^(-1), 
  where a is the upper triangle of the full matrix a;

  can only be called when from a MATLAB mex function;

  otherwise use the MATLAB engine
  
*/

{
  double*   dPtr;
  double   *inA,*inFullA,*inE,*inX;

  mxArray*  fullA;
  mxArray*  matE;
  mxArray*  matX;
  mxArray  *aux1,*aux2,*aux3;
  mxArray  *inputArray[2];

  int i,j;
  
  /* pack the matrix a */
  fullA = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(fullA);
  inFullA = dPtr;
  inA = a;
  for(i=0; i<n; i++){
    memcpy(inFullA+i,inA,(n-i)*sizeof(double));
    for(j=0;j<i;j++){
      *(inFullA+j) = *(dPtr + j*n + i);
    }
    inFullA += n;
    inA += n-i;
  }

  /*  mexCallMATLAB(0,&aux1,1,&fullA,"disp"); */
  
  /* get identity */
  matE = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(matE);
  inE  = dPtr;
  memset(dPtr,0,mxGetNumberOfElements(matE)*mxGetElementSize(matE));
  for(i=0; i<n; i++){
    *(inE+i) = (double)1.0;
    inE += n;
  }
  
  /*  mexCallMATLAB(0,&aux1,1,&matB,"disp"); */

  /* do computation */
  mexCallMATLAB(1,&aux1,1,&fullA,"chol");
  mexCallMATLAB(1,&aux2,1,&aux1,"transpose");

  inputArray[0] = aux2;
  inputArray[1] = matE;
  mexCallMATLAB(1,&aux3,2,inputArray,"mldivide");
  
  inputArray[0] = aux1;
  inputArray[1] = aux3;
  mexCallMATLAB(1,&matX,2,inputArray,"mldivide");

  /* unpack matrix X */
  dPtr =  mxGetPr(matX);
  inX = dPtr;
  inA = a;
  for(i = 0; i<n; i++){
    memcpy(inA,inX+i,(n-i)*sizeof(double));
    inX += n;
    inA += n-i;
  }
  
  /* mexCallMATLAB(0,&aux1,1,&matX,"disp");*/

  mxDestroyArray(aux1);
  mxDestroyArray(aux2);
  mxDestroyArray(aux3);
  mxDestroyArray(fullA);
  mxDestroyArray(matE);
  mxDestroyArray(matX);

  return;
}

void mexCallGaussInv(double* a,int n)
/*
  overwrites the positive definite nxn matrix a with a^(-1), 
  based on Gauss elimination, where a is the upper triangle of 
  the full matrix a; 

  can only be called when from a MATLAB mex function;

  otherwise use the MATLAB engine
  
*/

{
  double*   dPtr;
  double   *inA,*inFullA,*inE,*inX;

  mxArray*  fullA;
  mxArray*  matE;
  mxArray*  matX;
  mxArray  *inputArray[2];

  int i,j;
  
  /* pack the matrix a */
  fullA = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(fullA);
  inFullA = dPtr;
  inA = a;
  for(i=0; i<n; i++){
    memcpy(inFullA+i,inA,(n-i)*sizeof(double));
    for(j=0;j<i;j++){
      *(inFullA+j) = *(dPtr + j*n + i);
    }
    inFullA += n;
    inA += n-i;
  }

  /*  mexCallMATLAB(0,&aux1,1,&fullA,"disp"); */
  
  /* get identity */
  matE = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(matE);
  inE  = dPtr;
  memset(dPtr,0,mxGetNumberOfElements(matE)*mxGetElementSize(matE));
  for(i=0; i<n; i++){
    *(inE+i) = (double)1.0;
    inE += n;
  }
  
  /*  mexCallMATLAB(0,&aux1,1,&matB,"disp"); */

  /* do computation */
  inputArray[0] = fullA;
  inputArray[1] = matE;
  mexCallMATLAB(1,&matX,2,inputArray,"mldivide");

  /* unpack matrix X */
  dPtr =  mxGetPr(matX);
  inX = dPtr;
  inA = a;
  for(i = 0; i<n; i++){
    memcpy(inA,inX+i,(n-i)*sizeof(double));
    inX += n;
    inA += n-i;
  }
  
  /* mexCallMATLAB(0,&aux1,1,&matX,"disp");*/
  mxDestroyArray(fullA);
  mxDestroyArray(matE);
  mxDestroyArray(matX);

  return;
}

void mexCallSymEig(double* a,int n,double* diag,double* evecs,int tflag)
/*
  computes the eigenvectors of an nxn symmetric matrix A (only upper triangle
  stored) such that 
  
  A = evecs * diag * evecs'

  The return values are 
     - nxn rotation matrix evecs  if tflag = 0
                           evecs' if tflag = 1
     - n-dim vector diag 

  For both the vectors the memory has to be allocated outside the call
  
*/

{
  double*   dPtr;
  double   *inA,*inFullA;

  mxArray*  fullA;
  mxArray*  evT;            /* transposed (in MATLAB sense) 
			      eigenvector matrix, which is in C sense the 
			      proper evecs matrix */
  mxArray  *outputArray[2];
  int i,j;
  
  /* pack the matrix a */
  fullA = mxCreateDoubleMatrix(n,n,mxREAL);
  dPtr = mxGetPr(fullA);
  inFullA = dPtr;
  inA = a;
  for(i=0; i<n; i++){
    memcpy(inFullA+i,inA,(n-i)*sizeof(double));
    for(j=0;j<i;j++){
      *(inFullA+j) = *(dPtr + j*n + i);
    }
    inFullA += n;
    inA += n-i;
  }

  /*  mexCallMATLAB(0,&aux1,1,&fullA,"disp"); */
  
  /* do computation */
  mexCallMATLAB(2,outputArray,1,&fullA,"eig");

  /* unpack matrix evecs */
  if(tflag){
    dPtr = mxGetPr(outputArray[0]);
    memcpy(evecs,dPtr,n*n*sizeof(double));
  }
  else{
     mexCallMATLAB(1,&evT,1,&(outputArray[0]),"transpose");
     dPtr = mxGetPr(evT);
     memcpy(evecs,dPtr,n*n*sizeof(double));
     mxDestroyArray(evT);
  }

  /* unpack matrix diag */
  dPtr = mxGetPr(outputArray[1]);
  for(i = 0; i<n; i++){
    diag[i] = *(dPtr+i);
    dPtr+=n;
  }
    
  /* mexCallMATLAB(0,&aux1,1,&matX,"disp");*/
  mxDestroyArray(fullA);
  mxDestroyArray(outputArray[0]);
  mxDestroyArray(outputArray[1]);

  return;
}

/**************************************************************
  S E R V I C E  part
*************************************************************/

void mexCallDisp(double* a, int nRows, int nCols)
/*
 displays the  matrix a on the MATLAB screen
 can only be called from a MATLAB mex function;
 
*/

{
  mxArray *matA,*matAt,*dummy;
  
 
  /* pack matA */
  matA = dblMtx2mat(a,nCols,nRows);
  /* CAUTION: the definition of dblMtx is that the matrix comes 
     in row-wise storage. dblMtx will, however, fill the rows into columns,
     thus, before displying matA has to be transposed. */

  mexCallMATLAB(1,&matAt,1,&matA,"transpose");
  mexCallMATLAB(0,&dummy,1,&matAt,"disp");

  mxDestroyArray(matA);
  mxDestroyArray(matAt);

  return;
}
