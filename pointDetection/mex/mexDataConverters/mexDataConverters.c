/* service function for data conversion between MATLAB matrix/vector/scalars 
   and C/C++ data arrays */

#include <string.h>
#include <math.h>
#include "mexUtils.h"
#include "macros.h"

int mat2dblMtx(double** buf, int* nx, int* ny, const mxArray* mat)
{
  int  nDims;
  const int* dims;
  double*  dPtr;

  *buf = NULL;
  *nx = *ny = 0;
  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  nDims = mxGetNumberOfDimensions(mat);
  if(nDims != 2)
    return(0);
  dims = mxGetDimensions(mat);
  if( (dims[0] <= 0) || (dims[1] <=0))
    return(0);
  if(!mxIsDouble(mat))
    return(0);

  *nx = dims[0];
  *ny = dims[1];
  *buf = (double*)mxMalloc(dims[0]*dims[1]*sizeof(double));
  dPtr = mxGetPr(mat);
  memcpy((void*)*buf,(void*)dPtr,dims[0]*dims[1]*sizeof(double));

#ifdef VERBOSE
  mexPrintf("%d x %d double matrix converted\n",*nx,*ny);
#endif
  return(1);
}

mxArray* dblMtx2mat(double* buf, int nx, int ny)
{
  mxArray* mat;
  double *dPtr;

  if((!buf)||(!nx)||(!ny)){
    /* create empty matrix, which is a 0x0 matrix */
    mat=mxCreateDoubleMatrix(0,0,mxREAL);
  }
  else{
    mat=mxCreateDoubleMatrix(nx,ny,mxREAL);
    dPtr = mxGetPr(mat);
    memcpy((void*)dPtr,(void*)buf,nx*ny*sizeof(double));
  }

  return(mat);
}

int mat2uCharMtx(unsigned char** buf, int* nx, int* ny, const mxArray* mat)
{
  int  nDims;
  const int* dims;
  unsigned char*  dPtr;

  *buf = NULL;
  *nx=*ny=0;
  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  nDims = mxGetNumberOfDimensions(mat);
  if(nDims != 2)
    return(0);
  dims = mxGetDimensions(mat);
  if( (dims[0] <= 0) || (dims[1] <=0))
    return(0);
  if(!mxIsUint8(mat))
    return(0);

  *nx = dims[0];
  *ny = dims[1];
  *buf = (unsigned char*)mxMalloc(dims[0]*dims[1]*sizeof(unsigned char));
  dPtr = (unsigned char*)mxGetPr(mat);
  memcpy((void*)*buf,(void*)dPtr,dims[0]*dims[1]*sizeof(unsigned char));

#ifdef VERBOSE
  mexPrintf("%d x %d uint8 matrix converted\n",*nx,*ny);
#endif
  return(1);
}

mxArray* uCharMtx2mat(unsigned char* buf, int nx, int ny)
{
  mxArray* mat;
  unsigned char *dPtr;
  int  dims[2];

  if((!buf)||(!nx)||(!ny)){
    /* create empty matrix, which is a 0x0 matrix */
    dims[0] = 0;
    dims[1] = 0;
    mat=mxCreateNumericArray(2,dims,mxUINT8_CLASS,mxREAL);
  }
  else{
    dims[0] = nx;
    dims[1] = ny;
    mat=mxCreateNumericArray(2,dims,mxUINT8_CLASS,mxREAL);
    dPtr = (unsigned char*)mxGetPr(mat);
    memcpy((void*)dPtr,(void*)buf,nx*ny*sizeof(unsigned char));
  }
  return(mat);
}

int mat2intScalar(int* val,mxArray* mat)
{
  double* pData;

  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  if(mxGetNumberOfDimensions(mat) != 2)
    return(0);
  if((mxGetN(mat)*mxGetM(mat)) != 1)
    return(0);
  
  pData = mxGetPr(mat);
  *val = DROUND(*pData);
  return(1);
}

mxArray* intScalar2mat(int val)
{
  mxArray* mat;
  double*  pData;

  mat = mxCreateDoubleMatrix(1,1,mxREAL);
  pData = mxGetPr(mat);
  pData[0] = (double)val;

  return(mat);
}

int mat2dblScalar(double* val,const mxArray* mat)
{
  double* pData;

  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  if(mxGetNumberOfDimensions(mat) != 2)
    return(0);
  if((mxGetN(mat)*mxGetM(mat)) != 1)
    return(0);
  
  pData = mxGetPr(mat);
  *val = *pData;
  return(1);
}

mxArray* dblScalar2mat(double val)
{
  mxArray* mat;
  double*  pData;

  mat = mxCreateDoubleMatrix(1,1,mxREAL);
  pData = mxGetPr(mat);
  pData[0] = val;

  return(mat);
}

int mat2dblVec(double** vals, int* nVals, const mxArray* mat)
{
  int n,m;
  double* dPtr;

  *vals = NULL;
  *nVals = 0;

  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  if(mxGetNumberOfDimensions(mat) != 2)
    return(0);
  n = mxGetN(mat);
  m = mxGetM(mat);
  if((n-1)*(m-1))
    return(0);
  *nVals = n*m;
  *vals = (double*)mxMalloc(*nVals*sizeof(double));
  dPtr = mxGetPr(mat);
  memcpy(*vals,dPtr,*nVals*sizeof(double));
  return(1);
}

int mat2intVec(int** vals, int* nVals, mxArray* mat)
     /* the vals ptr IS ALLOCATED in the function */
{
  double* pData;
  int i,n,m;

  *vals = NULL;
  *nVals = 0;
  if(!mat)
    return(0);
  if(mxIsEmpty(mat))
    return(0);
  if(mxGetNumberOfDimensions(mat) != 2)
    return(0);
  n = mxGetN(mat);
  m = mxGetM(mat);
  if((n-1)*(m-1))
    return(0);
  *nVals = n*m;
  *vals = (int*)mxMalloc(*nVals * sizeof(int));
  pData = mxGetPr(mat);
  for(i = 0; i<*nVals; i++){
    (*vals)[i] = DROUND(pData[i]);
  }
  return(1);
}

void string2charBuf(char** buf,int* pBufSze,mxArray* pString)
{
  if(mxIsEmpty(pString)){
    *pBufSze = 0;
    *buf = NULL;
#ifdef VERBOSE
    mexPrintf("Empty string entered\n");
#endif
    return;
  }
  
  *pBufSze = mxGetM(pString)*mxGetN(pString) + 1;
  *buf = (char*)mxCalloc(*pBufSze,sizeof(char));
  mxGetString(pString,*buf,*pBufSze);
#ifdef VERBOSE
  mexPrintf("String entered: %s\n",*buf);
#endif
  
}


