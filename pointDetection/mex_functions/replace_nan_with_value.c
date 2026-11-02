/*
 * replace_nan_with_value.c - replace nan in a matrix
 *
 * Implemented by Xiongtao Ruan.
 * 
 * This is a MEX-file for MATLAB.
 *
 * mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' replace_nan_with_value.c
 */ 

#include <math.h>
#include <stdint.h>
#include <omp.h>
#include "mex.h"

/*
#ifndef mwSize
#define mwSize int
#endif 
*/

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    void *nan_mat, *nonnan_mat, *val;
    uint64_t nDims, n_element, i;
    mwSize *pDims;

    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
    if(nrhs!=2)
        mexErrMsgIdAndTxt("data:inputError", "Two input required.");
    if(nlhs!=1)
        mexErrMsgIdAndTxt("data:outputError", "One output required.");

    /*  create a pointer to the input matrix x */
    nan_mat = mxGetPr(prhs[0]);
    val = mxGetPr(prhs[1]);

    /*  get the dimensions of the matrix input x */
    nDims = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    pDims = (mwSize*) mxGetDimensions(prhs[0]);
    n_element = mxGetNumberOfElements(prhs[0]);
    mxClassID mDType = mxGetClassID(prhs[0]);
    
	switch (mDType){
	    case mxSINGLE_CLASS:
		    break;
	    case mxDOUBLE_CLASS:
		    break;
	    default:
		    mexErrMsgIdAndTxt("data:TypeError", "Only single and double types are supported!\n");
		    return;
	}
    /* a sparse array has only its nonzero values in memory; the value is read as an
       element of the input's class, so it needs at least that many bytes */
    if(mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("data:inputError", "The input must be a full array.");
    if(!(mxIsNumeric(prhs[1]) || mxIsLogical(prhs[1])) || mxIsSparse(prhs[1])
       || mxGetNumberOfElements(prhs[1]) * mxGetElementSize(prhs[1]) < mxGetElementSize(prhs[0]))
        mexErrMsgIdAndTxt("data:inputError", "The value must be a full numeric scalar of the class of the input or double.");

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateUninitNumericArray(nDims, pDims, mDType, mxREAL);
    nonnan_mat = mxGetPr(plhs[0]);

    #pragma omp parallel for
    for (uint64_t i=0;i<n_element;i++) {
	    switch (mDType){
	        case mxSINGLE_CLASS:
		        if (isnan(*((float*)(nan_mat) + i)))
		            *((float*)nonnan_mat + i) = *((float*)val);
		        else
		            *((float*)nonnan_mat + i) = *((float*)(nan_mat) + i);
		        break;
	        case mxDOUBLE_CLASS:
		        if (isnan(*((double*)(nan_mat) + i)))
		            *((double*)nonnan_mat + i) = *((double*)val);
		        else
		            *((double*)nonnan_mat + i) = *((double*)(nan_mat) + i);
		        break;
	    }
    }

    // memcpy(mxGetPr(plhs[0]), nan_mat, n_element * sizeof(double) ) ;

    /*  create a C pointer to a copy of the output matrix */
}
