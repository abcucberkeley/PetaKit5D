/*
 * compute_gradient_hessian_mex.c - compute gradient and hessian for gaussian curve fitting
 *
 * Implemented by Xiongtao Ruan.
 * 
 * This is a MEX-file for MATLAB.
 */ 

#include "math.h"
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
    double *nan_mat, *nonnan_mat, *val;
    int n_element,i;
    mwSize *pDims;

    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
    if(nrhs!=2) 
        mexErrMsgTxt("Two input required.");
    if(nlhs!=1) 
        mexErrMsgTxt("One output required.");

    /*  create a pointer to the input matrix x */
    nan_mat = mxGetPr(prhs[0]);
    val = mxGetPr(prhs[1]);

    /*  get the dimensions of the matrix input x */
    pDims = mxGetDimensions(prhs[0]);
    n_element = mxGetNumberOfElements(prhs[0]);

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateUninitNumericArray(3, pDims, mxDOUBLE_CLASS, mxREAL);

    nonnan_mat = mxGetPr(plhs[0]);

    for (i=0;i<n_element;i++) {
        if (isnan(*(nan_mat + i)))
            *(nonnan_mat + i) = *val;
        else
            *(nonnan_mat + i) = *(nan_mat + i);
    }
    
    // memcpy(mxGetPr(plhs[0]), nan_mat, n_element * sizeof(double) ) ;

    /*  create a C pointer to a copy of the output matrix */
}
