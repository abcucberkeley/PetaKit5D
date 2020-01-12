/*
 * compute_gradient_hessian_mex.c - compute gradient and hessian for gaussian curve fitting
 *
 * Implemented by Xiongtao Ruan.
 * 
 * This is a MEX-file for MATLAB.
 */ 

#include "math.h"
#include "mex.h"

#ifndef mwSize
#define mwSize int
#endif 


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *x_mat,*y_mat,*z_mat,*beta_xyz,*sigma_xy,*sigma_z,*h_mat;
  double beta_x, beta_y, beta_z;
  double x, y, z;
  int mrows,ncols,row;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=6) 
    mexErrMsgTxt("Six input required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /*  create a pointer to the input matrix x */
  x_mat = mxGetPr(prhs[0]);
  y_mat = mxGetPr(prhs[1]);
  z_mat = mxGetPr(prhs[2]);
  beta_xyz = mxGetPr(prhs[3]);
  sigma_xy = mxGetPr(prhs[4]);
  sigma_z = mxGetPr(prhs[5]);
  
  /*  get the dimensions of the matrix input x */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,1, mxREAL);
    
  h_mat = mxGetPr(plhs[0]);
  
  beta_x = *beta_xyz;
  beta_y = *(beta_xyz + 1);
  beta_z = *(beta_xyz + 2);
  
  for (row=0;row<mrows;row++) {
      x = *(x_mat + row);
      y = *(y_mat + row);
      z = *(z_mat + row);
        
      // update for h_mat
      *(h_mat + row) = exp(-0.5 * (((x - beta_x) * (x - beta_x) + (y - beta_y) * (y - beta_y)) / (*sigma_xy * *sigma_xy) + (z - beta_z) * (z - beta_z) / (*sigma_z * *sigma_z))); 
  }
    
  /*  create a C pointer to a copy of the output matrix */
}
