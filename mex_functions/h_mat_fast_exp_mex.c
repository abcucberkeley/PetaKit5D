/*
 * h_mat_fast_exp_mex.c - compute h_mat for gaussian curve fitting efficiently
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

 /* fast exp function is from 
  * https://stackoverflow.com/questions/10552280/fast-exp-calculation-possible-to-improve-accuracy-without-losing-too-much-perfo
  */
 /* max. rel. error <= 1.73e-3 on [-87,88] */
 float fast_exp (float x)
 {
   volatile union {
     float f;
     unsigned int i;
   } cvt;

   /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
   float t = x * 1.442695041f;
   float fi = floorf (t);
   float f = t - fi;
   int i = (int)fi;
   cvt.f = (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
   cvt.i += (i << 23);                                          /* scale by 2^i */
   return cvt.f;
 }



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
      *(h_mat + row) = fast_exp(-0.5 * (((x - beta_x) * (x - beta_x) + (y - beta_y) * (y - beta_y)) / (*sigma_xy * *sigma_xy) + (z - beta_z) * (z - beta_z) / (*sigma_z * *sigma_z))); 
  }
    
}
