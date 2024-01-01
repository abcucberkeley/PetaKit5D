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
  double *x_c_mat,*y_c_mat,*z_c_mat,*h_A_mat,*g_mat,*sigma_xy,*sigma_z,*grad,*hessian;
  double x_c, y_c, z_c, h_A, g, w_1, w_2, w_3;
  int mrows,ncols,row;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=7) 
    mexErrMsgTxt("Seven input required.");
  if(nlhs!=2) 
    mexErrMsgTxt("Two output required.");
  
  /*  create a pointer to the input matrix x */
  x_c_mat = mxGetPr(prhs[0]);
  y_c_mat = mxGetPr(prhs[1]);
  z_c_mat = mxGetPr(prhs[2]);
  h_A_mat = mxGetPr(prhs[3]);
  g_mat = mxGetPr(prhs[4]);
  sigma_xy = mxGetPr(prhs[5]);
  sigma_z = mxGetPr(prhs[6]);
  
  /*  get the dimensions of the matrix input x */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(3,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
    
  grad = mxGetPr(plhs[0]);
  hessian = mxGetPr(plhs[1]);
  
  for (row=0;row<mrows;row++) {
      x_c = *(x_c_mat + row);
      y_c = *(y_c_mat + row);
      z_c = *(z_c_mat + row);
      h_A = *(h_A_mat + row);
      g = *(g_mat + row);
      
      // intermediate temperate variables
      w_1 = 2 * h_A * h_A;
      w_2 = 2 * h_A * g;
      w_3 = w_1 - w_2;
      
      // update for gradient      
      *grad += -w_2 * x_c;
      *(grad + 1) += -w_2 * y_c;
      *(grad + 2) += -w_2 * z_c;
  
      // update for hessian
      *hessian += w_3 * x_c * x_c + w_2 / (*sigma_xy * *sigma_xy);
      *(hessian + 1) += w_3 * x_c * y_c;
      *(hessian + 2) += w_3 * x_c * z_c;
      *(hessian + 4) += w_3 * y_c * y_c + w_2 / (*sigma_xy * *sigma_xy);
      *(hessian + 5) += w_3 * y_c * z_c;
      *(hessian + 8) += w_3 * z_c * z_c + w_2 / (*sigma_z * *sigma_z);
  }
  
  // get the transposed terms  
  *(hessian + 3) = *(hessian + 1);
  *(hessian + 6) = *(hessian + 2);
  *(hessian + 7) = *(hessian + 5);
  
  /*  create a C pointer to a copy of the output matrix */
}
