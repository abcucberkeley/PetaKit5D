/*
 * compute_gradient_hessian_mex.c - calculate the cross terms between different Gaussians in the Hessian
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
    double *x_c_mat,*y_c_mat,*z_c_mat,*A_h_mat,*hessian;
    double x_c_i, y_c_i, z_c_i, h_A_i, x_c_j, y_c_j, z_c_j, h_A_ij;
    int mrows, ncols, n, s, s_t, row, i, j, i_row, j_row;

    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
    if(nrhs!=4) 
        mexErrMsgTxt("Four input required.");
    if(nlhs!=1) 
        mexErrMsgTxt("One output required.");

    /*  create a pointer to the input matrix x */
    x_c_mat = mxGetPr(prhs[0]);
    y_c_mat = mxGetPr(prhs[1]);
    z_c_mat = mxGetPr(prhs[2]);
    A_h_mat = mxGetPr(prhs[3]);

    /*  get the dimensions of the matrix input x */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    n = 3 * ncols;

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);

    hessian = mxGetPr(plhs[0]);

    for (row=0;row<mrows;row++) {
        for (i=0; i<ncols-1; i++) {
            i_row = i * mrows + row;
            x_c_i = *(x_c_mat + i_row);
            y_c_i = *(y_c_mat + i_row);
            z_c_i = *(z_c_mat + i_row);
            h_A_i = *(A_h_mat + i_row);
            
            for (j=i+1; j<ncols; j++) {
                j_row = j * mrows + row;
                x_c_j = *(x_c_mat + j_row);
                y_c_j = *(y_c_mat + j_row);
                z_c_j = *(z_c_mat + j_row);
                
                h_A_ij = 2 * h_A_i * *(A_h_mat + j_row);
                
                s = 3 * j * n + 3 * i;
                *(hessian + s) += h_A_ij * x_c_i * x_c_j;
                *(hessian + s + 1) += h_A_ij * y_c_i * x_c_j;
                *(hessian + s + 2) += h_A_ij * z_c_i * x_c_j;
                *(hessian + s + n) += h_A_ij * x_c_i * y_c_j;
                *(hessian + s + n + 1) += h_A_ij * y_c_i * y_c_j;
                *(hessian + s + n + 2) += h_A_ij * z_c_i * y_c_j;
                *(hessian + s + 2 * n) += h_A_ij * x_c_i * z_c_j;
                *(hessian + s + 2 * n + 1) += h_A_ij * y_c_i * z_c_j;
                *(hessian + s + 2 * n + 2) += h_A_ij * z_c_i * z_c_j;
            }
        }
    }
    
    // get the transposed parts
    for (i=0; i<ncols-1; i++) {
        i_row = i * mrows + row;

        for (j=i+1; j<ncols; j++) {
            j_row = j * mrows + row;

            s = 3 * j * n + 3 * i;
            s_t = 3 * i * n + 3 * j;
            
            *(hessian + s_t) += *(hessian + s);
            *(hessian + s_t + n) += *(hessian + s + 1);
            *(hessian + s_t + 2 * n) += *(hessian + s + 2);
            *(hessian + s_t + 1) += *(hessian + s + n);
            *(hessian + s_t + n + 1) += *(hessian + s + n + 1);
            *(hessian + s_t + 2 * n + 1) += *(hessian + s + n + 2);
            *(hessian + s_t + 2) += *(hessian + s + 2 * n);
            *(hessian + s_t + n + 2) += *(hessian + s + 2 * n + 1);
            *(hessian + s_t + 2 * n + 2) += *(hessian + s + 2 * n + 2);
        }
    }
    
}
