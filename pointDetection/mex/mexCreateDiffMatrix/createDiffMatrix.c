/* MATLAB C-MEX
 *
 * createDistanceMatrix.c *
 *
 * First version: Aaron Ponti - 02/11/26
 *
 * See createDiffMatrix.m for detailed help.
 *
 * Compilation:
 * Mac/Linux: mex  createDiffMatrix.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -output createDiffMatrix createDiffMatrix.cpp
 */


#include "mex.h"

void calcDiffMatrix(double *dX, double *dY, double *M, double *N, int Mrows, int Nrows)
{
    double	mX=0, mY=0;
    double	nX=0, nY=0;
    double  *posM, *posN;
    int  i,j;
    
    for (i=0;i<Nrows;i++) {
        
        /* Get source position */
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
        
        for (j=0;j<Mrows;j++) {
        
            /* Get target position */
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);

			/* Calculate differences */
			*dX++=nX-mX;
			*dY++=nY-mY;

		}
    }
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    /* Initialize pointers for 2 inputs and 1 output */
    double *M, *N;
	double *dX, *dY;

	/* Initialize int variables to store matrix dimensions */
	int Mrows,Mcols;
	int Nrows,Ncols;

    /* Check that the number of input and output parameters is valid */
	if(nrhs != 2)
		mexErrMsgTxt("Two input parameters required.");
	if(nlhs != 2)
		mexErrMsgTxt("Two output parameters required.");

	/* Read input parameter dimensions */
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);

	
	/* Check input parameter dimension */
	if ((Mcols!=2) || (Ncols!=2))
		mexErrMsgTxt("Only 2D point coordinates are supported.");
	
	/* Create matrix for the return arguments dX and dY */
	plhs[0]=mxCreateDoubleMatrix(Mrows,Nrows, mxREAL);
    plhs[1]=mxCreateDoubleMatrix(Mrows,Nrows, mxREAL);
    
    /* Assign pointers to each input and output */
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);
	dX=mxGetPr(plhs[0]);
	dY=mxGetPr(plhs[1]);

   	/* Call the C function */
	calcDiffMatrix(dX,dY,M,N,Mrows,Nrows);
	
}
