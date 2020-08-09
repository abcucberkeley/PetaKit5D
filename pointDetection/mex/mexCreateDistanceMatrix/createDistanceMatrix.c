/*
 *  createDistanceMatrix [MATLAB C-MEX]
 *
 *
 * Files:
 *
 * createDistanceMatrix.c - MEX interface
 * distmat.h              - function prototypes
 * distmat.c              - function definitions
 *
 * See createDistanceMatrix.m for detailed help.
 *
 * First version: Aaron Ponti - 02/08/28
 *
 * Compilation:
 * Mac/Linux: mex distmat.c  -output createDistanceMatrix createDistanceMatrix.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" distmat.c -output createDistanceMatrix createDistanceMatrix.c
 */

#include "mex.h"
#include "distmat.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    /* Initialize pointers for 2 inputs and 1 output */
    double *M, *N;
	double *D;

	/* Initialize int variables to store matrix dimensions */
	int Mrows,Mcols;
	int Nrows,Ncols;

    /* Check that the number of input and output parameters is valid */
	if(nrhs != 2)
		mexErrMsgTxt("Two input parameters required.");
	if(nlhs > 1)
		mexErrMsgTxt("One output parameter required.");	
	
	/* Read input parameter dimensions */
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);
	
	/* Check input parameter dimension */
	if ((Mcols>3) || (Ncols>3))
		mexErrMsgTxt("Point coordinates in more than 3 dimensions are not supported.");
	
	if (Mcols!=Ncols)
		mexErrMsgTxt("The points in the coordinate matrices have different number of dimensions.");

	/* Create matrix for the return argument D */
	plhs[0]=mxCreateDoubleMatrix(Mrows,Nrows, mxREAL);
    
    /* Assign pointers to each input and output */
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);
	D=mxGetPr(plhs[0]);

   	/* Call the correcponding C function */
	if (Mcols==1) { calcDistMatrix1D(D,M,N,Mrows,Nrows); }
	if (Mcols==2) { calcDistMatrix2D(D,M,N,Mrows,Nrows); }
	if (Mcols==3) { calcDistMatrix3D(D,M,N,Mrows,Nrows); }
	
}
