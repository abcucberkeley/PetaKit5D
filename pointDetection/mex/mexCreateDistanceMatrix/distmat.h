/* ------------------------------------------------------- */
/*                                                         */
/* distmat.h [distance matrix - function prototypes]       */
/*                                                         */
/* ------------------------------------------------------- */
/*                                                         */
/* First version: Aaron Ponti - 02/08/28                   */
/*                                                         */
/* ------------------------------------------------------- */


#ifndef DISTMAT_H_

	#define DISTMAT_H_

	/* Include math.h */
	#include <math.h>

	/* Function prototypes */

	/* calcDistMatrix1D */
	void calcDistMatrix1D(double *D, double *M, double *N, int Mrows, int Nrows);
	
	/* calcDistMatrix2D */
	void calcDistMatrix2D(double *D, double *M, double *N, int Mrows, int Nrows);
	
	/* calcDistMatrix3D */
	void calcDistMatrix3D(double *D, double *M, double *N, int Mrows, int Nrows);

#endif

