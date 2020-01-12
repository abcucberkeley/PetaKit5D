/* ------------------------------------------------------- */
/*                                                         */
/* distmat.c [distance matrix - function definitions]      */
/*                                                         */
/* ------------------------------------------------------- */
/*                                                         */
/* First version: Aaron Ponti - 02/08/28                   */
/*                                                         */
/* ------------------------------------------------------- */

#include "distmat.h"

/* #ifndef DISTMAT_H_
//     
//     #define DISTMAT_H_
// 
//     #include "distmat.h"
//    
// #endif
*/
/* calcDistMatrix1D */

void calcDistMatrix1D(double *D, double *M, double *N, int Mrows, int Nrows)
{
    int i,j;
    
    for (i=0;i<Nrows;i++) {
        
        for (j=0;j<Mrows;j++) {
        
            /* Store calculated distance in D */
            *D++=*(N+i)-*(M+j);            
        }
    }
}

/* calcDistMatrix2D */

void calcDistMatrix2D(double *D, double *M, double *N, int Mrows, int Nrows)
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
            
             /* Store calculated distance in D */
            *D++=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX));           
        }
    }
}

/* calcDistMatrix3D */

void calcDistMatrix3D(double *D, double *M, double *N, int Mrows, int Nrows)
{
    double	mX=0, mY=0, mZ=0;
    double	nX=0, nY=0, nZ=0;
    double  *posM, *posN;
    int  i,j;
    
    for (i=0;i<Nrows;i++) {
        
        /* Get source position */
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
		nZ=*(posN+2*Nrows);
        
        for (j=0;j<Mrows;j++) {
        
            /* Get target position */
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
			mZ=*(posM+2*Mrows);
            
            /* Store calculated distance in D */
            *D++=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX)+(nZ-mZ)*(nZ-mZ));
            
        }
    }
}
