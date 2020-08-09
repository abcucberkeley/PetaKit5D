/*
 * adtest_mex.c - perform adtest
 *
 * The testing function is copied and modified from Francis Agust's stat.h.
 * The gateway function is implemented by Xiongtao Ruan.
 * 
 * This is a MEX-file for MATLAB.
 */ 

#include <string.h> 
#include <math.h>
#include <mex.h>

#define NA 7 // number of allowed alpha values

int compDouble(const void *x, const void *y) {

    double dx = *(double *)x;
    double dy = *(double *)y;

    if (dx < dy) {
        return -1;
    } else if (dx > dy) {
        return +1;
    }
    return 0;
}

unsigned char adtest(const double *data, const int N, const int adCase, const double mu, const double sigma, const double alpha) {

    if (adCase<0 || adCase>4) {
        mexErrMsgTxt("'adCase' must be an integer between 0 and 4.");
    }

    double *sdata = (double*)malloc(sizeof(double)*N);
    memcpy(sdata, data, N*sizeof(double));
    qsort(sdata, N, sizeof(double), compDouble);

    double A2 = 0.0;
    double r2 = sqrt(2.0);

    int i;

    /* Normal CDF */
    double *z = (double*)malloc(sizeof(double)*N);
    for (i=0;i<N;++i) {
        z[i] = 0.5 * (1 + erf((sdata[i]-mu)/(r2*sigma)));
    }

    /* A-D test statistic */
    for (i=0;i<N;++i) {
        A2 += (2.0*i+1.0)*(log(z[i]) + log(1.0-z[N-1-i]));
    }
    A2 = (-1.0*N) - A2/(1.0*N);

    if (adCase==3) {
        A2 *= 1.0 + 0.75/N + 2.25/(N*N);
    }
    if (adCase==4) {
        A2 *= 1.0 + 0.6/N;
    }

    free(z);
    free(sdata);

    /* Look-up table for critical values */
    static const double alphaVec[NA] = {0.25,  0.15,  0.10,  0.05,  0.025, 0.01, 0.005};
    static const double ctable[35] =  {1.248, 1.610, 1.933, 2.492, 3.070, 3.857, 4.500,  /* case 0: normal, mu/sigma known */
                                       0.644, 0.782, 0.894, 1.087, 1.285, 1.551, 1.756,  /* case 1: normal, mu unknown, sigma known */
                                       1.072, 1.430, 1.743, 2.308, 2.898, 3.702, 4.324,  /* case 2: normal, mu known, sigma unknown */
                                       0.470, 0.561, 0.631, 0.752, 0.873, 1.035, 1.159,  /* case 3: normal, mu/sigma unknown */
                                       0.736, 0.916, 1.062, 1.321, 1.591, 1.959, 2.244}; /* case 4: exponential, µ unknown */

    /* get critical value from lookup table */
    int ai = NA; // number of values in alpha vector + 1
    for (i=0;i<NA;++i) {
        if ( fabs(alpha-alphaVec[i]) < 1e-10 ) {
            ai = i;
            break;
        }
    }
    if (ai==NA) {
        mexErrMsgTxt("Admissible 'alpha' values: 0.25,  0.15,  0.10,  0.05,  0.025, 0.01, 0.005");
    }
    return A2>ctable[ai + NA*adCase];
}


unsigned char adtest_sorted_data(const double *data, const int N, const int adCase, const double mu, const double sigma, const double alpha) {

    if (adCase<0 || adCase>4) {
        mexErrMsgTxt("'adCase' must be an integer between 0 and 4.");
    }

    double A2 = 0.0;
    double r2 = sqrt(2.0);

    int i;

    /* Normal CDF */
    double *z = (double*)malloc(sizeof(double)*N);
    for (i=0;i<N;++i) {
        z[i] = 0.5 * (1 + erf((data[i]-mu)/(r2*sigma)));
    }

    /* A-D test statistic */
    for (i=0;i<N;++i) {
        A2 += (2.0*i+1.0)*(log(z[i]) + log(1.0-z[N-1-i]));
    }
    A2 = (-1.0*N) - A2/(1.0*N);

    if (adCase==3) {
        A2 *= 1.0 + 0.75/N + 2.25/(N*N);
    }
    if (adCase==4) {
        A2 *= 1.0 + 0.6/N;
    }

    // free(z);
    // free(sdata);

    /* Look-up table for critical values */
    static const double alphaVec[NA] = {0.25,  0.15,  0.10,  0.05,  0.025, 0.01, 0.005};
    static const double ctable[35] =  {1.248, 1.610, 1.933, 2.492, 3.070, 3.857, 4.500,  /* case 0: normal, mu/sigma known */
                                       0.644, 0.782, 0.894, 1.087, 1.285, 1.551, 1.756,  /* case 1: normal, mu unknown, sigma known */
                                       1.072, 1.430, 1.743, 2.308, 2.898, 3.702, 4.324,  /* case 2: normal, mu known, sigma unknown */
                                       0.470, 0.561, 0.631, 0.752, 0.873, 1.035, 1.159,  /* case 3: normal, mu/sigma unknown */
                                       0.736, 0.916, 1.062, 1.321, 1.591, 1.959, 2.244}; /* case 4: exponential, µ unknown */

    /* get critical value from lookup table */
    int ai = NA; // number of values in alpha vector + 1
    for (i=0;i<NA;++i) {
        if ( fabs(alpha-alphaVec[i]) < 1e-10 ) {
            ai = i;
            break;
        }
    }
    if (ai==NA) {
        mexErrMsgTxt("Admissible 'alpha' values: 0.25,  0.15,  0.10,  0.05,  0.025, 0.01, 0.005");
    }
    return A2>ctable[ai + NA*adCase];
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *x_mat, *mu, *sigma, *alpha, *mode, *h;
  int N;
  
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=5) 
    mexErrMsgTxt("Five input required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /*  create a pointer to the input matrix x */
  x_mat = mxGetPr(prhs[0]);
  mu = mxGetPr(prhs[1]);
  sigma = mxGetPr(prhs[2]);
  alpha = mxGetPr(prhs[3]);
  mode = mxGetPr(prhs[4]);
  
  /*  get the dimensions of the matrix input x */
  N = mxGetNumberOfElements(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleScalar(0);
    
  h = mxGetPr(plhs[0]);
  *h = adtest_sorted_data(x_mat, N, *mode, *mu, *sigma, *alpha);
    
  /*  create a C pointer to a copy of the output matrix */
}

