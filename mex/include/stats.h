/* Collection of functions for computing statistical tests
 *
 * References:
 * [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions",
 *     National Bureau of Standards, 1972
 *
 * (c) Francois Aguet, 2011 (last modified May 23, 2012)
 * */

#ifndef STATS_H
#define STATS_H

#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#if defined(_WIN32) || defined(_WIN64)
    #define fmax max
    #define fmin min
    #define erf gsl_sf_erf
    #include <minmax.h>
    #include <gsl/gsl_sf_erf.h>
#endif

#define NA 7 // number of allowed alpha values


// CDF of the F distribution        
// [1]; see also fcdf.m (Matlab built-in)
double fcdf(const double x, const double v1, const double v2) {        
    double xx, p;
    if (v2 <= x*v1) {
        xx = v2/(v2+x*v1);
        // Upper incomplete beta
        p = 1.0 - gsl_sf_beta_inc(0.5*v2, 0.5*v1, xx);
    } else {
        double num = v1*x;
        xx = num/(num+v2);
        // Lower incomplete beta
        p = gsl_sf_beta_inc(0.5*v1, 0.5*v2, xx);
    }
    return p;
}


/* Code from:
 * Marsaglia, G., W.W. Tsang, and J. Wang (2003), "Evaluating Kolmogorov's
 * Distribution", Journal of Statistical Software, vol. 8, issue 18.
 * This is the implementation used in kstest.m in Matlab.
 */
void mMultiply(double *A, double *B, double *C, int m)
{ 
    int i,j,k;
    double s;
    for(i=0;i<m;i++)
        for(j=0; j<m; j++) {
        s=0.0; 
        for(k=0;k<m;k++) 
            s+=A[i*m+k]*B[k*m+j];
        C[i*m+j]=s;
        }
}

void mPower(double *A, int eA, double *V, int *eV, int m, int n)
{ 
    double *B;
    int eB,i;
    if(n==1) {
        for (i=0;i<m*m;i++)
            V[i]=A[i];
        *eV=eA; 
        return;
    } 
    mPower(A, eA, V, eV, m, n/2);
    B = (double*)malloc((m*m)*sizeof(double));
    mMultiply(V, V, B, m);
    eB = 2*(*eV);
    if (n%2==0) {
        for (i=0;i<m*m;i++)
            V[i]=B[i];
        *eV=eB;
    } else {
        mMultiply(A,B,V,m);
        *eV=eA+eB;
    } 
    if (V[(m/2)*m+(m/2)]>1e140) {
        for(i=0;i<m*m;i++)
            V[i] = V[i]*1e-140;
        *eV+=140;
    }
    free(B);
}

double K(int n, double d) {
    int k, m, i, j, g, eH, eQ;
    double h, s, *H, *Q;
    /* OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL */
    s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt((double)n)+1.409/n)*s);
    k = (int)(n*d)+1;
    m=2*k-1;
    h=k-n*d;
    H = (double*)malloc((m*m)*sizeof(double));
    Q = (double*)malloc((m*m)*sizeof(double));
    for (i=0;i<m;i++)
        for(j=0;j<m;j++)
            if(i-j+1<0)
                H[i*m+j]=0;
            else
                H[i*m+j]=1;
    for(i=0;i<m;i++) {
        H[i*m] -= pow(h, i+1);
        H[(m-1)*m+i] -= pow(h, (m-i));
    }
    H[(m-1)*m] += (2*h-1>0 ? pow(2*h-1, m) : 0);
    for (i=0;i<m;i++)
        for (j=0;j<m;j++)
            if (i-j+1>0)
                for (g=1;g<=i-j+1;g++)
                    H[i*m+j]/=g;
    eH=0;
    mPower(H, eH, Q, &eQ, m, n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++) {
        s = s*i/n;
        if (s<1e-140) {
            s*=1e140;
            eQ-=140;
        }
    }
    s*=pow(10., eQ);
    free(H);
    free(Q);
    return s;
}



/* Adapted from "Numerical Recipes in C" 3rd ed. pp. 334-335 & 736-738 */
double pks(double z) {
    if (z < 0.0)
        mexErrMsgTxt("z value for KS dist. must be positive.");
    if (z == 0.0)
        return 0.0;
    double z2 = z*z;
    if (z < 1.18) {
        double y = exp(-1.23370055013616983/z2); /* exp(-pi^2/(8*z^2)) */
        return 2.506628274631 / z * (y + pow(y,9) + pow(y,25) + pow(y,49)); /* sqrt(2*pi) */
    } else {
        double x = exp(-2*z2);
        return 1.0 - 2.0*(x - pow(x,4) - pow(x,9));
    }
}

/* Q_ks = 1 - P_ks */
double qks(double z) {
    if (z < 0.0)
        mexErrMsgTxt("z value for KS dist. must be positive.");
    if (z == 0.0)
        return 1.0;
    double z2 = z*z;
    if (z < 1.18) {
        double y = exp(-1.23370055013616983/z2); /* exp(-pi^2/(8*z^2)) */
        return 1.0 - 2.506628274631 / z * (y + pow(y,9) + pow(y,25) + pow(y,49)); /* sqrt(2*pi) */
    } else {
        double x = exp(-2*z2);
        return 2.0*(x - pow(x,4) - pow(x,9));
    }
}


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


/* adapted from Numerical Recipes in C */
double ksone(double *data, int N, double mu, double sigma) {
    
    double *sdata = (double*)malloc(sizeof(double)*N);
    memcpy(sdata, data, N*sizeof(double));
    
    qsort(sdata, N, sizeof(double), compDouble);
    
    double normCDF;
    int j;
    double D = 0.0;
    double dt, fn, fo = 0.0;
    
    for (j=0; j<N; ++j) {
        fn = (j+1.0)/N; /* CDF */
        normCDF = 0.5 * (1.0 + erf((sdata[j]-mu)/(sqrt(2.0)*sigma)));
        
        dt = fmax(fabs(fo - normCDF), fabs(fn - normCDF));
        if (dt > D)
            D = dt;
        fo = fn;
    }
    free(sdata);
    /* return 1.0 - K(N,D); this is a much more accurate solution, but also much slower */
    double en = sqrt((double)N);
    return qks((en + 0.12 + 0.11/en)*D);
}





/* Anderson-Darling test
 *
 * For the test and its derivation, see
 * [1] Anderson & Darling, Ann. Math. Stat. 23, 1952
 * [2] Anderson & Darling, J. Am. Stat. Assoc. 49, 1954
 *
 * Critical values taken from 
 * [3] R.B. D'Agostino and M.A. Stephens, Goodness-of-Fit Techniques, ed. M. Dekker, Inc., 1986
 */
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
#endif // STATS_H
