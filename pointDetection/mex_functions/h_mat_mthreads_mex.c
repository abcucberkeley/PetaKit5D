/*
 * compute_gradient_hessian_mex.c - compute gradient and hessian for gaussian curve fitting
 * This is a multi-threading version of the computing
 *
 * Implemented by Xiongtao Ruan. Dec 26, 2019
 *
 * Mutli-threading mechanism is adapted from Yair Altman's Max-in-place package. 
 * http://UndocumentedMatlab.com/blog/multi-threaded-mex
 * 
 */ 

#include "math.h"
#include "mex.h"

/* undef needed for LCC compiler */
#undef EXTERN_C
#ifdef _WIN32
    #include <windows.h>
    #include <process.h>
#else
    #include <pthread.h>
#endif

/* Macros */
#if !defined(MAX)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif



/* Main processing loop function */
void main_loop(const double beta_x, const double beta_y, const double beta_z, 
        const double sigma_xy_sq, const double sigma_z_sq, 
        const double *x_mat, const double *y_mat, const double *z_mat, 
        double *h_mat, int startIdx, int endIdx)
{
    double x, y, z;
    
    /* Loop through all matrix coordinates */
    for (int idx=startIdx; idx<=endIdx; idx++)
    {
        x = *(x_mat + idx);
        y = *(y_mat + idx);
        z = *(z_mat + idx);

        // update for h_mat
        *(h_mat + idx) = exp(-0.5 * (((x - beta_x) * (x - beta_x) + 
                (y - beta_y) * (y - beta_y)) / sigma_xy_sq + 
                (z - beta_z) * (z - beta_z) / sigma_z_sq));         
    }
}

/* Computation function in threads */
#ifdef _WIN32
  unsigned __stdcall thread_func(void *ThreadArgs_) {
#else
  void thread_func(void *ThreadArgs_) {
#endif
    double **ThreadArgs = ThreadArgs_;  /* void* => double** */
    const mxArray** prhs = (const mxArray**) ThreadArgs[0];
    
    double beta_x = ThreadArgs[0][0];
    double beta_y = ThreadArgs[1][0];    
    double beta_z = ThreadArgs[2][0];        
    double sigma_xy_sq = ThreadArgs[3][0];
    double sigma_z_sq = ThreadArgs[4][0];    
    double *x_mat = ThreadArgs[5];
    double *y_mat = ThreadArgs[6];
    double *z_mat = ThreadArgs[7];
    double *h_mat = ThreadArgs[8];
    int ThreadID = (int) ThreadArgs[9][0];
    int startIdx = (int) ThreadArgs[10][0];
    int endIdx   = (int) ThreadArgs[11][0];
    /*mexPrintf("Starting thread #%d: idx=%d:%d\n", ThreadID, startIdx, endIdx); */

    /* Run the main processing function */
    main_loop(beta_x, beta_y, beta_z, sigma_xy_sq, sigma_z_sq, x_mat, y_mat, z_mat, h_mat, startIdx, endIdx);

    /* Explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
        _endthreadex( 0 );
        return 0;
    #else
        pthread_exit(NULL);
    #endif
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *x_mat,*y_mat,*z_mat,*beta_xyz,*sigma_xy,*sigma_z,*h_mat,*Nthreadsd; 
    double beta_x, beta_y, beta_z;
    double x, y, z;
    int mrows,ncols,row;

    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
    if(nrhs!=7) 
    mexErrMsgTxt("Seven input required.");
    if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

    /*  create a pointer to the input matrix x */
    x_mat = mxGetPr(prhs[0]);
    y_mat = mxGetPr(prhs[1]);
    z_mat = mxGetPr(prhs[2]);
    beta_xyz = mxGetPr(prhs[3]);
    sigma_xy = mxGetPr(prhs[4]);
    sigma_z = mxGetPr(prhs[5]);
    Nthreadsd = mxGetPr(prhs[6]);
    
    /* save time for sigma_xy and sigma_z */
    double sigma_xy_sq = *sigma_xy * *sigma_xy;
    double sigma_z_sq = *sigma_z * *sigma_z;

    /*  get the dimensions of the matrix input x */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(mrows,1, mxREAL);

    h_mat = mxGetPr(plhs[0]);

    /* Get the number of threads from the Matlab engine (maxNumCompThreads) */
//     mxArray *matlabCallOut[1] = {0};
//     mxArray *matlabCallIn[1]  = {0};
//     mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
//     double *Nthreadsd = mxGetPr(matlabCallOut[0]);
    int Nthreads = (int) Nthreadsd[0];
    
    int n1 = mrows;
    
    beta_x = *beta_xyz;
    beta_y = *(beta_xyz + 1);
    beta_z = *(beta_xyz + 2);

//     for (row=0;row<mrows;row++) {
//       x = *(x_mat + row);
//       y = *(y_mat + row);
//       z = *(z_mat + row);
// 
//       // update for h_mat
//       *(h_mat + row) = exp(-0.5 * (((x - beta_x) * (x - beta_x) + (y - beta_y) * (y - beta_y)) / (*sigma_xy * *sigma_xy) + (z - beta_z) * (z - beta_z) / (*sigma_z * *sigma_z))); 
//     }
    
    
    if (Nthreads == 1) {

        /* Process the inputs directly (not via a thread) */
        main_loop(beta_x, beta_y, beta_z, sigma_xy_sq, sigma_z_sq, x_mat, y_mat, z_mat, h_mat, 0, n1-1);

    } else {  /* multi-threaded */

        /* Allocate memory for handles of worker threads */
        #ifdef _WIN32
            HANDLE    *ThreadList = (HANDLE*)   malloc(Nthreads*sizeof(HANDLE));
        #else
            pthread_t *ThreadList = (pthread_t*)malloc(Nthreads*sizeof(pthread_t));
        #endif

        /* Allocate memory for the thread arguments (attributes) */
        double **ThreadID, **ThreadStartIdx, **ThreadEndIdx, ***ThreadArgs;
        double *ThreadID1, *ThreadStartIdx1, *ThreadEndIdx1, **ThreadArgs1;

        ThreadID       = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadStartIdx = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadEndIdx   = (double **) malloc( Nthreads* sizeof(double *) );
        ThreadArgs     = (double ***)malloc( Nthreads* sizeof(double **) );

        /* Launch the requested number of threads */
        int i;
        int threadBlockSize = ceil( ((double)n1) / Nthreads );
        for (i=0; i<Nthreads; i++)
        {
            /* Create thread ID */
            ThreadID1 = (double *)malloc( 1* sizeof(double) );
            ThreadID1[0] = i;
            ThreadID[i] = ThreadID1;

            /* Compute start/end indexes for this thread */
            ThreadStartIdx1 = (double *) malloc( sizeof(double) );
            ThreadStartIdx1[0] = i * threadBlockSize;
            ThreadStartIdx[i] = ThreadStartIdx1;

            ThreadEndIdx1 = (double *) malloc( sizeof(double) );
            ThreadEndIdx1[0] = MIN((i+1)*threadBlockSize, n1) - 1;
            ThreadEndIdx[i] = ThreadEndIdx1;

            /* Prepare thread input args */
            ThreadArgs1 = (double **) malloc( 12 * sizeof(double*) );
            ThreadArgs1[0] = & beta_x;
            ThreadArgs1[1] = & beta_y;
            ThreadArgs1[2] = & beta_z;
            ThreadArgs1[3] = & sigma_xy_sq;
            ThreadArgs1[4] = & sigma_z_sq;
            ThreadArgs1[5] = x_mat;
            ThreadArgs1[6] = y_mat;
            ThreadArgs1[7] = z_mat;
            ThreadArgs1[8] = h_mat;
            ThreadArgs1[9] = ThreadID[i];
            ThreadArgs1[10] = ThreadStartIdx[i];
            ThreadArgs1[11] = ThreadEndIdx[i];

            ThreadArgs[i] = ThreadArgs1;

            /* Launch the thread with its associated args */
            #ifdef _WIN32
                ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &thread_func, ThreadArgs[i], 0, NULL);
            #else
                pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &thread_func, ThreadArgs[i]);
            #endif
        }

        /* Wait for all the treads to finish working */
        #ifdef _WIN32
            for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
            for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
        #else
            for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
        #endif

        /* Free the memory resources allocated for the threads */
//         for (i=0; i<Nthreads; i++)
//         {
//             free(ThreadArgs[i]);
//             free(ThreadID[i]);
//             free(ThreadStartIdx[i]);
//             free(ThreadEndIdx[i]);
//         }
// 
//         free(ThreadArgs);
//         free(ThreadID );
//         free(ThreadStartIdx);
//         free(ThreadEndIdx);
//         free(ThreadList);
    }

    return;    
    
}
