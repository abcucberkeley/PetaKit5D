#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    float* im = (float*)mxGetPr(prhs[0]);
    double xstep = (double)*mxGetPr(prhs[1]);
    uint64_t nint = (uint64_t)*mxGetPr(prhs[2]);
    uint8_t Reverse = (uint8_t)mxIsLogicalScalarTrue(prhs[3]);
    
    uint64_t* sz = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dim[3];
    dim[0] = sz[0];
    dim[1] = sz[1];
    dim[2] = (sz[2] - 1) * nint + 1;
    plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
    
    float* im_int = (float*)mxGetPr(plhs[0]);
    
    uint64_t xa = (uint64_t)ceil(xstep);
    double* s_mat = (double*)malloc((nint-1)*sizeof(double));
    double* t_mat = (double*)malloc((nint-1)*sizeof(double));
    double* sw_mat = (double*)malloc((nint-1)*sizeof(double));
    double* tw_mat = (double*)malloc((nint-1)*sizeof(double));
    
    
    for(uint64_t i = 1; i < nint; i++){
        double s;
        double t;
        if(Reverse){
            s = ceil(xstep) - xstep * (double)i / (double)nint;
            t = xstep * (double)(nint - i) / (double)nint;
        }
        else{
            s = xstep * (double)(nint - i) / (double)nint;
            t = ceil(xstep) - xstep * (double)i / (double)nint;
        }
        
        //distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);
        
        s_mat[i - 1] = s;
        t_mat[i - 1] = t;
        sw_mat[i - 1] = sw;
        tw_mat[i - 1] = tw;
    }
    
    
    
    #pragma omp parallel for
    for(uint64_t z = 0; z < sz[2]; z++){
        if(z == sz[2]-1){
            uint64_t zIndex = z * nint;
            for(uint64_t i = 0; i < sz[1]; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    im_int[j+(i*sz[0])+(zIndex*sz[1]*sz[0])] = im[j+(i*sz[0])+(z*sz[1]*sz[0])];
                }
            }
            continue;
        }
        
        float* im_s = (float*)malloc(sz[0]*(sz[1]+xa)*sizeof(float));
        float* im_t = (float*)malloc(sz[0]*(sz[1]+xa)*sizeof(float));
        
        if(Reverse){
            for(uint64_t i = 0; i < sz[1]+xa; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    //im_s[j+((i+xa)*sz[0])] = im[j+(i*sz[0])+(z*sz[1]*sz[0])];
                    //im_t[j+(i*sz[0])] = im[j+(i*sz[0])+((z+1)*sz[1]*sz[0])];
                    if(i < xa) im_s[j+(i*sz[0])] = 0;
                    else im_s[j+(i*sz[0])] = im[j+((i-xa)*sz[0])+(z*sz[1]*sz[0])];
                    if(i>=sz[1]) im_t[j+(i*sz[0])] = 0;
                    else im_t[j+(i*sz[0])] = im[j+(i*sz[0])+((z+1)*sz[1]*sz[0])];
                }
            }
        }
        else{
            for(uint64_t i = 0; i < sz[1]+xa; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    //im_s[j+(i*sz[0])] = im[j+(i*sz[0])+(z*sz[1]*sz[0])];
                    //im_t[j+((i+xa)*sz[0])] = im[j+(i*sz[0])+((z+1)*sz[1]*sz[0])];
                    if(i>=sz[1]) im_s[j+(i*sz[0])] = 0;
                    else im_s[j+(i*sz[0])] = im[j+(i*sz[0])+((z+1)*sz[1]*sz[0])];
                    if(i < xa) im_t[j+(i*sz[0])] = 0;
                    else im_t[j+(i*sz[0])] = im[j+((i-xa)*sz[0])+(z*sz[1]*sz[0])];
                }
            }
        }
        
        
        for(uint64_t ind = 0; ind < nint; ind++){
            uint64_t zint = z * nint + ind;
            
            if(!ind){
                for(uint64_t i = 0; i < sz[1]; i++){
                    for(uint64_t j = 0; j < sz[0]; j++){
                        im_int[j+(i*sz[0])+(zint*sz[1]*sz[0])] = im[j+(i*sz[0])+(z*sz[1]*sz[0])];
                    }
                }
                continue;
            }
            
            
            uint64_t s = (uint64_t)s_mat[ind - 1];
            uint64_t t = (uint64_t)t_mat[ind - 1];
            float sw = sw_mat[ind - 1];
            float tw = tw_mat[ind - 1];
            
            for(uint64_t i = 0; i < sz[1]; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    im_int[j+(i*sz[0])+(zint*sz[1]*sz[0])] =
                            ((im_s[j+((i+s)*sz[0])] * (1.0-sw) + im_s[j+((i+s+1)*sz[0])] * sw)/2.0)
                            + ((im_t[j+((i+t)*sz[0])] * (1.0-tw) + im_t[j+((i+t+1)*sz[0])] * tw)/2.0);
                }
            }
        }
        
        free(im_s);
        free(im_t);
    }
    free(s_mat);
    free(t_mat);
    free(sw_mat);
    free(tw_mat);
}