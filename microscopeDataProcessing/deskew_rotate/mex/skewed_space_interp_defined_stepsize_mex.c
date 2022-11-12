#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType != mxSINGLE_CLASS){
        mexErrMsgIdAndTxt("skewed_space_interp:dataTypeError","Only Data of type Single is supported");
    }
    float* im = (float*)mxGetPr(prhs[0]);
    const double xstep = (double)*mxGetPr(prhs[1]);
    const double stepsize = (double)*mxGetPr(prhs[2]);
    const uint8_t Reverse = (uint8_t)mxIsLogicalScalarTrue(prhs[3]);

    const uint64_t* sz = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t nint = (uint64_t) (floor(round((double)(sz[2] - 1) / stepsize * 10000) / 10000) + 1);
    uint64_t dim[3];
    dim[0] = sz[0];
    dim[1] = sz[1];
    dim[2] = nint;
    plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);

    float* im_int = (float*)mxGetPr(plhs[0]);

    const uint64_t xa = (uint64_t)ceil(xstep);
    double* s_mat = (double*)malloc((nint)*sizeof(double));
    double* t_mat = (double*)malloc((nint)*sizeof(double));
    double* sw_mat = (double*)malloc((nint)*sizeof(double));
    double* tw_mat = (double*)malloc((nint)*sizeof(double));
    double* zw_mat = (double*)malloc((nint)*sizeof(double));
    // z start and end
    uint64_t* zs_ind_mat = (uint64_t*)malloc((sz[2])*sizeof(uint64_t));
    uint64_t* zt_ind_mat = (uint64_t*)malloc((sz[2])*sizeof(uint64_t));
    
    uint64_t zind = 0;
    zs_ind_mat[0] = 0;
    zs_ind_mat[sz[2] - 1] = nint;
    zt_ind_mat[sz[2] - 1] = nint;
    for(uint64_t i = 0; i < nint; i++){
        double s;
        double t;
        double zout;
        double zw;
        // zout = ((double)(sz[2] - 1)) / ((double)(nint - 1)) * ((double)i);
        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        if(Reverse){
            s = ceil(xstep) - xstep * zw;
            t = xstep * (1 - zw);
        }
        else{
            s = xstep * (1 - zw);
            t = ceil(xstep) - xstep * zw;
        }

        //distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        s_mat[i] = s;
        t_mat[i] = t;
        sw_mat[i] = sw;
        tw_mat[i] = tw;
        zw_mat[i] = zw;
        
        // get start and end indices mapping
        if(floor(zout) > zind){
            zs_ind_mat[zind + 1] = i;
            zt_ind_mat[zind + 1] = i;
            zt_ind_mat[zind] = i - 1;
            zind += 1;
        }
        else{
            zt_ind_mat[zind] = i;
        }
        // printf("%llu %f %llu %f %f %f %f %f\n", i, zout, zind, s, t, sw, tw, zw);
    }
    /*
    for(uint64_t z = 0; z < sz[2]; z++){
        printf("%llu %llu %llu\n", z, zs_ind_mat[z], zt_ind_mat[z]);
    }
    */

    #pragma omp parallel for
    for(uint64_t z = 0; z < sz[2]; z++){
        if(z == sz[2]-1){
            if(zs_ind_mat[z] >= nint) 
                continue;
            uint64_t zIndex = nint-1;
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

        
        for(uint64_t ind = zs_ind_mat[z]; ind <= zt_ind_mat[z]; ind++){
            // printf("%llu %llu\n", z, ind);
            uint64_t zint = ind;

            uint64_t s = (uint64_t)s_mat[zint];
            uint64_t t = (uint64_t)t_mat[zint];
            float sw = sw_mat[zint];
            float tw = tw_mat[zint];
            float zw = zw_mat[zint];

            if(!sw){
                for(uint64_t i = 0; i < sz[1]; i++){
                    for(uint64_t j = 0; j < sz[0]; j++){
                        im_int[j+(i*sz[0])+(zint*sz[1]*sz[0])] = im[j+(i*sz[0])+(z*sz[1]*sz[0])];
                    }
                }
                continue;
            }
            for(uint64_t i = 0; i < sz[1]; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    im_int[j+(i*sz[0])+(zint*sz[1]*sz[0])] =
                            ((im_s[j+((i+s)*sz[0])] * (1.0-sw) + im_s[j+((i+s+1)*sz[0])] * sw) * (1.0 - zw))
                            + ((im_t[j+((i+t)*sz[0])] * (1.0-tw) + im_t[j+((i+t+1)*sz[0])] * tw) * zw);
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
    free(zw_mat);
    free(zs_ind_mat);
    free(zt_ind_mat);
}
