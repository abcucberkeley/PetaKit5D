#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include <string.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c


void skewed_space_interp_defined_stepsize(const float* im, float* im_int, const float xstep, const float stepsize, const uint8_t Reverse,
                                          const uint64_t sz[3], const uint64_t outDim[3]){

    const uint64_t nint = outDim[2];
    const uint64_t nxy = sz[0] * sz[1];

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
    
    #pragma omp parallel for
    for(uint64_t z = 0; z < sz[2]; z++){
        if(z == sz[2]-1){
            if(zs_ind_mat[z] >= nint) continue;

            uint64_t zint = nint-1;
            memcpy(im_int+zint*nxy, im+z*nxy, nxy*sizeof(float));
            continue;
        }

        float* im_s = (float*)malloc(sz[0]*(sz[1]+xa)*sizeof(float));
        float* im_t = (float*)malloc(sz[0]*(sz[1]+xa)*sizeof(float));

        if(Reverse){
            memset(im_s, 0, xa*sz[0]*sizeof(float));
            memcpy(im_s+xa*sz[0], im+z*nxy, nxy*sizeof(float));

            memcpy(im_t, im+(z+1)*nxy, nxy*sizeof(float));            
            memset(im_t+nxy, 0, xa*sz[0]*sizeof(float));
        }
        else{
            memcpy(im_s, im+z*nxy, nxy*sizeof(float));
            memset(im_s+nxy, 0, xa*sz[0]*sizeof(float));

            memset(im_t, 0, xa*sz[0]*sizeof(float));            
            memcpy(im_t+xa*sz[0], im+(z+1)*nxy, nxy*sizeof(float));
        }
        
        for(uint64_t ind = zs_ind_mat[z]; ind <= zt_ind_mat[z]; ind++){
            // printf("%llu %llu\n", z, ind);
            uint64_t zint = ind;

            uint64_t s = (uint64_t)s_mat[zint];
            uint64_t t = (uint64_t)t_mat[zint];
            float sw = sw_mat[zint];
            float tw = tw_mat[zint];
            float zw = zw_mat[zint];

            if((Reverse && !sw) || (!Reverse && !tw)){
                memcpy(im_int+zint*nxy, im+z*nxy, nxy*sizeof(float));
                continue;
            }
    
            // printf("%llu %llu %f %f %f\n", s, t, sw, tw, zw);
            #pragma omp parallel for 
            for(uint64_t i = 0; i < sz[1]; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    im_int[j+(i*sz[0])+(zint*nxy)] =
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


void skewed_space_interp_defined_stepsize_16bit(const uint16_t* im, uint16_t* im_int, const float xstep, const float stepsize, const uint8_t Reverse,
                                          const uint64_t sz[3], const uint64_t outDim[3]){

    const uint64_t nint = outDim[2];
    const uint64_t nxy = sz[0] * sz[1];

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

    #pragma omp parallel for
    for(uint64_t z = 0; z < sz[2]; z++){
        if(z == sz[2]-1){
            if(zs_ind_mat[z] >= nint) continue;

            uint64_t zint = nint-1;
            memcpy(im_int+zint*nxy, im+z*nxy, nxy*sizeof(uint16_t));
            continue;
        }

        uint16_t* im_s = (uint16_t*)malloc(sz[0]*(sz[1]+xa)*sizeof(uint16_t));
        uint16_t* im_t = (uint16_t*)malloc(sz[0]*(sz[1]+xa)*sizeof(uint16_t));

        if(Reverse){
            memset(im_s, 0, xa*sz[0]*sizeof(uint16_t));
            memcpy(im_s+xa*sz[0], im+z*nxy, nxy*sizeof(uint16_t));

            memcpy(im_t, im+(z+1)*nxy, nxy*sizeof(uint16_t));
            memset(im_t+nxy, 0, xa*sz[0]*sizeof(uint16_t));
        }
        else{
            memcpy(im_s, im+z*nxy, nxy*sizeof(uint16_t));
            memset(im_s+nxy, 0, xa*sz[0]*sizeof(uint16_t));

            memset(im_t, 0, xa*sz[0]*sizeof(uint16_t));
            memcpy(im_t+xa*sz[0], im+(z+1)*nxy, nxy*sizeof(uint16_t));
        }

        for(uint64_t ind = zs_ind_mat[z]; ind <= zt_ind_mat[z]; ind++){
            // printf("%llu %llu\n", z, ind);
            uint64_t zint = ind;

            uint64_t s = (uint64_t)s_mat[zint];
            uint64_t t = (uint64_t)t_mat[zint];
            float sw = sw_mat[zint];
            float tw = tw_mat[zint];
            float zw = zw_mat[zint];

            if((Reverse && !sw) || (!Reverse && !tw)){
                memcpy(im_int+zint*nxy, im+z*nxy, nxy*sizeof(uint16_t));
                continue;
            }

            // printf("%llu %llu %f %f %f\n", s, t, sw, tw, zw);
            #pragma omp parallel for
            for(uint64_t i = 0; i < sz[1]; i++){
                for(uint64_t j = 0; j < sz[0]; j++){
                    im_int[j+(i*sz[0])+(zint*nxy)] = (uint16_t)
                        (((im_s[j+((i+s)*sz[0])] * (1.0-sw) + im_s[j+((i+s+1)*sz[0])] * sw) * (1.0 - zw))
                        + ((im_t[j+((i+t)*sz[0])] * (1.0-tw) + im_t[j+((i+t+1)*sz[0])] * tw) * zw) + 0.5);
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


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

    if (nrhs != 4) mexErrMsgIdAndTxt("volume_warp:inputError","Number of input arguments must be 3");
    if (nlhs != 1) mexErrMsgIdAndTxt("volume_warp:outputError","Number of output arguments must be 1");

    const float xstep = (float)*mxGetPr(prhs[1]);
    const float stepsize = (float)*mxGetPr(prhs[2]);
    const uint8_t Reverse = (uint8_t)mxIsLogicalScalarTrue(prhs[3]);

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t sz[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        sz[t] = dimsA[t];
    }

    uint64_t nint = (uint64_t) (floor(round((float)(sz[2] - 1) / stepsize * 10000) / 10000) + 1);
    uint64_t outDim[3];
    outDim[0] = sz[0];
    outDim[1] = sz[1];
    outDim[2] = nint;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxSINGLE_CLASS){
        float* im = (float*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxSINGLE_CLASS, mxREAL);
        float* im_int = (float*)mxGetPr(plhs[0]);

        skewed_space_interp_defined_stepsize(im, im_int, xstep, stepsize, Reverse, sz, outDim);
    }
    else if(mDType == mxUINT16_CLASS){
        uint16_t* im = (uint16_t*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxUINT16_CLASS, mxREAL);
        uint16_t* im_int = (uint16_t*)mxGetPr(plhs[0]);

        skewed_space_interp_defined_stepsize_16bit(im, im_int, xstep, stepsize, Reverse, sz, outDim);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}
