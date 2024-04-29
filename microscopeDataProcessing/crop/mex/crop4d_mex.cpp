#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop4d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop4d_mex.c

template <typename T>
void crop4d_mex(const T* const &orig, T* &region, const uint64_t* const &dims, const uint64_t* const &bbox, const uint64_t &bytes){
    const uint64_t sizeX = (bbox[4] - bbox[0]) * bytes;
    const uint64_t regionShapeX = bbox[4]-bbox[0];
    const uint64_t regionShapeXY = regionShapeX*(bbox[5]-bbox[1]);
    const uint64_t regionShapeXYZ = regionShapeXY*(bbox[6]-bbox[2]);    
    const uint64_t origShapeX = dims[0];
    const uint64_t origShapeXY = origShapeX*dims[1];
    const uint64_t origShapeXYZ = origShapeXY*dims[2];
    const uint16_t nthread=omp_get_max_threads();
    
    // check memory continuity
    const bool xyContinuous = bbox[0]==0 && bbox[4]==dims[0];
    const bool xyzContinuous = xyContinuous && bbox[1]==0 && bbox[5]==dims[1];
    const bool xyztContinuous = xyzContinuous && bbox[2]==0 && bbox[6]==dims[2];
    
    // printf("%llu %llu %llu %llu %llu %llu %llu\n", sizeX, regionShapeX, regionShapeXY, regionShapeXYZ, origShapeX, origShapeXY, origShapeXYZ);
    // printf("%llu %llu %llu %d\n", xyContinuous, xyzContinuous, xyztContinuous, nthread);
    
    if(nthread == 1){
        if (xyztContinuous){
            memcpy(region, orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY+bbox[3]*origShapeXYZ, regionShapeXYZ*(bbox[7] - bbox[3])*bytes);
            return;
        }

        if (xyzContinuous){
            T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY;
            T* region_ptr = region;
            for(uint64_t t = bbox[3]; t < bbox[7]; t++){
                memcpy(region_ptr+(t-bbox[3])*regionShapeXYZ, orig_ptr+t*origShapeXYZ, regionShapeXYZ*bytes);
            }
            return;
        }

        if (xyContinuous){
            T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX;
            T* region_ptr = region;
            for(uint64_t t = bbox[3]; t < bbox[7]; t++){
                for(uint64_t z = bbox[2]; z < bbox[6]; z++){
                    memcpy(region_ptr+((z-bbox[2])*regionShapeXY+(t-bbox[3])*regionShapeXYZ), orig_ptr+(z*origShapeXY+t*origShapeXYZ), regionShapeXY*bytes);
                }
            }
            return;
        }

        // only continuous on x
        T* orig_ptr = (T*)orig+bbox[0];
        T* region_ptr = region;
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
	        for(uint64_t z = bbox[2]; z < bbox[6]; z++){
		        for(uint64_t y = bbox[1]; y < bbox[5]; y++){
		            memcpy(region_ptr+((y-bbox[1])*regionShapeX+(z-bbox[2])*regionShapeXY+(t-bbox[3])*regionShapeXYZ), orig_ptr+(y*origShapeX+z*origShapeXY+t*origShapeXYZ), sizeX);
		        }
	        }
        }
        return;
    }

    // multithread conditions
    if (xyzContinuous && bbox[7] - bbox[3] >= nthread){
        T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY;
        T* region_ptr = region;
        #pragma omp parallel for
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
            memcpy(region_ptr+(t-bbox[3])*regionShapeXYZ, orig_ptr+t*origShapeXYZ, regionShapeXYZ*bytes);
        }
        return;
    }
    
    if (xyContinuous && (bbox[7] - bbox[3]) * (bbox[6] - bbox[2]) >= nthread){
        T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX;
        T* region_ptr = region;
        #pragma omp parallel for collapse(2) 
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
            for(uint64_t z = bbox[2]; z < bbox[6]; z++){
                memcpy(region_ptr+((z-bbox[2])*regionShapeXY+(t-bbox[3])*regionShapeXYZ), orig_ptr+(z*origShapeXY+t*origShapeXYZ), regionShapeXY*bytes);
            }
        }
        return;
    }

    T* orig_ptr = (T*)orig+bbox[0];
    T* region_ptr = region;
    #pragma omp parallel for collapse(3) 
    for(uint64_t t = bbox[3]; t < bbox[7]; t++){
        for(uint64_t z = bbox[2]; z < bbox[6]; z++){
	        for(uint64_t y = bbox[1]; y < bbox[5]; y++){
	            memcpy(region_ptr+((y-bbox[1])*regionShapeX+(z-bbox[2])*regionShapeXY+(t-bbox[3])*regionShapeXYZ), orig_ptr+(y*origShapeX+z*origShapeXY+t*origShapeXYZ), sizeX);
	        }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("crop:inputError","Number of input arguments must be 2");
    
    if (nlhs != 1) mexErrMsgIdAndTxt("crop:outputError","Number of output arguments must be 1");

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dims[4] = {1, 1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++) dims[t] = dimsA[t];

    if(mxGetN(prhs[1]) != 8) mexErrMsgIdAndTxt("crop:inputError","Input range for bbox is not 8");

    uint64_t bbox[8];
    for(uint64_t t=0; t<8; t++){
        if(t < 4) bbox[t] = (uint64_t)*(mxGetPr(prhs[1])+t)-1;
        else  bbox[t] = (uint64_t)*(mxGetPr(prhs[1])+t);
    }

    // check if bbox is out of bound
    for(uint64_t t=0; t<4; t++){
        if(bbox[t+4] <= bbox[t]) mexErrMsgIdAndTxt("crop:inputError","Upper bound for bbox is smaller than lower bound");
        if(bbox[t+4] > dims[t]) mexErrMsgIdAndTxt("crop:inputError","Upper bound for bbox is invalid");
    }
    
    // printf("%llu %llu %llu %llu %llu %llu %llu %llu\n", bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5], bbox[6], bbox[7]);
    uint64_t rdims[4] = {1, 1, 1, 1};
    for(uint64_t t=0; t<4; t++){
        rdims[t] = bbox[t+4] - bbox[t];
    }    

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxUINT8_CLASS){
        uint64_t bytes = 1;
        uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateNumericArray(4,(mwSize*)rdims,mxUINT8_CLASS, mxREAL);        
        uint8_t* region = (uint8_t*)mxGetPr(plhs[0]);
        crop4d_mex<uint8_t>(orig, region, dims, bbox, bytes);
    }
    else if(mDType == mxUINT16_CLASS){
        uint64_t bytes = 2;
        uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateNumericArray(4,(mwSize*)rdims,mxUINT16_CLASS, mxREAL);        
        uint16_t* region = (uint16_t*)mxGetPr(plhs[0]);
        crop4d_mex<uint16_t>(orig, region, dims, bbox, bytes);
    }
    else if(mDType == mxSINGLE_CLASS){
        uint64_t bytes = 4;
        float* orig = (float*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateNumericArray(4,(mwSize*)rdims,mxSINGLE_CLASS, mxREAL);        
        float* region = (float*)mxGetPr(plhs[0]);
        crop4d_mex<float>(orig, region, dims, bbox, bytes);
    }
    else if(mDType == mxDOUBLE_CLASS){
        uint64_t bytes = 8;
        double* orig = (double*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateNumericArray(4,(mwSize*)rdims,mxDOUBLE_CLASS, mxREAL);        
        double* region = (double*)mxGetPr(plhs[0]);
        crop4d_mex<double>(orig, region, dims, bbox, bytes);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
}
