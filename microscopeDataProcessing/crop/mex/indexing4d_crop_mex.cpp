#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_crop_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_crop_mex.c

template <typename T>
void indexing4d_mex(T* &orig, const T* const &region, const uint64_t* const &dims, const uint64_t* const &rdims, const uint64_t* const &bbox, const uint64_t* const &rbbox, const uint64_t &bytes){
    const uint64_t sizeX = (bbox[4] - bbox[0]) * bytes;
    const uint64_t regionShapeX = rdims[0];
    const uint64_t regionShapeXY = rdims[0]*rdims[1];
    const uint64_t regionShapeXYZ = rdims[0]*rdims[1]*rdims[2];
    const uint64_t origShapeX = dims[0];
    const uint64_t origShapeXY = dims[0]*dims[1];
    const uint64_t origShapeXYZ = dims[0]*dims[1]*dims[2];
    const uint16_t nthread=omp_get_max_threads();
    
    // check memory continuity
    const bool xyContinuous = bbox[0]==0 && bbox[4]==dims[0] && rbbox[0]==0 && rbbox[4]==rdims[0];
    const bool xyzContinuous = xyContinuous && bbox[1]==0 && bbox[5]==dims[1] && rbbox[1]==0 && rbbox[5]==rdims[1];
    const bool xyztContinuous = xyzContinuous && bbox[2]==0 && bbox[6]==dims[2] && rbbox[2]==0 && rbbox[6]==rdims[2];
    
    // printf("%llu %llu %llu %llu %llu %llu %llu\n", sizeX, regionShapeX, regionShapeXY, regionShapeXYZ, origShapeX, origShapeXY, origShapeXYZ);
    // printf("%llu %llu %llu %d\n", xyContinuous, xyzContinuous, xyztContinuous, nthread);
    
    if(nthread == 1){
        if (xyztContinuous){
            memcpy(orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY+bbox[3]*origShapeXYZ, region+rbbox[0]+rbbox[1]*regionShapeX+rbbox[2]*regionShapeXY+rbbox[3]*regionShapeXYZ, origShapeXYZ*(bbox[7] - bbox[3])*bytes);
            return;   
        }

        if (xyzContinuous){
            T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY;
            T* region_ptr = (T*)region+rbbox[0]+rbbox[1]*regionShapeX+rbbox[2]*regionShapeXY;
            for(uint64_t t = bbox[3]; t < bbox[7]; t++){
                memcpy(orig_ptr+t*origShapeXYZ, region_ptr+(t-bbox[3]+rbbox[3])*regionShapeXYZ, origShapeXY*(bbox[6]-bbox[2])*bytes);
            }
            return;
        }

        if (xyContinuous){
            T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX;
            T* region_ptr = (T*)region+rbbox[0]+rbbox[1]*regionShapeX;
            for(uint64_t t = bbox[3]; t < bbox[7]; t++){
                for(uint64_t z = bbox[2]; z < bbox[6]; z++){
                    memcpy(orig_ptr+(z*origShapeXY+t*origShapeXYZ), region_ptr+((z-bbox[2]+rbbox[2])*regionShapeXY+(t-bbox[3]+rbbox[3])*regionShapeXYZ), sizeX*(bbox[5]-bbox[1]));
                }
            }
            return;
        }
        
        // only continuous on x
        T* orig_ptr = (T*)orig+bbox[0];
        T* region_ptr = (T*)region+rbbox[0];
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
	        for(uint64_t z = bbox[2]; z < bbox[6]; z++){
		        for(uint64_t y = bbox[1]; y < bbox[5]; y++){
		            memcpy(orig_ptr+(y*origShapeX+z*origShapeXY+t*origShapeXYZ), region_ptr+((y-bbox[1]+rbbox[1])*regionShapeX+(z-bbox[2]+rbbox[2])*regionShapeXY+(t-bbox[3]+rbbox[3])*regionShapeXYZ), sizeX);
		        }
	        }
        }
        return;
    }
    
    // multithread conditions
    if (xyzContinuous && bbox[7] - bbox[3] >= nthread){
        T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX+bbox[2]*origShapeXY;
        T* region_ptr = (T*)region+rbbox[0]+rbbox[1]*regionShapeX+rbbox[2]*regionShapeXY;        
        #pragma omp parallel for
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
            memcpy(orig_ptr+t*origShapeXYZ, region_ptr+(t-bbox[3]+rbbox[3])*regionShapeXYZ, origShapeXY*(bbox[6]-bbox[2])*bytes);
        }
        return;
    }
    
    if (xyContinuous && (bbox[7] - bbox[3]) * (bbox[6] - bbox[2]) >= nthread){
        T* orig_ptr = (T*)orig+bbox[0]+bbox[1]*origShapeX;
        T* region_ptr = (T*)region+rbbox[0]+rbbox[1]*regionShapeX;        
        #pragma omp parallel for collapse(2) 
        for(uint64_t t = bbox[3]; t < bbox[7]; t++){
            for(uint64_t z = bbox[2]; z < bbox[6]; z++){
                memcpy(orig_ptr+(z*origShapeXY+t*origShapeXYZ), region_ptr+((z-bbox[2]+rbbox[2])*regionShapeXY+(t-bbox[3]+rbbox[3])*regionShapeXYZ), sizeX*(bbox[5]-bbox[1]));
            }
        }
        return;
    }

    T* orig_ptr = (T*)orig+bbox[0];
    T* region_ptr = (T*)region+rbbox[0];
    #pragma omp parallel for collapse(3) 
    for(uint64_t t = bbox[3]; t < bbox[7]; t++){
        for(uint64_t z = bbox[2]; z < bbox[6]; z++){
	        for(uint64_t y = bbox[1]; y < bbox[5]; y++){
	            memcpy(orig_ptr+(y*origShapeX+z*origShapeXY+t*origShapeXYZ), region_ptr+((y-bbox[1]+rbbox[1])*regionShapeX+(z-bbox[2]+rbbox[2])*regionShapeXY+(t-bbox[3]+rbbox[3])*regionShapeXYZ), sizeX);
	        }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 && nrhs !=4) mexErrMsgIdAndTxt("indexing:inputError","Number of input arguments must be 3 or 4");
    
    if (nlhs != 0) mexErrMsgIdAndTxt("indexing:outputError","Number of output arguments must be 0");

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dims[4] = {1, 1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++) dims[t] = dimsA[t];

    if(mxGetN(prhs[1]) != 8) mexErrMsgIdAndTxt("indexing:inputError","Input range for bbox is not 8");

    uint64_t bbox[8];
    for(uint64_t t=0; t<8; t++){
        if(t < 4) bbox[t] = (uint64_t)*(mxGetPr(prhs[1])+t)-1;
        else  bbox[t] = (uint64_t)*(mxGetPr(prhs[1])+t);
    }

    // check if bbox is out of bound
    for(uint64_t t=0; t<4; t++){
        if(bbox[t+4] <= bbox[t]) mexErrMsgIdAndTxt("indexing:inputError","Upper bound for bbox is smaller than lower bound");
        if(bbox[t+4] > dims[t]) mexErrMsgIdAndTxt("indexing:inputError","Upper bound for bbox is invalid");
    }
    
    uint64_t* rdimsT = (uint64_t*)mxGetDimensions(prhs[2]);
    uint64_t rdims[4] = {1, 1, 1, 1};
    uint64_t nrdim = (uint64_t) mxGetNumberOfDimensions(prhs[2]);
    for(uint64_t t=0; t < nrdim; t++) rdims[t] = rdimsT[t];

    uint64_t rbbox[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    if (nrhs==3){
        for(uint64_t t=0; t < 4; t++) rbbox[t+4]= rdims[t];
    }else if (nrhs==4){
        if(mxGetN(prhs[3]) != 8) mexErrMsgIdAndTxt("indexing:inputError","Input range for region bbox is not 8");
        for(uint64_t t=0; t<8; t++){
            if(t < 4) rbbox[t] = (uint64_t)*(mxGetPr(prhs[3])+t)-1;
            else rbbox[t] = (uint64_t)*(mxGetPr(prhs[3])+t);
        }        
    }

    // printf("%llu %llu %llu %llu %llu %llu %llu %llu\n", rbbox[0], rbbox[1], rbbox[2], rbbox[3], rbbox[4], rbbox[5], rbbox[6], rbbox[7]);
    // printf("%llu %llu %llu %llu %llu %llu %llu %llu\n", bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5], bbox[6], bbox[7]);

    // check if region bbox is out of bound and if bbox and region bbox has the same size
    for(uint64_t t=0; t<4; t++){
        if(rbbox[t+4] <= rbbox[t]) mexErrMsgIdAndTxt("indexing:inputError","Upper bound for region bbox is smaller than lower bound");
        if(rbbox[t+4] > rdims[t]) mexErrMsgIdAndTxt("indexing:inputError","Upper bound for region bbox is invalid");
        if(bbox[t+4]-bbox[t] != rbbox[t+4]-rbbox[t]) mexErrMsgIdAndTxt("indexing:inputError","The range of region bbox is different from that of bbox");
    }
           
    mxClassID mDType = mxGetClassID(prhs[0]);
    mxClassID mDType_region = mxGetClassID(prhs[2]);
    if(mDType != mDType_region) mexErrMsgIdAndTxt("indexing:inputError","The data type of the region does not match that of the data!");

    if(mDType == mxLOGICAL_CLASS){
        uint64_t bytes = 1;
        bool* orig = (bool*)mxGetPr(prhs[0]);
        bool* region = (bool*)mxGetPr(prhs[2]);
        indexing4d_mex<bool>(orig, region, dims, rdims, bbox, rbbox, bytes);
    }
    else if(mDType == mxUINT8_CLASS){
        uint64_t bytes = 1;
        uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        uint8_t* region = (uint8_t*)mxGetPr(prhs[2]);
        indexing4d_mex<uint8_t>(orig, region, dims, rdims, bbox, rbbox, bytes);
    }
    else if(mDType == mxUINT16_CLASS){
        uint64_t bytes = 2;
        uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        uint16_t* region = (uint16_t*)mxGetPr(prhs[2]);
        indexing4d_mex<uint16_t>(orig, region, dims, rdims, bbox, rbbox, bytes);
    }
    else if(mDType == mxSINGLE_CLASS){
        uint64_t bytes = 4;
        float* orig = (float*)mxGetPr(prhs[0]);
        float* region = (float*)mxGetPr(prhs[2]);
        indexing4d_mex<float>(orig, region, dims, rdims, bbox, rbbox, bytes);
    }
    else if(mDType == mxDOUBLE_CLASS){
        uint64_t bytes = 8;
        double* orig = (double*)mxGetPr(prhs[0]);
        double* region = (double*)mxGetPr(prhs[2]);
        indexing4d_mex<double>(orig, region, dims, rdims, bbox, rbbox, bytes);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
}
