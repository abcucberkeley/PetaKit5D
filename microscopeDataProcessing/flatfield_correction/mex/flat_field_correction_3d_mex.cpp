#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <type_traits>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' flat_field_correction_3d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c

// Template function to check if the type is int, uint16_t, or uint8_t
template<typename T>
constexpr bool is_integral_type() {
    return std::is_same<T, int>::value || std::is_same<T, uint16_t>::value || std::is_same<T, uint8_t>::value;
}

// Template function to round the output if it's an integral type
template<typename T>
T round_output(float value) {
    if constexpr (is_integral_type<T>()) {
        return static_cast<T>(value + 0.5); // Round to nearest integer
    } else {
        return static_cast<T>(value); // No rounding for other types
    }
}

template <typename T, typename U>
void flat_field_correction_3d_mex(const T* const &imat, const float* const &fmat, const float* const &bmat, U* &nmat, const bool &constOffset, const float &constOffsetValue, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ) {
    const uint64_t shapeXY = shapeX * shapeY;

    // get the number of thread
    const uint16_t nthread=omp_get_max_threads();

    if (nthread > 1){
        #pragma omp parallel for collapse(3)
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    float f = 0;
                    uint64_t ind_xy = x + y * shapeX;
                    uint64_t ind_xyz = ind_xy + z * shapeXY;
                    
                    // Calculate fmat and dmat pointers once
                    const T* i_ptr = imat + ind_xyz;
                    const float* f_ptr = fmat + ind_xy;
                    const float* b_ptr = bmat + ind_xy;
                    U* n_ptr = nmat + ind_xyz;
                    
                    f = (static_cast<float>(*i_ptr) - (*b_ptr)) / (*f_ptr);
                    f = f >= 0 ? f : 0;
                    if(constOffset)
                        *n_ptr = round_output<U>(f + constOffsetValue);
                    else
                        *n_ptr = round_output<U>(f + (*b_ptr));
                }
            }
        }
    }
    else{
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    float f = 0;
                    uint64_t ind_xy = x + y * shapeX;
                    uint64_t ind_xyz = ind_xy + z * shapeXY;
                    
                    // Calculate fmat and dmat pointers once
                    const T* i_ptr = imat + ind_xyz;
                    const float* f_ptr = fmat + ind_xy;
                    const float* b_ptr = bmat + ind_xy;
                    U* n_ptr = nmat + ind_xyz;
                    
                    f = (static_cast<float>(*i_ptr) - (*b_ptr)) / (*f_ptr);
                    f = f >= 0 ? f : 0;
                    if(constOffset)
                        *n_ptr = round_output<U>(f + constOffsetValue);
                    else
                        *n_ptr = round_output<U>(f + (*b_ptr));
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4 && nrhs != 5) mexErrMsgIdAndTxt("flat_field_correction_3d_mex:inputError","Number of input arguments must be 3 or 4");
    if (nlhs != 1) mexErrMsgIdAndTxt("flat_field_correction_3d_mex:outputError","Number of output arguments must be 1");

    uint64_t* dimsI = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t* dimsF = (uint64_t*)mxGetDimensions(prhs[1]);
    uint64_t* dimsB = (uint64_t*)mxGetDimensions(prhs[2]);
    uint64_t dims[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        dims[t] = dimsI[t];
    }

    uint64_t ndimF = (uint64_t) mxGetNumberOfDimensions(prhs[1]);
    if (ndimF != 1 && ndimF != 2) mexErrMsgIdAndTxt("flat_field_correction_3d_mex:inputError","The flat field image must be 1d or 2d");    
    for(uint64_t t=0; t < ndimF; t++){
        if(dimsI[t] != dimsF[t]){
            mexErrMsgIdAndTxt("flat_field_correction_3d_mex:inputError","The size of flat field image does not matach that of the data");
        }
    }

    uint64_t ndimB = (uint64_t) mxGetNumberOfDimensions(prhs[2]);
    if (ndimB != 1 && ndimB != 2) mexErrMsgIdAndTxt("flat field correction:inputError","The background image must be 1d or 2d");
    for(uint64_t t=0; t < ndimB; t++){
        if(dimsI[t] != dimsB[t]){
            mexErrMsgIdAndTxt("flat_field_correction_3d_mex:inputError","The size of background image does not matach that of the data");
        }
    }

    uint64_t shapeX = dims[0];
    uint64_t shapeY = dims[1];
    uint64_t shapeZ = dims[2];
    
    bool castDataType = (bool) mxGetScalar(prhs[3]);

    bool constOffset = false;
    float constOffsetValue = 0;
    if (nrhs == 5) {
        constOffset = true;
        constOffsetValue = (float) mxGetScalar(prhs[4]);
    }

    mxClassID mDType = mxGetClassID(prhs[0]);
    mxClassID mDType_f = mxGetClassID(prhs[1]);
    mxClassID mDType_b = mxGetClassID(prhs[2]);
    if(mDType_f != mxSINGLE_CLASS || mDType_b != mxSINGLE_CLASS) mexErrMsgIdAndTxt("flat_field_correction_3d_mex:inputError","The data type of the flat field and background images must be single!");
    if(mDType == mxUINT16_CLASS){
        uint16_t* imat = (uint16_t*)mxGetPr(prhs[0]);
        float* fmat = (float*)mxGetPr(prhs[1]);
        float* bmat = (float*)mxGetPr(prhs[2]);
        
        if(castDataType){
            plhs[0] = mxCreateNumericArray(3, (mwSize*)dims, mxUINT16_CLASS, mxREAL);
            uint16_t* nmat = (uint16_t*)mxGetPr(plhs[0]);

            flat_field_correction_3d_mex(imat, fmat, bmat, nmat, constOffset, constOffsetValue, shapeX, shapeY, shapeZ);
        } else {
            plhs[0] = mxCreateNumericArray(3, (mwSize*)dims, mxSINGLE_CLASS, mxREAL);
            float* nmat = (float*)mxGetPr(plhs[0]);

            flat_field_correction_3d_mex(imat, fmat, bmat, nmat, constOffset, constOffsetValue, shapeX, shapeY, shapeZ);
        }
    }    
    else if(mDType == mxSINGLE_CLASS){
        float* imat = (float*)mxGetPr(prhs[0]);
        float* fmat = (float*)mxGetPr(prhs[1]);
        float* bmat = (float*)mxGetPr(prhs[2]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)dims, mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);
        
        flat_field_correction_3d_mex(imat, fmat, bmat, nmat, constOffset, constOffsetValue, shapeX, shapeY, shapeZ);
    }
    else{
        mexErrMsgIdAndTxt("flat_field_correction_3d_mex:dataTypeError","Data type not suppported");
    }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
