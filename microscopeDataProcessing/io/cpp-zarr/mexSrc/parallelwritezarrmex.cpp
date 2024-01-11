#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <omp.h>
#include <stddef.h>
#ifdef _WIN32
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include <fstream>
#include <algorithm>
#include "blosc.h"
#include "mex.h"
#include "../src/parallelwritezarr.h"
#include "../src/helperfunctions.h"
#include "../src/zarr.h"
#include "../src/parallelreadzarr.h"
#include "zlib.h"

//compile
//mex -v COPTIMFLAGS="-DNDEBUG -O3" CFLAGS='$CFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 '-L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid parallelWriteZarr.c helperFunctions.c parallelReadZarr.c

//With zlib
//mex -v COPTIMFLAGS="-DNDEBUG -O3" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O3 -DNDEBUG" CFLAGS='$CFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 '-L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid -lz parallelWriteZarr.c helperFunctions.c parallelReadZarr.c

//mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 zarrMex.c
//
//Windows
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-IC:\Program Files (x86)\bloscZarr\include' '-LC:\Program Files (x86)\bloscZarr\lib' -lblosc '-IC:\Program Files (x86)\cJSON\include\' '-LC:\Program Files (x86)\cJSON\lib' -lcjson '-IC:\Program Files (x86)\blosc\include' '-LC:\Program Files (x86)\blosc\lib' -lblosc2 parallelWriteZarr.c parallelReadZarr.c helperFunctions.c

// TODO: FIX MEMORY LEAKS
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs < 2) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 2 arguments");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string");
    if(mxIsEmpty(prhs[1])) mexErrMsgIdAndTxt("zarr:inputError","All input data axes must be of at least size 1");

    std::vector<uint64_t> startCoords = {0,0,0};
    std::vector<uint64_t> endCoords = {0,0,0};
    bool crop = false;
    bool useUuid = true;
    uint8_t bboxIndex = 0;
    bool sparse = true;

    // Dims are 1 by default
    uint64_t iDims[3] = {1,1,1};

    std::string folderName(mxArrayToString(prhs[0]));
    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    folderName = expandTilde(folderName.c_str());
    #endif
    
    // Check if metadata exists that we can use or if we have to create new metadata
    zarr Zarr;
    if(!fileExists(folderName+"/.zarray")){
        Zarr.set_fileName(folderName);
    }
    else{
        try{
            Zarr = zarr(folderName);
        }
        catch(const std::string &e){
            if(e.find("metadataFileMissing") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
            }
            else if(e == "metadataIncomplete"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray file");
            }
            else mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred\n");
        }
    }

    for(int i = 2; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n",i+1);
        if(!mxIsChar(prhs[i])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        std::string currInput = mxArrayToString(prhs[i]);

        if(currInput == "uuid"){
            useUuid = (bool)*((mxGetPr(prhs[i+1])));
        }
        else if(currInput == "bbox"){
            // Skip bbox if it is empty
            if(!mxGetN(prhs[i+1])) continue;
            else if(mxGetN(prhs[i+1]) == 6){
                crop = true;
                startCoords[0] = (uint64_t)*(mxGetPr(prhs[i+1]))-1;
                startCoords[1] = (uint64_t)*((mxGetPr(prhs[i+1])+1))-1;
                startCoords[2] = (uint64_t)*((mxGetPr(prhs[i+1])+2))-1;
                endCoords[0] = (uint64_t)*((mxGetPr(prhs[i+1])+3));
                endCoords[1] = (uint64_t)*((mxGetPr(prhs[i+1])+4));
                endCoords[2] = (uint64_t)*((mxGetPr(prhs[i+1])+5));
    
    
                uint64_t* iDimsT = (uint64_t*)mxGetDimensions(prhs[1]);
                uint64_t niDims = (uint64_t) mxGetNumberOfDimensions(prhs[1]);
                for(uint64_t i = 0; i < niDims; i++) iDims[i] = iDimsT[i];
    
                if(startCoords[0]+1 < 1 ||
                   startCoords[1]+1 < 1 ||
                   startCoords[2]+1 < 1) mexErrMsgIdAndTxt("zarr:inputError","Lower bounds must be at least 1");
    
                if(endCoords[0]-startCoords[0] > iDims[0] ||
                   endCoords[1]-startCoords[1] > iDims[1] ||
                   endCoords[2]-startCoords[2] > iDims[2]) mexErrMsgIdAndTxt("zarr:inputError","Bounds are invalid for the input data size");
            }
            else if(mxGetN(prhs[i+1]) == 3) bboxIndex = i+1; 
            else mexErrMsgIdAndTxt("zarr:inputError","Input range is not 6 or 3");
        }
        else if(currInput == "cname"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            Zarr.set_cname(mxArrayToString(prhs[i+1]));
        }
        else if(currInput == "subfolders"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","subfolders must be an array of 3 numbers\n");
            Zarr.set_subfolders({(uint64_t)*(mxGetPr(prhs[i+1])),
                               (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                               (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "chunks"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunks must be an array of 3 numbers\n");
            Zarr.set_chunks({(uint64_t)*(mxGetPr(prhs[i+1])),
                             (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                             (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "chunk_shape"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunk_shape must be an array of 3 numbers\n");
            Zarr.set_shard(true);
            Zarr.set_chunk_shape({(uint64_t)*(mxGetPr(prhs[i+1])),
                               (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                               (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "sparse"){
            sparse = (bool)*((mxGetPr(prhs[i+1])));
        }
        else if(currInput == "dimension_separator"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a string\n");
            const std::string dimension_separator(mxArrayToString(prhs[i+1]));
            if(dimension_separator != "." && dimension_separator != "/") mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a . or /\n");
            Zarr.set_dimension_separator(dimension_separator);
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
            Currently Supported Names: uuid, bbox, cname, subfolders, chunks, chunk_shape, sparse\n",currInput.c_str());
        }
    }

    void* zarrC = NULL;

    mxClassID mDType = mxGetClassID(prhs[1]);
    if(mDType == mxUINT8_CLASS){
        Zarr.set_dtype("<u1");
    }
    else if(mDType == mxUINT16_CLASS){
        Zarr.set_dtype("<u2");
    }
    else if(mDType == mxSINGLE_CLASS){
        Zarr.set_dtype("<f4");
    }
    else if(mDType == mxDOUBLE_CLASS){
        Zarr.set_dtype("<f8");
    }

    if(!crop){
        uint64_t nDims = (uint64_t)mxGetNumberOfDimensions(prhs[1]);
        if(nDims < 2 || nDims > 3) mexErrMsgIdAndTxt("zarr:inputError","Input data must be 2D or 3D");

        uint64_t* dims = (uint64_t*)mxGetDimensions(prhs[1]);
        if(nDims == 3) Zarr.set_shape({dims[0],dims[1],dims[2]});
        else Zarr.set_shape({dims[0],dims[1],1});
        if(bboxIndex){
            Zarr.set_chunks({(uint64_t)*(mxGetPr(prhs[bboxIndex])),
                            (uint64_t)*((mxGetPr(prhs[bboxIndex])+1)),
                            (uint64_t)*((mxGetPr(prhs[bboxIndex])+2))});
        }
        try{
            Zarr.write_zarray();
        }
        catch(const std::string &e){
            if(e == "unsupportedCompressor"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",Zarr.get_cname().c_str());
            }
            else if(e.find("cannotOpenZarray") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
            }
            else mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred\n");
        }
    }
    else{
        Zarr.set_shape({endCoords[0],endCoords[1],endCoords[2]});

        if(fileExists(folderName+"/.zarray")){
            if(endCoords[0]-startCoords[0] != iDims[0] ||
               endCoords[1]-startCoords[1] != iDims[1] ||
               endCoords[2]-startCoords[2] != iDims[2]) mexErrMsgIdAndTxt("zarr:inputError","Bounding box size does not match the size of the input data");
        }
        else {
            try{
                Zarr.write_zarray();
            }
            catch(const std::string &e){
                if(e == "unsupportedCompressor"){
                    mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",Zarr.get_cname().c_str());
                }
                else if(e.find("cannotOpenZarray") != std::string::npos){
                    mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
                }
                else mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred\n");
            }
        }

        const std::string dtypeT(Zarr.get_dtype());
        try{
            Zarr = zarr(folderName);
        }
        catch(const std::string &e){
            if(e.find("metadataFileMissing") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
            }
            else if(e == "metadataIncomplete"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray file");
            }
            else mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred\n");
        }
        
        if(dtypeT != Zarr.get_dtype()){
            uint64_t size = (endCoords[0]-startCoords[0])*
                (endCoords[1]-startCoords[1])*
                (endCoords[2]-startCoords[2]);

            uint64_t bitsT = 0;
            if(dtypeT == "<u1") bitsT = 8;
            else if(dtypeT == "<u2") bitsT = 16;
            else if(dtypeT == "<f4") bitsT = 32;
            else if(dtypeT == "<f8") bitsT = 64;
            else mexErrMsgIdAndTxt("tiff:dataTypeError","Cannont convert to passed in data type. Data type not suppported");


            if(Zarr.get_dtype() == "<u1"){
                zarrC = malloc(size*sizeof(uint8_t));
                if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
                else if(bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
                else if(bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
            }
            else if(Zarr.get_dtype() == "<u2"){
                zarrC = malloc(size*sizeof(uint16_t));
                if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
                else if (bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
                else if (bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
            }
            else if(Zarr.get_dtype() == "<f4"){
                zarrC = malloc(size*sizeof(float));
                if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
                else if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
                else if(bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
            }
            else if(Zarr.get_dtype() == "<f8"){
                zarrC = malloc(size*sizeof(double));
                if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
                else if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
                else if(bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("zarr:dataTypeError","Cannont convert to passed in data type. Data type not suppported");
            }
        }
    }


    //endCoords[0]-startCoords[0]
    if(endCoords[0] > Zarr.get_shape(0) ||
       endCoords[1] > Zarr.get_shape(1) ||
       endCoords[2] > Zarr.get_shape(2)) mexErrMsgIdAndTxt("zarr:inputError","Upper bound is invalid");
    if(!crop){
        endCoords = {Zarr.get_shape(0),Zarr.get_shape(1),Zarr.get_shape(2)};
        startCoords = {0,0,0};
    }
    const std::vector<uint64_t> writeShape({endCoords[0]-startCoords[0],
                                      endCoords[1]-startCoords[1],
                                      endCoords[2]-startCoords[2]});
    //printf("%s startCoords[0]yz: %d %d %d endCoords[0]yz: %d %d %d chunkxyz: %d %d %d writeShape[0]yz: %d %d %d\n",Zarr.get_fileName().c_str(),startCoords[0],startCoords[1],startCoords[2],endCoords[0],endCoords[1],endCoords[2],Zarr.get_chunks(0),Zarr.get_chunks(1),Zarr.get_chunks(2),writeShape[0],writeShape[1],writeShape[2]);

    Zarr.set_chunkInfo(startCoords, endCoords);
    bool err = 0;
    if(Zarr.get_dtype() == "<u1"){
        uint64_t bits = 8;
        uint8_t* zarrArr;
        if(zarrC) zarrArr = (uint8_t*)zarrC;
        else zarrArr =  (uint8_t*)mxGetPr(prhs[1]);
        err = parallelWriteZarr(Zarr, (void*)zarrArr, startCoords, endCoords, writeShape, bits, useUuid, crop, sparse);
    }
    else if(Zarr.get_dtype() == "<u2"){
        uint64_t bits = 16;
        uint16_t* zarrArr;
        if(zarrC) zarrArr = (uint16_t*)zarrC;
        else zarrArr = (uint16_t*)mxGetPr(prhs[1]);
        err = parallelWriteZarr(Zarr, (void*)zarrArr, startCoords, endCoords, writeShape, bits, useUuid, crop, sparse);
    }
    else if(Zarr.get_dtype() == "<f4"){
        uint64_t bits = 32;
        float* zarrArr;
        if(zarrC) zarrArr = (float*)zarrC;
        else zarrArr = (float*)mxGetPr(prhs[1]);
        err = parallelWriteZarr(Zarr, (void*)zarrArr, startCoords, endCoords, writeShape, bits, useUuid, crop, sparse);
    }
    else if(Zarr.get_dtype() == "<f8"){
        uint64_t bits = 64;
        double* zarrArr;
        if(zarrC) zarrArr = (double*)zarrC;
        else zarrArr = (double*)mxGetPr(prhs[1]);
        err = parallelWriteZarr(Zarr, (void*)zarrArr, startCoords, endCoords, writeShape, bits, useUuid, crop, sparse);
    }
    else{
        free(zarrC);
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }

    // zarrC is either a copy for data conversion or NULL
    free(zarrC);

    if(err) mexErrMsgIdAndTxt("zarr:writeError",Zarr.get_errString().c_str());
    
}