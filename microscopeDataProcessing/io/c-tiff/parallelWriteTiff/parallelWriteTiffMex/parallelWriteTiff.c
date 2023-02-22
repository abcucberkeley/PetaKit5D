#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "lzwEncode.h"
#include "tiffio.h"
#include "omp.h"
#include "mex.h"
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' '-L/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' -ltiff /clusterfs/fiona/matthewmueller/parallelTiffTesting/main.c
//mex COMPFLAGS='$COMPFLAGS /openmp' '-IC:\Program Files (x86)\tiff\include\' '-LC:\Program Files (x86)\tiff\lib\' -ltiffd.lib C:\Users\Matt\Documents\parallelTiff\main.cpp

//zlib
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' '-I/global/home/groups/consultsw/sl-7.x86_64/modules/zlib/1.2.11/include/' '-L/global/home/groups/consultsw/sl-7.x86_64/modules/zlib/1.2.11/lib' -lz '-L/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' -ltiff parallelWriteTiff.c

//lzw
//mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' '-L/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' -ltiff parallelWriteTiff.c lzw.c

//libtiff 4.4.0
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" LDOPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/clusterfs/fiona/matthewmueller/software/tiff-4.4.0/include' '-L/clusterfs/fiona/matthewmueller/software/tiff-4.4.0/lib' -ltiff parallelWriteTiff.c lzwEncode.c

// Handle the tilde character in filenames on Linux/Mac
#ifndef _WIN32
#include <wordexp.h>
char* expandTilde(char* path) {
    wordexp_t expPath;
    wordexp(path, &expPath, 0);
    return expPath.we_wordv[0];
}
#endif

static void mkdirRecursive(const char *dir) {
    char tmp[8192];
    char *p = NULL;
    size_t len;
    #ifdef _WIN32
    char fileSep = '\\';
    #else
    char fileSep = '/';
    #endif
    int status;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    //printf(tmp);
    if (tmp[len - 1] == fileSep)
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++){
        //printf("p: %s\n",p);
        //printf("*p: %c",*p);
        //printf("fileSep: %c\n",fileSep);
        if (*p == fileSep) {
            *p = 0;

            #ifdef _WIN32
            mkdir(tmp);
            #else
            mkdir(tmp, 0775);
            #endif

            chmod(tmp, 0775);
            *p = fileSep;
        }
    }
    #ifdef _WIN32
    mkdir(tmp);
    #else
    mkdir(tmp, 0775);
    #endif
    chmod(tmp, 0775);
}


void DummyHandler(const char* module, const char* fmt, va_list ap)
{
    // ignore errors and warnings
}

void writeTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, const void* tiffOld, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint64_t stripsPerDir, uint64_t* cSizes, const char* mode){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    int32_t w;
    int safeMode = 0;
    uint64_t extraBytes = 2000;
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        uint64_t len = 0;
        uint8_t* compr = (uint8_t*)malloc((((x*stripSize)*(bits/8))+(extraBytes*(bits/8)))*2+1);
        for(uint64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice || safeMode) break;

            for (uint64_t i = 0; i*stripSize < y; i++)
            {
                uint8_t* cArr = (uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8));
                if (stripSize*(i+1) > y){
                    len = (y-(stripSize*i))*x*(bits/8);
                }
                else{
                    len = stripSize*x*(bits/8);
                }
                if(dir == z+startSlice-1 && len == (y-(stripSize*i))*x*(bits/8)){
                    uint8_t* cArrL = (uint8_t*)malloc(len+(extraBytes*(bits/8)));
                    memcpy(cArrL,cArr,len);
                    cSizes[i+(dir*stripsPerDir)] = lzwEncode(cArrL,compr,len+(extraBytes*(bits/8)));
                    if(cSizes[i+(dir*stripsPerDir)] > len){
                        free(cArrL);
                        #pragma omp critical
                        {
                            safeMode = 1;
                        }
                        break;
                    }
                    memcpy(cArr,compr,cSizes[i+(dir*stripsPerDir)]);
                    free(cArrL);
                    continue;
                }
                cSizes[i+(dir*stripsPerDir)] = lzwEncode(cArr,compr,len+(extraBytes*(bits/8)));
                if(cSizes[i+(dir*stripsPerDir)] > len){
                    #pragma omp critical
                    {
                        safeMode = 1;
                    }
                    break;
                }
                memcpy(cArr,compr,cSizes[i+(dir*stripsPerDir)]);
            }
        }
        free(compr);
    }
    uint8_t** comprA = NULL;
    if(safeMode){
        // Restore the array as it may have changed
        #pragma omp parallel for collapse(3)
        for(uint64_t dir = 0; dir < z; dir++){
            for(uint64_t i = 0; i < x; i++){
                for(uint64_t j = 0; j < y; j++){
                    switch(bits){
                        case 8:
                            ((uint8_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint8_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                            break;
                        case 16:
                            ((uint16_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                            break;
                        case 32:
                            ((float*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((float*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                            break;
                        case 64:
                            ((double*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((double*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                            break;
                    }
                }
            }
        }
        uint64_t totalStrips = stripsPerDir*z;
        comprA = (uint8_t**)malloc(totalStrips*sizeof(uint8_t*));
        #pragma omp parallel for
        for(uint64_t i = 0; i < totalStrips; i++){
            comprA[i] = malloc((((x*stripSize)*(bits/8))+(extraBytes*(bits/8)))*2+1);
        }

        #pragma omp parallel for
        for(w = 0; w < numWorkers; w++){
            uint64_t len = 0;
            uint8_t* compr = (uint8_t*)malloc((((x*stripSize)*(bits/8))+(extraBytes*(bits/8)))*2+1);
            for(uint64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
                if(dir>=z+startSlice) break;

                for (uint64_t i = 0; i*stripSize < y; i++)
                {
                    uint8_t* cArr = comprA[i+(dir*stripsPerDir)];
                    if (stripSize*(i+1) > y){
                        len = (y-(stripSize*i))*x*(bits/8);
                    }
                    else{
                        len = stripSize*x*(bits/8);
                    }

                    if(dir == z+startSlice-1 && len == (y-(stripSize*i))*x*(bits/8)){
                        uint8_t* cArrL = (uint8_t*)malloc(len+(extraBytes*(bits/8)));
                        memcpy(cArrL,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)),len);
                        cSizes[i+(dir*stripsPerDir)] = lzwEncode(cArrL,compr,len+(extraBytes*(bits/8)));
                        memcpy(cArr,compr,cSizes[i+(dir*stripsPerDir)]);
                        free(cArrL);
                        continue;
                    }
                    cSizes[i+(dir*stripsPerDir)] = lzwEncode((uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)),compr,len+(extraBytes*(bits/8)));
                    memcpy(cArr,compr,cSizes[i+(dir*stripsPerDir)]);

                }
            }
            free(compr);
        }
    }

    uint64_t totalStrips = stripsPerDir*z;
    uint64_t cN = totalStrips;
    uint64_t cSize = 0;
    #pragma omp parallel for reduction (+:cSize)
    for (uint64_t i = 0; i < cN; i++){
        cSize=cSize+cSizes[i];
    }

    // Check if data actually compressed well
    uint64_t uncSize = x*y*z*(bits/8);
    double comprFactor = .9;
    uint8_t compress = 1;
    if(uncSize*comprFactor < cSize){
        cSize = uncSize;
        compress = 0;
        // Restore the array as it may have changed
        if(!safeMode){
            #pragma omp parallel for collapse(3)
            for(uint64_t dir = 0; dir < z; dir++){
                for(uint64_t i = 0; i < x; i++){
                    for(uint64_t j = 0; j < y; j++){
                        switch(bits){
                            case 8:
                                ((uint8_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint8_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                                break;
                            case 16:
                                ((uint16_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                                break;
                            case 32:
                                ((float*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((float*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                                break;
                            case 64:
                                ((double*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((double*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                                break;
                        }
                    }
                }
            }
        }
    }

    // Add bytes for some mandatory tags
    uint64_t manBytes = 24;
    cSize += manBytes*z;
    //printf("Compressed Size is: %llu bytes\n",cSize);

    //mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported. %d %d %d %d",x,y,z,bits);
    //Close to UINT32_MAX. Want some extra room for incorrect size calculation for now.
    TIFF* tif = NULL;
    if(!strcmp(mode,"w")){
        tif = TIFFOpen(fileName, (cSize < 3.8e9) ? "w" : "w8");
        if(!tif){
            mexErrMsgIdAndTxt("tiff:threadError","Error: File \"%s\" cannot be opened",fileName);
        }
    }
    else if(!strcmp(mode,"a")){
        tif = TIFFOpen(fileName, "r");
        if(!tif){
            mexErrMsgIdAndTxt("tiff:threadError","Error: File \"%s\" cannot be opened",fileName);
        }
        uint64_t xTemp = 1,yTemp = 1,zTemp = 1;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &xTemp);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &yTemp);
        uint16_t s = 0, m = 0, t = 1;
        while(TIFFSetDirectory(tif,t)){
            s = t;
            t *= 8;
            if(s > t){
                t = 65535;
                printf("Number of slices > 32768");
                break;
            }
        }
        while(s != t){
            m = (s+t+1)/2;
            if(TIFFSetDirectory(tif,m)){
                s = m;
            }
            else{
                if(m > 0) t = m-1;
                else t = m;
            }
        }
        zTemp = s+1;
        cSize += (xTemp*yTemp*zTemp)*(bits/8);
        TIFFClose(tif);
        tif = TIFFOpen(fileName, (cSize < 3.8e9) ? "a" : "a8");
        if(!tif){
            mexErrMsgIdAndTxt("tiff:threadError","Error: File \"%s\" cannot be opened",fileName);
        }
    }
    else{
        printf("Error: mode \"%s\" is not supported. Use w or a for mode type", mode);
        return;
    }

    uint8_t err = 0;
    //char errString[10000];
    uint64_t len = 0;
    for(uint64_t dir = startSlice; dir < z; dir++){
        if(dir>=z+startSlice || err) break;
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, x);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, y);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, stripSize);
        if(compress) TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
        else TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, 1);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

        if(bits >= 32){
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        }

        for (int64_t i = 0; i*stripSize < y; i++)
        {
            if(compress){
                if(!comprA){
                    TIFFWriteRawStrip(tif,i,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), cSizes[i+(dir*stripsPerDir)]);
                }
                else TIFFWriteRawStrip(tif,i,comprA[i+(dir*stripsPerDir)], cSizes[i+(dir*stripsPerDir)]);
            }
            else{
                if (stripSize*(i+1) > y){
                    len = (y-(stripSize*i))*x*(bits/8);
                }
                else{
                    len = stripSize*x*(bits/8);
                }
                TIFFWriteRawStrip(tif,i,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
            }

        }
        TIFFWriteDirectory(tif);
    }
    if(comprA){
        #pragma omp parallel for
        for(uint64_t i = 0; i < totalStrips; i++){
            free(comprA[i]);
        }
        free(comprA);
    }
    TIFFClose(tif);
    //if(err) mexErrMsgIdAndTxt("tiff:threadError",errString);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs < 2) mexErrMsgIdAndTxt("tiff:inputError","This function requires at least 2 arguments");

    // Check if the fileName is a char array or matlab style
    char* fileName = NULL;
    if(!mxIsClass(prhs[0], "string")){
        if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("tiff:inputError","The first argument must be a string");
        fileName = mxArrayToString(prhs[0]);
    }
    else{ 
        mxArray* mString[1];
        mxArray* mCharA[1];

        // Convert string to char array
        mString[0] = mxDuplicateArray(prhs[0]);
        mexCallMATLAB(1, mCharA, 1, mString, "char");
        fileName = mxArrayToString(mCharA[0]);
    }

    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    if(strchr(fileName,'~')) fileName = expandTilde(fileName);
    #endif

    // Check if folder exists, if not then make it (recursive if needed)
    char* folderName = strdup(fileName);
    char *lastSlash = NULL;
    #ifdef _WIN32
    lastSlash = strrchr(folderName, '\\');
    #else
    lastSlash = strrchr(folderName, '/');
    #endif
    if(lastSlash){
        *lastSlash = '\0';
        FILE* f = fopen(folderName,"r");
        if(f){
            fclose(f);
        }
        else{
            mkdirRecursive(folderName);
        }
    }
    free(folderName);

    TIFFSetWarningHandler(DummyHandler);
    const char* mode;
    if(nrhs > 2){
         mode = mxArrayToString(prhs[2]);
    }
    else{
        mode = "w";
    }
    int nDims = (int) mxGetNumberOfDimensions(prhs[1]);
    if(nDims < 2 || nDims > 3) mexErrMsgIdAndTxt("tiff:inputError","Data must be 2D or 3D");
    uint64_t* dims = (uint64_t*) mxGetDimensions(prhs[1]);



    uint64_t x = dims[1],y = dims[0],z = dims[2],bits = 0, startSlice = 0;

    // For 2D images MATLAB passes in the 3rd dim as 0 so we set it to 1;
    if(!z){
        z = 1;
    }

    mxClassID mDType = mxGetClassID(prhs[1]);
    if(mDType == mxUINT8_CLASS){
        bits = 8;
    }
    else if(mDType == mxUINT16_CLASS){
        bits = 16;
    }
    else if(mDType == mxSINGLE_CLASS){
        bits = 32;
    }
    else if(mDType == mxDOUBLE_CLASS){
        bits = 64;
    }

    //mexErrMsgIdAndTxt("tiff:inputError","TESTING");

    uint64_t stripSize = 512;
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    uint64_t totalStrips = stripsPerDir*z;
    uint64_t* cSizes = (uint64_t*)malloc(totalStrips*sizeof(uint64_t));

    uint64_t dim[3];
    dim[0] = y;
    dim[1] = x;
    dim[2] = z;

    if(bits == 8){
        uint8_t* tiffOld = (uint8_t*)mxGetPr(prhs[1]);
        uint8_t* tiff = (uint8_t*)malloc(x*y*z*(bits/8));

        #pragma omp parallel for collapse(3)
        for(uint64_t dir = 0; dir < z; dir++){
            for(uint64_t j = 0; j < y; j++){
                for(uint64_t i = 0; i < x; i++){
                    ((uint8_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint8_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                }
            }
        }
        writeTiffParallel(x,y,z,fileName, (void*)tiff, (void*)tiffOld, bits, startSlice, stripSize, stripsPerDir, cSizes, mode);
        free(tiff);
    }
    else if(bits == 16){

        uint16_t* tiffOld = (uint16_t*)mxGetPr(prhs[1]);
        uint16_t* tiff = (uint16_t*)malloc(x*y*z*(bits/8));

        #pragma omp parallel for collapse(3)
        for(uint64_t dir = 0; dir < z; dir++){
            for(uint64_t j = 0; j < y; j++){
                for(uint64_t i = 0; i < x; i++){
                    ((uint16_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                }
            }
        }
        writeTiffParallel(x,y,z,fileName, (void*)tiff, (void*)tiffOld, bits, startSlice, stripSize, stripsPerDir, cSizes, mode);
        free(tiff);
    }
    else if(bits == 32){
        float* tiffOld = (float*)mxGetPr(prhs[1]);
        float* tiff = (float*)malloc(x*y*z*(bits/8));

        #pragma omp parallel for collapse(3)
        for(uint64_t dir = 0; dir < z; dir++){
            for(uint64_t j = 0; j < y; j++){
                for(uint64_t i = 0; i < x; i++){
                    ((float*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((float*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                }
            }
        }
        writeTiffParallel(x,y,z,fileName, (void*)tiff, (void*)tiffOld, bits, startSlice, stripSize, stripsPerDir, cSizes, mode);
        free(tiff);
    }
    else if(bits == 64){
        double* tiffOld = (double*)mxGetPr(prhs[1]);
        double* tiff = (double*)malloc(x*y*z*(bits/8));

        #pragma omp parallel for collapse(3)
        for(uint64_t dir = 0; dir < z; dir++){
            for(uint64_t j = 0; j < y; j++){
                for(uint64_t i = 0; i < x; i++){
                    ((double*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((double*)tiffOld)[j+(i*y)+((dir-startSlice)*(x*y))];
                }
            }
        }
        writeTiffParallel(x,y,z,fileName, (void*)tiff, (void*)tiffOld, bits, startSlice, stripSize, stripsPerDir, cSizes, mode);
        free(tiff);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
    free(cSizes);
}
