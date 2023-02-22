#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>

#include "tiffio.h"
#include "omp.h"
#include "mex.h"
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' '-L/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' -ltiff parallelReadTiff.c
//mex COMPFLAGS='$COMPFLAGS /openmp' '-IC:\Program Files (x86)\tiff\include\' '-LC:\Program Files (x86)\tiff\lib\' -ltiffd.lib C:\Users\Matt\Documents\parallelTiff\main.cpp

//libtiff 4.4.0
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" LDOPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/clusterfs/fiona/matthewmueller/software/tiff-4.4.0/include' '-L/clusterfs/fiona/matthewmueller/software/tiff-4.4.0/lib' -ltiff parallelReadTiff.c

// Handle the tilde character in filenames on Linux/Mac
#ifndef _WIN32
#include <wordexp.h>
char* expandTilde(char* path) {
    wordexp_t expPath;
    wordexp(path, &expPath, 0);
    return expPath.we_wordv[0];
}
#endif

void DummyHandler(const char* module, const char* fmt, va_list ap)
{
    // ignore errors and warnings
}

// Backup method in case there are errors reading strips
void readTiffParallelBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif) mexErrMsgIdAndTxt("tiff:threadError","Thread %d: File \"%s\" cannot be opened\n",w,fileName);

        void* buffer = malloc(x*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice) break;

            int counter = 0;
            while(!TIFFSetDirectory(tif, (uint64_t)dir) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,dir,counter+1);
                counter++;
            }

            for (int64_t i = 0; i < y; i++)
            {
                TIFFReadScanline(tif, buffer, i, 0);
                if(!flipXY){
                    memcpy(tiff+((i*x)*bytes),buffer,x*bytes);
                    continue;
                }
                //loading the data into a buffer
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j];
                        }
                            break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j];
                        }
                            break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((float*)buffer)[j];
                        }
                            break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((double*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((double*)buffer)[j];
                        }
                            break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
}

void readTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    uint16_t compressed = 1;
    TIFF* tif = TIFFOpen(fileName, "r");
    TIFFGetField(tif, TIFFTAG_COMPRESSION, &compressed);
    

    

    int32_t w;
    uint8_t errBak = 0;
    uint8_t err = 0;
    char errString[10000];
    if(compressed > 1 || z < 32768){
        TIFFClose(tif);
        #pragma omp parallel for
        for(w = 0; w < numWorkers; w++){

            uint8_t outCounter = 0;
            TIFF* tif = TIFFOpen(fileName, "r");
            while(!tif){
                tif = TIFFOpen(fileName, "r");
                if(outCounter == 3){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                    }
                    continue;
                }
                outCounter++;
            }

            void* buffer = malloc(x*stripSize*bytes);
            for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
                if(dir>=z+startSlice || err) break;

                uint8_t counter = 0;
                while(!TIFFSetDirectory(tif, (uint64_t)dir) && counter<3){
                    counter++;
                    if(counter == 3){
                        #pragma omp critical
                        {
                            err = 1;
                            sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                        }
                    }
                }
                if(err) break;
                for (int64_t i = 0; i*stripSize < y; i++)
                {

                    //loading the data into a buffer
                    int64_t cBytes = TIFFReadEncodedStrip(tif, i, buffer, stripSize*x*bytes);
                    if(cBytes < 0){
                        #pragma omp critical
                        {
                            errBak = 1;
                            err = 1;
                            sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,i);
                        }
                        break;
                    }
                    if(!flipXY){
                        memcpy(tiff+((i*stripSize*x)*bytes),buffer,cBytes);
                        continue;
                    }
                    switch(bits){
                        case 8:
                            // Map Values to flip x and y for MATLAB
                            for(int64_t k = 0; k < stripSize; k++){
                                if((k+(i*stripSize)) >= y) break;
                                for(int64_t j = 0; j < x; j++){
                                    ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j+(k*x)];
                                }
                            }
                                    break;
                        case 16:
                            // Map Values to flip x and y for MATLAB
                            for(int64_t k = 0; k < stripSize; k++){
                                if((k+(i*stripSize)) >= y) break;
                                for(int64_t j = 0; j < x; j++){
                                    ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j+(k*x)];
                                }
                            }
                                    break;
                        case 32:
                            // Map Values to flip x and y for MATLAB
                            for(int64_t k = 0; k < stripSize; k++){
                                if((k+(i*stripSize)) >= y) break;
                                for(int64_t j = 0; j < x; j++){
                                    ((float*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((float*)buffer)[j+(k*x)];
                                }
                            }
                                    break;
                        case 64:
                            // Map Values to flip x and y for MATLAB
                            for(int64_t k = 0; k < stripSize; k++){
                                if((k+(i*stripSize)) >= y) break;
                                for(int64_t j = 0; j < x; j++){
                                    ((double*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((double*)buffer)[j+(k*x)];
                                }
                            }
                                    break;
                    }
                }
            }
            free(buffer);
            TIFFClose(tif);
        }
    }
    else{
        uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
        #ifdef _WIN32
        int fd = open(fileName,O_RDONLY | O_BINARY);
        #else
        int fd = open(fileName,O_RDONLY);
        #endif
        if(fd == -1) mexErrMsgIdAndTxt("disk:threadError","File \"%s\" cannot be opened from Disk\n",fileName);

        if(!tif) mexErrMsgIdAndTxt("tiff:threadError","File \"%s\" cannot be opened\n",fileName);
        uint64_t offset = 0;
        uint64_t* offsets = NULL;
        TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
        uint64_t* byteCounts = NULL;
        TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &byteCounts);
        if(!offsets || !byteCounts) mexErrMsgIdAndTxt("tiff:threadError","Could not get offsets or byte counts from the tiff file\n");
        offset = offsets[0];
        uint64_t fOffset = offsets[stripsPerDir-1]+byteCounts[stripsPerDir-1];
        uint64_t zSize = fOffset-offset;
        TIFFSetDirectory(tif,1);
        TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
        uint64_t gap = offsets[0]-fOffset;
    
        lseek(fd, offset, SEEK_SET);


        TIFFClose(tif);
        uint64_t curr = 0;
        uint64_t bytesRead = 0;
        // TESTING
        // Not sure if we will need to read in chunks like for ImageJ
        for(uint64_t i = 0; i < z; i++){
            bytesRead = read(fd,tiff+curr,zSize);
            curr += bytesRead;
            lseek(fd,gap,SEEK_CUR);
        }
        close(fd);
        uint64_t size = x*y*z*(bits/8);
        void* tiffC = malloc(size);
        memcpy(tiffC,tiff,size);
        #pragma omp parallel for
        for(uint64_t k = 0; k < z; k++){
            for(uint64_t j = 0; j < x; j++){
                for(uint64_t i = 0; i < y; i++){
                    switch(bits){
                        case 8:
                            ((uint8_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint8_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 16:
                            ((uint16_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint16_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 32:
                            ((float*)tiff)[i+(j*y)+(k*x*y)] = ((float*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 64:
                            ((double*)tiff)[i+(j*y)+(k*x*y)] = ((double*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                    }
                }
            }
        }
        free(tiffC);
    }
    if(err){
        if(errBak) readTiffParallelBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else mexErrMsgIdAndTxt("tiff:threadError",errString);
    }
}

// Backup method in case there are errors reading strips
void readTiffParallel2DBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (y-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif) mexErrMsgIdAndTxt("tiff:threadError","Thread %d: File \"%s\" cannot be opened\n",w,fileName);

        void* buffer = malloc(x*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice) break;

            int counter = 0;
            while(!TIFFSetDirectory(tif, (uint64_t)0) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,dir,counter+1);
                counter++;
            }

            for (int64_t i = (w*batchSize); i < ((w+1)*batchSize); i++)
            {
                if(i >= y) break;
                TIFFReadScanline(tif, buffer, i, 0);
                if(!flipXY){
                    memcpy(tiff+((i*x)*bytes),buffer,x*bytes);
                    continue;
                }
                //loading the data into a buffer
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j];
                        }
                            break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j];
                        }
                            break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((float*)buffer)[j];
                        }
                            break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        for(int64_t j = 0; j < x; j++){
                            ((double*)tiff)[((j*y)+i)+((dir-startSlice)*(x*y))] = ((double*)buffer)[j];
                        }
                            break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
}

void readTiffParallel2D(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    int32_t batchSize = (stripsPerDir-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    uint8_t err = 0;
    uint8_t errBak = 0;
    char errString[10000];


    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        uint8_t outCounter = 0;
        TIFF* tif = TIFFOpen(fileName, "r");
        while(!tif){
            tif = TIFFOpen(fileName, "r");
            if(outCounter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                }
                continue;
            }
            outCounter++;
        }

        void* buffer = malloc(x*stripSize*bytes);


        uint8_t counter = 0;
        while(!TIFFSetDirectory(tif, 0) && counter<3){
            printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,0,counter+1);
            counter++;
            if(counter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                }
            }
        }
        for (int64_t i = (w*batchSize); i < (w+1)*batchSize; i++)
        {
            if(i*stripSize >= y || err) break;
            //loading the data into a buffer
            int64_t cBytes = TIFFReadEncodedStrip(tif, i, buffer, stripSize*x*bytes);
            if(cBytes < 0){
                #pragma omp critical
                {
                    errBak = 1;
                    err = 1;
                    sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,i);
                }
                break;
            }
            if(!flipXY){
                memcpy(tiff+((i*stripSize*x)*bytes),buffer,cBytes);
                continue;
            }
            switch(bits){
                case 8:
                    // Map Values to flip x and y for MATLAB
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint8_t*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 16:
                    // Map Values to flip x and y for MATLAB
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint16_t*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 32:
                    // Map Values to flip x and y for MATLAB
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+(k+(i*stripSize)))] = ((float*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 64:
                    // Map Values to flip x and y for MATLAB
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((double*)tiff)[((j*y)+(k+(i*stripSize)))] = ((double*)buffer)[j+(k*x)];
                        }
                    }
                            break;
            }
        }
        free(buffer);
        TIFFClose(tif);
    }

    if(err) {
        if(errBak) readTiffParallel2DBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else mexErrMsgIdAndTxt("tiff:threadError",errString);
    }
}

// Reading images saved by ImageJ
void readTiffParallelImageJ(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    #ifdef _WIN32
    int fd = open(fileName,O_RDONLY | O_BINARY);
    #else
    int fd = open(fileName,O_RDONLY);
    #endif
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) mexErrMsgIdAndTxt("tiff:threadError","File \"%s\" cannot be opened\n",fileName);
    uint64_t offset = 0;
    uint64_t* offsets = NULL;
    TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
    if(offsets) offset = offsets[0];

    TIFFClose(tif);
    lseek(fd, offset, SEEK_SET);
    uint64_t bytes = bits/8;
    //#pragma omp parallel for
    /*
    for(uint64_t i = 0; i < z; i++){
    uint64_t cOffset = x*y*bytes*i;
    //pread(fd,tiff+cOffset,x*y*bytes,offset+cOffset);
    read(fd,tiff+cOffset,x*y*bytes);
    }*/
    uint64_t chunk = 0;
    uint64_t tBytes = x*y*z*bytes;
    uint64_t bytesRead;
    uint64_t rBytes = tBytes;
    if(tBytes < INT_MAX) bytesRead = read(fd,tiff,tBytes);
    else{
        while(chunk < tBytes){
            rBytes = tBytes-chunk;
            if(rBytes > INT_MAX) bytesRead = read(fd,tiff+chunk,INT_MAX);
            else bytesRead = read(fd,tiff+chunk,rBytes);
            chunk += bytesRead;
        }
    }
    close(fd);
    // Swap endianess for types greater than 8 bits
    // TODO: May need to change later because we may not always need to swap
    if(bits > 8){
        #pragma omp parallel for
        for(uint64_t i = 0; i < x*y*z; i++){
            switch(bits){
                case 16:
                    //((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] & 0xff) >> 8) | (((uint16_t*)tiff)[i] << 8);
                    //((uint16_t*)tiff)[i] = bswap_16(((uint16_t*)tiff)[i]);
                    ((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] << 8) & 0xff00) | ((((uint16_t*)tiff)[i] >> 8) & 0x00ff);
                    break;
                case 32:
                    //((num & 0xff000000) >> 24) | ((num & 0x00ff0000) >> 8) | ((num & 0x0000ff00) << 8) | (num << 24)
                    //((float*)tiff)[i] = bswap_32(((float*)tiff)[i]);
                    ((uint32_t*)tiff)[i] = ((((uint32_t*)tiff)[i] << 24) & 0xff000000 ) |
                        ((((uint32_t*)tiff)[i] <<  8) & 0x00ff0000 ) |
                        ((((uint32_t*)tiff)[i] >>  8) & 0x0000ff00 ) |
                        ((((uint32_t*)tiff)[i] >> 24) & 0x000000ff );
                    break;
                case 64:
                    //((double*)tiff)[i] = bswap_64(((double*)tiff)[i]);
                    ((uint64_t*)tiff)[i] = ( (((uint64_t*)tiff)[i] << 56) & 0xff00000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 40) & 0x00ff000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 24) & 0x0000ff0000000000UL ) |
                        ( (((uint64_t*)tiff)[i] <<  8) & 0x000000ff00000000UL ) |
                        ( (((uint64_t*)tiff)[i] >>  8) & 0x00000000ff000000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 24) & 0x0000000000ff0000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 40) & 0x000000000000ff00UL ) |
                        ( (((uint64_t*)tiff)[i] >> 56) & 0x00000000000000ffUL );
                    break;
            }

        }
    }
    // Find a way to do this in-place without making a copy
    if(flipXY){
        uint64_t size = x*y*z*(bits/8);
        void* tiffC = malloc(size);
        memcpy(tiffC,tiff,size);
        #pragma omp parallel for
        for(uint64_t k = 0; k < z; k++){
            for(uint64_t j = 0; j < x; j++){
                for(uint64_t i = 0; i < y; i++){
                    switch(bits){
                        case 8:
                            ((uint8_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint8_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 16:
                            ((uint16_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint16_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 32:
                            ((float*)tiff)[i+(j*y)+(k*x*y)] = ((float*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 64:
                            ((double*)tiff)[i+(j*y)+(k*x*y)] = ((double*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                    }
                }
            }
        }
        free(tiffC);
    }
}

uint8_t isImageJIm(const char* fileName){
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return 0;
    char* tiffDesc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &tiffDesc)){
        if(strstr(tiffDesc, "ImageJ")){
            return 1;
        }
    }
    return 0;
}

uint64_t imageJImGetZ(const char* fileName){
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return 0;
    char* tiffDesc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &tiffDesc)){
        if(strstr(tiffDesc, "ImageJ")){
            char* nZ = strstr(tiffDesc,"images=");
            if(nZ){
                nZ+=7;
                char* temp;
                return strtol(nZ,&temp,10);
            }
        }
    }
    return 0;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
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

    uint8_t flipXY = 1;
    //uint8_t flipXY = 0;


    //if(nrhs > 2){
    //    flipXY = (uint8_t)*(mxGetPr(prhs[2]));
    //}


    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) mexErrMsgIdAndTxt("tiff:inputError","File \"%s\" cannot be opened",fileName);

    uint64_t x = 1,y = 1,z = 1,bits = 1, startSlice = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);

    if(nrhs == 1){
        uint16_t s = 0, m = 0, t = 1;
        while(TIFFSetDirectory(tif,t)){
            s = t;
            t *= 8;
            if(s > t){
                t = 65535;
                printf("Number of slices > 32768\n");
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
        z = s+1;
    }
    else{
        if(mxGetN(prhs[1]) != 2){
            mexErrMsgIdAndTxt("tiff:inputError","Input range is not 2");
        }
        else{
            startSlice = (uint64_t)*(mxGetPr(prhs[1]))-1;
            z = (uint64_t)*((mxGetPr(prhs[1])+1))-startSlice;
            if (!TIFFSetDirectory(tif,startSlice+z-1) || !TIFFSetDirectory(tif,startSlice)){
                mexErrMsgIdAndTxt("tiff:rangeOutOfBound","Range is out of bounds");
            }
        }
    }

    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    uint64_t stripSize = 1;
    TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &stripSize);
    TIFFClose(tif);

    uint8_t imageJIm = 0;
    if(isImageJIm(fileName)){
        imageJIm = 1;
        uint64_t tempZ = imageJImGetZ(fileName);
        if(tempZ) z = tempZ;
    }

    uint64_t dim[3];
    dim[0] = y;
    dim[1] = x;
    dim[2] = z;



    // Case for ImageJ
    if(imageJIm){
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
    // Case for 2D
    else if(z <= 1){
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
    // Case for 3D
    else{
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
}
