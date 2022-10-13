#include <stdlib.h>
#include <math.h>
#include "parallelReadTiff.h"
//mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' '-L/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.1.0/libtiff/' -ltiff /clusterfs/fiona/matthewmueller/parallelTiffTesting/main.c


void DummyHandler(const char* module, const char* fmt, va_list ap)
{
    // ignore errors and warnings
}

void readTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;

    uint64_t bytes = bits/8;

    int32_t w;
    uint8_t err = 0;
    char errString[10000];
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif){
            #pragma omp critical
            {
            err = 1;
            sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
            }
        }
        void* buffer = malloc(x*stripSize*bytes);
        for(int64_t dir = startSlice+(w*batchSize); dir < startSlice+((w+1)*batchSize); dir++){
            if(dir>=z+startSlice || err) break;

            uint8_t counter = 0;
            while(!TIFFSetDirectory(tif, (uint64_t)dir) && counter<3){
                printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,dir,counter+1);
                counter++;
            }

            for (int64_t i = 0; i*stripSize < y; i++)
            {

                //loading the data into a buffer
                switch(bits){
                    case 8:
                        // Map Values to flip x and y for MATLAB
                        TIFFReadEncodedStrip(tif, i,(uint8_t*)buffer, stripSize*x*(bits/8));
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint8_t*)buffer)[j+(k*x)];
                            }
                        }
                        break;
                    case 16:
                        // Map Values to flip x and y for MATLAB
                        TIFFReadEncodedStrip(tif, i,(uint16_t*)buffer, stripSize*x*(bits/8));
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((uint16_t*)buffer)[j+(k*x)];
                            }
                        }
                        break;
                    case 32:
                        // Map Values to flip x and y for MATLAB
                        TIFFReadEncodedStrip(tif, i,(float*)buffer, stripSize*x*(bits/8));
                        for(int64_t k = 0; k < stripSize; k++){
                            if((k+(i*stripSize)) >= y) break;
                            for(int64_t j = 0; j < x; j++){
                                ((float*)tiff)[((j*y)+(k+(i*stripSize)))+((dir-startSlice)*(x*y))] = ((float*)buffer)[j+(k*x)];
                            }
                        }
                        break;
                    case 64:
                        // Map Values to flip x and y for MATLAB
                        TIFFReadEncodedStrip(tif, i,(double*)buffer, stripSize*x*(bits/8));
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
    if(err) printf("%s\n", errString);
}

void readTiffParallel2D(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize){
    int32_t numWorkers = omp_get_max_threads();
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    int32_t batchSize = (stripsPerDir-1)/numWorkers+1;

    uint64_t bytes = bits/8;

    int32_t w;
    uint8_t err = 0;
    char errString[10000];


    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){

        TIFF* tif = TIFFOpen(fileName, "r");
        if(!tif){
            #pragma omp critical
            {
                err = 1;
                sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
            }
        }

        void* buffer = malloc(x*stripSize*bytes);


        uint8_t counter = 0;
        while(!TIFFSetDirectory(tif, 0) && counter<3){
            printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,0,counter+1);
            counter++;
        }
        for (int64_t i = (w*batchSize); i < (w+1)*batchSize; i++)
        {
            if(i*stripSize >= y || err) break;
            //loading the data into a buffer
            switch(bits){
                case 8:
                    // Map Values to flip x and y for MATLAB
                    TIFFReadEncodedStrip(tif, i,(uint8_t*)buffer, stripSize*x*(bits/8));
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((uint8_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint8_t*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 16:
                    // Map Values to flip x and y for MATLAB
                    TIFFReadEncodedStrip(tif, i,(uint16_t*)buffer, stripSize*x*(bits/8));
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((uint16_t*)tiff)[((j*y)+(k+(i*stripSize)))] = ((uint16_t*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 32:
                    // Map Values to flip x and y for MATLAB
                    TIFFReadEncodedStrip(tif, i,(float*)buffer, stripSize*x*(bits/8));
                    for(int64_t k = 0; k < stripSize; k++){
                        if((k+(i*stripSize)) >= y) break;
                        for(int64_t j = 0; j < x; j++){
                            ((float*)tiff)[((j*y)+(k+(i*stripSize)))] = ((float*)buffer)[j+(k*x)];
                        }
                    }
                            break;
                case 64:
                    // Map Values to flip x and y for MATLAB
                    TIFFReadEncodedStrip(tif, i,(double*)buffer, stripSize*x*(bits/8));
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
    if(err) printf("%s\n", errString);
}

void* readTiffParallelWrapper(const char* fileName)
{
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return NULL;

    uint64_t x = 1,y = 1,z = 1,bits = 1, startSlice = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);

    uint64_t s = 0, m = 0, t = 1;
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
    z = s+1;

    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    uint64_t stripSize = 1;
    TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &stripSize);
	TIFFClose(tif);

    if(z <= 1){
        if(bits == 8){
            uint8_t* tiff = (uint8_t*)malloc(x*y*z*sizeof(uint8_t));
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 16){
            uint16_t* tiff = (uint16_t*)malloc(x*y*z*sizeof(uint16_t));
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 32){
            float* tiff = (float*)malloc(x*y*z*sizeof(float));
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 64){
            double* tiff = (double*)malloc(x*y*z*sizeof(double));
            readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else{
            return NULL;
        }
    
    }
    else{
        if(bits == 8){
            uint8_t* tiff = (uint8_t*)malloc(x*y*z*sizeof(uint8_t));
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 16){
            uint16_t* tiff = (uint16_t*)malloc(x*y*z*sizeof(uint16_t));
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 32){
            float* tiff = (float*)malloc(x*y*z*sizeof(float));
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else if(bits == 64){
            double* tiff = (double*)malloc(x*y*z*sizeof(double));
            readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize);
            return (void*)tiff;
        }
        else{
            return NULL;
        }
    }

    // Should never get here but return NULL if we do
    return NULL;
}

uint64_t* getImageSize(const char* fileName){

    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) printf("File \"%s\" cannot be opened",fileName);

    uint64_t x = 1,y = 1,z = 1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
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

    TIFFClose(tif);
    uint64_t* dims = (uint64_t*)malloc(3*sizeof(uint64_t));
    dims[0] = y;
    dims[1] = x;
    dims[2] = z;
    return dims;
}

// Returns number of bits the tiff file is.
uint64_t getDataType(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) printf("File \"%s\" cannot be opened",fileName);

    uint64_t bits = 1;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFClose(tif);

    return bits;


}
