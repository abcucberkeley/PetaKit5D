#include <cstdint>
#include <cmath>
#include <cstring>
#include <omp.h>
#include "tiffio.h"
#include "lzwencode.h"
#include "helperfunctions.h"

uint8_t writeTiffSingle(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose){
    TIFF* tif = NULL;
    if(!strcmp(mode,"w")){
        tif = TIFFOpen(fileName, "w8");
        //tif = TIFFOpen(fileName, (cSize < 3.8e9) ? "w" : "w8");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
    }
    else if(!strcmp(mode,"a")){
        tif = TIFFOpen(fileName, "r");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
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
        TIFFClose(tif);
        tif = TIFFOpen(fileName, "a8");
        //tif = TIFFOpen(fileName, (cSize < 3.8e9) ? "a" : "a8");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
    }
    else{
        printf("Error: mode \"%s\" is not supported. Use w or a for mode type", mode);
        return 1;
    }

    uint64_t len = 0;
    //uint64_t lenCounter = 0;
    for(uint64_t dir = startSlice; dir < z; dir++){
        if(dir>=z+startSlice) break;
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, x);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, y);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, stripSize);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
        //TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        //if(compress) TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
        //else TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, 1);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

        if(bits >= 32){
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        }

        for (int64_t i = 0; i*stripSize < y; i++)
        {   
                if (stripSize*(i+1) > y){
                    len = (y-(stripSize*i))*x*(bits/8);
                }
                else{
                    len = stripSize*x*(bits/8);
                }
                /*
                if(transpose){
                    lenCounter = 0;
                    for(uint64_t jT = i*stripSize; jT < y; jT++){
                        if(lenCounter == len) break;
                        for(uint64_t iT = 0; iT < x; iT++){
                            if(lenCounter == len) break;
                            switch(bits){
                                case 8:
                                    ((uint8_t*)tiff)[lenCounter] = ((uint8_t*)tiffOld)[jT+(iT*y)+((dir-startSlice)*(x*y))];
                                    break;
                                case 16:
                                    ((uint16_t*)tiff)[lenCounter] = ((uint16_t*)tiffOld)[jT+(iT*y)+((dir-startSlice)*(x*y))];
                                    break;
                                case 32:
                                    ((float*)tiff)[lenCounter] = ((float*)tiffOld)[jT+(iT*y)+((dir-startSlice)*(x*y))];
                                    break;
                                case 64:
                                    ((double*)tiff)[lenCounter] = ((double*)tiffOld)[jT+(iT*y)+((dir-startSlice)*(x*y))];
                                    break;
                            }
                            lenCounter++;
                        }
                    }
                    TIFFWriteEncodedStrip(tif, i, tiff, len);
                }
                */
                if(transpose) TIFFWriteEncodedStrip(tif, i, (uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
                else TIFFWriteEncodedStrip(tif,i,(uint8_t*)tiffOld+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
                //TIFFWriteRawStrip(tif,i,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
            

        }
        TIFFWriteDirectory(tif);
    }
    TIFFClose(tif);
    return 0;
}


uint8_t writeTiffParallel(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose){
    int32_t numWorkers = omp_get_max_threads();
    
    if(numWorkers == 1){
        return writeTiffSingle(x, y, z, fileName, tiff, tiffOld, bits, startSlice, stripSize, mode, transpose);
    }
    
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    uint64_t totalStrips = stripsPerDir*z;
    uint64_t* cSizes = (uint64_t*)malloc(totalStrips*sizeof(uint64_t));
    
    int32_t batchSize = (z-1)/numWorkers+1;
    uint8_t safeMode = 0;
    uint64_t extraBytes = 2000;

    #pragma omp parallel for
    for(int32_t w = 0; w < numWorkers; w++){

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

                // TESTING: May use this for the parallel version in the future
                // Maybe make a function for the transpose and return instead of breaking?
                /*
                uint64_t lenCounter = 0;
                for(uint64_t jT = i*stripSize; jT < y; jT++){
                    if(lenCounter == len) break;
                    for(uint64_t iT = 0; iT < x; iT++){
                        if(lenCounter == len) break;
                        ((uint16_t*)tiff)[iT+(jT*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)tiffOld)[jT+(iT*y)+((dir-startSlice)*(x*y))];
                        lenCounter++;
                    }
                }
                */

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
            for(uint64_t j = 0; j < y; j++){
                for(uint64_t i = 0; i < x; i++){
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
        comprA = (uint8_t**)malloc(totalStrips*sizeof(uint8_t*));
        #pragma omp parallel for
        for(uint64_t i = 0; i < totalStrips; i++){
            comprA[i] = (uint8_t*)malloc((((x*stripSize)*(bits/8))+(extraBytes*(bits/8)))*2+1);
        }

        #pragma omp parallel for
        for(int32_t w = 0; w < numWorkers; w++){
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
            if(transpose){
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
            else memcpy(tiff,tiffOld,x*y*z*(bits/8));
        }
    }

    // Add bytes for some mandatory tags
    uint64_t manBytes = 24;
    cSize += manBytes*z;

    // Close to UINT32_MAX. Want some extra room for incorrect size calculation for now.
    TIFF* tif = NULL;
    if(!strcmp(mode,"w")){
        tif = TIFFOpen(fileName, (cSize < 3.8e9) ? "w" : "w8");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
    }
    else if(!strcmp(mode,"a")){
        tif = TIFFOpen(fileName, "r");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
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
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
    }
    else{
        printf("Error: mode \"%s\" is not supported. Use w or a for mode type", mode);
        return 1;
    }

    uint64_t len = 0;
    for(uint64_t dir = startSlice; dir < z; dir++){
        if(dir>=z+startSlice) break;
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
    free(cSizes);
    TIFFClose(tif);
	return 0;
}

uint8_t writeTiffParallelWrapper(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* data, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose){
    int32_t numWorkers = omp_get_max_threads();
    void* tiff = (void*)malloc(x*y*z*(bits/8));
    if(transpose){
        // Only use omp if there is more than one thread
        if(numWorkers > 1){
            if(bits == 8){
                #pragma omp parallel for collapse(3)
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((uint8_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint8_t*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }
            }
            else if(bits == 16){
                #pragma omp parallel for collapse(3)
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((uint16_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else if(bits == 32){
                #pragma omp parallel for collapse(3)
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((float*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((float*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else if(bits == 64){
                #pragma omp parallel for collapse(3)
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((double*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((double*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else{
                free(tiff);
                return 1;
            }
        }
        else{
            if(bits == 8){
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((uint8_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint8_t*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }
            }
            else if(bits == 16){
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((uint16_t*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((uint16_t*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else if(bits == 32){
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((float*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((float*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else if(bits == 64){
                for(uint64_t dir = 0; dir < z; dir++){
                    for(uint64_t j = 0; j < y; j++){
                        for(uint64_t i = 0; i < x; i++){
                            ((double*)tiff)[i+(j*x)+((dir-startSlice)*(x*y))] = ((double*)data)[j+(i*y)+((dir-startSlice)*(x*y))];
                        }
                    }
                }        
            }
            else{
                free(tiff);
                return 1;
            }
        }
    }
    else memcpy(tiff,data,x*y*z*(bits/8));
    
    uint8_t ret = writeTiffParallel(x, y, z, fileName, tiff, data, bits, startSlice, stripSize, mode, transpose);
    free(tiff);
    return ret;
}

void writeTiffParallelHelper(const char* fileName, void* tiffOld, uint64_t bits, const char* mode, uint64_t x, uint64_t y, uint64_t z, uint64_t startSlice, uint8_t flipXY)
{
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

	if(!mode) mode = "w";

	// For 2D images MATLAB passes in the 3rd dim as 0 so we set it to 1;
	if(!z){
		z = 1;
	}

	uint64_t stripSize = 512;

	writeTiffParallelWrapper(x,y,z,fileName, tiffOld, bits, startSlice, stripSize, mode, flipXY);
}
