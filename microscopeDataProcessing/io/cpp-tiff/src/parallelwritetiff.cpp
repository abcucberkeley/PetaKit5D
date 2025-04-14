#include <cstdint>
#include <cmath>
#include <cstring>
#include <omp.h>
#include <chrono>
#include <thread>
#include <future>
#include "tiffio.h"
#include "lzwencode.h"
#include "helperfunctions.h"

uint8_t writeTiffSingle(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression){
    TIFF* tif = NULL;
    if(!strcmp(mode,"w")){
        tif = TIFFOpen(fileName, "w8");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
			return 1;
        }
    }
    else if(!strcmp(mode,"a")){
        tif = TIFFOpen(fileName, "a8");
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
    int compressionType = COMPRESSION_NONE;
    if(compression != "none"){
        compressionType = COMPRESSION_LZW;
    }
    for(uint64_t dir = startSlice; dir < z; dir++){
        if(dir>=z+startSlice) break;
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, x);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, y);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, stripSize);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, compressionType);
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
                TIFFWriteEncodedStrip(tif, i, (uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
                //if(transpose) TIFFWriteEncodedStrip(tif, i, (uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
                //else TIFFWriteEncodedStrip(tif,i,(uint8_t*)tiffOld+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
                //TIFFWriteRawStrip(tif,i,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)), len);
        }
        TIFFWriteDirectory(tif);
    }
    TIFFClose(tif);
    return 0;
}

uint8_t writeTiffThread(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* tiff, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const uint64_t stripsPerDir, uint8_t** &comprA, uint64_t* cSizes, const std::string &compression){
    TIFF* tif = NULL;
    if(!strcmp(mode,"w")){
        tif = TIFFOpen(fileName, "w8");
        if(!tif){
            printf("Error: File \"%s\" cannot be opened",fileName);
            return 1;
        }
    }
    else if(!strcmp(mode,"a")){
        tif = TIFFOpen(fileName, "a8");
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
    int compressionType = COMPRESSION_NONE;
    bool compress = false;
    if(compression != "none"){
        compress = true;
        compressionType = COMPRESSION_LZW;
    }
    for(uint64_t dir = startSlice; dir < z; dir++){
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, x);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, y);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, stripSize);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, compressionType);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, 1);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

        if(bits >= 32){
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
        }

        for (int64_t i = 0; i*stripSize < y; i++)
        {
            if(compress){
                while(!cSizes[i+(dir*stripsPerDir)]){
                    std::this_thread::sleep_for(std::chrono::microseconds(1));
                }
                TIFFWriteRawStrip(tif,i,comprA[i+(dir*stripsPerDir)], cSizes[i+(dir*stripsPerDir)]);
                free(comprA[i+(dir*stripsPerDir)]);
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
	free(comprA);
	free(cSizes);
    TIFFClose(tif);
    return 0;
}

uint8_t writeTiffParallel(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression){
    int32_t numWorkers = omp_get_max_threads();
    
    if(numWorkers == 1){
        return writeTiffSingle(x, y, z, fileName, tiff, tiffOld, bits, startSlice, stripSize, mode, transpose, compression);
    }
    
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    uint64_t totalStrips = stripsPerDir*z;
    uint64_t extraBytes = 2000;
    uint8_t** comprA = NULL;
    uint64_t* cSizes = (uint64_t*)calloc(totalStrips, sizeof(uint64_t));

    std::future<uint8_t> writerThreadResult = std::async(std::launch::async, writeTiffThread, x, y, z, fileName, tiff, bits, startSlice, stripSize, mode, stripsPerDir, std::ref(comprA), cSizes, std::ref(compression));
    if(compression != "none"){
        comprA = (uint8_t**)malloc(totalStrips*sizeof(uint8_t*));

        #pragma omp parallel for
        for(uint64_t dir = startSlice; dir < z; dir++){
            uint64_t len = 0;
            for (uint64_t i = 0; i*stripSize < y; i++)
            {
                comprA[i+(dir*stripsPerDir)] = (uint8_t*)malloc((((x*stripSize)*(bits/8))+(extraBytes*(bits/8)))*2+1);
                if (stripSize*(i+1) > y){
                    len = (y-(stripSize*i))*x*(bits/8);
                }
                else{
                    len = stripSize*x*(bits/8);
                }

                if(dir == z-1 && len == (y-(stripSize*i))*x*(bits/8)){
                    uint8_t* cArrL = (uint8_t*)malloc(len+(extraBytes*(bits/8)));
                    memcpy(cArrL,(uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)),len);
                    cSizes[i+(dir*stripsPerDir)] = lzwEncode(cArrL,comprA[i+(dir*stripsPerDir)],len+(extraBytes*(bits/8)));
                    free(cArrL);
                    continue;
                }
                cSizes[i+(dir*stripsPerDir)] = lzwEncode((uint8_t*)tiff+((((i*stripSize)*x)+((dir-startSlice)*(x*y)))*(bits/8)),comprA[i+(dir*stripsPerDir)],len+(extraBytes*(bits/8)));
            }
        }
    }
    return writerThreadResult.get();
}

uint8_t writeTiffParallelWrapper(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* data, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression){
    int32_t numWorkers = omp_get_max_threads();
    void* tiff = nullptr;

    if(transpose){
        tiff = (void*)malloc(x*y*z*(bits/8));
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
    else{
        return writeTiffParallel(x, y, z, fileName, data, data, bits, startSlice, stripSize, mode, transpose, compression);
    }
    
    uint8_t ret = writeTiffParallel(x, y, z, fileName, tiff, data, bits, startSlice, stripSize, mode, transpose, compression);
    free(tiff);
    return ret;
}

void writeTiffParallelHelper(const char* fileName, const void* tiffOld, uint64_t bits, const char* mode, uint64_t x, uint64_t y, uint64_t z, uint64_t startSlice, const bool transpose, const std::string &compression)
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

    writeTiffParallelWrapper(x, y, z, fileName, tiffOld, bits, startSlice, stripSize, mode, transpose, compression);
}
